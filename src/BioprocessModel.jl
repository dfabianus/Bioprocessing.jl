include("./KineticsLibrary.jl")
include("./BioprocessControllers.jl")
using .KineticsLibrary
import .BioprocessControllers as control
using ModelingToolkit, DifferentialEquations, Plots
using ModelingToolkit.IfElse
using Unitful
using Distributions

##### FINAL INPUTS #####
components = [:X, :S, :O2, :CO2]
balances = [:C, :DOR]
molar_masses = [26.5, 30, 32, 44]
elemental_matrix = [
    1 1 0 1;
    4.113 0 -4 4
]
catalyst = 1
inflows = [:R]
c_inflows = zeros(length(components))
c_inflows[2,1] = 100

@parameters t
Dt = Differential(t)
@variables V(t), [description="reactor volume [L]"]
@variables (m(t)[1:length(components)], [description="mass of component [g]"])
@variables (c(t)[1:length(components)], [description="concentration of component [g/l]"])
@variables (r(t)[1:length(components)], [description="turnover rate of component [g/h]"])
@variables (q(t)[1:length(components)], [description="specific rate of component [g/g/h]"])
@variables (Q(t)[1:length(components)], [description="Incoming massflow of component [g/h]"])
@variables (F(t)[1:length(inflows)], [description="Incoming volume flows [L/h]"])
@parameters (c_F[1:length(components)], [description="Incoming concentrations [g/L]"])
@parameters (M[1:length(components)], [description="Molar mass [g/mol]"])
@parameters (E[1:length(balances), 1:length(components)], [description="Elemental matrix"])

ceq = c .~ m ./ V
req = [r[1] ~ -0.55 * r[2],
       r[2] ~ -1.2 * c[2] / (0.1 + c[2]) * m[1],
       r[3] ~ -1.44 * r[2],
       r[4] ~ 1.2 * r[2],
]
deq = Dt.(m) ~ Q + r
Feq = F .~ 0.005
Qeq = Q .~ F .* c_F
qeq = q .~ r ./ m[catalyst]
Veq = Dt(V) ~ sum(F)
p = [c_F[n] => i for (n,i) in enumerate(c_inflows)]
Meq = [M[n] => i for (n,i) in enumerate(molar_masses)]

# for (m,row) in enumerate(eachrow(elemental_matrix))
#     Eeq[m] = [E[m,n] => i for (n,i) in enumerate(row)]
# end

# Eeq = [E[n,m] => i for (n,i) in enumerate(elemental_matrix)]
#p = c_F => c_inflows
#Meq = M => molar_masses
#Eeq = E => elemental_matrix

V0 = 1.5
x0 = [V => V0; [m[n] => i for (n,i) in enumerate([0.25, 25, 0 ,0])]]

@named sys = ODESystem([ceq; req; deq; Qeq; qeq; Feq; Veq])
sys_simp = structural_simplify(sys)
tspan = (0.0, 11)
prob = ODEProblem(sys_simp, x0, tspan, [p; Meq; Eeq])
sol = solve(prob, TRBDF2())
# Plot the solution
p1 = plot(sol, idxs=q, title = "Bioprocess")

Constants = @constants(
    x_OAIR = 0.2094,    [description="mole fraction of oxygen in air [-]"],
    x_CAIR = 0.0004,    [description="mole fraction of CO2 in air [-]"],
    V_m = 22.414,       [description="molar gas volume at normal conditions [L/mol]"],
)

@parameters(
    k_S,                [tunable = true, description="substrate affinity constant [g/L]"],
    k_O2,               [tunable = true, description="oxygen affinity constant [g/L]"],
    c_SR,               [tunable = false, description="substrate concentration in feed [g/L]"],
    q_Smax,             [tunable = true, description="maximum specific substrate uptake rate [g/g/h]"],
    q_O2max,            [tunable = true, description="maximum specific oxygen uptake rate [g/g/h]"],
    k_Lamax,            [tunable = true, description="maximum oxygen mass transfer coefficient [1/h]"],
    N_Stmax,            [tunable = true, description="maximum stirrer speed [rpm]"],
)

DynamicParameters = @variables( # learn by PINN?
    Y_XS(t),        [description="biomass substrate yield [g/g]"],
    Y_XO2(t),       [description="biomass oxygen yield [g/g]"],
    q_Scap(t),      [description="specific substrate uptake capacity [g/g/h]"],
    q_O2cap(t),     [description="specific oxygen uptake capacity [g/g/h]"],
    k_D(t),         [description="specific death rate [g/g/h]"],
    c_O2sat(t),     [description="oxygen saturation concentration [g/L]"],
    k_La(t),        [description="oxygen mass transfer coefficient [1/h]"],
)

@variables(
    V(t), [description="reactor volume [L]"],
    m_S_met(t), [description="metabolized substrate mass [g]"],
)

StatesMasses = @variables(
    m_X(t), [description="total biomass [g]"],
    m_Xv(t), [description="living biomass [g]"],
    m_Xd(t), [description="dead biomass [g]"],
    m_S(t), [description="substrate mass [g]"],
    m_P(t), [description="product mass [g]"],
    m_O2(t), [description="dissolved oxygen [g]"],
    m_CO2(t), [description="dissolved CO2 [g]"]
)

StatesConcentrations = @variables(
    c_X(t), [description="total biomass concentration [g/L]"],
    c_Xv(t), [description="living biomass concentration [g/L]"],
    c_Xd(t), [description="dead biomass concentration [g/L]"],
    c_S(t), [description="substrate concentration [g/L]"],
    c_P(t), [description="product concentration [g/L]"],
    c_O2(t), [description="dissolved oxygen [g/L]"],
    c_CO2(t), [description="dissolved CO2 [g/L]"]
)

VolumetricRates = @variables(
    r_X(t), [description="growth rate [g/h]"],
    r_Xv(t), [description="growth rate [g/h]"],
    r_Xd(t), [description="death rate [g/h]"],
    r_S(t), [description="substrate uptake rate [g/h]"],
    r_P(t), [description="substrate uptake rate [g/h]"],
    r_O2(t),[description="oxygen uptake rate [g/h]"],
    r_CO2(t),[description="carbon evolution rate [g/h]"],
)

SpecificRates = @variables(
    q_X(t), [description="specific growth rate [g/g/h]"],
    q_Xv(t), [description="specific growth rate [g/g/h]"],
    q_Xd(t), [description="specific death rate [g/g/h]"],
    q_S(t), [description="specific substrate uptake rate [g/g/h]"],
    q_P(t), [description="specific substrate uptake rate [g/g/h]"],
    q_O2(t),[description="specific oxygen uptake rate [g/g/h]"],
    q_CO2(t),[description="specific carbon evolution rate [g/g/h]"],
)

MeasuredOutputs = @variables(
    Q_O2(t), [description="Oxygen supply rate [mol/h]", output = true],
    Q_CO2(t), [description="Carbon exhaust rate [mol/h]", output = true],
    p_O2(t), [description="dissolved oxygen tension [%]", output = true],
    ϑ_L(t), [description="liquid temperature [mol/h]", output = true],
    pH(t), [description="pH value [-]", output = true],
)

# ControlInputs = @variables(
#     F_R(t), [description="feed rate [L/h]", input = true, bounds = (0, Inf)],
#     F_B(t), [description="base titration rate [L/h]", input = true, bounds = (0, Inf)],
#     F_A(t), [description="acid titration rate [L/h]", input = true, bounds = (0, Inf)],
#     N_St(t), [description="stirrer speed [rpm]", input = true, bounds = (0, Inf)],
#     F_AIR(t), [description="Air aeration rate [L/h]", input = true, bounds = (0, Inf)],
#     F_O2(t), [description="ogygen aeration rate [L/h]", input = true, bounds = (0, Inf)],
#     N_St(t), [description="stirrer speed [rpm]", input = true, bounds = (0, Inf)],
# )

ControlInputs = @variables(
    F_R(t), [description="feed rate [L/h]", bounds = (0, Inf)],
    F_B(t), [description="base titration rate [L/h]", bounds = (0, Inf)],
    F_A(t), [description="acid titration rate [L/h]", bounds = (0, Inf)],
    N_St(t), [description="stirrer speed [rpm]", bounds = (0, Inf)],
    F_AIR(t), [description="Air aeration rate [L/h]", bounds = (0, Inf)],
    F_O2(t), [description="ogygen aeration rate [L/h]", bounds = (0, Inf)],
    N_St(t), [description="stirrer speed [rpm]", bounds = (0, Inf)],
)

c = StatesConcentrations .~ StatesMasses ./ V
q = SpecificRates .~ VolumetricRates ./ m_X

M =[q_X q_S; q_O2 q_CO2]

c2 = c[[1,3,4]]
pvar = [
    q_Scap ~ 1.2,
    # q_O2cap ~ 2.4,
    Y_XS ~ 0.55,
    # Y_XO2 ~ 1.44,
    k_D ~ 0.004,
]

eqs = [
    # c_O2sat	~ 100,
    # k_La ~ 0,#k_Lamax * N_St / N_Stmax,
    q_S ~ monod(c_S,q_Scap, k_S),
    # q_O2 ~ monod(c_O2,q_O2cap, k_O2),
    q_X ~ Y_XS * q_S,#minimum([Y_XS * q_S, Y_XO2 * q_O2]),
    # OTR ~ k_La * (c_O2sat - c_O2),
    Dt(V) ~ F_R,
    Dt(m_X) ~ q_X * m_X - k_D * m_X,
    Dt(m_Xd) ~ k_D * m_X,
    Dt(m_S) ~ F_R * c_SR - q_S * m_X,
    Dt(m_S_met) ~ q_S * m_X,
    # Dt(m_O2) ~ OTR - q_O2 * m_X,
]

ctr = [
    F_R ~ control.qS(0.05, m_X, c_SR),
    # N_St ~ 500,
]

x0 = [
    V => 1.5,
    m_X => 0.25,
    m_Xd => 0,
    m_S => 25,
    m_S_met => 0,
    # m_O2 => 0,
]

p = [
    k_S => 0.1,
    # k_O2 => 0.01,
    c_SR => 400,
    # k_Lamax => 600,
    # N_Stmax => 1200,
]

@named sys = ODESystem([eqs; ctr;pvar;c2])
sys_simp = structural_simplify(sys)
tspan = (0.0, 11)
prob = ODEProblem(sys_simp, x0, tspan, p)
sol = solve(prob, TRBDF2())
# Plot the solution
p1 = plot(sol, title = "Bioprocess")
p2 = plot(sol, idxs = q_S, title = "substrate uptake rate")
p3 = plot(sol, idxs = F_R, title = "feed rate")

plot(p1, p2, p3, layout = (2, 2))

# The most simple form of a batch description
# No kinetics, therefore no coupling between substrate availability and metabolism
# just valid when qS = qSmax.
# BATCH_01 = [
#     Dt(m_X) ~ q_Smax * Y_XS * m_X,
#     Dt(m_S) ~ -q_Smax * m_X
# ]

# BATCH_02 = [
#     q_S ~ monod(m_S,q_Smax, k_S),
#     Dt(m_X) ~ q_S * Y_XS * m_X,
#     Dt(m_S) ~ -q_S * m_X
# ]

# FEDBATCH_03 = [
#     F_R ~ q_Sw * m_X / c_SR,
#     q_S ~ monod(m_S,q_Smax, k_S),
#     Dt(V) ~ F_R,
#     Dt(m_X) ~ q_S * Y_XS * m_X,
#     Dt(m_S) ~ F_R * c_SR - q_S * m_X
# ]

# open_loop_qS() = q_Sw * m_X / c_SR
# open_loop_qS(t_start) = IfElse.ifelse(t > t_start, F_R_ffqS(), 0)


# function BIOPROCESS(u=0; name=:model)
#     eqs = [
#         #F_R ~ u,
#         q_S ~ monod(m_S,q_Smax, k_S),
#         Dt(V) ~ F_R,
#         Dt(m_X) ~ q_S * Y_XS * m_X,
#         Dt(m_S) ~ F_R * c_SR - q_S * m_X
#     ]
#     @named sys = ODESystem(eqs)
#     return structural_simplify(sys), sys
# end

###############################################################
###############################################################
# @named sys = ODESystem(FEDBATCH_01)
# simp_model_1 = structural_simplify(sys)
# sys, sys_raw = BIOPROCESS(open_loop_qS(10))
# # Convert from a symbolic to a numerical problem to simulate
# tspan = (0.0, 24)
# prob = ODEProblem(sys, [0.25,25], tspan,[1.2,0.09,0.5])
# prob = ODEProblem(sys, [1,0.25,25], tspan,[1.2,0.2,0.5,400,0.01])

# # Solve the ODE
# sol = solve(prob, TRBDF2())

