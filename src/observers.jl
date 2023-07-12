module Observers
export x_OGin, x_CGin, EXH2O, INERT, Q_CO2, Q_O2, Q_S, r, RQ, θ, convert_units, kalman, kalman_vec, kalman_state_derivative, kalman_flow_rate

using DataFrames 
using LinearAlgebra
using Interpolations
using HampelFilter
using GLM
using CairoMakie

######## KALMAN FILTERING ######################################
function kalman(zₖ, xₖ = [0.0 0.0]', Pₖ = [10 0; 0 10]; Δt = 1.0, Q = [0 0; 0 0.01], R = 0.04)
    A = [1.0 Δt; 0.0 1.0]
    C = [1.0 0.0]
    # State prediction
    xₖ = A * xₖ
    Pₖ = A * Pₖ * A' + Q
    # Kalman Gain
    K = Pₖ * C'* inv(C * Pₖ * C' .+ R)

    # Output / estimate
    x_hat = xₖ .+ K * (zₖ .- C * xₖ)
        
    # Covariance error
    P_hat = Pₖ .- K * C * Pₖ

    return x_hat, P_hat
end
function kalman_vec(t, z, xₖ = [0 0]', Pₖ = [10 0; 0 10]; Q = [0 0; 0 0.01], R = 0.04) 
    if length(z) < 2
        return kalman(z, xₖ, Pₖ; Δt=t, Q=Q, R=R)
    else
        x = xₖ
        P = vec(Pₖ)
        Δt = diff(t,dims=1)
        for (Δtₖ,zₖ) in zip(Δt, z[2:end])
            xₖ, Pₖ = kalman(zₖ, xₖ, Pₖ; Δt=Δtₖ, Q=Q, R=R)
            x = [x xₖ]
            P = [P vec(Pₖ)]
        end
        return x,P
    end
end
function kalman_state_derivative(t, z, xₖ = [0 0]', Pₖ = [10 0; 0 10]; Q = [0 0; 0 0.01], R = 0.04) 
    x,P = kalman_vec(t, z, xₖ, Pₖ; Q = Q, R = R) 
    state = x[1,:]
    derivative = x[2,:]
    return state, derivative, x, P
end
function kalman_flow_rate(t, z, density = 1000, init_vol = 1.5, xₖ = [0 0]', Pₖ = [10 0; 0 10]; Q = [0 0; 0 0.01], R = 0.04) 
    state, derivative, x, P = kalman_state_derivative(t, z, xₖ, Pₖ; Q = Q, R = R) 
    return (
        volume = state ./ density .+ init_vol, 
        flow_rate = -derivative ./ density,
        x = x,
        P = P) 
end
function F_R(y::DataFrame, weight_signal_name::String="m_R"; density = 1000, init_vol = 1.5, xₖ = [0 0]', Pₖ = [10 0; 0 10], Q = [0 0; 0 0.01], R = 0.04) 
    # This is the feed rate calculation method using kalman filtering.
    # More methods can be implemented based on other approaches e.g. moving averages
    return kalman_flow_rate(y.time, y[!,weight_signal_name], density, init_vol, xₖ, Pₖ; Q = Q, R = R).flow_rate
end
function V_L(y::DataFrame, weight_signal_name::String="m_L", density = 1000, init_vol = 1.5, xₖ = [0 0]', Pₖ = [10 0; 0 10]; Q = [0 0; 0 0.01], R = 0.04) 
    # This is the feed rate calculation method using kalman filtering.
    # More methods can be implemented based on other approaches e.g. moving averages
    return kalman_flow_rate(y.time, y[!,weight_signal_name], density, init_vol, xₖ, Pₖ; Q = Q, R = R).volume
end

######## BIOPROCESS MODELS ####################################
x_OGin(F_AIR,F_O2;x_OAIR=0.2094) = (F_AIR.*x_OAIR+F_O2)./(F_AIR+F_O2)
x_OGin(y::DataFrame;x_OAIR=0.2094) = x_OGin(y.F_AIR,y.F_O2;x_OAIR=x_OAIR)

x_CGin(F_AIR,F_O2;x_CAIR=0.0004) = (F_AIR.*x_CAIR)./(F_AIR+F_O2)
x_CGin(y::DataFrame;x_CAIR=0.0004) = x_CGin(y.F_AIR,y.F_O2,x_CAIR=x_CAIR)

EXH2O(x_WET=0.2094;x_OAIR=0.2094)=1-(x_WET./x_OAIR)

INERT(x_OGin,x_CGin,x_O2,x_CO2;EXH2O=EXH2O()) = (1 .-x_OGin.-x_CGin)./(1 .-x_O2.-x_CO2.-EXH2O) 
INERT(y::DataFrame;EXH2O=EXH2O()) = INERT(x_OGin(y),x_CGin(y),y.x_O2,y.x_CO2;EXH2O=EXH2O) 

Q_CO2(F_AIR,F_O2,x_CO2,x_CGin,INERT;V_M=22.414) = (F_AIR+F_O2)./V_M.*(x_CO2.*INERT-x_CGin) # mol/h
Q_CO2(y::DataFrame;V_M=22.414,EXH2O=EXH2O()) = Q_CO2(y.F_AIR,y.F_O2,y.x_CO2,x_CGin(y),INERT(y,EXH2O=EXH2O)) # mol/h

Q_O2(F_AIR,F_O2,x_O2,x_OGin,INERT;V_M=22.414) = (F_AIR+F_O2)./V_M.*(x_O2.*INERT-x_OGin) # mol/h
Q_O2(y::DataFrame;V_M=22.414,EXH2O=EXH2O()) = Q_O2(y.F_AIR,y.F_O2,y.x_O2,x_OGin(y),INERT(y,EXH2O=EXH2O);V_M=22.414) # mol/h

Q_S(F_R,c_SR,M_S) = -F_R .* c_SR ./ M_S
Q_S(y::DataFrame,c_SR,M_S) = y.F_R .* c_SR ./ M_S

r(q,m_x) = q .* m_x

RQ(CER,OUR) = OUR ./ CER
RQ(y::DataFrame) = -y.Q_O2 ./ y.Q_CO2

function c2m(y_offline::DataFrame, y_online::DataFrame, c) 
    V_L_interp = LinearInterpolation(y_online.time, y_online.V_L,extrapolation_bc=Line())
    y_offline[!,"m_"*c] = y_offline[!,c].*V_L_interp(y_offline.time)
    return y_offline
end
m2c(y_offline::DataFrame, y_online::DataFrame, m) = y[!,m].*V_L

function θ(df,f,args...;kwargs...) 
    df[!,Symbol(f)] = f(df,args...;kwargs...) 
    return df
end

function convert_units(df,unit_ref_dict::Dict)
    [df[!,key] = df[!,key].*value for (key,value) in unit_ref_dict]
    return df
end

######## ELEMENTAL BALANCING: K2S1 ##################################

Q(y::DataFrame) = [
    y.Q_S, #substrate supply rate
    y.Q_O2, # oxygen supply rate
    y.Q_CO2 # carbon exhaust rate
]

function K2S1(rm, carbon, gamma, sigma)
    # Q is a vector of the supply rates and exhaust rates measured at the bioreactor
    # Units are mol/h
    # gamma includes [gamma_c, gamma_m]. 
    #r_m = Q

    Em = [carbon[2] carbon[3] carbon[4]
        gamma[2] gamma[3] gamma[4]] 

    Ec = [carbon[1]; gamma[1]]

    #Ec_star: Moore Penrose pseudo inverse of Ec
    Ec_star=(inv(Ec'*Ec))*Ec'

    # R: Redundancy matrix
    R=Em-Ec*Ec_star*Em

    # Rred: Reduced redundancy matrix (containing just the independent rows of R
    U,S,V=svd(R)
    Sconv=[1 0]
    C=Sconv*S
    K=C*S'*U'
    Rred=K*R

    # eps: residual vector
    eps = Rred * rm

    # P: Residual variance covariance matrix
    P = Rred * sigma *Rred'

    # Reconciliation of measured and calculated rates
    delta = (sigma*Rred'*inv(P) * Rred)* rm
    rm_best = rm-delta
    xc_best = -Ec_star*Em*rm_best

    # Sum of weighted squares of residuals
    h = eps' * inv(P) * eps

    # Calculate the function outputs
    return (r=vcat(xc_best, rm_best), h=h)
end

function K2S1_model(V_L, x, Qm, Qc=0, Mm=[30,44,32], Mc=[26.5], carbon=[1,1,1,0], gamma = [4.113, 4, 0, -4], sigma=Diagonal([0.03,0.03,0.03]); kinetics=monod_kinetics, dynamics=bioreactor_continuous)
    rm = kinetics(x,Qm.*Mm, V_L)
    r_hat, h = K2S1(rm./Mm, carbon, gamma, sigma)
    Q = vcat(Qc,Qm)
    M = vcat(Mc,Mm)
    dx = dynamics(x,Q.*M,r_hat.*M)
    return dx, r_hat, h
end

function K2S1_fitting_problem(V_L, x, Qm, Qc=0, Mm=[30,44,32], Mc=[26.5], carbon=[1,1,1,0], gamma = [4.113, 4, 0, -4], sigma=Diagonal([0.03,0.03,0.03]); kinetics=monod_kinetics_special, dynamics=bioreactor_continuous)
    dx = K2S1_model(V_L, x, Qm, Qc, Mm, Mc, carbon, gamma, sigma; kinetics, dynamics)
    prob = ODEProblem(f2,u0,tspan,p)
end

function K2S1_dynamical(t_S, V_L, x, Qm, Qc=0, Mm=[30,44,32], Mc=[26.5], carbon=[1,1,1,0], gamma = [4.113, 4, 0, -4], sigma=Diagonal([0.03,0.03,0.03]); kinetics=monod_kinetics, dynamics=bioreactor)
    rm = kinetics(x,Qm.*Mm, V_L)
    r_hat, h = K2S1(rm./Mm, carbon, gamma, sigma)
    if x[1] < 10
        #r_hat[1] = -0.37 * r_hat[2]
        r_hat[1] = -0.5 * r_hat[2]
    end
    Q = vcat(Qc,Qm)
    M = vcat(Mc,Mm)
    x = dynamics(x,Q.*M,r_hat.*M,t_S)
    return x, r_hat, h
end
function K2S1_dynamical_monod(t_S, V_L, x, Qm, Qc=0, Mm=[30,44,32], Mc=[26.5], carbon=[1,1,1,0], gamma = [4.113, 4, 0, -4], sigma=Diagonal([0.03,0.03,0.03]); dynamics=bioreactor, qSmax=1.25, qSmax2=0.25, kS=0.01)
    rm = monod_kinetics_special(x,Qm.*Mm, V_L; q_Smax=qSmax, qSmax_2 = qSmax2, k_S=kS)
    r_hat, h = K2S1(rm./Mm, carbon, gamma, sigma)
    if x[1] < 12
        #r_hat[1] = -0.37 * r_hat[2]
        #r_hat[1] = -0.7 * r_hat[2]
    end
    Q = vcat(Qc,Qm)
    M = vcat(Mc,Mm)
    x = dynamics(x,Q.*M,r_hat.*M,t_S)
    return x, r_hat, h
end
function monod_kinetics_special(x,Qm,V_L; q_Smax=1.25, qSmax_2 = 0.25, k_S=0.01)
    # from Qm build rm
    m_X, m_S, _, _ = x
    Q_S, Q_CO2, Q_O2 = Qm
    if m_X > 20
        q_Smax = qSmax_2
    end
    q_S = q_Smax .* (m_S./V_L) ./ (k_S .+ (m_S./V_L))
    r_S = -q_S .* m_X
    return [r_S, Q_CO2, Q_O2]
end

function K2S1_vec(t, V_L, xₖ, Qm, Qc=0, Mm=[30,44,32], Mc=[26.5], carbon=[1,1,1,0], gamma = [4.113, 4, 0, -4], sigma=Diagonal([0.03,0.03,0.03]); kinetics=monod_kinetics, dynamics=bioreactor)
    x = xₖ
    r = zeros(size(xₖ))
    h = 0
    Δt = diff(t,dims=1)
    for (Δtₖ,Qmₖ, V_Lk) in zip(Δt, eachrow(Qm[2:end,:]), V_L)
        xₖ, r_hatₖ, hₖ = K2S1_dynamical(Δtₖ, V_Lk, xₖ, Qmₖ, Qc, Mm, Mc, carbon, gamma, sigma; kinetics=kinetics, dynamics=dynamics)
        x = [x xₖ]
        r = [r r_hatₖ]
        h = [h hₖ]
    end
    return x,r,h
end
function K2S1_vec_monod(t, V_L, xₖ, Qm, Qc=0, Mm=[30,44,32], Mc=[26.5], carbon=[1,1,1,0], gamma = [4.113, 4, 0, -4], sigma=Diagonal([0.03,0.03,0.03]); dynamics=bioreactor, qSmax=1.25, qSmax2=0.25, kS=0.01)
    x = xₖ
    r = zeros(size(xₖ))
    h = 0
    Δt = diff(t,dims=1)
    for (Δtₖ,Qmₖ, V_Lk) in zip(Δt, eachrow(Qm[2:end,:]), V_L)
        xₖ, r_hatₖ, hₖ = K2S1_dynamical_monod(Δtₖ, V_Lk, xₖ, Qmₖ, Qc, Mm, Mc, carbon, gamma, sigma; dynamics=dynamics, qSmax=qSmax, qSmax2=qSmax2, kS=kS)
        x = [x xₖ]
        r = [r r_hatₖ]
        h = [h hₖ]
    end
    return x,r,h
end

function K2S1_vec(y::DataFrame, xₖ, Qc=0, Mm=[30,44,32], Mc=[26.5], carbon=[1,1,1,0], gamma =  [4.113, 4, 0, -4], sigma=Diagonal([0.03,0.03,0.03]); kinetics=monod_kinetics, dynamics=bioreactor)
    Qm = Matrix(y[:,[:Q_S, :Q_CO2, :Q_O2]])
    x,r,h = K2S1_vec(y.time,  y.V_L, xₖ, Qm, Qc, Mm, Mc, carbon, gamma, sigma; kinetics=kinetics, dynamics=dynamics)
    y.mX_K2S1 = [xi[1] for xi in eachcol(x)]
    y.mS_K2S1 = [xi[2] for xi in eachcol(x)]
    y.rX_K2S1 = [ri[1] for ri in eachcol(r)]
    y.rS_K2S1 = [ri[2] for ri in eachcol(r)]
    y.rC_K2S1 = [ri[3] for ri in eachcol(r)]
    y.rO_K2S1 = [ri[4] for ri in eachcol(r)]
    #y.h_K2S1 = h
    return y
end
function K2S1(y::DataFrame, carbon=[1,1,1,0], gamma = [4, -4, 0, 4.113], sigma=Diagonal([0.03,0.03,0.03]))
    Q = y[:,[:Q_S, :Q_CO2, :Q_O2]]
    y.rX_K2S1 = [K2S1([-S,O2,CO2],carbon,gamma,sigma).r[1] for (S,O2,CO2) in zip(Q[:,1], Q[:,2], Q[:,3])]
    return y
end
function monod_kinetics(x,Qm,V_L; q_Smax=1.2,k_S=0.4)
    # from Qm build rm
    m_X, m_S, _, _ = x
    Q_S, Q_CO2, Q_O2 = Qm
    q_S = q_Smax .* (m_S./V_L) ./ (k_S .+ (m_S./V_L))
    r_S = -q_S .* m_X
    return [r_S, Q_CO2, Q_O2]
end
function no_kinetics(x,Qm,V_L)
    Q_S, Q_CO2, Q_O2 = Qm
    r_S = -Q_S
    return [r_S, Q_CO2, Q_O2]
end
function bioreactor(x,Q,r,t_S,F_out=0)
    dx = bioreactor_continuous(x,Q,r,F_out)
    x = x .+ t_S .* dx
    return x
end
function bioreactor_continuous(x,Q,r,F_out=0)
    dx = [Qᵢ+rᵢ-xᵢ*F_out for (xᵢ,Qᵢ,rᵢ) in zip(x,Q,r)]
    return dx
end

##### ARCHIVE #######
mol2g(x, M) = x .* M


struct ReactorMeasurement
    # This struct can later be extended with measurements.jl
    # later push to another module (maybe ReactorMeasurements or so)
    F_AIR::Float64      # L/h
    F_O2::Float64       # L/h
    m_R::Float64        # L/h
    x_O2::Float64       # -
    x_CO2::Float64      # -
    m_L::Float64        # L
end
function create_reactorMeasurementVector(F_AIR,F_O2,m_R,x_O2,x_CO2,V_L)
    return [ReactorMeasurement(F_AIRi,F_O2i,m_Ri,x_O2i,x_CO2i,V_Li) 
    for (F_AIRi,F_O2i,m_Ri,x_O2i,x_CO2i,V_Li) in zip(F_AIR,F_O2,m_R,x_O2,x_CO2,V_L)]
end

function hampfil(y::DataFrame, signal_name="m_R",k=11,thr=3)
    return hampel(y[!,signal_name],k,thr)
end

function XOD(y::DataFrame; plot=true)
    XOD_model = lm(@formula(DCW ~ OD), y)
    y.XOD = predict(XOD_model)
    return y, XOD_model
end


##### K2S1 #####
struct K2S1_param
    E
    Sigma
end
end