using CSV
using DataFrames

function get_elementalMatrix(df,components,elements)
    return Matrix(df[in(elements).(df.Element),components])
end

# read in data/elemental_balancing.xlsx
df = CSV.read("data/ElementalMatrix.csv", DataFrame)

# define the components and elements
components  =           ["O2", "CO2", "C6H12O6", "EscherichiaColi"]
measuredComponents =    ["O2", "CO2", "C6H12O6"]
elements    =           ["C", "DoR"]

# get the elemental matrix
elementalMatrix = get_elementalMatrix(df,components,elements)
elementalMatrixMeasured = get_elementalMatrix(df,measuredComponents,elements)

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