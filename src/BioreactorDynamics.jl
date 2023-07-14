import YAML
C = YAML.load_file("data/naturalConstants.yml")

change(in, out, transfer, reaction) = in + out + transfer + reaction

x_OGin(F_AIR,F_O2) = (F_AIR.*C.xOAIR.value+F_O2)./(F_AIR+F_O2)

x_CGin(F_AIR,F_O2) = (F_AIR.*C.xCAIR.value)./(F_AIR+F_O2)

EXH2O(x_WET=0.2094)=1-(x_WET./C.xOAIR.value)

INERT(x_OGin,x_CGin,x_O2,x_CO2;EXH2O=EXH2O()) = (1 .-x_OGin.-x_CGin)./(1 .-x_O2.-x_CO2.-EXH2O) 

Q_CO2(F_AIR,F_O2,x_CO2,x_CGin,INERT) = (F_AIR+F_O2)./C.VnM.value.*(x_CO2.*INERT-x_CGin) # mol/h

Q_O2(F_AIR,F_O2,x_O2,x_OGin,INERT) = (F_AIR+F_O2)./C.VnM.value.*(x_O2.*INERT-x_OGin) # mol/h

Q_S(F_R,c_SR,M_S) = -F_R .* c_SR ./ M_S

# Add other dynamics here