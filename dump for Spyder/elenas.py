# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 12:16:28 2023

@author: Elena Burcea
"""

import ast
import pandas as pd
from cruise_conditions_calculator import cruise_conditions_calculator

cruise_speed, Re_cruise, cruise_air_density, a_cruise, T_cruise, cruise_dynamic_viscosity = cruise_conditions_calculator()

  
# Imports wing data
with open("aerodynamic_outputs.dat") as f:
    aerodynamic_outputs = f.read()
  
aerodynamic_outputs = ast.literal_eval(aerodynamic_outputs)

root_chord = aerodynamic_outputs["chord root"]
tip_chord = aerodynamic_outputs["chord tip"]
root_twist = aerodynamic_outputs["twist root"]
tip_twist = aerodynamic_outputs["twist tip"]
taper = aerodynamic_outputs["taper"]
wing_span = aerodynamic_outputs["span"]
wing_sweep = aerodynamic_outputs["wing sweep"]
wing_area = aerodynamic_outputs["S_ref"]
naca_series = aerodynamic_outputs["NACA"]
aver_chord = (root_chord+tip_chord)/2
wing_area = aerodynamic_outputs["S_ref"]
CL_wing = aerodynamic_outputs["CL"]
CD_wing = aerodynamic_outputs["CD"]
CM_wing = aerodynamic_outputs["CM"]

Cdalfa0= 0.41457972610255545


# Imports derivatives
derivatives = pd.read_csv("derivatives.dat")

CL_alpha, CD_alpha, CM_alpha = derivatives["wrt pitch"]
CL_psi, CD_psi, CM_psi = derivatives["wrt yaw"]
CL_phi, CD_phi, CM_phi = derivatives["wrt roll"]




#write equations
#compare values of derivatives
#determine if ok or not ok ....


#formulas and calculation
(dalpha_dt)=2; #based on A330 data 
V=238.07;
rho=1.2256;
S=80;
cmq=CM_alpha*aver_chord*(dalpha_dt)/2/V;
xac=aver_chord*0.25;
#xcg-xnp=cmalpha/clalpha
cmu=1.2/(0.5*rho*V*V*S**0.1*1*CM_alpha*aver_chord); # deltae as 0.1 and dcm/ddeltae as -1.2 and xcg-xac as 0

#longitudinal dimensional stability derivatives
#IY-moment of inertia to be received from TANIKA
#the values for qo taken from constant excel, S and c to be taken from tan
S=30;#ranodm
q0=31470.08465;
chord=10.58; #andre
Iy=101;# completely random
CL_deltaE=0.4 #completely random
m=240000;
CD0=0.024135;
CL0=0.683 #calculated by me
C_mq=1 #random
CDu=1 #random value
CD_deltaE=1 #random value
CMT_alpha=1 #random value
CM0=-199.0893 #andre value
CMT_alpha=0; #steady state flight
CL_alphadot=1; #random value
CM_deltaE=1 #random value
Cxtu=1; #random
Cxt0=1; #random
CL_u=1; #random
CL_q=1; #random
Cmtu=1; #random
CmT0=1; #ranodm
CM_alphadot=1; #random

M_alpha=(q0*S*CM_alpha*chord/Iy);
Z_alpha=(-q0*S*(CL_alpha+CD0)/m);
Z_deltaE=(-q0*S*CL_deltaE/m);
X_alpha=(-q0*S*(CD_alpha-CL0))/m;
Mq=(q0*S*(chord*chord)*C_mq)/(2*Iy*V);



Xu= (-q0*S*(CDu+2*CD0))/(m*V);
X_deltaE=(-q0*S*CD_deltaE/m);
M_u=(q0*S*chord*(cmu+2*CM0));
MT_alpha=(q0*S*chord*CMT_alpha)/Iy;
Z_alphadot=(-q0*S*chord*CL_alphadot)/(2*m*V);


M_deltaE=(q0*S*chord*CM_deltaE)/Iy;

XTu=(q0*S*(Cxtu+2*Cxt0))/(m*V);
Z_u=(-q0*S*(CL_u+2*CL0))/(m*V);
Z_q=(q0*S*chord*CL_q)/(2*m*V);
MT_u=(q0*S*chord*(Cmtu+2*CmT0))/(Iy*V);
M_alphadot=(q0*S*(chord*chord)*CM_alphadot)/(2*Iy*V);


#lateral-directional stabitliy derivatives
CYb=1; #random
CY_deltaR=1; #random
Ix=100; #random, to get from Tan
CI_deltaA=1; #random 
b=40; #half-span, assumed to be 80/2
Cnp=1; #ranodm
Iz=100; #random, need from Tan
CYp=1; #random
CI_beta=1; #random
CI_deltaR=1; #random
Cnr=1; #random
Cyr=1; #random
CIp=-1; #random
Cn_beta=1; #random
Cn_deltaA=1; #random
CY_deltaA=1; #random
CIr=1; #random
CnT_beta=1; #random
Cn_deltaR=1; #random

Y_beta=(q0*S*CYb)/m;
Y_deltaR=(q0*S*CY_deltaR)/m;
L_deltaA=(q0*S*b*CI_deltaA)/Ix
Np=(q0*S*b*b*Cnp)/(2*Iz*V);

Yp=(q0*S*b*CYp)/(2*m*V);
L_beta=(q0*S*b*CI_beta)/Ix;
L_deltaR=(q0*S*b*CI_deltaR)/Ix;
Nr=(q0*S*b*b*Cnr)/(Iz*2*V);

Yr=(q0*S*b*Cyr)/(2*m*V);
Lp=(q0*S*b*b*CIp)/(2*Ix*V);
N_beta=(q0*S*b*Cn_beta)/Iz;
N_deltaA=(q0*S*b*Cn_deltaA)/Iz;

Y_deltaA=(q0*S*CY_deltaA)/m;
Lr=(q0*S*b*b*CIr)/(2*Ix*V);
NT_beta=(q0*S*b*CnT_beta)/Iz;
N_deltaR=(q0*S*b*Cn_deltaR)/Iz

CTxu=1; #random value

#criteria for static stability of aircraft

if CM_alpha<0: print("CM_alpha is ok")
if CM_alpha>0: print("CM_alpha is not ok")
if cmq<0: print ("cmq is ok");
if cmq>0: print ("cmq is not ok");
if cmu>0: print ("cmu is ok");
if cmu<0: print ("cmu is not ok");
if CL_alpha>0: print ("CL_aplha is ok");
if CL_alpha<0: print ("CL_aplha is not ok");
if(Cn_beta)>0: print ("Cn_beta is ok");
if(Cn_beta)<0: print ("Cn_beta is not ok");       
if Cnr<0: print ("Cnr is ok");
if Cnr>0: print ("Cnr is not ok");                
if CYb<0: print ("CYb is ok");
if CYb>0: print ("CYb is not ok");   
if (CTxu-CDu)<0: print ("CTxu-CDu is ok");
if (CTxu-CDu)>0: print ("CTxu-CDu is not ok"); 
if (CI_beta)<0: print ("CI_beta is ok");
if CI_beta>0: print ("CI_beta is not ok"); 
if (CIp)<0: print ("CIp is ok");
if CIp>0: print ("CIp is not ok"); 