#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#import necessary package and collect initial variables
import openmdao.api as om
import ast
import pandas as pd
from cruise_conditions_calculator import cruise_conditions_calculator
from aerosandbox import *
import aerosandbox.numpy as np
import all_derivatives_calculator

print('stability loaded')
# In[ ]:


#define discipline group and call analyses 
class stab(om.Group):
    def initialize(self):
        
        #read all excel constants        
        initial_vars = pd.read_excel('constants.xlsx')
        
        #sort for just the stability constants
        all_stab_vars = np.where(initial_vars['stab'] == True, [initial_vars['variable name'], initial_vars['value'], initial_vars['description']], None)
        
        #create a DataFrame for the systems variables to pull from in individual analyses
        all_stab_vars = pd.DataFrame(all_stab_vars).dropna(axis = 1).transpose()
        
        #pass variables into the correct openMDAO format, passing name, description and value
        i=0
        for variable in all_stab_vars[0]:
            value = all_stab_vars.iloc[i][1]
            description = all_stab_vars.iloc[i][2]
            self.options.declare(variable, default = value, desc = description)
            i+=1
        i=0
        
    def setup(self):
        #set analysis (subsystem) names and what inputs they take and give
        self.add_subsystem('derivatives_calc', derivatives_calc())
    
    def configure(self):
        #promote all variables (lazy option, they can be connected individually)
        self.promotes('derivatives_calc',any=['*'])


# In[ ]:


class transfer_func(om.ExplicitComponent):
    def initialize(self):
        
        #create 2 empty files to write to 
        #pd.write('sthsrth.xls')
        
        #take constants from group a level above
        self.options.declare('v',prob.model.stab.options['V'])
        self.options.declare('rho',prob.model.stab.options['rho'])
        self.options.declare('q0',prob.model.stab.options['q0'])
        self.options.declare('mach',prob.model.stab.options['mach'])
    
    def setup(self):
        self.add_input('chord_root')
        self.add_input('chord_tip')
        self.add_input('twist_root')
        self.add_input('twist_tip')
        self.add_input('taper')
        self.add_input('b') #span
        self.add_input('wing_sweep')
        self.add_input('S_ref')
        
        self.add_input('Cma')
        self.add_input('CL_DE')
        self.add_input('m')
        self.add_input('CD_DE')
        self.add_input('CM0')
        self.add_input('CM_DE')
        self.add_input('Ixz')
        
        self.add_input('Ix')
        self.add_input('Iy')
        self.add_input('Iz')
    
    def compute(self,inputs,outputs):
        df = pd.read_csv("derivatives.dat")

        df = df.to_dict()

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
        derivatives['CD_DE'][0]

        V=238.07;
        rho=1.2256;
        S=80;
        cmq=derivatives['Cmq'][0];
        xac=aver_chord*0.25;
        #xcg-xnp=cmalpha/clalpha
        cmu=1.2/(0.5*rho*V*V*S**0.1*1*derivatives['Cma'][0]*aver_chord); # deltae as 0.1 and dcm/ddeltae as -1.2 and xcg-xac as 0

        #longitudinal dimensional stability derivatives
        #IY-moment of inertia to be received from TANIKA
        #the values for qo taken from constant excel, S and c to be taken from tan
        S=30;#ranodm
        g=9.81; #gravitational acceleration
        q0=31470.08465;
        chord=10.58; #andre
        Iy=0.00959;# Tan
        CL_deltaE=0.4 #completely random
        m=240000;
        CD0=0.024135;
        CL0=0.683 #calculated by me
        # C_mq=1 #random
        CDu=1 #random value
        CD_deltaE=1 #random value
        # CMT_alpha=1 #random value
        CM0=-199.0893 #andre value
        CMT_alpha=0; #steady state flight
        CL_alphadot=1; #random value
        CM_deltaE=1 #random value
        Mach=0.7 #constant
        Cxt0=CD0; #flight mechanics formula
        Cxtu=-2*Cxt0;
        CL_u=1; #random
        # CL_q=1; #random
        Cmtu=1; #random
        CmT0=1; #ranodm
        CM_alphadot=1; #random

        M_alpha=(q0*S*derivatives['Cma'][0]*chord/Iy);
        Z_alpha=(-q0*S*(derivatives['CLa'][0]+CD0)/m);
        Z_deltaE=(-q0*S*CL_deltaE/m); #bad
        X_alpha=(-q0*S*(derivatives['CDa'][0]-CL0))/m;
        Mq=(q0*S*(chord*chord)*derivatives['Cmq'][0])/(2*Iy*V);



        Xu= (-q0*S*(CDu+2*CD0))/(m*V); #bad?
        X_deltaE=(-q0*S*derivatives['CD_DE'][0]/m);
        M_u=(q0*S*chord*(cmu+2*CM0)); 
        MT_alpha=(q0*S*chord*CMT_alpha)/Iy;
        Z_alphadot=(-q0*S*chord*CL_alphadot)/(2*m*V); #to be updated with andres new script


        M_deltaE=(q0*S*chord*CM_deltaE)/Iy;

        XTu=(q0*S*(Cxtu+2*Cxt0))/(m*V);
        Z_u=(-q0*S*(CL_u+2*CL0))/(m*V); #bad
        Z_q=(q0*S*chord*derivatives['CLq'][0])/(2*m*V);
        MT_u=(q0*S*chord*(Cmtu+2*CmT0))/(Iy*V); #bad
        M_alphadot=(q0*S*(chord*chord)*CM_alphadot)/(2*Iy*V);


        M_deltaE=(q0*S*chord*CM_deltaE)/Iy;

        #start trasnfer functions
        #LONGITUDINAL
        #dyanmic longitudinal eq of motions
        theta=5*3.14 /180 #assume theta is 5 degrees and convert to radians- pitch angle
        alpha=5*3.14 /180 #assume alpha is 5 degrees and convert to radians
        deltaE=5*3.14 /180 #assume deltaE is 5 degrees and convert to radians
        q=1; #random???
        alphadot=1; #random
        u=V; #completely random
        u_dot=-g*theta*1+Xu*u+XTu*u+X_alpha*alpha+X_deltaE*deltaE;
        w_dot=V*q-g*theta*0+Z_u*u+Z_alpha*alpha+Z_alphadot*alphadot+Z_q*q+Z_deltaE*deltaE;



        #forward velocity to deflection angle transfer function
        sin_theta=0; #pitch
        cos_theta=1 # pitch assumed as 0

        A_u=X_deltaE*(V-Z_alphadot);
        B_u=-X_deltaE*((V-Z_alphadot)*Mq+Z_alpha+M_alphadot*(V+Z_q))+Z_deltaE*X_alpha;
        C_u=X_deltaE*(Mq*Z_alpha+M_alphadot*g*sin_theta-(M_alpha+MT_alpha)*(V+Z_q))+Z_deltaE*(-M_alphadot*g*cos_theta-X_alpha*Mq)+M_deltaE*(X_alpha*(V+Z_q)-(V-Z_alphadot)*g*cos_theta);
        D_u=X_deltaE*(M_alpha*MT_alpha)*g*sin_theta-Z_deltaE*M_alpha*g*cos_theta+M_deltaE*(Z_alpha*g*cos_theta-Z_alpha*g*sin_theta);
        s=1 # from Laplace transformation, now random value
        N_u=A_u*s*s*s+B_u*s*s+C_u*s+D_u;


        #determining the determinant D1

        A=V-Z_alphadot;
        B=-(V-Z_alphadot)*(Xu+XTu+Mq)-Z_alpha-M_alphadot*(V+Z_q);
        C=(Xu+XTu)*(Mq*(V-Z_alphadot)+Z_alpha+M_alphadot*(V+Z_q))+Mq*Z_alpha-Z_u*X_alpha+M_alphadot*g*sin_theta-(M_alpha+MT_u)*(V+Z_q);
        D=g*sin_theta*(M_alpha+MT_alpha-M_alphadot*(Xu+XTu))+g*cos_theta*(Z_u*M_alphadot+(M_u+MT_u)*(V-Z_alphadot))+(M_u+MT_u)*(-X_alpha*(V+Z_q))+Z_u*X_alpha*Mq+(Xu+XTu)*((M_alpha+MT_alpha)*(V+Z_q)-Mq*Z_alpha);
        E=g*cos_theta*((M_alpha+MT_alpha)*Z_u-Z_alpha*(M_u+MT_u))+g*sin_theta*((M_u+MT_u)*X_alpha-(Xu+XTu)*(M_alpha+MT_alpha));

        poly=np.array([A, B, C, D, E])
        roots=np.roots(poly)
        print("Roots for D1 are", roots)



        #short period
        from math import sqrt
        Wn_SP=sqrt(Z_alpha*Mq/V-M_alpha);
        Zeta_SP=(-(Mq+Z_alpha/V+M_alphadot)/(2*Wn_SP))/10000;
        eta_SP=-Zeta_SP*Wn_SP;
        W_SP=Wn_SP*sqrt(1-Zeta_SP*Zeta_SP);
        period=2*3.14/W_SP;
        t=0.693/(abs(eta_SP));

        #phugoid approximation
        W_nP=sqrt(-Z_u*g/V);
        Zeta_P=-Xu/(2*W_nP);
        W_P=W_nP*sqrt(1-Zeta_P*Zeta_P);
        eta_P=-Zeta_P*W_nP;
        period=2*3.14/W_P*10000;
        t=0.693/(abs(eta_P));

        #plt.plot(Zeta_SP,W_nP, marker='o')
        #plt.show()
        #write to 2 files

        #dynamic lateral stability derivatives
        # CYb=1; #random
        # CY_deltaR=1; #random
        Ix=0.00803; # Tan
        # CI_deltaA=1; #random 
        b=40; #half-span, assumed to be 80/2
        # Cnp=1; #ranodm
        Iz=0.00163; #Tan
        # CYp=1; #random
        # CI_beta=1; #random
        # CI_deltaR=1; #random
        # Cnr=1; #random
        # Cyr=1; #random
        # CIp=-1; #random
        # Cn_beta=1; #random
        # Cn_deltaA=1; #random
        # CY_deltaA=1; #random
        # CIr=1; #random
        CnT_beta=1; #random
        # Cn_deltaR=1; #random

        Y_beta=(q0*S*derivatives['CYb'][0])/m;
        Y_deltaR=(q0*S*derivatives['CYr'][0])/m;
        L_deltaA=(q0*S*b*derivatives['Cla'][0])/Ix 
        Np=(q0*S*b*b*derivatives['Cnp'][0])/(2*Iz*V);

        Yp=(q0*S*b*derivatives['CYp'][0])/(2*m*V);
        L_beta=(q0*S*b*derivatives['Clb'][0])/Ix;
        L_deltaR=(q0*S*b*derivatives['Clr'][0])/Ix; 
        Nr=(q0*S*b*b*derivatives['Cnr'][0])/(Iz*2*V);

        Yr=(q0*S*b*derivatives['CYr'][0])/(2*m*V);
        Lp=(q0*S*b*b*derivatives['Clp'][0])/(2*Ix*V);
        N_beta=(q0*S*b*derivatives['Cnb'][0])/Iz;
        N_deltaA=(q0*S*b*derivatives['Cna'][0])/Iz;

        Y_deltaA=(q0*S*derivatives['CYa'][0])/m;
        Lr=(q0*S*b*b*derivatives['Clr'][0])/(2*Ix*V); 
        NT_beta=(q0*S*b*CnT_beta)/Iz; #bad
        N_deltaR=(q0*S*b*derivatives['Cnr'][0])/Iz

        CTxu=1; #random value


# In[ ]:


class derivatives_calc(om.ExplicitComponent):
    def initialize(self):
        #take constants from group a level above
        self.options.declare('v',prob.model.stab.options['v'])
        self.options.declare('rho',prob.model.stab.options['rho'])
        self.options.declare('q0',prob.model.stab.options['q0'])
        self.options.declare('mach',prob.model.stab.options['mach'])
        
        #carry this on to declare all constants used from constants file
        
    def setup(self):
        self.add_input('chord_root')
        self.add_input('chord_tip')
        self.add_input('twist_root')
        self.add_input('twist_tip')
        self.add_input('taper')
        self.add_input('b') #span
        self.add_input('wing_sweep')
        self.add_input('S_ref')
        
        self.add_input('Cma')
        self.add_input('CL_DE')
        self.add_input('m')
        self.add_input('CD_DE')
        self.add_input('CM0')
        self.add_input('CM_DE')
        self.add_input('Ixz')
        
        self.add_input('Ix')
        self.add_input('Iy')
        self.add_input('Iz')
        
        #add all inputs that will come from other disciplines 
        
        #set the constraints as outputs and then add the constraint on them as follows:
        #self.add_output('Cma')                    ##not sure it can be an input and an output
        #prob.model.add_constraint('Cma', upper = 0)
        
        self.add_output('Cmq')
        #prob.model.add_constraint('Cmq', upper = 0)
        
        self.add_output('cmu')
        #prob.model.add_constraint('cmu', lower = 0)
        
        self.add_output('CLa')
        #prob.model.add_constraint('CLa', lower = 0)
        
        self.add_output('Cnb')
        #prob.model.add_constraint('Cnb', lower = 0)
        
        self.add_output('Cnr')
        #prob.model.add_constraint('Cnr', upper = 0)
        
        self.add_output('CYb')
        #prob.model.add_constraint('CYb', upper = 0)
        
        self.add_output('CTCD')
        #prob.model.add_constraint('CTCD', upper = 0)
        
        self.add_output('Clb')
        #prob.model.add_constraint('Clb', upper = 0)
        
        self.add_output('Clp')
        #prob.model.add_constraint('Clp', upper = 0)
        
        self.add_discrete_input('aerodynamic_outputs',val = {"alpha": 0.0, 
                                                             "chord tip": 8.400740631102945, 
                                                             "chord root": 12.630388118952684, 
                                                             "twist tip": 6.135423187343126, 
                                                             "twist root": 7.557428613406701, 
                                                             "wing sweep": 15.0, 
                                                             "t over c tip": 0.12003453278116566, 
                                                             "t over c root": 0.12003453278116566, 
                                                             "CD": 0.024134617289702367, 
                                                             "CL": 0.6276310000000003, 
                                                             "CM": -199.0892756189604, 
                                                             "span": 80, 
                                                             "taper": 0.5, 
                                                             "S_ref": 672.0592504882356, 
                                                             "NACA": "23012"})
        
    def compute(self,inputs,outputs,discrete_inputs, discrete_outputs):
        
        #start by redefining the constants and inputs if you want to use them in a simpler form:
        V = self.options['v']
        rho = self.options['rho']
        q0 = self.options['q0']
        
        
        root_chord = inputs['chord_root']
        tip_chord = inputs['chord_tip']
        aver_chord = (root_chord+tip_chord)/2
        
        Cdalfa0= 0.41457972610255545
        
        elevator_deflection_angle = 20
        elevator_hinge_position = 0.70
        span_flap_start = 5
        span_flap_end = 15
        flap_deflection_angle = 20
        flap_hinge_position = 0.70
        span_aileron_start = 33
        span_aileron_end = 38
        aileron_deflection_angle = 20
        aileron_hinge_position = 0.7
        
        
        aerodynamic_outputs = discrete_inputs['aerodynamic_outputs']
        
        # Calculates all stability coefficients
        derivatives = all_derivatives_calculator.calculate(aerodynamic_outputs, 
                                                         elevator_deflection_angle, 
                                                         elevator_hinge_position,    
                                                         span_flap_start, 
                                                         span_flap_end,
                                                         flap_deflection_angle,
                                                         flap_hinge_position, 
                                                         span_aileron_start, 
                                                         span_aileron_end,
                                                         aileron_deflection_angle,
                                                         aileron_hinge_position)

        CL_wing = derivatives["CL"]
        CD_wing = derivatives["CD"]
        CM_wing = derivatives["Cm"]
        # # Saves to a text file all the derivatives
        # all_derivatives = pd.DataFrame.from_dict(all_derivatives)
        # all_derivatives.to_csv("derivatives.dat", index=False)



        #write equations
        #compare values of derivatives
        #determine if ok or not ok ....


        #formulas and calculation
        (dalpha_dt)=2; #based on A330 data 
        S=inputs['S_ref']
        
        Cmq=derivatives['Cmq'];
        
        xac=aver_chord*0.25;
        #xcg-xnp=cmalpha/clalpha
        cmu=1.2/(0.5*rho*V*V*S**0.1*1*derivatives['Cma']*aver_chord); # deltae as 0.1 and dcm/ddeltae as -1.2 and xcg-xac as 0

        #longitudinal dimensional stability derivatives
        #IY-moment of inertia to be received from TANIKA
        #the values for qo taken from constant excel, S and c to be taken from tan

        chord=aver_chord; #andre
        Iy = inputs['Iy'];# Tan
        CL_DE=inputs['CL_DE'] #andre
        m=inputs['m'];
        CD0=0.024135;
        CL0=0.683 #calculated by me
        # C_mq=1 #random
        CDu=1 #random value
        CD_DE=inputs['CD_DE']; #andre
        # CMT_alpha=1 #random value
        CM0=inputs['CM0'] #andre value
        CMT_alpha=0; #steady state flight
        CL_alphadot=1; #random value
        CM_DE=inputs['CM_DE']; #andre
        Mach=self.options['mach']
        Cxt0=CD0; #flight mechanics formula
        Cxtu=-2*Cxt0;
        CL_u=1; #random
        # CL_q=1; #random
        Cmtu=1; #random
        CmT0=1; #ranodm
        CM_alphadot=1; #random

        M_alpha=(q0*S*derivatives['Cma']*chord/Iy);
        Z_alpha=(-q0*S*(derivatives['CLa']+CD0)/m);
        Z_deltaE=(-q0*S*CL_DE/m); #bad
        X_alpha=(-q0*S*(derivatives['CDa']-CL0))/m;
        Mq=(q0*S*(chord*chord)*derivatives['Cmq'])/(2*Iy*V);



        Xu= (-q0*S*(CDu+2*CD0))/(m*V); #bad?
        X_deltaE=(-q0*S*CD_DE)/m;
        M_u=(q0*S*chord*(cmu+2*CM0)); 
        MT_alpha=(q0*S*chord*CMT_alpha)/Iy;
        Z_alphadot=(-q0*S*chord*CL_alphadot)/(2*m*V); #to be updated with andres new script


        M_deltaE=(q0*S*chord*CM_DE)/Iy;

        XTu=(q0*S*(Cxtu+2*Cxt0))/(m*V);
        Z_u=(-q0*S*(CL_u+2*CL0))/(m*V); #bad
        Z_q=(q0*S*chord*derivatives['CLq'])/(2*m*V);
        MT_u=(q0*S*chord*(Cmtu+2*CmT0))/(Iy*V); #bad
        M_alphadot=(q0*S*(chord*chord)*CM_alphadot)/(2*Iy*V);


        #lateral-directional stabitliy derivatives
        # CYb=1; #random
        # CY_deltaR=1; #random
        Ix=0.00803; # Tan
        # CI_deltaA=1; #random 
        b=40; #half-span, assumed to be 80/2
        # Cnp=1; #ranodm
        Iz=0.00163; #Tan
        # CYp=1; #random
        # CI_beta=1; #random
        # CI_deltaR=1; #random
        # Cnr=1; #random
        # Cyr=1; #random
        # CIp=-1; #random
        # Cn_beta=1; #random
        # Cn_deltaA=1; #random
        # CY_deltaA=1; #random
        # CIr=1; #random
        CnT_beta=1; #random
        # Cn_deltaR=1; #random

        Y_beta=(q0*S*derivatives['CYb'])/m;
        Y_deltaR=(q0*S*derivatives['CYr'])/m;
        L_deltaA=(q0*S*b*derivatives['Cla'])/Ix 
        Np=(q0*S*b*b*derivatives['Cnp'])/(2*Iz*V);

        Yp=(q0*S*b*derivatives['CYp'])/(2*m*V);
        L_beta=(q0*S*b*derivatives['Clb'])/Ix;
        L_deltaR=(q0*S*b*derivatives['Clr'])/Ix; 
        Nr=(q0*S*b*b*derivatives['Cnr'])/(Iz*2*V);

        Yr=(q0*S*b*derivatives['CYr'])/(2*m*V);
        Lp=(q0*S*b*b*derivatives['Clp'])/(2*Ix*V);
        N_beta=(q0*S*b*derivatives['Cnb'])/Iz;
        N_deltaA=(q0*S*b*derivatives['Cna'])/Iz;
        
        Y_deltaA=(q0*S*derivatives['CYa'])/m;
        Lr=(q0*S*b*b*derivatives['Clr'])/(2*Ix*V); 
        NT_beta=(q0*S*b*CnT_beta)/Iz; #bad
        N_deltaR=(q0*S*b*derivatives['Cnr'])/Iz

        CTxu=1; #random value
        
        
        cmq = 1
        
        #outputs['CM_alpha'] = Cma
        #outputs['Cmq'] = Cma
        #outputs['cmu'] = cmq
        #outputs['CLa'] = cmq
        #outputs['Cnb'] = cmq
        #outputs['Cnr'] = cmq
        #outputs['CYb'] = cmq
        #outputs['CTCD'] = cmq
        #outputs['Clb'] = cmq
        #outputs['Clp'] = cmq


# In[ ]:




