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
        self.add_input('sweep')
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
        
        self.add_output('CM_alpha')
        #
        
        self.add_output('cmu')
        #
        
        self.add_output('CLa')
        #
        
        self.add_output('Cnb')
        #
        
        self.add_output('Cnr')
        #
        
        self.add_output('CYb')
        #
        
        self.add_output('CTCD')
        #
        
        self.add_output('Clb')
        #
        
        self.add_output('Clp')
        #
        
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
        CL_alpha = derivatives["CL_alpha_dot"] #, CD_alpha, CM_alpha
        CM_alpha = derivatives["CM_alpha_dot"]

        
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
        
        outputs['CM_alpha'] = CM_alpha
        outputs['Cmq'] = Cmq
        outputs['cmu'] = cmu
        outputs['CLa'] = CL_alpha
        outputs['Cnb'] = derivatives['Cnb']
        outputs['Cnr'] = derivatives['Clr']
        outputs['CYb'] = derivatives['CYb']
        outputs['CTCD'] = CTxu-CDu
        outputs['Clb'] = derivatives['Clb']
        outputs['Clp'] = derivatives['Clp']


# In[ ]:




