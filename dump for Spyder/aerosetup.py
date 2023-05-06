#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#import necessary package and collect initial variables
import math
import json
import numpy as np
import pandas as pd
import prelim_wing_design
import openmdao.api as om
import matplotlib.pyplot as plt
from openaerostruct.utils import plot_wing
from openaerostruct.geometry.utils import generate_mesh
from openaerostruct.geometry.geometry_group import Geometry 
from wing_pressure_loads_calculator import wing_pressure_loads
from openaerostruct. aerodynamics.aero_groups import AeroPoint 
import naca_five_digit_aerofoil_coordinates_calculator as naca_calc
from cruise_conditions_calculator import cruise_conditions_calculator
import ast
# import all_derivatives_calculator

print('aero loaded')
# In[ ]:


#define discipline group and call analyses 
class aero(om.Group):
    def initialize(self):
        
        #read all excel constants        
        initial_vars = pd.read_excel('constants.xlsx')
        
        #sort for just the systems constants
        all_aero_vars = np.where(initial_vars['aero'] == True, [initial_vars['variable name'], initial_vars['value'], initial_vars['description']], None)
        
        #create a DataFrame for the systems variables to pull from in individual analyses
        all_aero_vars = pd.DataFrame(all_aero_vars).dropna(axis = 1).transpose()
        
        #pass variables into the correct openMDAO format, passing name, description and value
        i=0
        for variable in all_aero_vars[0]:
            value = all_aero_vars.iloc[i][1]
            description = all_aero_vars.iloc[i][2]
            self.options.declare(variable, default = value, desc = description)
            i+=1
        i=0
        
    def setup(self):
        
        self.add_subsystem("aero_IDO", aero_IDO(), promotes=["*"])
        #self.add_subsystem('derivatives_calc',derivatives_calc(), promotes = ['*'])
        
    def configure(self):
        #promote all variables (lazy option, they can be connected individually), will connect automatically if the variable names are the same 
        self.promotes('aero_IDO',any=['*'])
        #self.promotes('derivatives_calc',any=['*'])
        


# In[ ]:


class aero_IDO(om.ExplicitComponent):
    def initialize(self):
        #take constants from group a level above
        self.options.declare('naca_series',prob.model.aero.options['naca_series'])
        self.options.declare('max_b',prob.model.aero.options['max_allowable_wing_span'])
        
    def setup(self):
        self.add_input('m')

        self.add_output('taper')
        ##self.add_output('b')
        self.add_output('chord_tip')
        self.add_output('chord_root')
        self.add_output('sweep')
        self.add_output('wing_area')
        self.add_output('aspect_ratio')
        
        self.add_discrete_output('aerodynamic_outputs', val = {"alpha": 0.0, 
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
                                                             "NACA": "23012"})# typical values used for initialisation
        
    def compute(self,inputs,outputs, discrete_inputs, discrete_outputs):
        #start by redefining your variables using inputs and constants (called self.options)
        mass = inputs['m'][0] # Q4J: is this line right or should it follow the same format as the two below?  ***
        naca_series = self.options['naca_series']
        max_allowable_wing_span = self.options['max_b']
        
        #write to file, aero IDO will read them and run and output to a different file
        keyaeronumbers = pd.DataFrame(list([mass,naca_series,max_allowable_wing_span]))
        keyaeronumbers.to_csv("keyAero.dat", index=False) ###
        #in your IDO use the keyAero file
        
        # Run Aero IDO with different problem and setup names so as to not interfere with global MDO problem
        get_ipython().run_line_magic('run', 'aero_ido_vf8.py')
        
        # Reads Aero IDO output file 
        with open("aerodynamic_outputs.dat") as f:
            aerodynamic_outputs = f.read()
        aerodynamic_outputs = ast.literal_eval(aerodynamic_outputs)
        
        discrete_outputs['aerodynamic_outputs'] = aerodynamic_outputs
        
        #add outputs to connect to other disciplines
        outputs['chord_tip'] = aerodynamic_outputs["chord tip"]
        outputs['chord_root'] = aerodynamic_outputs["chord root"]
        ##outputs["b"] = aerodynamic_outputs["span"]
        outputs["sweep"] = aerodynamic_outputs["wing sweep"]
        outputs["wing_area"] = aerodynamic_outputs["S_ref"]
        outputs["aspect_ratio"] = aerodynamic_outputs["S_ref"]/aerodynamic_outputs["span"]
        outputs["taper"] = 0.45


# In[ ]:


'''class derivatives_calc(om.ExplicitComponent):
    def setup(self):
        self.add_output('derivatives')
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
        
    def compute(self,inputs,outputs, discrete_inputs, discrete_outputs):
        
        # Default control surfaces definition simply used to compute derivatives
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
        
        # Runs the derivatives calculator 
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
        outputs['derivatives'] = all_derivatives

        # Saves to a text file all the derivatives in case post-processing is required
        all_derivatives = pd.DataFrame.from_dict(all_derivatives)
        all_derivatives.to_csv("derivatives.dat", index=False)'''

