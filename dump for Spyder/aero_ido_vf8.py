# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 11:49:57 2023

@author: Andre Santos
"""

import json
import numpy as np
import pandas as pd
import prelim_wing_design
import openmdao.api as om
from openaerostruct.geometry.utils import generate_mesh
from openaerostruct.geometry.geometry_group import Geometry 
from openaerostruct. aerodynamics.aero_groups import AeroPoint 
import naca_five_digit_aerofoil_coordinates_calculator as naca_calc
from cruise_conditions_calculator import cruise_conditions_calculator

# EXTERNAL INPUTS
keyaeronumbers = pd.read_csv("keyAero.dat")
naca_series = str(int(keyaeronumbers.iloc[1]))
max_allowable_wing_span = float(keyaeronumbers.iloc[2])
mass = float(keyaeronumbers.iloc[0])

# Acquires flight conditions at cruise
cruise_speed, Re_cruise, cruise_air_density, a_cruise, T_cruise, eta_cruise = cruise_conditions_calculator()

# Calculates basic wing characteristics based on the preliminary wing design approach
aver_chord_0 = prelim_wing_design.avg_chord(mass, max_allowable_wing_span)     # [m]
wing_area_target = prelim_wing_design.S_Ref(mass, max_allowable_wing_span)     # [m^2]
cl_cruise_target = prelim_wing_design.CL_cruise(mass, max_allowable_wing_span) # []

# Generates unit-aerofoil coordinates
x_aerofoil, y_aerofoil =  naca_calc.naca_five_digit_aerofoil_coordinates_calculator(naca_series)
upper_x, lower_x = x_aerofoil[round(len(x_aerofoil)/2):], x_aerofoil[:round(len(x_aerofoil)/2)]
upper_y, lower_y = y_aerofoil[round(len(y_aerofoil)/2):], y_aerofoil[:round(len(y_aerofoil)/2)]
upper_x, lower_x = upper_x.astype('float64'), lower_x.astype('float64')
upper_y, lower_y = upper_y.astype('float64'), lower_y.astype('float64')
del x_aerofoil, y_aerofoil

# Calculate's the airfoil's maximum thickness to chord ratio
t_over_c_max = max(upper_y - np.flip(lower_y))/1
c_max_t = (upper_x[upper_y - np.flip(lower_y) == t_over_c_max][0] + lower_x[np.flip(upper_y)-lower_y==t_over_c_max][0])/2

# Create a dictionary to store options about the mesh
# Though it's not necessary, including the initial span and chord speeds up the solution
mesh_dict_aero={
"num_y": 15,
"num_x": 9,
"wing_type": "rect",
"symmetry" : True,
"span": max_allowable_wing_span,
"chord": aver_chord_0,
"chord_cos_spacing": 1,
"span_cos_spacing": 1,
}

# Generates the mesh based on the settings defined above

mesh_aero = generate_mesh(mesh_dict_aero)

# Creates a dictionary with info and options about the wing
surface = {
"name": "wing",
"symmetry": True,
"S_ref_type": "projected",
"fem_model_type": "wingbox",
"data_x_upper": upper_x,
"data_x_lower": lower_x,
"data_y_upper": upper_y,
"data_y_lower": lower_y,
"twist_cp": np.zeros([2]),               # Initializes the twist angles of the wing
"mesh": mesh_aero,
"CM0": 0.0,
"CL0": 0.0,                              # CL of the surface at alpha = 0
"CD0": 0.00632,                          # CD of the surface at alpha = 0  
"k_lam": 0.05,                           # Percentage of chord with laminar flow used for viscous analysis
"t_over_c_cp": np.array([t_over_c_max,
                t_over_c_max]),          # Thickness-to-chord ratio at control points (i.e. root and tip)
"c_max_t": c_max_t,                      # Chordwise location of maximum thickness
"with_viscous" :True,                    # Defines the analysis and viscous/inviscous
"with_wave": False,                      # Defines if the analysis comprises wave drag or not
"sweep": 15,                             # Leading edge sweep angle
"chord_cp": [aver_chord_0,
             aver_chord_0],              # Initializes the chord length at control points (i.e. root and tip)
"taper":np.array([0.45]),                 # Initializes the taper ratio
"span":max_allowable_wing_span           # Defines wingspan, wingtip to wingtip
}

# Creates the OpenMDAO problem
prob_aero = om.Problem()
indep_var_comp = om.IndepVarComp()
indep_var_comp.add_output("v", val=-cruise_speed, units="m/s")                  # Cruise speed relative to outside air
indep_var_comp.add_output ("alpha", val=0, units="deg")                         # Assumed stead-state angle of attack
indep_var_comp.add_output("Mach_number", val=cruise_speed/a_cruise)             # 'True' Mach number at cruise conditions
indep_var_comp.add_output("re", val=Re_cruise, units="1/m")                     # Reynolds number at cruise conditions
indep_var_comp.add_output("rho", val=cruise_air_density, units="kg/m**3")       # Air density at cruise altitude
indep_var_comp.add_output("cg", val=np.zeros(3), units="m")                     # Initialization of the centre of gravity vector
indep_var_comp.add_output("CM", val=0, )

# Adds the IndepVarComp to the problem model
prob_aero.model.add_subsystem("prob_aero_vars", indep_var_comp, promotes=["*"])

# Creates and adds a group that handles the geometry for the use of the aerodynamic lifting surface
geom_group = Geometry(surface=surface)
prob_aero.model.add_subsystem(surface["name"], geom_group)

# Creates the aero point group for the cruise analysis; an instance of `AeroPoint` should be created for each flight condition
aero_group = AeroPoint(surfaces=[surface],user_specified_Sref=True)
point_name = "aero_point_0"
prob_aero.model.add_subsystem(point_name, aero_group, promotes_inputs=["v", "alpha", "Mach_number", "re", "rho", "cg"])

# Connects the mesh from the geometry component to the analysis point
prob_aero.model.connect(surface["name"] + ".mesh", point_name + "." + surface["name"] + ".def_mesh")

# Performs the connections with the modified names within the 'aero_states' group
prob_aero.model.connect(surface["name"] + ".mesh", point_name + ".aero_states." + surface["name"] + "_def_mesh")


# Imports the Scipy Optimizer and set the driver of the problem to use it, which defaults to an SLSQP optimisation method
prob_aero.driver = om.ScipyOptimizeDriver()
prob_aero.driver.options["tol"] = 1e-9
prob_aero.driver.options["debug_print"] = ["nl_cons", "objs", "desvars"]

# Sets up and adds the design variables, the constraints and the objective
prob_aero.model.add_design_var(point_name + ".wing_perf.t_over_c", lower=0, upper=2) # although not exactly a design variable, oas doesn't work without this
prob_aero.model.add_design_var("wing.twist_cp", lower=0, upper=15) 
prob_aero.model.add_design_var("wing.chord_cp", lower=6.5, upper=15, indices=[-1])

# The objective is: MinDrag While Area and CL_cruise requirements are met
prob_aero.model.add_constraint(point_name + ".wing.S_ref", equals=wing_area_target)
prob_aero.model.add_constraint(point_name + ".wing_perf.CL", equals=cl_cruise_target)
prob_aero.model.add_objective(point_name + ".wing_perf.CD", scaler=1e4)

# Sets up and runs the optimisation problem
prob_aero.setup()
prob_aero.run_driver()

aerodynamic_outputs = {
    "alpha": prob_aero["alpha"][0],
    "chord tip": prob_aero["wing.chord_cp"][0],
    "chord root": prob_aero["wing.chord_cp"][1],
    "twist tip": prob_aero["wing.twist_cp"][0],
    "twist root": prob_aero["wing.twist_cp"][1],
    "wing sweep": prob_aero["wing.sweep"][0],
    "t over c tip": prob_aero["wing.t_over_c_cp"][0],
    "t over c root": prob_aero["wing.t_over_c_cp"][1],
    "CD": prob_aero[point_name + ".wing_perf.CD"][0],  
    "CL": prob_aero[point_name + ".wing_perf.CL"][0],
    "CM": prob_aero[point_name + ".total_perf.CM"][1],
    "span":max_allowable_wing_span,
    "taper": prob_aero["wing.taper"][0],
    "S_ref": prob_aero[point_name + ".wing.S_ref"][0],
    "NACA": naca_series
    }

# Saves to a text file the wing characteristics
with open('aerodynamic_outputs.dat', 'w') as convert_file:
	convert_file.write(json.dumps(aerodynamic_outputs))
