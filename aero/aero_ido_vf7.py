# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 11:49:57 2023

@author: Andre Santos
"""

import math
import json
import numpy as np
import pandas as pd
import prelim_wing_design
import openmdao.api as om
import derivatives_calculator
import wing_section_properties
import matplotlib.pyplot as plt
import all_derivatives_calculator
from openaerostruct.utils import plot_wing
from openaerostruct.geometry.utils import generate_mesh
from openaerostruct.geometry.geometry_group import Geometry 
from wing_pressure_loads_calculator import wing_pressure_loads
from openaerostruct. aerodynamics.aero_groups import AeroPoint 
import naca_five_digit_aerofoil_coordinates_calculator as naca_calc
from cruise_conditions_calculator import cruise_conditions_calculator

# EXTERNAL INPUTS
naca_series = '23012'
max_allowable_wing_span = 80 # [m]
mass = 155600 # [kg]

constants_file_path = "D:/University/Year 4/MECH5080M Team Project/1 - Aerodynamics Works/Aero_IDO/constants.xlsx"
mass = pd.read_excel(constants_file_path)
mass = mass[mass["variable name"]=="M_initial"]["value"][12]

# House Keeping
to_plot_wing = False # If False avoids openning oas plotting window for wing
to_plot_aerofoils = False # If False does not plot tip and root aerofoils

# Calculates basic wing characteristics based on the preliminary wing design approach
aver_chord_0 = prelim_wing_design.avg_chord(mass, max_allowable_wing_span)     # [m]
wing_area_target = prelim_wing_design.S_Ref(mass, max_allowable_wing_span)     # [m^2]
cl_cruise_target = prelim_wing_design.CL_cruise(mass, max_allowable_wing_span) # []

# Acquires flight conditions at cruise
cruise_speed, Re_cruise, cruise_air_density, a_cruise, T_cruise, eta_cruise = cruise_conditions_calculator()

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
mesh_dict={
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
mesh = generate_mesh(mesh_dict)


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
"mesh": mesh,
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
"taper":np.array([0.5]),                 # Initializes the taper ratio
"span":max_allowable_wing_span           # Defines wingspan, wingtip to wingtip
}

# Creates the OpenMDAO problem
prob = om.Problem()
indep_var_comp = om.IndepVarComp()
indep_var_comp.add_output("v", val=-cruise_speed, units="m/s")                  # Cruise speed relative to outside air
indep_var_comp.add_output ("alpha", val=0, units="deg")                         # Assumed stead-state angle of attack
indep_var_comp.add_output("Mach_number", val=cruise_speed/a_cruise)             # 'True' Mach number at cruise conditions
indep_var_comp.add_output("re", val=Re_cruise, units="1/m")                     # Reynolds number at cruise conditions
indep_var_comp.add_output("rho", val=cruise_air_density, units="kg/m**3")       # Air density at cruise altitude
indep_var_comp.add_output("cg", val=np.zeros(3), units="m")                     # Initialization of the centre of gravity vector
indep_var_comp.add_output("CM", val=0, )

# Adds the IndepVarComp to the problem model
prob.model.add_subsystem("prob_vars", indep_var_comp, promotes=["*"])

# Creates and adds a group that handles the geometry for the use of the aerodynamic lifting surface
geom_group = Geometry(surface=surface)
prob.model.add_subsystem(surface["name"], geom_group)

# Creates the aero point group for the cruise analysis; an instance of `AeroPoint` should be created for each flight condition
aero_group = AeroPoint(surfaces=[surface],user_specified_Sref=True)
point_name = "aero_point_0"
prob.model.add_subsystem(point_name, aero_group, promotes_inputs=["v", "alpha", "Mach_number", "re", "rho", "cg"])

# Connects the mesh from the geometry component to the analysis point
prob.model.connect(surface["name"] + ".mesh", point_name + "." + surface["name"] + ".def_mesh")

# Performs the connections with the modified names within the 'aero_states' group
prob.model.connect(surface["name"] + ".mesh", point_name + ".aero_states." + surface["name"] + "_def_mesh")


# Imports the Scipy Optimizer and set the driver of the problem to use it, which defaults to an SLSQP optimisation method
prob.driver = om.ScipyOptimizeDriver()
prob.driver.options["tol"] = 1e-9
prob.driver.options["debug_print"] = ["nl_cons", "objs", "desvars"]

if to_plot_wing:
    recorder = om.SqliteRecorder("aero_ido.db")
    prob.driver.add_recorder(recorder)
    prob.driver.recording_options["record_derivatives"] = True
    prob.driver.recording_options["includes"] = ["*"]


# Sets up and adds the design variables, the constraints and the objective
# prob.model.add_design_var("wing.sweep", lower=0, upper=15) # Wing sweep was a design variable, but because it always tended to 15deg, it was fixed to 15deg
# prob.model.add_design_var("wing.span", lower=20, upper=80) # Span would be a design variable, but for the same reason as sweep, it was fixed to the max allowable
# prob.model.add_design_var("alpha", lower=-10, upper=10)  # The only reason this line of code was left here is that it is a legacy from the version of 2021/22
prob.model.add_design_var(point_name + ".wing_perf.t_over_c", lower=0, upper=2) # although not exactly a design variable, oas doesn't work without this
prob.model.add_design_var("wing.twist_cp", lower=0, upper=15) 
prob.model.add_design_var("wing.chord_cp", lower=6.5, upper=15, indices=[-1])
# The objective is: MinDrag While Area and CL_cruise requirements are met
prob.model.add_constraint(point_name + ".wing.S_ref", equals=wing_area_target)
prob.model.add_constraint(point_name + ".wing_perf.CL", equals=cl_cruise_target)
prob.model.add_objective(point_name + ".wing_perf.CD", scaler=1e4)

# Sets up and runs the optimisation problem
prob.setup()

prob.run_driver()

print("CD = " + str(prob[point_name + ".wing_perf.CD"][0]))
print("CL = " + str(prob[point_name + ".wing_perf.CL"][0]))

aerodynamic_outputs = {
    "alpha": prob["alpha"][0],
    "chord tip": prob["wing.chord_cp"][0],
    "chord root": prob["wing.chord_cp"][1],
    "twist tip": prob["wing.twist_cp"][0],
    "twist root": prob["wing.twist_cp"][1],
    "wing sweep": prob["wing.sweep"][0],
    "t over c tip": prob["wing.t_over_c_cp"][0],
    "t over c root": prob["wing.t_over_c_cp"][1],
    "CD": prob[point_name + ".wing_perf.CD"][0],  
    "CL": prob[point_name + ".wing_perf.CL"][0],
    "CM": prob[point_name + ".total_perf.CM"][1],
    "span":max_allowable_wing_span,
    "taper": prob["wing.taper"][0],
    "S_ref": prob[point_name + ".wing.S_ref"][0],
    "NACA": naca_series
    }

# Saves to a text file the wing characteristics
with open('aerodynamic_outputs.dat', 'w') as convert_file:
	convert_file.write(json.dumps(aerodynamic_outputs))


# Displays final wing
if to_plot_wing:
    print("Displaying wing plot...")
    args = [[],[]]
    args[1] = "aero_ido.db"
    plot_wing.disp_plot(args=args)


# Creates tip aerofoil
twist_tip = math.radians(aerodynamic_outputs["twist tip"])

# Re-scales the aerofoil
x_tip_scaled_upper = (aerodynamic_outputs["chord tip"]*upper_x)
y_tip_scaled_upper = (aerodynamic_outputs["chord tip"]*upper_y)
x_tip_scaled_lower = (aerodynamic_outputs["chord tip"]*lower_x)
y_tip_scaled_lower = (aerodynamic_outputs["chord tip"]*lower_y)

# Rotates the aerofoil to match the twist angle at the tip
x_tip_upper=(lambda x, y: x*math.cos(twist_tip) - y*math.sin(twist_tip))(x_tip_scaled_upper,y_tip_scaled_upper)
y_tip_upper=(lambda x, y: x*math.sin(twist_tip) + y*math.cos(twist_tip))(x_tip_scaled_upper,y_tip_scaled_upper)
x_tip_lower=(lambda x, y: x*math.cos(twist_tip) - y*math.sin(twist_tip))(x_tip_scaled_lower,y_tip_scaled_lower)
y_tip_lower=(lambda x, y: x*math.sin(twist_tip) + y*math.cos(twist_tip))(x_tip_scaled_lower,y_tip_scaled_lower)

# Plots the tip aerofoil
if to_plot_aerofoils:
    plt.plot(x_tip_upper, y_tip_upper)
    plt.plot(x_tip_lower, y_tip_lower)
    plt.title('Tip Aerofoil')
    plt.axis('equal')
    plt.show()


# Creates root aerofoil
twist_root = math.radians(aerodynamic_outputs["twist root"])

# Re-scales the aerofoil
x_root_scaled_upper = (aerodynamic_outputs["chord root"]*upper_x)
y_root_scaled_upper = (aerodynamic_outputs["chord root"]*upper_y)
x_root_scaled_lower = (aerodynamic_outputs["chord root"]*lower_x)
y_root_scaled_lower = (aerodynamic_outputs["chord root"]*lower_y)

# Rotates the aerofoil to match the twist angle at the tip
x_root_upper=(lambda x, y: x*math.cos(twist_root) - y*math.sin(twist_root))(x_root_scaled_upper,y_root_scaled_upper)
y_root_upper=(lambda x, y: x*math.sin(twist_root) + y*math.cos(twist_root))(x_root_scaled_upper,y_root_scaled_upper)
x_root_lower=(lambda x, y: x*math.cos(twist_root) - y*math.sin(twist_root))(x_root_scaled_lower,y_root_scaled_lower)
y_root_lower=(lambda x, y: x*math.sin(twist_root) + y*math.cos(twist_root))(x_root_scaled_lower,y_root_scaled_lower)

# Plots the tip aerofoil
if to_plot_aerofoils:
    plt.plot(x_root_upper, y_root_upper)
    plt.plot(x_root_lower, y_root_lower)
    plt.title('Root Aerofoil')
    plt.axis('equal')
    plt.show()

    
# Calculates the 2nd moments of area for the tip and root aerofoils
Ixx_tip, Iyy_tip, Ixy_tip = wing_section_properties.InertiaMomentCalculator(x_tip_upper, y_tip_upper, x_tip_lower, y_tip_lower)
Ixx_root, Iyy_root, Ixy_root = wing_section_properties.InertiaMomentCalculator(x_root_upper, y_root_upper, x_root_lower, y_root_lower)

# Exports the 2nd moments of area to a text file
pd.DataFrame({"Aerofoil":["Tip","Root"], "Ixx m^4":[Ixx_tip,Ixx_root], "Iyy m^4":[Iyy_tip,Iyy_root], "Ixy m^4":[Ixy_tip,Ixy_root]}).to_csv("Ixx_Iyy_Ixy.dat", index=False)


# Calculates the "pressures" exerted on the tip and root aerofoils
upper_tip_pressure, lower_tip_pressure = wing_pressure_loads(aerodynamic_outputs["CL"], aerodynamic_outputs["CD"],  x_tip_upper, y_tip_upper, x_tip_lower, y_tip_lower)
upper_root_pressure, lower_root_pressure = wing_pressure_loads(aerodynamic_outputs["CL"], aerodynamic_outputs["CD"], x_root_upper, y_root_upper, x_root_lower, y_root_lower)

# Calculates the global pressures exerted on the wing's upper and lower surfaces
upper_surface_pressure = (upper_tip_pressure + upper_root_pressure)/2
lower_surface_pressure = (lower_tip_pressure + lower_root_pressure)/2

# Export these pressures to a .dat file
pd.DataFrame({"upper_surface (Pa)":[upper_surface_pressure], "lower_surface (Pa)":[lower_surface_pressure]}).to_csv("wing_loads.dat", index=False)


# Exports the final wing's geometry
x_coords = prob["aero_point_0.wing.def_mesh"][:,:,0].flatten()
y_coords = prob["aero_point_0.wing.def_mesh"][:,:,1].flatten()
wing_coords = pd.DataFrame({"x":x_coords, "y":y_coords})
# Changes the dataset to solely include the vertices of the wing
wing_coords = pd.DataFrame({"x":[min(wing_coords[wing_coords['y']==min(wing_coords['y'])]['x']), 
                                    max(wing_coords[wing_coords['y']==min(wing_coords['y'])]['x']),
                                    min(wing_coords[wing_coords['y']==max(wing_coords['y'])]['x']),
                                    max(wing_coords[wing_coords['y']==min(wing_coords['y'])]['x'])],
                           "y":[min(wing_coords['y']),
                                    min(wing_coords['y']),
                                    max(wing_coords['y']),
                                    max(wing_coords['y'])]})
wing_coords = wing_coords.apply(lambda coord: round(coord / 1e-6) * 1e-6)
wing_coords.to_csv('wing_coordinates.dat', index=False)


# Selects information on control surfaces
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

# Calculates all stability coefficients
all_derivatives = all_derivatives_calculator.calculate(aerodynamic_outputs, 
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

# Saves to a text file all the derivatives
all_derivatives = pd.DataFrame.from_dict(all_derivatives)
all_derivatives.to_csv("derivatives.dat", index=False)


print("END: Aero IDO finished running !")