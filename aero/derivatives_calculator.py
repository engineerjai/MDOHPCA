# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 13:35:04 2023

@author: Andre Santos
"""

"""
This script includes all the functions used to calcualte aerodynamic stability 
derivatives for the stability anasyles.


"""

def CL_CD_CM_alpha(aerodynamic_outputs):
    """
       Calculates the stability derivivates that relate the lift, drag and 
       pitching moment coefficients to the angle of attack of the wing; 
       .i.e. d(CL)/dAoA, d(CL)/dAoA and d(CM)/dAoA.
    
    Args:
        aerodynamic_outputs       (list)   []
            aerodynamic_outputs["chord root"]     (float)    [m]
            aerodynamic_outputs["chord tip"]      (float)    [m]
            aerodynamic_outputs["twist root"]     (float)    [deg]
            aerodynamic_outputs["twist tip"]      (float)    [deg]
            aerodynamic_outputs["taper"]          (float)    []
            aerodynamic_outputs["span"]           (float)    [m]
            aerodynamic_outputs["winf sweep"]     (float)    [deg]
            aerodynamic_outputs["CL"]             (float)    []
            aerodynamic_outputs["CD"]             (float)    []
            aerodynamic_outputs["CM"]             (float)    []
            aerodynamic_outputs["S_ref"]          (float)    [m^2]
            aerodynamic_outputs["NACA"]           (string)   []
            
    Returns:
        CL_alpha                  (float)  []
        CD_alpha                  (float)  []
        CM_alpha                  (float)  []

    """
        
    import numpy as np
    import openmdao.api as om
    from openaerostruct.geometry.utils import generate_mesh
    from openaerostruct.geometry.geometry_group import Geometry 
    from openaerostruct. aerodynamics.aero_groups import AeroPoint 
    import naca_five_digit_aerofoil_coordinates_calculator as naca_calc
    from cruise_conditions_calculator import cruise_conditions_calculator
    
    # Creates variables with wing characteristics
    root_chord = aerodynamic_outputs["chord root"]
    tip_chord = aerodynamic_outputs["chord tip"]
    root_twist = aerodynamic_outputs["twist root"]
    tip_twist = aerodynamic_outputs["twist tip"]
    taper = aerodynamic_outputs["taper"]
    wing_span = aerodynamic_outputs["span"]
    wing_sweep = aerodynamic_outputs["wing sweep"]
    CL_a0 = aerodynamic_outputs["CL"]
    CD_a0 = aerodynamic_outputs["CD"]
    CM_a0 = aerodynamic_outputs["CM"]
    wing_area = aerodynamic_outputs["S_ref"]
    naca_series = aerodynamic_outputs["NACA"]
    aver_chord = (root_chord+tip_chord)/2
    
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
    
    # Creates a dictionary to store options about the mesh
    # Though it's not necessary, including the initial span and chord speeds up the solution
    mesh_dict={
    "num_y": 15,
    "num_x": 9,
    "wing_type": "rect",
    "symmetry" : True,
    "span": wing_span,
    "chord": aver_chord,
    "chord_cos_spacing": 1,
    "span_cos_spacing": 1,
    }
    
    # Generates the mesh based on the settings defined above
    mesh = generate_mesh(mesh_dict)
    
    
    # Creates a dictionary with info and options about the wing
    surface = {
    "name":"wing",
    "symmetry":True,
    "S_ref_type":"projected",
    "S_ref":wing_area,
    "fem_model_type":"wingbox",
    "data_x_upper":upper_x,
    "data_x_lower":lower_x,
    "data_y_upper":upper_y,
    "data_y_lower":lower_y,
    "twist_cp":np.array([tip_twist,
                          root_twist]),      # Initializes the twist angles of the wing
    "mesh":mesh,
    "CL0":CL_a0,                              # CL of the surface at alpha = 0
    "CD0":CD_a0,                          # CD of the surface at alpha = 0  
    "k_lam":0.05,                           # Percentage of chord with laminar flow used for viscous analysis
    "t_over_c_cp":np.array([t_over_c_max,
                    t_over_c_max]),          # Thickness-to-chord ratio at control points (i.e. root and tip)
    "c_max_t":c_max_t,                      # Chordwise location of maximum thickness
    "with_viscous":True,                    # Defines the analysis and viscous/inviscous
    "with_wave":False,                      # Defines if the analysis comprises wave drag or not
    "sweep":wing_sweep,                     # Leading edge sweep angle
    "chord_cp":[tip_chord,
                 root_chord],                # Initializes the chord length at control points (i.e. root and tip)
    "taper":np.array([taper]),                 # Initializes the taper ratio
    "span":wing_span                         # Defines wingspan, wingtip to wingtip
    }

    
    # Creates the OpenMDAO problem
    prob_der = om.Problem()
    indep_var_comp = om.IndepVarComp()
    indep_var_comp.add_output("v", val=-cruise_speed, units="m/s")                  # Cruise speed relative to outside air
    indep_var_comp.add_output ("alpha", val=-10, units="deg")                      # In this case, negative mean positively upwards
    indep_var_comp.add_output("Mach_number", val=cruise_speed/a_cruise)             # 'True' Mach number at cruise conditions
    indep_var_comp.add_output("re", val=Re_cruise, units="1/m")                     # Reynolds number at cruise conditions
    indep_var_comp.add_output("rho", val=cruise_air_density, units="kg/m**3")       # Air density at cruise altitude
    indep_var_comp.add_output("cg", val=np.zeros(3), units="m")                     # Initialization of the centre of gravity vector
    
    # Adds the IndepVarComp to the prob_derlem model
    prob_der.model.add_subsystem("prob_der_vars", indep_var_comp, promotes=["*"])
    
    # Creates and adds a group that handles the geometry for the use of the aerodynamic lifting surface
    geom_group = Geometry(surface=surface)
    prob_der.model.add_subsystem(surface["name"], geom_group)
    
    # Creates the aero point group for the cruise analysis; an instance of `AeroPoint` should be created for each flight condition
    aero_group = AeroPoint(surfaces=[surface],user_specified_Sref=True)
    point_name = "aero_point_0"
    prob_der.model.add_subsystem(point_name, aero_group, promotes_inputs=["v", "alpha", "Mach_number", "re", "rho", "cg"])
    
    # Connects the mesh from the geometry component to the analysis point
    prob_der.model.connect(surface["name"] + ".mesh", point_name + "." + surface["name"] + ".def_mesh")
    
    # Performs the connections with the modified names within the 'aero_states' group
    prob_der.model.connect(surface["name"] + ".mesh", point_name + ".aero_states." + surface["name"] + "_def_mesh")
    
    # Imports the Scipy Optimizer and set the driver of the prob_derlem to use it, which defaults to an SLSQP optimisation method
    prob_der.driver = om.ScipyOptimizeDriver()
    prob_der.driver.options["tol"] = 1e-9
    prob_der.driver.options["debug_print"] = ["nl_cons", "objs", "desvars"]
    
    # Sets up and adds the design variables, the constraints and the objective
    # Sets up a 'phathom' prob_derlem whereby the variable result will always be the upper end of the range
    # In this case, it's computing CL for AoA=10deg.
    prob_der.model.add_design_var(point_name + ".wing_perf.t_over_c", lower=0, upper=2) # although not exactly a design variable, oas doesn't work without this
    prob_der.model.add_constraint(point_name + ".wing.S_ref", lower=0, upper=wing_area, units="m**2")
    prob_der.model.add_design_var("alpha", lower=10, upper=30, units="deg")
    prob_der.model.add_objective(point_name + ".wing_perf.CD", scaler=1e4)

    # Sets up and runs the optimisation prob_derlem
    prob_der.setup()
    prob_der.run_driver()
    
    ### Uncomment for debugging purposes
    # print("CL = " + str(prob_der[point_name + ".wing_perf.CL"][0]))
    # print("CD = " + str(prob_der[point_name + ".wing_perf.CD"][0]))
    # print("CM = " + str(prob_der[point_name + ".total_perf.CM"][0]))
    # print("CM = " + str(prob_der[point_name + ".total_perf.CM"][1]))
    # print("CM = " + str(prob_der[point_name + ".total_perf.CM"][2]))
    # print("AoA = " + str(prob_der["alpha"][0]))
    
    CL_a10 = prob_der[point_name + ".wing_perf.CL"][0]
    CD_a10 = prob_der[point_name + ".wing_perf.CD"][0]
    CM_a10 = prob_der[point_name + ".total_perf.CM"][1]
    
    
    CL_alpha = (CL_a10 - CL_a0)/(10/180*np.pi)           # [per radian]
    CD_alpha = (CD_a10 - CD_a0)/(10/180*np.pi)           # [per radian]
    CM_alpha = (CM_a10 - CM_a0)/(10/180*np.pi)           # [per radian]

    
    return CL_alpha, CD_alpha, CM_alpha
    
    
    
    
    
    