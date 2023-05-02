import ast
from aerosandbox import *
import aerosandbox.numpy as np
from cruise_conditions_calculator import cruise_conditions_calculator


def calculate(aerodynamic_outputs, 
              elevator_deflection_angle, 
              elevator_hinge_position, 
              span_flap_start, 
              span_flap_end,
              flap_deflection_angle,
              flap_hinge_position, 
              span_aileron_start, 
              span_aileron_end,
              aileron_deflection_angle,
              aileron_hinge_position,
              cruise_altitude=18288):
    
    
    SLF_output_data = SLF_aero_analysis(aerodynamic_outputs,cruise_altitude)
    CL_alpha_dot, CM_alpha_dot = CL_CM_alpha_dot_calculator(aerodynamic_outputs, SLF_output_data, cruise_altitude)
    CLu, CDu, CMu = CLu_CDu_CMu_calculator(aerodynamic_outputs, SLF_output_data, cruise_altitude)
    CL_DE, CD_DE, CM_DE = elevator_derivatives(aerodynamic_outputs, elevator_deflection_angle, elevator_hinge_position, SLF_output_data, cruise_altitude)
    CL_DF, CD_DF, CM_DF = flaps_derivatives(aerodynamic_outputs, span_flap_start, span_flap_end, flap_deflection_angle, flap_hinge_position, SLF_output_data, cruise_altitude)
    CL_DA, CD_DA, CM_DA = flaps_derivatives(aerodynamic_outputs, span_aileron_start, span_aileron_end, aileron_deflection_angle, aileron_hinge_position, SLF_output_data, cruise_altitude)
    all_derivatives = {**SLF_output_data,
                       **{"CL_alpha_dot":CL_alpha_dot,
                          "CM_alpha_dot":CM_alpha_dot,
                          "CLu":CLu,
                          "CDu":CDu,
                          "CMu":CMu,
                          "CL_DE":CL_DE, 
                          "CD_DE":CD_DE, 
                          "CM_DE":CM_DE,
                          "CL_DF":CL_DF, 
                          "CD_DF":CD_DF, 
                          "CM_DF":CM_DF,
                          "CL_DA":CL_DA, 
                          "CD_DA":CD_DA, 
                          "CM_DA":CM_DA}}
    
    return all_derivatives


def SLF_aero_analysis(aerodynamic_outputs, cruise_altitude=18288):
    
    # Calculates cruise conditions  
    cruise_speed, Re_cruise, cruise_air_density, a_cruise, T_cruise, eta_cruise = cruise_conditions_calculator(cruise_altitude=18288)
    
    # Uncomment the following block of code to read directly from directory
    # # Imports wing data
    # with open("aerodynamic_outputs.dat") as f:
    #     aerodynamic_outputs = f.read()
    # aerodynamic_outputs = ast.literal_eval(aerodynamic_outputs)
    
    root_chord = aerodynamic_outputs["chord root"]
    tip_chord = aerodynamic_outputs["chord tip"]
    root_twist = aerodynamic_outputs["twist root"]
    tip_twist = aerodynamic_outputs["twist tip"]
    wing_span = aerodynamic_outputs["span"]
    wing_sweep = aerodynamic_outputs["wing sweep"]
    naca_series = aerodynamic_outputs["NACA"]
     
    Jeta_1 = Airplane(
        name="JETA",
        xyz_ref=[0, 0, 0], # CG location
        wings=[
            Wing(
                name="Main Wing",
                xyz_le=[0, 0, 0], # Coordinates of the wing's leading edge
                symmetric=True,
                xsecs=[ # The wing's cross ("X") sections
                    WingXSec(  # Root
                        xyz_le=[0, 0, 0], # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
                        chord=root_chord,
                        twist=root_twist, # degrees
                        airfoil=Airfoil(name="naca"+naca_series),
                    ),
                    WingXSec(  # Tip
                        xyz_le=[np.tan(np.radians(wing_sweep))*wing_span/2, wing_span/2, 0],
                        chord=tip_chord,
                        twist=tip_twist,
                        airfoil=Airfoil(name="naca"+naca_series),
                    )
                ]
            ),
            Wing(
                name="Horizontal Stabilizer",
                xyz_le=[40, 0, 2],
                symmetric=True,
                xsecs=[
                    WingXSec(  # root
                        xyz_le=[40, 0, 2],
                        chord=9,
                        twist=5,
                        airfoil=Airfoil(name="naca0012"),
                        control_surface_type='symmetric',  # Elevator
                        control_surface_deflection=0,
                        control_surface_hinge_point=0.75
                    ),
                    WingXSec(  # tip
                        xyz_le=[44, 11, 2],
                        chord=5,
                        twist=5,
                        airfoil=Airfoil(name="naca0012")
                    )
                ]
            ),
            Wing(
                name="Vertical Stabilizer",
                xyz_le=[40, 0, 2],
                symmetric=False,
                xsecs=[
                    WingXSec(
                        xyz_le=[40, 0, 2],
                        chord=9,
                        twist=0,
                        airfoil=Airfoil(name="naca0012"),
                        control_surface_type='symmetric',  # Rudder
                        control_surface_deflection=0,
                        control_surface_hinge_point=0.75
                    ),
                    WingXSec(
                        xyz_le=[45, 0, 13],
                        chord=4,
                        twist=0,
                        airfoil=Airfoil(name="naca0012")
                    )
                ]
            )
        ]
    )
    
    aero_problem = VortexLatticeMethod(
        airplane=Jeta_1,
        op_point=OperatingPoint(
            velocity=cruise_speed,
            alpha=0,
            beta=0,
            p=0,
            q=0,
            r=0,
            atmosphere=Atmosphere(altitude=cruise_altitude)
        ),
    )
    
    SLF_output_data = aero_problem.run_with_stability_derivatives(alpha=True, beta=True, p=True, q=True, r=True) # Runs and prints results to console
    
    return SLF_output_data
    

def flaps_derivatives(aerodynamic_outputs, span_flap_start, span_flap_end, flap_deflection_angle, flap_hinge_position, SLF_output_data, cruise_altitude=18288): 
    # Calculates cruise conditions  
    cruise_speed, Re_cruise, cruise_air_density, a_cruise, T_cruise, eta_cruise = cruise_conditions_calculator()
    
   # Uncomment the following block of code to read directly from directory
   # # Imports wing data
   # with open("aerodynamic_outputs.dat") as f:
   #     aerodynamic_outputs = f.read()
   # aerodynamic_outputs = ast.literal_eval(aerodynamic_outputs)
   
    
    root_chord = aerodynamic_outputs["chord root"]
    tip_chord = aerodynamic_outputs["chord tip"]
    root_twist = aerodynamic_outputs["twist root"]
    tip_twist = aerodynamic_outputs["twist tip"]
    wing_span = aerodynamic_outputs["span"]
    wing_sweep = aerodynamic_outputs["wing sweep"]
    naca_series = aerodynamic_outputs["NACA"]
     
    # Control surfaces position calculation
    a_W = wing_span*np.tan(np.radians(wing_sweep))/2+tip_chord # [m] distance used in intermediary calcultions
    alpha_W = np.arctan((a_W-root_chord)/wing_span*2) # [rad] angle between wing's traling edge and centerline
    chord_flap_start = a_W - np.tan(alpha_W)*(wing_span/2-span_flap_start)-np.tan(np.radians(wing_sweep))*span_flap_start
    chord_flap_end = a_W - np.tan(alpha_W)*(wing_span/2-span_flap_end)-np.tan(np.radians(wing_sweep))*span_flap_end
    
    Jeta_1F = Airplane(
        name="JETA",
        xyz_ref=[0, 0, 0], # CG location
        wings=[
            Wing(
                name="Main Wing",
                xyz_le=[0, 0, 0], # Coordinates of the wing's leading edge
                symmetric=True,
                xsecs=[ # The wing's cross ("X") sections
                    WingXSec(  # Root
                        xyz_le=[0, 0, 0], # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
                        chord=root_chord,
                        twist=root_twist, # degrees
                        airfoil=Airfoil(name="naca"+naca_series),
                    ),
                    WingXSec(  # Flaps start
                        xyz_le=[np.tan(np.radians(wing_sweep))*span_flap_start, span_flap_start, 0],
                        chord=chord_flap_start,
                        twist=root_twist-(root_twist-tip_twist)*(span_flap_start/wing_span),
                        airfoil=Airfoil(name="naca"+naca_series).add_control_surface(deflection=flap_deflection_angle, 
                                                                 hinge_point_x=flap_hinge_position),
                        control_surface_type='asymmetric'
                    ),
                    WingXSec(  # Flaps end
                        xyz_le=[np.tan(np.radians(wing_sweep))*span_flap_end, span_flap_end, 0],
                        chord=chord_flap_end,
                        twist=root_twist-(root_twist-tip_twist)*(span_flap_end/wing_span),
                        airfoil=Airfoil(name="naca"+naca_series).add_control_surface(deflection=flap_deflection_angle, 
                                                                 hinge_point_x=flap_hinge_position),
                        control_surface_type='asymmetric'
                    ),
                    WingXSec(  # Tip
                        xyz_le=[np.tan(np.radians(wing_sweep))*wing_span/2, wing_span/2, 0],
                        chord=tip_chord,
                        twist=tip_twist,
                        airfoil=Airfoil(name="naca"+naca_series),
                    )
                ]
            ),
            Wing(
                name="Horizontal Stabilizer",
                xyz_le=[40, 0, 2],
                symmetric=True,
                xsecs=[
                    WingXSec(  # root
                        xyz_le=[40, 0, 2],
                        chord=9,
                        twist=5,
                        airfoil=Airfoil(name="naca0012"),
                        control_surface_type='symmetric',  # Elevator
                        control_surface_deflection=0,
                        control_surface_hinge_point=0.75
                    ),
                    WingXSec(  # tip
                        xyz_le=[44, 11, 2],
                        chord=5,
                        twist=5,
                        airfoil=Airfoil(name="naca0012")
                    )
                ]
            ),
            Wing(
                name="Vertical Stabilizer",
                xyz_le=[40, 0, 2],
                symmetric=False,
                xsecs=[
                    WingXSec(
                        xyz_le=[40, 0, 2],
                        chord=9,
                        twist=0,
                        airfoil=Airfoil(name="naca0012"),
                        control_surface_type='symmetric',  # Rudder
                        control_surface_deflection=0,
                        control_surface_hinge_point=0.75
                    ),
                    WingXSec(
                        xyz_le=[45, 0, 13],
                        chord=4,
                        twist=0,
                        airfoil=Airfoil(name="naca0012")
                    )
                ]
            )
        ]
    )
    
    aero_problem = VortexLatticeMethod(
        airplane=Jeta_1F,
        op_point=OperatingPoint(
            velocity=cruise_speed,
            alpha=0,
            beta=0,
            p=0,
            q=0,
            r=0,
            atmosphere=Atmosphere(altitude=cruise_altitude)
        ),
    )
    
    flap_analysis_output_data = aero_problem.run()
    CL_DF = (flap_analysis_output_data["CL"]-SLF_output_data["CL"])/flap_deflection_angle
    CD_DF = (flap_analysis_output_data["CD"]-SLF_output_data["CD"])/flap_deflection_angle
    CM_DF = (flap_analysis_output_data["Cm"]-SLF_output_data["Cm"])/flap_deflection_angle
    
    return CL_DF, CD_DF, CM_DF
    
def elevator_derivatives(aerodynamic_outputs, elevator_deflection_angle, elevator_hinge_position, SLF_output_data, cruise_altitude=18288):
     
    # Calculates cruise conditions  
    cruise_speed, Re_cruise, cruise_air_density, a_cruise, T_cruise, eta_cruise = cruise_conditions_calculator()
    
   # Uncomment the following block of code to read directly from directory
   # # Imports wing data
   # with open("aerodynamic_outputs.dat") as f:
   #     aerodynamic_outputs = f.read()
   # aerodynamic_outputs = ast.literal_eval(aerodynamic_outputs)
   
    
    root_chord = aerodynamic_outputs["chord root"]
    tip_chord = aerodynamic_outputs["chord tip"]
    root_twist = aerodynamic_outputs["twist root"]
    tip_twist = aerodynamic_outputs["twist tip"]
    wing_span = aerodynamic_outputs["span"]
    wing_sweep = aerodynamic_outputs["wing sweep"]
    naca_series = aerodynamic_outputs["NACA"]
       
    Jeta_1E = Airplane(
        name="JETA",
        xyz_ref=[0, 0, 0], # CG location
        wings=[
            Wing(
                name="Main Wing",
                xyz_le=[0, 0, 0], # Coordinates of the wing's leading edge
                symmetric=True,
                xsecs=[ # The wing's cross ("X") sections
                    WingXSec(  # Root
                        xyz_le=[0, 0, 0], # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
                        chord=root_chord,
                        twist=root_twist, # degrees
                        airfoil=Airfoil(name="naca"+naca_series),
                    ),
                    WingXSec(  # Tip
                        xyz_le=[np.tan(np.radians(wing_sweep))*wing_span/2, wing_span/2, 0],
                        chord=tip_chord,
                        twist=tip_twist,
                        airfoil=Airfoil(name="naca"+naca_series),
                    )
                ]
            ),
            Wing(
                name="Horizontal Stabilizer",
                xyz_le=[40, 0, 2],
                symmetric=True,
                xsecs=[
                    WingXSec(  # root
                        xyz_le=[40, 0, 2],
                        chord=9,
                        twist=5,
                        airfoil = Airfoil(name="naca0012").add_control_surface(deflection=elevator_deflection_angle, hinge_point_x=elevator_hinge_position),
                        control_surface_type='symmetric',  # Elevator
                    ),
                    WingXSec(  # tip
                        xyz_le=[44, 11, 2],
                        chord=5,
                        twist=5,
                        airfoil= Airfoil(name="naca0012").add_control_surface(deflection=elevator_deflection_angle, hinge_point_x=elevator_hinge_position),
                        control_surface_type='symmetric',  # Elevator
                    )
                ]
            ),
            Wing(
                name="Vertical Stabilizer",
                xyz_le=[40, 0, 2],
                symmetric=False,
                xsecs=[
                    WingXSec(
                        xyz_le=[40, 0, 2],
                        chord=9,
                        twist=0,
                        airfoil=Airfoil(name="naca0012"),
                        control_surface_type='symmetric',  # Rudder
                        control_surface_deflection=0,
                        control_surface_hinge_point=0.75
                    ),
                    WingXSec(
                        xyz_le=[45, 0, 13],
                        chord=4,
                        twist=0,
                        airfoil=Airfoil(name="naca0012"),
                        control_surface_type='symmetric',  # Rudder
                        control_surface_deflection=0,
                        control_surface_hinge_point=0.75
                    )
                ]
            )
        ]
    )

    aero_problem = VortexLatticeMethod( # Analysis type: Vortex Lattice Method, version 3
        airplane=Jeta_1E,
        op_point=OperatingPoint(
            velocity=cruise_speed,
            alpha=0,
            beta=0,
            p=0,
            q=0,
            r=0,
            atmosphere=Atmosphere(altitude=cruise_altitude)
        ),
        
    )

    elevator_analysis_output_data = aero_problem.run()
    CL_DE = (elevator_analysis_output_data["CL"]-SLF_output_data["CL"])/elevator_deflection_angle
    CD_DE = (elevator_analysis_output_data["CD"]-SLF_output_data["CD"])/elevator_deflection_angle
    CM_DE = (elevator_analysis_output_data["Cm"]-SLF_output_data["Cm"])/elevator_deflection_angle
    
    return CL_DE, CD_DE, CM_DE

    



def CL_CM_alpha_dot_calculator(aerodynamic_outputs, SLF_output_data, cruise_altitude=18288):
     
    # Calculates cruise conditions  
    cruise_speed, Re_cruise, cruise_air_density, a_cruise, T_cruise, eta_cruise = cruise_conditions_calculator()
    
   # Uncomment the following block of code to read directly from directory
   # # Imports wing data
   # with open("aerodynamic_outputs.dat") as f:
   #     aerodynamic_outputs = f.read()
   # aerodynamic_outputs = ast.literal_eval(aerodynamic_outputs)
   
    
    root_chord = aerodynamic_outputs["chord root"]
    tip_chord = aerodynamic_outputs["chord tip"]
    root_twist = aerodynamic_outputs["twist root"]
    tip_twist = aerodynamic_outputs["twist tip"]
    wing_span = aerodynamic_outputs["span"]
    wing_sweep = aerodynamic_outputs["wing sweep"]
    naca_series = aerodynamic_outputs["NACA"]
       
    Jeta_1AA = Airplane(
        name="JETA",
        xyz_ref=[0, 0, 0], # CG location
        wings=[
            Wing(
                name="Main Wing",
                xyz_le=[0, 0, 0], # Coordinates of the wing's leading edge
                symmetric=True,
                xsecs=[ # The wing's cross ("X") sections
                    WingXSec(  # Root
                        xyz_le=[0, 0, 0], # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
                        chord=root_chord,
                        twist=root_twist, # degrees
                        airfoil=Airfoil(name="naca"+naca_series),
                    ),
                    WingXSec(  # Tip
                        xyz_le=[np.tan(np.radians(wing_sweep))*wing_span/2, wing_span/2, 0],
                        chord=tip_chord,
                        twist=tip_twist,
                        airfoil=Airfoil(name="naca"+naca_series),
                    )
                ]
            ),
            Wing(
                name="Horizontal Stabilizer",
                xyz_le=[40, 0, 2],
                symmetric=True,
                xsecs=[
                    WingXSec(  # root
                        xyz_le=[40, 0, 2],
                        chord=9,
                        twist=5,
                        airfoil=Airfoil(name="naca0012"),
                        control_surface_type='symmetric',  # Elevator
                        control_surface_deflection=0,
                        control_surface_hinge_point=0.75
                    ),
                    WingXSec(  # tip
                        xyz_le=[44, 11, 2],
                        chord=5,
                        twist=5,
                        airfoil=Airfoil(name="naca0012"),
                        control_surface_type='symmetric',  # Elevator
                        control_surface_deflection=0,
                        control_surface_hinge_point=0.75
                    )
                ]
            ),
            Wing(
                name="Vertical Stabilizer",
                xyz_le=[40, 0, 2],
                symmetric=False,
                xsecs=[
                    WingXSec(
                        xyz_le=[40, 0, 2],
                        chord=9,
                        twist=0,
                        airfoil=Airfoil(name="naca0012"),
                        control_surface_type='symmetric',  # Rudder
                        control_surface_deflection=0,
                        control_surface_hinge_point=0.75
                    ),
                    WingXSec(
                        xyz_le=[45, 0, 13],
                        chord=4,
                        twist=0,
                        airfoil=Airfoil(name="naca0012")
                    )
                ]
            )
        ]
    )
    
    aero_problem = VortexLatticeMethod(
        airplane=Jeta_1AA,
        op_point=OperatingPoint(
            velocity=cruise_speed,
            alpha=5,
            beta=0,
            p=0,
            q=0,
            r=0,
            atmosphere=Atmosphere(altitude=cruise_altitude)
        ),
    )

    output_A5 = aero_problem.run()
    CL_A5 = output_A5["CL"]
    CM_A5 = output_A5["Cm"]
    
    aero_problem = VortexLatticeMethod(
        airplane=Jeta_1AA,
        op_point=OperatingPoint(
            velocity=cruise_speed,
            alpha=10,
            beta=0,
            p=0,
            q=0,
            r=0,
            atmosphere=Atmosphere(altitude=cruise_altitude)
        ),
    )
    
    output_A10 = aero_problem.run()
    CL_A10 = output_A10["CL"]
    CM_A10 = output_A10["Cm"]
    
    
    CL_A0 = SLF_output_data["CL"]
    CL_alpha_dot = ((CL_A10-CL_A5)/5 + (CL_A5-CL_A0)/5)/5
    
    CM_A0 = SLF_output_data["Cm"]
    CM_alpha_dot = ((CM_A10-CM_A5)/5 + (CM_A5-CM_A0)/5)/5
    
    return CL_alpha_dot, CM_alpha_dot





def CLu_CDu_CMu_calculator(aerodynamic_outputs, SLF_output_data, cruise_altitude=18288):
     
    # Calculates cruise conditions  
    cruise_speed, Re_cruise, cruise_air_density, a_cruise, T_cruise, eta_cruise = cruise_conditions_calculator()
    
   # Uncomment the following block of code to read directly from directory
   # # Imports wing data
   # with open("aerodynamic_outputs.dat") as f:
   #     aerodynamic_outputs = f.read()
   # aerodynamic_outputs = ast.literal_eval(aerodynamic_outputs)
   
    
    root_chord = aerodynamic_outputs["chord root"]
    tip_chord = aerodynamic_outputs["chord tip"]
    root_twist = aerodynamic_outputs["twist root"]
    tip_twist = aerodynamic_outputs["twist tip"]
    wing_span = aerodynamic_outputs["span"]
    wing_sweep = aerodynamic_outputs["wing sweep"]
    naca_series = aerodynamic_outputs["NACA"]
       
    Jeta_1AB = Airplane(
        name="JETA",
        xyz_ref=[0, 0, 0], # CG location
        wings=[
            Wing(
                name="Main Wing",
                xyz_le=[0, 0, 0], # Coordinates of the wing's leading edge
                symmetric=True,
                xsecs=[ # The wing's cross ("X") sections
                    WingXSec(  # Root
                        xyz_le=[0, 0, 0], # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
                        chord=root_chord,
                        twist=root_twist, # degrees
                        airfoil=Airfoil(name="naca"+naca_series),
                    ),
                    WingXSec(  # Tip
                        xyz_le=[np.tan(np.radians(wing_sweep))*wing_span/2, wing_span/2, 0],
                        chord=tip_chord,
                        twist=tip_twist,
                        airfoil=Airfoil(name="naca"+naca_series),
                    )
                ]
            ),
            Wing(
                name="Horizontal Stabilizer",
                xyz_le=[40, 0, 2],
                symmetric=True,
                xsecs=[
                    WingXSec(  # root
                        xyz_le=[40, 0, 2],
                        chord=9,
                        twist=5,
                        airfoil=Airfoil(name="naca0012"),
                        control_surface_type='symmetric',  # Elevator
                        control_surface_deflection=0,
                        control_surface_hinge_point=0.75
                    ),
                    WingXSec(  # tip
                        xyz_le=[44, 11, 2],
                        chord=5,
                        twist=5,
                        airfoil=Airfoil(name="naca0012"),
                        control_surface_type='symmetric',  # Elevator
                        control_surface_deflection=0,
                        control_surface_hinge_point=0.75
                    )
                ]
            ),
            Wing(
                name="Vertical Stabilizer",
                xyz_le=[40, 0, 2],
                symmetric=False,
                xsecs=[
                    WingXSec(
                        xyz_le=[40, 0, 2],
                        chord=9,
                        twist=0,
                        airfoil=Airfoil(name="naca0012"),
                        control_surface_type='symmetric',  # Rudder
                        control_surface_deflection=0,
                        control_surface_hinge_point=0.75
                    ),
                    WingXSec(
                        xyz_le=[45, 0, 13],
                        chord=4,
                        twist=0,
                        airfoil=Airfoil(name="naca0012")
                    )
                ]
            )
        ]
    )
    
    aero_problem = VortexLatticeMethod(
        airplane=Jeta_1AB,
        op_point=OperatingPoint(
            velocity=cruise_speed/1.2,
            alpha=0,
            beta=0,
            p=0,
            q=0,
            r=0,
            atmosphere=Atmosphere(altitude=cruise_altitude)
        ),
    )

    output_U12 = aero_problem.run()
    
    CL_U12 = output_U12["CL"]
    CD_U12 = output_U12["CL"]
    CM_U12 = output_U12["Cm"]
    CL_U0 = SLF_output_data["CL"]
    CD_U0 = SLF_output_data["CD"]
    CM_U0 = SLF_output_data["Cm"]
    
    CLu = (CL_U12-CL_U0)/(cruise_speed/1.2-cruise_speed)
    CDu = (CD_U12-CD_U0)/(cruise_speed/1.2-cruise_speed)
    CMu = (CM_U12-CM_U0)/(cruise_speed/1.2-cruise_speed)
    
    return CLu, CDu, CMu

