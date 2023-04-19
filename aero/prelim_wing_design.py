# -*- coding: utf-8 -*-
"""
   Calculates the initial wing area and lift coefficient required.

Args:
    mass           (float)   [kg]:    Aircraft's total mass.
    span           (float)   [m]:     Wing span (Default = 80m)
Returns:
    wing_area      (float)   [m^2]:   Projected planform wing area.
    CL_cruise      (float)   []:      Cruise lift coefficient
"""

def innit(mass, output, span=80):
    # Imports weather conditions at cruise level
    from cruise_conditions_calculator import cruise_conditions_calculator
    cruise_speed, Re_cruise, cruise_air_density, a_cruise, T_cruise, eta_cruise = cruise_conditions_calculator()
    
    # Caculates area required assuming that Lift = Weight and that CL~=0.627631
    CL_cruise = 0.627631
    weight = 9.80665 * mass
    q = 0.5*cruise_air_density*cruise_speed**2
    wing_area = weight/(q*CL_cruise)
    
    # Calculate the average chord length based on a flat rectangular wing
    aver_chord = wing_area/span
    
    if output=='avg_chord':
        return aver_chord
    elif output=='S_Ref':
        return wing_area
    elif output=='CL_cruise':
        return CL_cruise
    

def S_Ref(mass, span=80):
    """ Calculates preliminary wing area
        Args: mass, span (Default=80m)"""
    wing_area = innit(mass, 'S_Ref', span,)
    return wing_area
    
def avg_chord(mass, span=80):
    """ Calculates preliminary average wing chord length
        Args: mass, span (Default=80m)"""
    aver_chord = innit(mass, 'avg_chord', span)
    return aver_chord

def CL_cruise(mass, span=80):
    """ Calculates preliminary cruise lift coefficient
        Args: mass, span (Default=80m)"""
    aver_chord = innit(mass, 'CL_cruise', span)
    return aver_chord

