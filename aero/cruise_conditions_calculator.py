# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 11:25:45 2023

@author: Andre Santos
"""

import warnings
import numpy as np
warnings.filterwarnings("ignore")

def cruise_conditions_calculator(cruise_mach_number=0.7, cruise_altitude=18288, Re_characteristic_length=1):
    
    """
       Calculates air properties and cruise flight conditions for a given cruise altitude.
       All inputs are optional; if not specified, default values are used.
    
    Opt Args:
        cruise_mach_number (float): IAS Mach number
        cruise_altitude (float): Altitude of cruise [m]
        Re_characteristic_length (float): Parameter used to calculate dimensional Reynolds number [m]
        
     Default Values:       
        cruise_mach_number        (float) []   (Default = 0.7)
        cruise_altitude           (float) [m]  (Default = 18288m = 60k ft)
        Re_characteristic_length  (float) [m]  (Default = 1)
    
    Returns:
        cruise_speed              (float) [m/s]
        Re_cruise                 (float) []
        cruise_air_density        (float) [kg/m^3]
        a_cruise                  (float) [m/s]
        T_cruise                  (float) [K]
        cruise_dynamic_viscosity  (float) [Pa.s] 
        
    """
    
    # Constants ISA
    T0 = 288.15             # Air temperature at sea level
    gama = 1.4              # Ratio of specific heats
    R = 287.26              # Universal Gas Constant [J/kg/K]
    g = 9.80665             # Gravitational acceleration [m/s^2]
    lapse_rate = 0.0065     # [°C/m]
    air_density_0 = 1.2256  # Air density at sea level [kg/m^3]
    a0 = 340.1              # Speed of sound at sea level [m/s^2]
    
    
    # Cruise altitude ambient conditions
    T_cruise = T0 - cruise_altitude*lapse_rate
    a_cruise = np.sqrt(R*T_cruise*gama)
    cruise_air_density = air_density_0*(T_cruise/T0)**(g/(lapse_rate*R)-1)
    cruise_dynamic_viscosity = 0.000001458*T_cruise**1.5/(T_cruise+110.4) # Formula source: https://www.omnicalculator.com/physics/kinematic-viscosity-of-air
    
    # Cruise speed calculation
    cruise_speed = cruise_mach_number*a0 # Indicated airspeed at cruise altitude
    
    # Reynolds number calculation
    try:
        # If the characteristic length is provided, calculates dimensional Reynolds number
        Re_cruise = cruise_speed*Re_characteristic_length*cruise_air_density/cruise_dynamic_viscosity
    except:
        # If no characteristic length is provided, calculates dimensionless Reynolds number
        Re_cruise = cruise_speed*cruise_air_density/cruise_dynamic_viscosity

    return cruise_speed, Re_cruise, cruise_air_density, a_cruise, T_cruise, cruise_dynamic_viscosity









