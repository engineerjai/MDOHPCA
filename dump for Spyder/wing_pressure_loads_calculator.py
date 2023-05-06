# -*- coding: utf-8 -*-
"""
Created on Sat Apr 15 16:17:18 2023

@author: André Santos
"""

import numpy as np

#def wing_pressure_loads(CL, CD, x_upper, y_upper, x_lower, y_lower, air_speed=cruise_speed, air_density=cruise_air_density):
def wing_pressure_loads(CL, CD, x_upper, y_upper, x_lower, y_lower, air_speed='Default', air_density='Default'):
    
    """Calculates the approximated pressure exerted on the top and bottom surfaces
       of an aerofoil with the given coordinates.
    
    Args:
        CL                 (float)           []:         Wing's lift coefficient
        CD                 (float)           []:         Wing's drag coefficient
        x_upper            (array of float)  [m]:        x-axis coordinates of aerofoil's upper surface
        y_upper            (array of float)  [m]:        y-axis coordinates of aerofoil's upper surface
        x_lower            (array of float)  [m]:        x-axis coordinates of aerofoil's bottom surface
        y_lower            (array of float)  [m]:        y-axis coordinates of aerofoil's bottom surface
        cruise_speed       (float)           [m/s]:      air speed (Default = cruise speed)
        air_density        (float)           [kg/m^3]:   outside air density (Default = air density at cruise altitude)
        
    Returns:
        upper_pressure     (float)           [Pa]:       Global pressure applied to the wing's top surface
        lower_pressure     (float)           [Pa]:       Global pressure applied to the wing's bottom surface
        
    NOTE: The method used in this script is only an approximation to the reality.
    For more accurate results please consult ESDU73012 or use a higher fidelity 
    CFD solver. This method only works only for 'quasi-symmetric' aerofoils.
    """

    # If the last two arguments were not provided, the default cruise options are given to them
    if (air_speed == 'Default') | (air_density == 'Default'):
        from cruise_conditions_calculator import cruise_conditions_calculator
        cruise_speed, Re_cruise, cruise_air_density, a_cruise, T_cruise, eta_cruise = cruise_conditions_calculator()
        if (air_speed == 'Default'): air_speed = cruise_speed
        if (air_density == 'Default'): air_density = cruise_air_density
    
    # If True, activates plotting for debugging purposes
    to_plot = False
    
    # Calculates dynamic pressure
    q = 0.5*air_density*air_speed**2
    
    # Pressure = Force/Area
    # Drag = 0.5*rho*CD*v^2*Area
    # Drag/Area = 0.5*rho*v^2*CD
    # Forward_Pressure = Drag/Area
    # Forward_Pressure = 0.5*rho*CD*v^2
    # Forward_Pressure = CD*q
    forward_pressure = q * CD
    
    # Pressure = Force/Area
    # Lift = 0.5*rho*CL*v^2*Area
    # Lift/Area = 0.5*rho*v^2*CL
    # Forward_Pressure = Lift/Area
    # Forward_Pressure = 0.5*rho*CL*v^2
    # Forward_Pressure = CL*q
    vertical_pressure = q * CL
    
    # Calculates the effective wetted area to which the drag force is applied to in the upper and lower surfaces
    L_seg = (lambda x_i0, x_i1, y_i0, y_i1: ((x_i1-x_i0)**2 + (y_i1-y_i0)**2)**0.5)     \
                                                (np.append(x_upper,x_lower),   \
                                                 np.append(np.append(x_upper[1:],x_lower),x_upper[0]), \
                                                 np.append(y_upper,y_lower)  , \
                                                 np.append(np.append(y_upper[1:],y_lower),y_upper[0]))
                                                    
    # Defines the limiting points and segments of the wetted upper surface
    y_max = max(y_upper)
    y_min = min(y_lower)
    A_upper_indices = np.zeros([x_upper.size*2], dtype='bool')
    A_upper_indices[x_upper.size:] = (y_upper<=y_max) & (x_upper<=x_upper[y_upper==y_max])
    A_lower_indices = np.zeros([x_lower.size*2], dtype='bool')
    A_lower_indices[:x_lower.size] = (y_lower>=y_min) & (x_lower<=x_lower[y_lower==y_min])
    
    # Plotting tools for debugging purposes
    if to_plot:
        import matplotlib.pyplot as plt
        plt.plot(x_upper[(y_upper<=y_max) & (x_upper<=x_upper[y_upper==y_max])], y_upper[(y_upper<=y_max) & (x_upper<=x_upper[y_upper==y_max])])
        plt.plot(x_lower[(y_lower>=y_min) & (x_lower<=x_lower[y_lower==y_min])], y_lower[(y_lower>=y_min) & (x_lower<=x_lower[y_lower==y_min])])
        plt.title('Aerofoil wetted effective area')
        plt.axis('equal')
        plt.show()
    
    # Calculates the effective wetted area for upper and lower surfaces
    A_upper = sum(L_seg[A_upper_indices])
    A_lower = sum(L_seg[A_lower_indices])
    
    # Destributes the vertical and forward pressures between the upper and lower surfaces of the wing
    upper_pressure = forward_pressure * (A_upper/(A_upper+A_lower))
    lower_pressure = vertical_pressure + (vertical_pressure-upper_pressure)
    
    
    return upper_pressure, lower_pressure