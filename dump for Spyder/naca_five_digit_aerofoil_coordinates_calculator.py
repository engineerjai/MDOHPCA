# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 11:39:45 2023

@author: Andre Santos
"""

import numpy as np
from k1_calculator import k1_calculator


def naca_five_digit_aerofoil_coordinates_calculator(naca_series, chord_length=1, number_of_points=100):

    """Generates x, y coordinates for a NACA aerofoil of five-digit series.
    
    Args:
        naca_series (str): NACA series, e.g. '23013'
        chord_length (float): Chord length of the aerofoil [m] (Default = 1m)
        number_of_points (int): Number of coordinates to create (Default = 100)
    
    Returns:
        (list of float, list of float): x and y coordinates of the airfoil
        
        NTS: Only non-reflex aerfoils (NACA""0"") work for now
        NTS: Expand function to include convex aerofoils (NACA""1"")
        NTS: The first two digits (NACALP) have finite combinations (221, 22 23, 24 and 25)
        NTS: Comment the code
        NTS = Note to Self
        
    """
 
    t = float(naca_series[3:5])/100
    
    x = np.arange(0, chord_length + chord_length/number_of_points, chord_length/number_of_points)
    y_t = 5*t*(0.2969*x**0.5-0.1260*x-0.3516*x**2+0.2843*x**3-0.1015*x**4)
    r = 1.1019*t**2
    k1 = k1_calculator(naca_series)
    
    y_c = np.zeros(number_of_points+1)
    y_c[x < r] = k1/6*(x[x < r]**3-3*r*x[x < r]**2+r**2*(3-r)*x[x < r])
    y_c[x > r] = k1/6*r**3*(1-x[x > r])
    theta = np.zeros(number_of_points+1)
    theta[x < r] = np.arctan(k1/6*(3*x[x < r]**2-6*r*x[x < r]+r**2*(3-r)))
    theta[x > r] = np.arctan(-k1/6*r**3)
    
    x_upper = x - y_t*np.sin(theta)
    x_lower = x + y_t*np.sin(theta)
    y_upper = y_c + y_t*np.cos(theta)
    y_lower = y_c - y_t*np.cos(theta)
    
    x = np.append(x_lower, np.flip(x_upper))
    y = np.append(y_lower, np.flip(y_upper))
    
    return x, y



