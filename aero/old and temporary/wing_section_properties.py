# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 12:17:54 2023

@author: Andre Santos
"""

import numpy as np

def InertiaMomentCalculator(x_upper, y_upper, x_lower, y_lower, t=0.001):

    """Calculates the 2nd moment of area (i.e. moment of inertia) of an aerofoil
       from its coordinates, using thin-walled section assumptions and parallel
       axis theorem.
    
    Args:
        x_upper (array of float)  [m]:   x-axis coordinates of aerofoil's upper surface
        y_upper (array of float)  [m]:   y-axis coordinates of aerofoil's upper surface
        x_lower (array of float)  [m]:   x-axis coordinates of aerofoil's bottom surface
        y_lower (array of float)  [m]:   y-axis coordinates of aerofoil's bottom surface
        t       (float)           [m]:   thickness of aerofoil skin; Default = 1 mm
        
    Returns:
        I_xx        (float)  [m^4]:   2nd moment of area about the x-axis
        I_yy        (float)  [m^4]:   2nd moment of area about the y-axis
        I_xy        (float)  [m^4]:   Product of inertia [m^4]
        x_centroid  (float)    [m]:   Aerofoil's centroid's abcissa
        y_centroid  (flaot)    [m]:   Aerofoil's centroid's ordinate
    """
    
    
    # STEP 1 - CALCULATE AEROFOIL'S CENTROID
    # Calculates each aerofoil segment's length and centroid coordinates
    L_seg = (lambda x_i0, x_i1, y_i0, y_i1: ((x_i1-x_i0)**2 + (y_i1-y_i0)**2)**0.5)     \
                                                (np.append(x_upper,np.flip(x_lower)),   \
                                                 np.append(np.append(x_upper[1:],np.flip(x_lower)),x_upper[0]), \
                                                 np.append(y_upper,np.flip(y_lower))  , \
                                                 np.append(np.append(y_upper[1:],np.flip(y_lower)),y_upper[0]))[:-1]
    x_seg_bar = (lambda x_i0, x_i1: (x_i0+x_i1)/2)(np.append(x_upper,np.flip(x_lower)), \
                                                   np.append(np.append(x_upper[1:],np.flip(x_lower)),x_upper[0]))[:-1]
    y_seg_bar = (lambda y_i0, y_i1: (y_i0+y_i1)/2)(np.append(y_upper,np.flip(y_lower)), \
                                                   np.append(np.append(y_upper[1:],np.flip(y_lower)),y_upper[0]))[:-1]
    
    # Calculates aerofoil's centroid coordinates
    x_centroid = sum(x_seg_bar*L_seg*t)/sum(L_seg*t)
    y_centroid = sum(y_seg_bar*L_seg*t)/sum(L_seg*t)
    
    
    # STEP 2 - CALCULATE THE MOMENTS OF INERTIA
    # Computes aerofoil segments' coordinates relative to aerofoil's centroid
    x_seg_bar = x_seg_bar - x_centroid
    y_seg_bar = y_seg_bar - y_centroid
    
    # Computes each segment's angle with the x-axis, in radians
    theta_seg = (lambda x_i0, x_i1, y_i0, y_i1: np.arctan((y_i1-y_i0)/(x_i1-x_i0)))     \
                                                (np.append(x_upper,np.flip(x_lower)),   \
                                                 np.append(np.append(x_upper[1:],np.flip(x_lower)),x_upper[0]), \
                                                 np.append(y_upper,np.flip(y_lower))  , \
                                                 np.append(np.append(y_upper[1:],np.flip(y_lower)),y_upper[0]))[:-1]
    
    # Computes the local moments of inertia for each segment
    I_xx_local = (L_seg**3)*t*(np.sin(theta_seg)**2)/12
    I_yy_local = (L_seg**3)*t*(np.cos(theta_seg)**2)/12
    I_xy_local = (L_seg**3)*t*(np.sin(2*theta_seg))/24
    
    # Computes the global moments of inertia for the aerofoil
    I_xx = sum(I_xx_local + L_seg*t*y_seg_bar**2)    
    I_yy = sum(I_yy_local + L_seg*t*x_seg_bar**2)    
    I_xy = sum(I_xy_local + L_seg*t*x_seg_bar*y_seg_bar)    
    
    
    return I_xx, I_yy, I_xy



def centroidCalculator(x_upper, y_upper, x_lower, y_lower, t=0.001):

    """Calculates the centroid of an aerofoil from its coordinates.
    
    Args:
        x_upper (array of float)  [m]:   x-axis coordinates of aerofoil's upper surface
        y_upper (array of float)  [m]:   y-axis coordinates of aerofoil's upper surface
        x_lower (array of float)  [m]:   x-axis coordinates of aerofoil's bottom surface
        y_lower (array of float)  [m]:   y-axis coordinates of aerofoil's bottom surface
        t       (float)           [m]:   thickness of aerofoil skin; Default = 1 mm
        
    Returns:
        x_centroid  (float)    [m]:   Aerofoil's centroid's abcissa
        y_centroid  (flaot)    [m]:   Aerofoil's centroid's ordinate
    """
    
    # Calculates each aerofoil segment's length and centroid coordinates
    L_seg = (lambda x_i0, x_i1, y_i0, y_i1: ((x_i1-x_i0)**2 + (y_i1-y_i0)**2)**0.5)     \
                                                (np.append(x_upper,np.flip(x_lower)),   \
                                                 np.append(np.append(x_upper[1:],np.flip(x_lower)),x_upper[0]), \
                                                 np.append(y_upper,np.flip(y_lower))  , \
                                                 np.append(np.append(y_upper[1:],np.flip(y_lower)),y_upper[0]))[:-1]
    x_seg_bar = (lambda x_i0, x_i1: (x_i0+x_i1)/2)(np.append(x_upper,np.flip(x_lower)), \
                                                   np.append(np.append(x_upper[1:],np.flip(x_lower)),x_upper[0]))[:-1]
    y_seg_bar = (lambda y_i0, y_i1: (y_i0+y_i1)/2)(np.append(y_upper,np.flip(y_lower)), \
                                                   np.append(np.append(y_upper[1:],np.flip(y_lower)),y_upper[0]))[:-1]
    
    # Calculates aerofoil's centroid coordinates
    x_centroid = sum(x_seg_bar*L_seg*t)/sum(L_seg*t)
    y_centroid = sum(y_seg_bar*L_seg*t)/sum(L_seg*t)
    
    return x_centroid, y_centroid    