# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 17:06:32 2023

@author: Andre Santos
"""

import numpy as np
from sympy import symbols, solve

def k2_k1_ratio_calculator(naca_series):
    
    k2_k1_table = {221:0.000764, 231:0.00677, 241:0.0303, 251:0.1355}
    
    camber_line = int(naca_series[:3])
    
    try:
        k2_k1 = k2_k1_table[camber_line]
    except:
        # Calculates approx. k2/k1 from https://archive.aoe.vt.edu/mason/Mason_f/CAtxtAppA.pdf
        # NTS: The following chunk of code does not work because line 24 does not apply for reflexed airfoils;
        # NTS: Ask for assistance regarding finding m so that Cmc/4 = 0 from thin airfoil theory
        xf = 0.05*float(naca_series[1])
        m = symbols('m')
        expr = m*(1-(m/3)**0.5)-xf
        m = solve(expr)
        m = [x.as_two_terms()[0] for x in m]
        m = float(m[0])
        k2_k1 = (3*(m-xf)**2-m**3)/(1-m)**3
        
    return k2_k1