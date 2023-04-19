# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 12:13:30 2023

@author: Nuno Santos
"""

def k1_calculator(naca_series):
    camber_line = float(naca_series[:3])
    
    k1_table = {210:361.40, 220:51.640, 230:15.957, 
          240:6.643, 250:3.230}
        
    k1 = k1_table[camber_line]
    return k1