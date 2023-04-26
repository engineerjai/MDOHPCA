# -*- coding: utf-8 -*-
"""
Landing Gear

calculate the tension on the landing gear, and determine if it is acceptable
landing gear is a constraint

Created on Wed Apr 26 13:20:24 2023

@author: Tanika Dodhia
"""
import math

#constants
gravity=9.81

#variables
aircraft_mass=20000
acceleration=43
angle=45
theta=math.radians(angle)

#limits
tension_max=900000

tension= (aircraft_mass*acceleration)/(2*math.sin(theta))

if tension>tension_max:
    print ('Not Ok')
if tension <tension_max:
    print ('Ok')