# -*- coding: utf-8 -*-
"""
Bending Moment Script

Created on Wed Apr 26 12:57:51 2023

@author: Tanika Dodhia
"""

import pandas as pd

#constants defined in literature
Ka=float(0.60)
Ki=float(0.036)
E= float(900000)

#chord length
chord_length=float(4)
#span
half_span=float(30)

upper_naca = pd.read_excel('naca23012_upper_coordinates.xlsx')
lower_naca= pd.read_excel('naca23012_lower_coordinates.xlsx')

zu=upper_naca.iloc[:,1]
zl=lower_naca.iloc[:,1]

for i in range(zu.size+1):
    t_general=zu-zl
    h_general=(zu-zl)*0.5
    
t=max(t_general)
h=max(h_general)

tau=float(t/chord_length)
epsilon=float(h/chord_length)

#calculations of area and bending inertia
A=Ka*chord_length*chord_length*tau
I= (Ki*(chord_length**4)*tau)*((tau**2)+(epsilon**2))
#bending stiffness
ei=E*I 

#bending moment
p=1.118
v=13
cl=0.6
wl=0.5*p*v*cl
bending_moment=ei*wl
