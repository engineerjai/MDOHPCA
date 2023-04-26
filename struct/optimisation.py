# -*- coding: utf-8 -*-
"""
Optimisation Script

This script looks at the outputs from abaqus, and reads them in from an '.xlsx'
It looks at the stress and deflection requirements, then returns the row with the minimum weight


Created on Tue Apr 25 23:36:15 2023

@author: Tanika Dodhia
"""

import pandas as pd

#import data
data_raw = pd.read_excel('DummyData.xlsx')
 
#constraints
max_stress=100000000
max_deflection=3

#filtering data using the constraints
stress_constraint=data_raw[data_raw.Stress.lt(max_stress).groupby(data_raw.Weight).transform('all')]
deflection_constraint=stress_constraint[stress_constraint.Deflection.lt(max_deflection).groupby(stress_constraint.Weight).transform('all')]
data_modified=deflection_constraint

#determining the chord length and span with the lowest weight for optimisation purposes
min_weight=data_modified[data_modified["Weight"]==data_modified["Weight"].min()]
span=min_weight.iloc[0,0]
chord_length=min_weight.iloc[0,1]
