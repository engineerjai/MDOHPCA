# -*- coding: utf-8 -*-
"""
Spars Locations Scripts

Script to determine the location of  the various spars for each chord length

Created on Wed Apr 26 12:28:14 2023

@author: Tanika Dodhia
"""

import numpy as np
import math

taper=0.45
half_span=30
root_chord= 7
tip_chord=root_chord*taper
sweep=15
sweep_rad=math.tan(math.radians(sweep))

#first spar 1/4 of the way through
#sweep angle is considered to be 0

spars=np.zeros((2,2))
spars[0]=np.array([0.25*root_chord, (0.25*tip_chord)+half_span*sweep_rad])
spars[1]=np.array([0.5*root_chord, (0.5*tip_chord)+half_span*sweep_rad])

#average spacing for spars determined from literature
#spar 25% and 50% along the wing - values determined from literature and the abaqus models
