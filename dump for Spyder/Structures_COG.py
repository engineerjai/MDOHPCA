# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Centre of Gravity 
    Structures Component
Author: Tanika Dodhia
Date: 12/03/2023

"""

import math
import matplotlib.pyplot as plt

#VARIABLES

#changeable
sweep= 0.015 #sweep angle wing
wing_x=30 #wing span

#fixed variables
mach =0.6 #mach number
passengers= 220 #number of passengers
crew= 12 #number of crew
engines= 2 #number of engines
nacelles=engines #number of nacelles
hyd_pressure = 344 #hydraulic pressure in bar
max_range = 8334*1000 #max range



#weights
take_off = 10000 #take-off weight
empty_weight = 5000 #empty weight
fuel_weight = 3000 # equivalent fuel weight

#fuselage
f_taper= math.radians(12) #fuselage taper angle
f_diameter = 7 #diameter fuselage
f_wet_area= 442 #fuselage wetted area from A350XWB forum (JELKS, 2021)

#interior
length_cock= 5 # cockpit length
length_cabin = 26.6 #cabin length as per A310 max (JELKS)
length_cargo = 1.534*4 #cargo length (4 pallets LD3 - need 8) wikipedia (JELKS)
forward_bulk = f_diameter/5 #x axis 

#wings
average_chord = 3 #average chord
s = 200 # wing area
b= 80 #wing span
aspect_ratio = 12 #aspect ratio
flap_area= s*0.172385620 #wing flap area
sweep= math.radians(sweep)
y_wing = f_diameter-0.2 #wing position
incident_angle = math.radians(0) # wing incidence angle
thic_chord_rat = 0.127 # thickness to chord ratio for the wing root



#landing gear
x_landing = length_cock #distance to front of landing gear

#tail
tail_hc = 5.5 #horizontal tail chord
tail_hs= 16 #horizontal tail span
tail_h_ar= tail_hs/tail_hc #horizontal aspect ratio
tail_h_angle = 30 #average tail angle horizontal A350 (JELKS)

tail_vc = 6 #vertical tail chord
tail_vh= 6 #vertical tail height
tail_v_ar= tail_vh/tail_vc #vertical aspect ratio
tail_v_angle = 32.5 #average tail angle vertical A350 (JELKS)

#realtive positions

x_spar_fwd=wing_x +0.2376*average_chord #distance to front spar
x_spar_aft = wing_x +0.56937*average_chord #distance to aft spar
x_trail_Edge = wing_x + average_chord # distance to trailing edge

length_fuselage=30 #fueslage length
tail_x=length_fuselage-tail_hc-0.5



length_main_lg= 0.75*length_fuselage*0.0254 #length of main landing gear
length_nose_lg = 0.7*length_fuselage*0.0254 #length of nose landing gear

#WEIGHTS


design_gw = 2.20462262*take_off #gross design weight pounds
load_limit=1.5*2 #1.5*load factor for limit
l_fuse_ft= length_fuselage*3.281 # fueslage length in feet
f_wet_area_ft2= f_wet_area*10.764 #fueslage wetted area in ft2
taper_ratio=0.2 #taper ratio
b_ft=b*3.281 #wing span in ft
constant_b= 0.75*((1+2*taper_ratio)/(1+taper_ratio))*((b_ft*math.tan(sweep))/1) #constant for wing span


#wing
s_ft = s*3.28084 #wing area in ft
control_surface = flap_area*10.7639 #area of control surface in ft2

weight_wing_lb=0.0051*((design_gw*load_limit)**0.557)*(s_ft**0.649)*(aspect_ratio**0.5)*((thic_chord_rat)**-0.4)*((1+taper_ratio)**0.1)*((math.cos(sweep))**-1)*control_surface*0.1
weight_wing=weight_wing_lb*0.453592 #lbs to kg


#tail horizontal
constant_uht = 1 #not everything moves
f_width = 1*3.28084 # fuselage width in ft
tail_hs_ft = tail_hs*3.28084 # tail span horizontal in ft
tail_hc_ft = tail_hc**3.28084 # tail length horizontal in ft
tail_area_h= tail_hs*tail_hc *10.7639 #tail area horizontal in ft2
gyration=0.3*tail_hc_ft #radius of gyration
elevator_area= 0.25*tail_area_h

weight_tailh_lb = 0.0379*constant_uht*((1+f_width/tail_hs_ft)**0.25)*(design_gw*0.639)*(load_limit**0.1)*(tail_area_h**0.75)*(tail_hc_ft*-1)*(gyration**0.704)*((math.cos(sweep))**-1)*(tail_h_ar**0.166)*(((1+elevator_area)/tail_area_h)**0.1)
weight_tailh = weight_tailh_lb*0.453592

#tail vertical

tail_vc_ft = tail_vc*3.28084 #vertical tail in feet
tail_vh_ft= tail_vh*3.28084 #vertical tail height in ft
tail_area_v= tail_vh*tail_vc *10.7639 #tail area vertical in ft2
gyration_v=tail_vc_ft #yaw gyration radius
root_taper_chord_ratio = 0.09 #root taper chord ratio


weight_tailv_lb = 0.0026*(2**0.225)*(design_gw*0.536)*(tail_vc_ft**-0.5)*(tail_area_v**0.5)*(gyration_v**0.875)*((math.cos(sweep))**-1)*(tail_v_ar**0.35)*(root_taper_chord_ratio**-0.5)
weight_tailv = weight_tailv_lb**0.453592




#coordinates


#wings
weight_wingxy=(wing_x+0.4*average_chord+0.25*b*math.tan(sweep), y_wing)
cg_wing =[(weight_wingxy[0]*weight_wing),(weight_wingxy[1]*weight_wing)]

#tail
#horizontal
tail_h_angle=math.radians(tail_h_angle) #deg to rad
weight_tailhxy=(tail_x+0.3*tail_hc+0.25*tail_hs*math.tan(tail_h_angle),f_diameter)
cg3h=[(weight_tailhxy[0]*weight_tailh),(weight_tailhxy[1]*weight_tailh)]

#vertical
tail_v_angle=math.radians(tail_v_angle) #deg to rad
weight_tailvxy=(tail_x+0.3*tail_vc+0.25*tail_vh*math.tan(tail_v_angle),f_diameter+0.5*tail_vh)
cg3v=[(weight_tailvxy[0]*weight_tailv),(weight_tailvxy[1]*weight_tailv)]


