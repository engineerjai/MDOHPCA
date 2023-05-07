#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#import necessary package and collect initial variables
import openmdao.api as om
import pandas as pd
import numpy as np
import math

print('structures loaded')
# In[ ]:


class struct(om.Group):
    def initialize(self):
        
        #read all excel constants        
        initial_vars = pd.read_excel('constants.xlsx')
        
        #sort for just the struct constants
        all_struct_vars = np.where(initial_vars['struct'] == True, [initial_vars['variable name'], initial_vars['value'], initial_vars['description']], None)
        
        #create a DataFrame for the struct variables to pull from in individual analyses
        all_struct_vars = pd.DataFrame(all_struct_vars).dropna(axis = 1).transpose()
        
        #pass variables into the correct openMDAO format, passing name, description and value
        i=0
        for variable in all_struct_vars[0]:
            value = all_struct_vars.iloc[i][1]
            description = all_struct_vars.iloc[i][2]
            self.options.declare(variable, default = value, desc = description)
            i+=1
        i=0
        
        
        WingData = pd.read_excel('WingData.xlsx')
        self.options.declare('WingData', WingData)
        
    def setup(self):
        #set analysis (subsystem) names and what inputs they take and give
        self.add_subsystem('bendingmoments', bendingmoments())
        self.add_subsystem('landing_gear_calcs', landing_gear_calcs())
        self.add_subsystem('spars', spars())
        self.add_subsystem('optimisation', optimisation())
        self.add_subsystem('structures_cog', structures_cog())
        #can optionally promote the variables at this step
        
    def configure(self):
        #promote all variables (lazy option, they can be connected individually)
        self.promotes('bendingmoments',any=['*'])
        self.promotes('landing_gear_calcs',any=['*'])
        self.promotes('spars',any=['*'])
        self.promotes('optimisation',any=['*'])
        self.promotes('structures_cog',any=['*'])


# In[ ]:


class bendingmoments(om.ExplicitComponent):
    def initialize(self):
        #take constants from group a level above
        self.options.declare('Ka',prob.model.struct.options['Ka'])
        self.options.declare('Ki',prob.model.struct.options['Ki'])
        self.options.declare('E',prob.model.struct.options['E'])
        self.options.declare('rho',prob.model.struct.options['rho'])
        self.options.declare('v',prob.model.struct.options['v'])
        
        self.options.declare('upper_naca', pd.read_excel('naca23012_upper_coordinates.xlsx'))
        self.options.declare('lower_naca', pd.read_excel('naca23012_lower_coordinates.xlsx'))
        
    def setup(self):
        #define inputs (variables to come from other disciplines) and default values
        self.add_input('CL', val = 0.6)
        self.add_input('c', val = 3)
        self.add_input('b', val = 60)
        #prob.model.add_design_var
        self.add_output('bending_moment')
        
    def compute(self,inputs,outputs):
        
        Ka = self.options['Ka']
        Ki = self.options['Ki']
        E = self.options['E']
        c = inputs['c']
        b = inputs['b']
        upper_naca = self.options['upper_naca']
        lower_naca = self.options['lower_naca']
        
        half_span = b

        zu=upper_naca.iloc[:,1]
        zl=lower_naca.iloc[:,1]

        for i in range(zu.size+1):
            t_general=zu-zl
            h_general=(zu-zl)*0.5

        t=max(t_general)
        h=max(h_general)

        tau=float(t/c)
        epsilon=float(h/c)

        #calculations of area and bending inertia
        A=Ka*c*c*tau
        I= (Ki*(c**4)*tau)*((tau**2)+(epsilon**2))
        #bending stiffness
        ei=E*I 

        #bending moment
        p=self.options['rho']
        v1=self.options ['v']
        cl = inputs['CL']
        wl=0.5*p*v1*cl
        bending_moment=ei*wl
        
        outputs['bending_moment'] = bending_moment


# In[ ]:


class spars(om.ExplicitComponent):
    def initialize(self):
        #take constants from group a level above
        pass
    def setup(self):
        #define inputs (variables to come from other disciplines) and default values
        self.add_input('taper', val = 0.45)
        self.add_input('b', val = 30) #design input?
        self.add_input('c', val = 7)#design input
        self.add_input('sweep', val = 15)
        
        self.add_output('spars', shape = (2,2))
                
    def compute(self,inputs,outputs):
        sweep_rad=math.tan(math.radians(inputs['sweep']))
        tip_chord=inputs['c']*inputs['taper']
        #first spar 1/4 of the way through
        #sweep angle is considered to be 0

        spars=np.zeros((2,2))
        
        spars[0][0] = (0.25*tip_chord)+inputs['b']*sweep_rad
        spars[0][1] = 0.25*inputs['c']
        spars[1][0] = (0.5*tip_chord)+inputs['b']*sweep_rad
        spars[1][1] = 0.5*inputs['c']
        outputs['spars'] = spars


# In[ ]:


class landing_gear_calcs(om.ExplicitComponent):
    def initialize(self):
        #take constants from group a level above
        self.options.declare('g',prob.model.struct.options['g'])
        self.options.declare('a',prob.model.struct.options['landing_acceleration'])
        self.options.declare('tension_max',prob.model.struct.options['tension_max'])
        
    def setup(self):
        #define inputs (variables to come from other disciplines) and default values
        self.add_input('m', val = 240000)
        self.add_output('landing_tension')
        prob.model.add_constraint('landing_tension', upper = self.options['tension_max'])
        
    def compute(self,inputs,outputs):
        tension = (inputs['m']*self.options['a'])/(2*math.sin(45))
        
        outputs['landing_tension'] = tension


# In[ ]:


class optimisation(om.ExplicitComponent):
    def initialize(self):
        #read in table
        
        #self.options.declare('DummyData', DummyData)
        self.options.declare('WingData', prob.model.struct.options['WingData'])
        
        #take constants from group a level above
        self.options.declare('tensile_stress',prob.model.struct.options['tensile_stress'])
        self.options.declare('compressive_stress',prob.model.struct.options['compressive_stress'])
        self.options.declare('max_deflection',prob.model.struct.options['max_deflection'])
        
    def setup(self):
        #define inputs (variables to come from other disciplines) and default values
        
        self.add_output('b')
        self.add_output('min_weight',shape = (1,13))
        self.add_output('Ix')
        self.add_output('Iy')
        self.add_output('Iz')
        self.add_output('Ixz')
        
    def compute(self,inputs,outputs):
        
        tensile_stress = self.options['tensile_stress']
        compressive_stress = self.options['compressive_stress']
        max_deflection = self.options['max_deflection']
        data_raw = self.options['WingData']
        
        #filtering data using the constraints
        tensile_stress_constraint=data_raw[data_raw.TensileStress.lt(tensile_stress).groupby(data_raw.Weight).transform('all')]
        compressive_constraint=tensile_stress_constraint[tensile_stress_constraint.CompressiveStress.lt(compressive_stress).groupby(tensile_stress_constraint.Weight).transform('all')]
        deflection_constraint=compressive_constraint[
            compressive_constraint.Deflection.lt(max_deflection).groupby(compressive_constraint.Weight).transform('all')]
        data_modified=deflection_constraint
        
        dv= len(data_modified)
        if dv==0:
            min_weight=3622
            half_span=80    #actually full_span 
            root_chord=12.6
            tip_chord=8.4
            Ix=2420160
            Iy=2589780
            Iz=177178
            Ixz=129704
            
        else: 
        
        #determining the chord length and span with the lowest weight for optimisation purposes
            min_weight= data_modified[data_modified["Weight"]==data_modified["Weight"].min()]
        #####################################
           # min_weight = pd.DataFrame([[60,5,4]])
            b=min_weight.iloc[0,0]
            c=min_weight.iloc[0,1]
            tip_chord=min_weight.iloc[0,2]
            Ix=min_weight.iloc[0,7]
            Iy=min_weight.iloc[0,8]
            Iz=min_weight.iloc[0,9]
            Ixz=min_weight.iloc[0,12]
        
        outputs['b'] = b
        outputs['min_weight'] = min_weight
        outputs['Ix']= Ix
        outputs['Iy']= Iy
        outputs['Iz']= Iz
        outputs['Ixz']=Ixz


# In[ ]:


class structures_cog(om.ExplicitComponent):
    def initialize(self):
        #take constants from group a level above
        self.options.declare('mach',prob.model.struct.options['mach'])
        self.options.declare('n_passengers',prob.model.struct.options['n_passengers'])
        self.options.declare('range',prob.model.struct.options['range'])
        
        self.options.declare('M_initial',prob.model.struct.options['M_initial'])
        self.options.declare('pp_mass',prob.model.struct.options['pp_mass'])
        
        self.options.declare('L_cock',prob.model.struct.options['L_cock'])
        self.options.declare('L_ce',prob.model.struct.options['L_ce'])
        self.options.declare('L_c1',prob.model.struct.options['L_c1'])
        
        self.options.declare('tail_hc',prob.model.struct.options['tail_hc'])
        self.options.declare('tail_hs',prob.model.struct.options['tail_hs'])
        self.options.declare('tail_h_angle', prob.model.struct.options['tail_h_angle'])
        
        self.options.declare('tail_vc',prob.model.struct.options['tail_vc'])
        self.options.declare('tail_vh',prob.model.struct.options['tail_vh'])
        self.options.declare('tail_v_angle',prob.model.struct.options['tail_v_angle'])
                
        self.options.declare('root_taper_chord_ratio',prob.model.struct.options['root_taper_chord_ratio'])
        
        self.options.declare('constant_uht',prob.model.struct.options['constant_uht'])
        self.options.declare('fuselage_diameter',prob.model.struct.options['fuselage_diameter'])
        self.options.declare('fuselage_taper',prob.model.struct.options['fuselage_taper'])
        self.options.declare('f_wet_area',prob.model.struct.options['f_wet_area'])
        self.options.declare('thic_chord_rat',prob.model.struct.options['thic_chord_rat'])

        
    def setup(self):
        #define inputs (variables to come from other disciplines) and default values
        
        self.add_input('sweep', val = 0.015)
        self.add_input('c', val=30)
        self.add_input('b', val=80)
        self.add_input('s', val=200) #wing area
        self.add_input('aspect_ratio', val = 12)
        self.add_input('taper_ratio', val = 0.45)
        self.add_input('length_fuselage', val = 30)
        self.add_input('root_x')
        
        #currently trying these outputs - unsure
        self.add_output('cog_structures', shape = (2))
        self.add_output('weight_structures')
        
    def compute(self,inputs,outputs):
        
        mach = self.options['mach']
        n_passengers=self.options['n_passengers']
        
        M_initial=self.options['M_initial']
        pp_mass=self.options['pp_mass']
        
        L_cock= self.options['L_cock']
        L_ce=self.options['L_ce']
        L_c1=self.options['L_c1']
        
        tail_hc=self.options['tail_hc']
        tail_hs=self.options['tail_hs']
        tail_h_angle=self.options['tail_h_angle']
        
        tail_vc=self.options['tail_vc']
        tail_vh=self.options['tail_vh']
        tail_v_angle=self.options['tail_v_angle']
        
        thic_chord_rat = self.options['thic_chord_rat']
        
        length_fuselage=inputs['length_fuselage']
        s = inputs['s']
        b = inputs['b']
        c = inputs['c']
        taper_ratio = inputs['taper_ratio']
        aspect_ratio = inputs['aspect_ratio']
        
        sweep = inputs['sweep']
        root_taper_chord_ratio= self.options['root_taper_chord_ratio']
        
        constant_uht=self.options['constant_uht']
        f_diameter = self.options['fuselage_diameter']
        f_wet_area = self.options['f_wet_area']
        fuselage_taper = self.options['fuselage_taper']
        
        wing_x = inputs['root_x']
        
        forward_bulk = f_diameter/5 #x axis 
        #wings
        flap_area= s*0.172385620 #wing flap area
        sweep= math.radians(sweep)
        y_wing = f_diameter-0.2 #wing position
        incident_angle = math.radians(0) # wing incidence angle
        
        #landing gear
        x_landing = L_cock #distance to front of landing gear
        
        #tail
        tail_h_ar= tail_hs/tail_hc #horizontal aspect ratio
        tail_v_ar= tail_vh/tail_vc #vertical aspect ratio
        
        #realtive positions
        x_spar_fwd=b +0.2376*c #distance to front spar
        x_spar_aft = b +0.56937*c #distance to aft spar
        x_trail_Edge = b + c # distance to trailing edge
        tail_x=length_fuselage-tail_hc-0.5
        
        length_main_lg= 0.75*length_fuselage*0.0254 #length of main landing gear
        length_nose_lg = 0.7*length_fuselage*0.0254 #length of nose landing gear
        
        #WEIGHTS
        design_gw = 2.20462262*M_initial #gross design weight pounds
        load_limit=1.5*2 #1.5*load factor for limit
        l_fuse_ft= length_fuselage*3.281 # fueslage length in feet
        f_wet_area_ft2= f_wet_area*10.764 #fueslage wetted area in ft2
        b_ft=b*3.281 #wing span in ft
        constant_b= 0.75*((1+2*taper_ratio)/(1+taper_ratio))*((b_ft*math.tan(sweep))/1) #constant for wing span
        
        #wing
        s_ft = s*3.28084 #wing area in ft
        control_surface = flap_area*10.7639 #area of control surface in ft2
        
        weight_wing_lb=0.0051*((design_gw*load_limit)**0.557)*(s_ft**0.649)*(aspect_ratio**0.5)*((thic_chord_rat)**-0.4)*((1+taper_ratio)**0.1)
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
        
        
        weight_tailv_lb = 0.0026*(2**0.225)*(design_gw*0.536)*(tail_vc_ft**-0.5)*(tail_area_v**0.5)*(gyration_v**0.875)*((math.cos(sweep))**-1)*(tail_v_ar**0.35)*(root_taper_chord_ratio**-0.5)
        weight_tailv = weight_tailv_lb**0.453592
        
        #coordinates
        
        
        #wings
        average_chord = c
        
        weight_wingxy=(wing_x+0.4*average_chord+0.25*b*math.tan(sweep), y_wing)
        cg_wing =[(weight_wingxy[0]*weight_wing),(weight_wingxy[1]*weight_wing)]
        
        #tail
        #horizontal
        tail_h_angle=math.radians(tail_h_angle) #deg to rad
        weight_tailhxy=[tail_x+0.3*tail_hc+0.25*tail_hs*math.tan(tail_h_angle),f_diameter]
        weight_tailhxy0 = weight_tailhxy[0]
        weight_tailhxy1 = weight_tailhxy[1]
        cgh=[(weight_tailhxy0*weight_tailv),(weight_tailhxy1*weight_tailv)]
        
        #vertical
        tail_v_angle=math.radians(tail_v_angle) #deg to rad
        weight_tailvxy=[tail_x+0.3*tail_vc+0.25*tail_vh*math.tan(tail_v_angle),f_diameter+0.5*tail_vh]
        weight_tailvxy0 = weight_tailvxy[0]
        weight_tailvxy1 = weight_tailvxy[1]
        cgv=[(weight_tailvxy0*weight_tailv),(weight_tailvxy1*weight_tailv)]
        
       # cog_structures = np.zeros(2)
       # cog_structures[0] = cg_wing[0]+cgh[0]+cgv[0]
       # cog_structures[1] = cg_wing[1]+cgh[1]+cgv[1]
        

       # weight_structures = weight_wingxy[0] + weight_tailhxy[0] + weight_tailvxy[0]
        
        weight_structures=weight_wing+weight_tailv
        cg_total=[(cg_wing[0]+cgh[0]+cgv[0]),(cg_wing[1]+cgh[1]+cgv[1])]
        cog_structures=[(cg_total[0]/weight_structures),(cg_total[1]/weight_structures)]
        
        outputs['cog_structures'] = cog_structures
        outputs['weight_structures'] = weight_structures
        

