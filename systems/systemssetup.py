#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#import necessary package and collect initial variables
import openmdao.api as om
import pandas as pd
import numpy as np


# In[ ]:


#define discipline group, get constants and call analyses 
class sys(om.Group):
    def initialize(self):
        
        #read all excel constants        
        initial_vars = pd.read_excel('constants.xlsx')
        
        #sort for just the systems constants
        all_sys_vars = np.where(initial_vars['sys'] == True, [initial_vars['variable name'], initial_vars['value'], initial_vars['description']], None)
        
        #create a DataFrame for the systems variables to pull from in individual analyses
        all_sys_vars = pd.DataFrame(all_sys_vars).dropna(axis = 1).transpose()
        
        #pass variables into the correct openMDAO format, passing name, description and value
        i=0
        for variable in all_sys_vars[0]:
            value = all_sys_vars.iloc[i][1]
            description = all_sys_vars.iloc[i][2]
            self.options.declare(variable, default = value, desc = description)
            i+=1
        i=0
        
    def setup(self):
        #set analysis (subsystem) names and what inputs they take and give
        self.add_subsystem('fuel_mass_calc', fuel_mass_calc())
        self.add_subsystem('tanks', tanks())
        self.add_subsystem('pipes', pipes())
        self.add_subsystem('engine', engine())
        self.add_subsystem('actuator', actuator())
        self.add_subsystem('systems_roundup', systems_roundup())
        #can optionally promote the variables at this step
        
        
    def configure(self):
        #promote all variables (lazy option, they can be connected individually)
        self.promotes('fuel_mass_calc',any=['*'])
        self.promotes('tanks',any=['*'])
        self.promotes('pipes',any=['*'])
        self.promotes('engine',any=['*'])
        self.promotes('actuator',any=['*'])
        self.promotes('systems_roundup',any=['*'])


# In[ ]:


class fuel_mass_calc(om.ExplicitComponent):
    def initialize(self):
        #take constants from group a level above
        self.options.declare('eta',prob.model.sys.options['eta'])
        self.options.declare('R',prob.model.sys.options['range'])
        self.options.declare('SFC',prob.model.sys.options['TSFC'])
        self.options.declare('M_initial',prob.model.sys.options['M_initial'])
        self.options.declare('v',prob.model.sys.options['v'])
        self.options.declare('g',prob.model.sys.options['g'])
        self.options.declare('S',prob.model.sys.options['S'])
        self.options.declare('rho',prob.model.sys.options['rho'])
        
    def setup(self):
        #define inputs (variables to come from other disciplines) and default values
        self.add_input('CD', val = 0.02)
        self.add_input('CL', val = 0.6)
        self.add_input('m')
        
        self.add_output('fuel_mass')        
        prob.model.add_constraint('fuel_mass', lower = 20000)
    def setup_partials(self):
        self.declare_partials('*', '*', method = 'fd')
    
    def compute(self,inputs,outputs):
        #define variables used in equations and where they come from (constant or input)
        eta = self.options['eta']
        R = self.options['R']*1000
        SFC = self.options['SFC']
        v = self.options['v']
        S = self.options['S']
        rho = self.options['rho']
        
        W0 = inputs['m']
        CD = inputs['CD']
        CL = inputs['CL']
        
        #fuel mass fraction derived from breguet range equation
        #Wf = np.exp(-R*SFC/(eta*CL/CD))
        
        #using flight mech equation for jet aircraft:
        root_terms = (2/(rho*S))*CL
        rt_W1 = W0**0.5 - ((R*SFC*CD)/(2*(root_terms**0.5)))
        
        W1 = rt_W1**2
        
        outputs['fuel_mass'] = W0 - W1


# In[ ]:


class tanks(om.ExplicitComponent):
    def initialize(self):
        #take constants from group a level above
        self.options.declare('LH2_rho',prob.model.sys.options['LH2_rho'])
        self.options.declare('tank_rho',prob.model.sys.options['tank_rho'])
        self.options.declare('insul_rho',prob.model.sys.options['insul_rho'])
        self.options.declare('tank_t',prob.model.sys.options['tank_t'])
        self.options.declare('flux',prob.model.sys.options['max_tank_flux'])
        self.options.declare('L_cock',prob.model.sys.options['L_cock'])
        self.options.declare('L_ce',prob.model.sys.options['L_ce'])
        self.options.declare('L_c1',prob.model.sys.options['L_c1'])
        self.options.declare('pp_mass',prob.model.sys.options['pp_mass'])
        self.options.declare('insul_kappa',prob.model.sys.options['insul_kappa'])
        self.options.declare('To',prob.model.sys.options['To'])
        self.options.declare('LH2_T',prob.model.sys.options['LH2_T'])
        self.options.declare('fuselage_diameter',prob.model.sys.options['fuselage_diameter'])

    def setup(self):
        #define inputs (variables to come from other disciplines) and default values
        self.add_input('fuel_mass', val = 30000)
        self.add_input('total_CG', shape = (1,2))

        self.add_input('tank_ratio') #ratio of tank diameter to fuselage diameter
        prob.model.add_constraint('tank_ratio',lower = 0.3, upper = 1)
        
        self.add_output('tank1_mass_full')
        self.add_output('tank1_mass_empty')
        self.add_output('tank2_mass_full')
        self.add_output('tank2_mass_empty')
        self.add_output('CGtp')
        self.add_output('tank1_x')
        self.add_output('tank2_x')
        self.add_output('tank_config_n')
        self.add_output('tank_i_t')
        self.add_output('length_fuselage')

    def compute(self,inputs,outputs):
        LH2_rho = self.options['LH2_rho']
        tank_rho = self.options['tank_rho']
        insul_rho = self.options['insul_rho']
        t = self.options['tank_t']
        flux = self.options['flux']
        L_cock = self.options['L_cock']
        L_ce = self.options['L_ce']
        L_c1 = self.options['L_c1']
        pp_mass = self.options['pp_mass']
        k = self.options['insul_kappa']
        To = self.options['To']
        LH2_T = self.options['LH2_T']
        
        f_mass = inputs['fuel_mass']
        R = inputs['tank_ratio']*self.options['fuselage_diameter']/2 #outer radius of tank
        
        
        f_vol = f_mass*(1.072/LH2_rho) #calculate fuel volume, uses allowance factor from book
        
        tV1 = f_vol*3/5 #ratios of tank volume using book values of roughly 3:2 forward to aft
        tV2 = f_vol*2/5
        
        #1D heat relationship through tank walls using Fourier equation, with a 2cm buffer:
        insul_t = (k*(To-LH2_T)/flux)+0.02
        outputs['tank_i_t'] = insul_t
        R_i = R-insul_t #inner radius of tanks
        
        tL1_c = (tV1-((4/3)*(np.pi*R_i**3)/1.6))/(np.pi*R_i**2) #length of cylindrical section of tank 1
        tL2_c = (tV2-((4/3)*(np.pi*R_i**3)/1.6))/(np.pi*R_i**2) #ellipse a/b ratio of 1.6 is used from book for end sections
        
        tL1 = tL1_c+2*(R_i/1.6) #total lengths of tanks 
        tL2 = tL2_c+2*(R_i/1.6)
        
        A1 = (np.pi*2*R_i*tL1_c)+4*np.pi*(((R_i**2)**1.6+2*((R_i**2/1.6)**1.6))/3)**(1/1.6) #surface area of tanks
        A2 = (np.pi*2*R_i*tL2_c)+4*np.pi*(((R_i**2)**1.6+2*((R_i**2/1.6)**1.6))/3)**(1/1.6) #uses surface area of ellipsoid for end parts
        
        tm1 = A1*(t*tank_rho+insul_t*insul_rho) #empty tank masses
        tm2 = A2*(t*tank_rho+insul_t*insul_rho)
        
        outputs['tank1_mass_empty'] = tm1
        outputs['tank2_mass_empty'] = tm2
        fulltank1 = tm1+f_mass*0.6
        fulltank2 = tm2+f_mass*0.4
        outputs['tank1_mass_full'] = fulltank1
        outputs['tank2_mass_full'] = fulltank2
        
        t_c_m = fulltank1 + fulltank2 + pp_mass #total full mass of tanks and payload
        
        #centre of gravity of tanks and payload using different configurations
        #configuration 1
        CGtp1 = ((L_cock+0.5*tL1)*fulltank1+(0.2*pp_mass*(L_cock+tL1+(L_c1/2)))+fulltank2*(L_cock+tL1+L_c1+0.5*tL2)+(0.8*pp_mass*(L_cock+tL1+L_c1+tL2+0.5*L_ce)))/t_c_m 
        tank1_x1 = L_cock + 0.5*tL1
        tank2_x1 = L_cock + tL1 + L_c1 + 0.5*tL2
        #configuration 2
        CGtp2 = ((L_cock+0.5*tL1)*fulltank1+(0.8*pp_mass*(L_cock+tL1+(L_ce/2)))+fulltank2*(L_cock+tL1+L_ce+0.5*tL2)+(0.2*pp_mass*(L_cock+tL1+L_ce+tL2+0.5*L_c1)))/t_c_m
        tank1_x2 = L_cock + 0.5*tL1
        tank2_x2 = L_cock + tL1 + L_ce + 0.5*tL2
        #configuration 3
        CGtp3 = ((L_cock+0.5*L_c1)*0.2*pp_mass+(fulltank1*(L_cock+L_c1+(tL1/2)))+fulltank2*(L_cock+tL1+L_c1+0.5*tL2)+(0.8*pp_mass*(L_cock+tL1+L_c1+tL2+0.5*L_ce)))/t_c_m
        tank1_x3 = L_cock + L_c1 + 0.5*tL1
        tank2_x3 = L_cock + tL1 + L_c1 + 0.5*tL2
        #configuration 4
        CGtp4 = ((L_cock+0.5*L_ce)*0.8*pp_mass+(fulltank1*(L_cock+L_ce+(tL1/2)))+fulltank2*(L_cock+tL1+L_ce+0.5*tL2)+(0.2*pp_mass*(L_cock+tL1+L_ce+tL2+0.5*L_c1)))/t_c_m
        tank1_x4 = L_cock + L_ce + 0.5*tL1
        tank2_x4 = L_cock + tL1 + L_ce + 0.5*tL2
        #configuration 5
        CGtp5 = ((L_cock+0.5*tL1)*fulltank1+(fulltank2*(L_cock+tL1+(tL2/2)))+pp_mass*(L_cock+tL1+0.5*L_c1+tL2+0.5*L_ce))/t_c_m
        tank1_x5 = L_cock+0.5*tL1
        tank2_x5 = L_cock + tL1 + 0.5*tL2
        #configuration 6
        CGtp6 = ((L_cock+0.5*L_c1)*0.2*pp_mass+(fulltank1*(L_cock+L_c1+(tL1/2)))+fulltank2*(L_cock+tL1+L_ce+L_c1+0.5*tL2)+(0.8*pp_mass*(L_cock+tL1+L_c1+0.5*L_ce)))/t_c_m
        tank1_x6 = L_cock + L_c1 + 0.5*tL1
        tank2_x6 = L_cock + L_c1 + tL1 + L_ce + 0.5*tL2
        #configuration 7
        CGtp7 = ((L_cock+0.5*tL1)*fulltank1+(pp_mass*(L_cock+tL1+0.5*L_c1+0.5*L_ce))+fulltank2*(L_cock+tL1+L_c1+L_ce+0.5*tL2))/t_c_m 
        tank1_x7 = L_cock + 0.5*tL1
        tank2_x7 = L_cock + L_c1 + tL1 + L_ce + 0.5*tL2
        #configuration 8
        CGtp8 = ((L_cock+0.5*L_c1+0.5*L_ce)*pp_mass+(fulltank1*(L_cock+L_c1+L_ce+(tL1/2)))+fulltank2*(L_cock+tL1+L_ce+L_c1+0.5*tL2))/t_c_m
        tank1_x8 = L_cock + L_c1 + L_ce + 0.5*tL1
        tank2_x8 = L_cock + L_c1 + tL1 + L_ce + 0.5*tL2
        
        total_CG = inputs['total_CG'][0]
        
        CGtp_array = [CGtp1, CGtp2, CGtp3, CGtp4, CGtp5, CGtp6, CGtp7, CGtp8]
        tank1_x_array = [tank1_x1, tank1_x2, tank1_x3, tank1_x4, tank1_x5, tank1_x6, tank1_x7, tank1_x8]
        tank2_x_array = [tank2_x1, tank2_x2, tank2_x3, tank2_x4, tank2_x5, tank2_x6, tank2_x7, tank2_x8]
        
        differences = abs(CGtp_array - (np.ones((8,1))*total_CG[0]))
        config = np.argmin(differences) 
        
        outputs['tank_config_n'] = config + 1
        #optimal configuration output
        outputs['CGtp'] = CGtp_array[config]
        outputs['tank1_x'] = tank1_x_array[config]
        outputs['tank2_x'] = tank2_x_array[config]
        
        outputs['length_fuselage'] =  tL1 + tL2 + L_cock + L_c1 + L_ce #[tl1, tl2, L_cock, L_c1, L_ce]
        #add conditions for each tank being empty to see how the CG changes
        #could optionally add tank insulation width study


# In[ ]:


class pipes(om.ExplicitComponent):
    def initialize(self):
        self.options.declare('insul_kappa',prob.model.sys.options['insul_kappa'])
        self.options.declare('insul_rho',prob.model.sys.options['insul_rho'])
        self.options.declare('di',prob.model.sys.options['di'])
        self.options.declare('L_H',prob.model.sys.options['L_H'])
        self.options.declare('mol_H',prob.model.sys.options['mol_H'])
        self.options.declare('LH2_Cp',prob.model.sys.options['LH2_Cp'])
        self.options.declare('LH2_max_T',prob.model.sys.options['LH2_max_T'])
        self.options.declare('LH2_T',prob.model.sys.options['LH2_T'])
        self.options.declare('boost_eta',prob.model.sys.options['boost_eta'])
        self.options.declare('boost_m_eta',prob.model.sys.options['boost_m_eta'])
        self.options.declare('boost_P',prob.model.sys.options['boost_P'])
        self.options.declare('boost_power_max',prob.model.sys.options['boost_power_max'])
        self.options.declare('m_dot_max',prob.model.sys.options['m_dot_max'])
        self.options.declare('m_dot',prob.model.sys.options['m_dot'])
        self.options.declare('To',prob.model.sys.options['To'])
        
    def setup(self):
        self.add_input('engine_x', val = 30)
        self.add_input('engine_y', val = 10)
        self.add_input('tank1_x', val = 8)
        self.add_input('tank2_x', val = 40)        
        self.add_output('pipe_CG')
        self.add_output('pipe_mass')
        self.add_output('pipe_i_t')
        
    def compute(self, inputs, outputs):
        #find maximum length of pipe (from tank to engine) and determine the heat gain with different insulation thicknesses
        #determine total pipe length and mass to determine pipe CG
        tank1_x = inputs['tank1_x']   #tank 1 x position
        tank2_x = inputs['tank2_x']
        engine_x = inputs['engine_x'] #both engines x pos
        engine_y = inputs['engine_y'] #both engines distance y (spanwise) from centreline
        
        
        maxpipe1_2 = abs(tank1_x - engine_x)+abs(engine_y) #longest route will be through cross-feed, tank to alternate engine
        maxpipe2_1 = abs(tank2_x - engine_x)+abs(engine_y)
        max_p = max(maxpipe1_2,maxpipe2_1) #maximum pipe length
        
        #thermal relationships using takeoff conditions for largest efficiency losses and cruise for largest pipeline losses
        T0 = self.options['LH2_T']
        #maximum allowed heat gain in J/kg
        Qmax = (self.options['LH2_Cp']*(self.options['LH2_max_T']-T0))+(self.options['L_H']/self.options['mol_H'])
        
        #enthalpy rise from boost pump motor (from [3] pg 129)
        del_h_m_imp = (1-self.options['boost_m_eta'])*(550/778)*self.options['boost_power_max']/(self.options['m_dot_max']*2.20462262)
        del_h_m = 2324.446*del_h_m_imp
        
        #enthalpy rise from boost pump pump
        del_h_p_imp = self.options['boost_P']*(0.0175 + 0.0412*((1/self.options['boost_eta'])-1))
        del_h_p = 2324.446*del_h_p_imp
        
        #remaining allowed heat gain, uses cruise for maximum losses
        h_L = Qmax - del_h_m - del_h_p
        pipe_heat_losses = h_L * self.options['m_dot'] #in watts
        
        d_i = self.options['di'] #inner diameter
        r_i = d_i/2 #inner radius
        k = self.options['insul_kappa']
        To = self.options['To']
        
        #minimum insulation thickness required + 1cm buffer
        r_o = r_i * np.exp((2*np.pi*k*max_p*(To-T0))/pipe_heat_losses) + 0.01        
        d_o = 2*r_o #outside diameter of pipes
        al_t = 0.0041   #radiation and structural aluminium jacket thickness
        inner_t = al_t   #inner pipe thickness, made of stainless steel

        pipe_insul_t = r_o - al_t - inner_t - r_i    #pipe insulation thickness based on thermal relationship
        pipe_L = 2*(maxpipe1_2+maxpipe2_1)+abs(tank1_x - tank2_x)     #total length of pipe estimate
        
        #components from jacket, insulation and inner
        insul_area = (((d_o/2)-al_t)**2-(inner_t+(d_i/2))**2)*np.pi     #area of inuslation 
        insul_kg_per_m = insul_area*self.options['insul_rho']
        
        jacket_area = np.pi*d_o*al_t
        jacket_kg_per_m = jacket_area*2823
        
        inner_area = np.pi*d_i*inner_t
        inner_kg_per_m = inner_area*7750

        
        #pipes CG, assuming symmetrical so CG_y = 0, and masses are proportional to length so only lengths and x_pos are needed
        #lengths of different sections 
        #section from tank 1 to engine 1
        L_1_1_x = abs(tank1_x - engine_x)    #length of section from tank 1 to engine 1 in x direction only
        x_1_1_x = (tank1_x + engine_x)/2   #x position of midpoint of section from tank 1 to engine 1
        L_1_1_y = engine_y
        x_1_1_y = engine_x
        
        #tank 1 to engine 2 has the same absolute properties as engine positions are symetrical
        L_1_2_x = L_1_1_x
        x_1_2_x = x_1_1_x
        L_1_2_y = L_1_1_y
        x_1_2_y = x_1_1_y
        
        #tank 2 to engine 1
        L_2_1_x = abs(tank2_x - engine_x)    #length of section from tank 1 to engine 1 in x direction only
        x_2_1_x = (tank2_x + engine_x)/2   #x position of midpoint of section from tank 1 to engine 1
        L_2_1_y = engine_y
        x_2_1_y = engine_x
        
        #tank 2 to engine 2
        L_2_2_x = L_2_1_x
        x_2_2_x = x_2_1_x
        L_2_2_y = L_2_1_y
        x_2_2_y = x_2_1_y
        
        #cross feed between tanks 
        L_x = abs(tank1_x - tank2_x)
        x_x = (tank1_x + tank2_x)/2  
        
        #pipes CG
        pipes_CG_comp = (L_1_1_x*x_1_1_x)+(L_1_1_y*x_1_1_y)+(L_1_2_x*x_1_2_x)+(L_1_2_y*x_1_2_y)+(L_2_1_x*x_2_1_x)+(L_2_1_y*x_2_1_y)+(L_2_2_x*x_2_2_x)+(L_2_2_y*x_2_2_x)
        pipe_L = L_1_1_x + L_1_1_y + L_1_2_x + L_1_2_y + L_2_1_x + L_2_1_y + L_2_2_x + L_2_2_y
        outputs['pipe_CG'] = pipes_CG_comp/pipe_L
            
        #total pipe mass
        outputs['pipe_mass'] = pipe_L*(inner_kg_per_m+jacket_kg_per_m+insul_kg_per_m)
        outputs['pipe_i_t'] = d_o


# In[ ]:


class engine(om.ExplicitComponent):
    def initialize(self):
        self.options.declare('T_W',prob.model.sys.options['T_W'])
        self.options.declare('T_m',prob.model.sys.options['T_m'])
    
    def setup(self):
        self.add_input('fuselage_diameter', val = 8)
        self.add_input('wingspan', val = 80)
        self.add_input('sweep', val = 30) #sweep angle (degrees)
        self.add_input('root_x', val = 15) #x-coord of root of wing
        self.add_input('m')
        
        self.add_output('engine_mass')
        self.add_output('engine_x')
        self.add_output('engine_y')
        
    def compute(self,inputs,outputs):
        #engine size based off the original constraint diagram calculation
        T_W = self.options['T_W']
        T_m = self.options['T_m']
        engine_mass = T_m*T_W*inputs['m']*9.81/2000  #mass per engine
        outputs['engine_mass'] = engine_mass
        
        #decide on engine position, will be an implicit function with the pipes calculation in previous cell as well as the bending moment from structures
        #ensure to add bounds so the engines are not too far or too close
        engine_y = ((inputs['wingspan']/2)-(inputs['fuselage_diameter']/2))/2 #estimate of halfway between fuselage and wing tip, will need to work with structures to make implicit function
        outputs['engine_y'] = engine_y
        sweep_r = inputs['sweep']*np.pi/180  #sweep in radians
        outputs['engine_x'] = engine_y*np.tan(sweep_r)+inputs['root_x']


# In[ ]:


class actuator(om.ExplicitComponent):
    #box for control surface requirement calculations 
    #calculate mass of systems based on load requirements 
    #use positions and mass to calculate a CG
    def initialize(self):
        self.options.declare('kgpN',prob.model.sys.options['kgpN']) #hydraulic actuator mass per newton of force available
    
    def setup(self):
        #flap and aileron locations and lift requirements along the span (y direction)
        self.add_input('flap_pos_y', val = 12) #y-distance from centreline of fuselage, can be passed through as an array if more than 2 are used
        self.add_input('flap_req', val = 50000) #total lift requirement from flaps (N)
        self.add_input('flap_n', val = 2) #number of flaps
        self.add_input('flap_L', val = 1) #flap length used for torque
        self.add_input('aileron_pos_y', val = 35) 
        self.add_input('aileron_req', val = 5000)
        self.add_input('aileron_n', val = 2) #number of ailerons
        self.add_input('aileron_L', val = 0.5) #aileron length
        self.add_input('wing_t', val = 1) #wing thickness also used for torque
        self.add_input('sweep', val = 30) #sweep angle (degrees)
        self.add_input('root_x', val = 10) #x-coord of root of wing
        
        
        self.add_output('each_flap_actuator_mass') #will be passed to structures and systems roundup
        self.add_output('each_aileron_actuator_mass')
        self.add_output('actuators_CG_x')
        
    def compute(self,inputs,outputs):
        moment_F = (inputs['flap_req']/inputs['flap_n'])*(inputs['flap_L']/2)   #torque required assuming the force on flap acts halfway along the flap
        arm_L = inputs['wing_t']/2
        flap_hydraulic_force = moment_F/arm_L  #force required from flap actuator assuming the distance from the hinge to the actuator line is half the width of the wing
        
        moment_a = (inputs['aileron_req']/inputs['aileron_n'])*(inputs['aileron_L']/2)   #torque required assuming the force on flap acts halfway along the flap
        aileron_hydraulic_force = moment_a/arm_L
        
        kgpN = self.options['kgpN']
        
        outputs['each_flap_actuator_mass'] = flap_hydraulic_force*kgpN
        outputs['each_aileron_actuator_mass'] = aileron_hydraulic_force*kgpN
        
        sweep_r = inputs['sweep']*np.pi/180  #sweep in radians
        flap_pos_x = inputs['root_x']+(inputs['flap_pos_y']*np.tan(sweep_r))
        aileron_pos_x = inputs['root_x']+(inputs['aileron_pos_y']*np.tan(sweep_r))
        
        #component terms are the position, which can be an array of multiple x positions but only for 1 wing and the associated mass
        outputs['actuators_CG_x'] = (flap_pos_x*(flap_hydraulic_force*kgpN)*2 + aileron_pos_x*(aileron_hydraulic_force*kgpN)*2)/(flap_hydraulic_force*kgpN*inputs['flap_n'] +aileron_hydraulic_force*kgpN*inputs['aileron_n'])
        


# In[ ]:


class systems_roundup(om.ExplicitComponent):
    #add additional components that are a fixed mass (and probably position too)
    #collate other systems analyses to obtain one single systems mass and CG to pass to structures
    
    def setup(self):
        #inputs from other systems calcs
        self.add_input('each_flap_actuator_mass', val = 60) 
        self.add_input('each_aileron_actuator_mass', val = 10) 
        self.add_input('flap_n', val = 2) 
        self.add_input('aileron_n', val = 2) 
        self.add_input('actuators_CG_x', val = 35) 
        
        self.add_input('engine_mass', val = 2000) 
        self.add_input('engine_x', val = 35) 
        
        self.add_input('pipe_mass', val = 500) 
        self.add_input('pipe_CG', val = 35) 
        
        self.add_input('CGtp', val = 35) 
        self.add_input('tank1_mass_full', val = 3000)
        self.add_input('tank2_mass_full', val = 2000) 
        self.add_input('tank1_x', val = 15) 
        self.add_input('tank2_x', val = 40) 
        
        
        self.add_input('fuselage_length', val = 80) 
        
        
        self.add_output('systems_mass') 
        self.add_output('systems_CG') 
        
    def compute(self,inputs,outputs):
        #actuators
        actuators_mass = inputs['each_flap_actuator_mass']*inputs['flap_n']+inputs['each_aileron_actuator_mass']*inputs['aileron_n']
        actuator_CG_x = inputs['actuators_CG_x']
        
        #total engine
        engines_mass = inputs['engine_mass']*2 + 4*18.7 #engines plus high pressure engine-mounted pumps
        engine_CG_x = inputs['engine_x']
        
        #full tanks
        tanks_mass = inputs['tank1_mass_full']+inputs['tank2_mass_full']
        full_tanks_CG_x = inputs['CGtp']
        
        #empty tank 1, empty tank 2, both tanks empty
        #####
        
        #pipes
        pipes_mass = inputs['pipe_mass']
        pipes_CG_x = inputs['pipe_CG']
        
        #additional systems
        ## Vent boil off system is placed at the back 
        vent_mass = 60
        vent_x = inputs['fuselage_length']-2
        
        ## Boost pumps are placed at the engines
        boost_pump_mass = 6*28.8  #multistage centrifugal pump for great efficiency and immersion in LH2, with an electric motor driving it, 3 per tank
        boost_pump_CG_x = (inputs['tank1_x']+inputs['tank1_x'])/2
        
        #APU also at the back
        APU_mass = 150
        APU_x = inputs['fuselage_length']-5
        
        #fuselage shell mass 
        fuselage_mass = inputs['fuselage_length']*200
        
        #totals
        ##product of masses and x-positions
        sys_CG_comp = (actuators_mass*actuator_CG_x) + (engines_mass*engine_CG_x) + (tanks_mass*full_tanks_CG_x) + (pipes_mass*pipes_CG_x) + (vent_mass*vent_x) + (boost_pump_mass*boost_pump_CG_x) + (APU_mass*APU_x)
        ##total mass
        sys_mass = actuators_mass + engines_mass + tanks_mass + pipes_mass + vent_mass + boost_pump_mass + APU_mass + fuselage_mass
        
        outputs['systems_mass'] = sys_mass
        outputs['systems_CG'] = sys_CG_comp/sys_mass

