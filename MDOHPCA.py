#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#    Multi-Disciplinary Design Optimisation of a Hydrogen Powered Commercial Aircraft (MDOHPCA)
#    this is the master file that runs each of the individual disciplines 
#    It sets up each of the analyses, then converts it into an openMDAO problem with a certain structure, then solves for minimum mass



# In[ ]:


import pandas as pd
import numpy as np
import openmdao.api as om


# In[ ]:


class MDOHPCA(om.Group):
    """
    Top level group containing the MDO.
    """           
    def setup(self):
        """
        set up each person's disciplines that can contain analyses
        """
        self.add_subsystem('aero', aero())
        self.add_subsystem('sys', sys())
        self.add_subsystem('stab', stab())
        self.add_subsystem('struct', struct())
        
        self.add_subsystem('mass', mass())
        
    def configure(self):
        #promote all variables (lazy option, they can be connected individually)
        self.promotes('sys',any=['*'])
        self.promotes('struct',any=['*'])
        
        self.promotes('mass',any=['*'])


# In[ ]:


class mass(om.ExplicitComponent):
    def setup(self):
        self.add_input('systems_mass')
        self.add_input('weight_structures')
        self.add_input('cog_structures', shape = (2,))
        self.add_input('systems_CG')
        
        self.add_output('m')
        prob.model.add_objective('m')
        
        self.add_output('total_CG', shape = (1,2))
        
    def setup_partials(self):
        self.declare_partials('*', '*', method = 'fd')
        
    def compute(self,inputs,outputs):
        #calculation of mass from systems and structures
        mass = inputs['systems_mass'] + inputs['weight_structures']
        
        outputs['m'] = mass
        
        #calculation of CG from systems and structures
        struct_CG = inputs['cog_structures']
        sys_CG = inputs['systems_CG']
        
        struct_comp_x = struct_CG[0]*inputs['weight_structures']
        struct_comp_y = struct_CG[1]*inputs['weight_structures']
        
        sys_comp_x = sys_CG*inputs['systems_mass']
        
        total_CG_x = (struct_comp_x + sys_comp_x)/mass
        total_CG_y = struct_comp_y/mass
        
        total_CG = [total_CG_x, total_CG_y]
        outputs['total_CG'] = total_CG


# In[ ]:


get_ipython().run_line_magic('run', 'struct/structsetup.py')
get_ipython().run_line_magic('run', 'stab/stabsetup.py')
get_ipython().run_line_magic('run', 'systems/systemssetup.py')
get_ipython().run_line_magic('run', 'aero/aerosetup_2.py')


# In[ ]:


prob = om.Problem(model = MDOHPCA())
#prob.model.add_subsystem('aero', om.Group())


# In[ ]:


#setting inputs to reduce ambiguities where inputs differ in different disciplines
prob.model.set_input_defaults('CL', val = 0.6)
prob.model.set_input_defaults('sweep', val = 15)
prob.model.set_input_defaults('c', val = 4)
prob.model.set_input_defaults('b', val = 60)
prob.model.set_input_defaults('root_x', val = 10)
prob.model.set_input_defaults('m', val = 240000)


prob.model.add_design_var('tank_ratio', lower = 0.3, upper = 1)
prob.model.set_input_defaults('tank_ratio', val = 1)


# In[ ]:


#declare problem options (driver, optimiser)

#prob.nonlinear_solver = om.NonlinearBlockGS()
#prob.nonlinear_solver = om.NewtonSolver(solve_subsystems=False)
prob.model.nonlinear_solver = om.NonlinearBlockGS()
prob.model.nonlinear_solver.options['iprint'] = 2
prob.model.nonlinear_solver.options['maxiter'] = 200

#prob.nonlinear_solver.options['maxiter'] = 100

prob.driver = om.ScipyOptimizeDriver()
prob.driver.options['optimizer'] = 'SLSQP'
prob.driver.options['tol'] = 1e-5

#setup problem to be run
prob.setup()


# In[ ]:


prob.run_model()


# In[ ]:


#prob.run_driver()


# In[ ]:


prob.model.get_design_vars()


# In[ ]:


prob.model.get_val('tank_ratio')


# In[ ]:


#get properties of specific constants in disciplines
#prob.model.sys.options._dict['pp_mass']

prob.model.list_outputs()


# In[ ]:


mass = 6
naca_series = '23014'
max_allowable_wing_span = 90
keyaeronumbers = pd.DataFrame(list([mass,naca_series,max_allowable_wing_span]))
keyaeronumbers


# In[ ]:


keyaeronumbers = pd.read_csv("keyAero.dat")


# In[ ]:


f = pd.DataFrame(keyaeronumbers)


# In[ ]:


str(int(f.iloc[1]))


# In[ ]:




