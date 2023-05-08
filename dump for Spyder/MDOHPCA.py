# -*- coding: utf-8 -*-
"""
Created on Sun May  7 17:18:17 2023

@author: Jai
"""
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
        
        self.add_subsystem('sys', sys())
        self.add_subsystem('stab', stab())
        self.add_subsystem('struct', struct())
        self.add_subsystem('aero', aero())
        
        self.add_subsystem('mass', mass())
        self.set_order(['mass','aero','stab','struct','sys']) #aero inputs necessary first before others
        
    def configure(self):
        #promote all variables (lazy option, they can be connected individually)
        self.promotes('sys',any=['*'])
        self.promotes('struct',any=['*'])
        self.promotes('aero',any=['*'])
        self.promotes('stab',any=['*'])
        
        self.promotes('mass',any=['*'])

        
# In[ ]:


class mass(om.ExplicitComponent):
    def setup(self):
        self.add_input('systems_mass')
        self.add_input('weight_structures')
        self.add_input('cog_structures', shape = (2,))
        self.add_input('systems_CG')
        
        self.add_output('m')
        
        
        self.add_output('total_CG', shape = (1,2))
        
    def setup_partials(self):
        self.declare_partials('*', '*', method = 'fd')
        
    def compute(self,inputs,outputs):
        #calculation of mass from systems and structures
        mass = inputs['systems_mass'] + inputs['weight_structures']
        if mass < 10000:
            mass = 50000 #initialize + prevent divergence
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
from subprocess import run
#run('python struct\structsetup.py')
with open('systemssetup.py', 'r') as sy:
    exec(sy.read())
with open('aerosetup.py', 'r') as ae:
    exec(ae.read())
with open('stabsetup.py', 'r') as sta:
    exec(sta.read())
with open('structsetup.py', 'r') as st:
    exec(st.read())

 #get_ipython().run_line_magic('run', 'struct/structsetup.py')
#get_ipython().run_line_magic('run', 'stab/stabsetup.py')
#get_ipython().run_line_magic('run', 'systems/systemssetup.py')
#get_ipython().run_line_magic('run', 'aero/aerosetup.py')

# In[ ]:


prob = om.Problem(model = MDOHPCA(), reports = True)
#prob.model.add_subsystem('aero', om.Group())


# In[ ]:


#setting inputs to reduce ambiguities where inputs differ in different disciplines
prob.model.set_input_defaults('CL', val = 0.6)
prob.model.set_input_defaults('CD', val = 0.02)

prob.model.set_input_defaults('sweep', val = 15)
prob.model.set_input_defaults('taper', val = 0.4)
prob.model.set_input_defaults('total_CG', val = 10)

prob.model.set_input_defaults('c', val = 4)
prob.model.set_input_defaults('b', val = 60)
prob.model.set_input_defaults('root_x', val = 10)
prob.model.set_input_defaults('m', val = 240000)


# In[ ]:
#constraints
prob.model.add_constraint('fuel_mass', lower = 5000, upper = 100000, linear = False)
prob.model.add_constraint('landing_tension', upper = 90000000)

prob.model.add_constraint('Cmq', upper = 0)
prob.model.add_constraint('CM_alpha', upper = 0)
prob.model.add_constraint('cmu', lower = 0)
prob.model.add_constraint('CLa', lower = 0)
prob.model.add_constraint('Cnb', lower = 0)
prob.model.add_constraint('Cnr', upper = 0)
prob.model.add_constraint('CYb', upper = 0)
prob.model.add_constraint('CTCD', upper = 0.01)
prob.model.add_constraint('Clb', upper = 0)
prob.model.add_constraint('Clp', upper = 0)

prob.model.add_objective('m', scaler = 1e-4)

# In[ ]:
#declare problem options (of driver and optimiser)

#prob.nonlinear_solver = om.NonlinearBlockGS()
#prob.nonlinear_solver = om.NewtonSolver(solve_subsystems=False)
prob.model.nonlinear_solver = om.NonlinearBlockGS()
prob.model.nonlinear_solver.options['iprint'] = 2
prob.model.nonlinear_solver.options['rtol'] = 1e-10
prob.model.nonlinear_solver.options['maxiter'] = 30
prob.model.nonlinear_solver.options['use_aitken'] = True
prob.model.nonlinear_solver.options['restart_from_successful'] = True


prob.driver = om.ScipyOptimizeDriver()
prob.driver.options['optimizer'] = 'SLSQP'
prob.driver.options['tol'] = 1e-5
prob.driver.options['disp'] = True

#prob.driver.options["debug_print"] = ["nl_cons", "objs", "desvars"]
#setup problem to be run
prob.setup()
prob.final_setup()


#write XDSM model
#from omxdsm import write_xdsm
#write_xdsm(prob, filename='MDOHPCA_trial', out_format='tex', show_browser=True,
#           quiet=False, output_side='left',include_solver=True)#
# In[ ]:


prob.run_model()


# In[ ]:


#prob.run_driver()



# In[ ]:


#get properties of specific constants in disciplines
#prob.model.sys.options._dict['pp_mass']

prob.model.list_outputs()
prob.model.list_inputs()


# In[ ]:







