#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 14:11:07 2018

@author: brh
"""

def main(argv):
    
    import SimulationRun as sr
    import numpy as np
    
    
    #temp2 = temp['head']
    
    #
    #print(temp2[0,:])
    #print(s.np_array('x'))
    
    #print(s.dom_params)
    #print(s.dataframe())
    #solutionhead = s.np_array('x')
    #solutionX = s.np_array('x')
    
    #print(solutionhead)
    
    #sys.exit()

    #def custom(self):
    #    pass
    
    #--------------------------------------------------------------------------                       
    # RUN-TIME OUTPUT CONTROLS
    # useful for compact printout during debugging
    np.set_printoptions(precision=10)


    #--------------------------------------------------------------------------
    # SETUP  
    # Provide a dictionary that replaces the methods of setting.###
    #   with the values in the dictionary. See SettingDefault.py for values
    # At a minimum, the 'user_case_name' must be provided that
    #   matches one of the predefined cases. 
    #   See SimulationRun.predefined_case_setup()               
    
    # Waller Creek baseline - starting from time 0
#    custom = {'user_case_name'     : 'Waller_Creek_baseline', 
#              'plot_element_limits': [190,250],
#              'plot_types'         : ['free_surface','flowrate']}
     
    custom = {'user_case_name'     : 'flow_over_a_bump_128', 
              'flow_over_a_bump_NX': 256,
              'steps_max': 30000,
              'heightBC': 0.262,
              'time_total': 300,
              'plot_time_interval' : 1,
              'print_time_header_step_interval':10,
              'cfl_max':  0.7,
              'cfl_increase_dt': 0.6,
              'mannings_n_global_default': 0.0,
              'mannings_n_for_buffer': 0.03,
              'method_momentum' :  'T00',
              'method_hydjump_face': 'momentum_match',
              'method_Q_adjust_inconsistent': False,
              'method_Q_damp_oscillation_coef': 0.01,
              'method_Q_damp_oscillation_type': None, 
              'method_Q_adjust_Vshape' : True,
              'method_Q_adjust_Vshape_coef' :1.0,
              'method_use_cosangle' : True,
              'plot_element_limits': [20, 75],
              'plot_types'         : ['free_surface','flowrate'],
              'swashes_bin'        : '/Users/brh/Box Sync/AA_Sync/CODE/SWASHES-1.03.00/bin/swashes'}
    
#              'IC_type'              :'restart_file',        
#              'restart_time'         : 145,
#              'restart_time_units'   :'seconds',
#              'restart_volumefile'   :   'volume_output_20180423_0610.txt',
#              'restart_flowratefile' :   'flowrate_output_20180423_0610.txt',
#              'restart_directory'    : 'restart',

              #'flow_over_a_bump_NX': 128,
              #'plot_element_limits': [35,65],
              #'flow_over_a_bump_NX': 256,
              #'plot_element_limits': [70,130],
    

    #--------------------------------------------------------------------------
    # run the simulation
    [trun, setting] = sr.SimulationRun().simulation_toplevel(custom)
    

if __name__ == '__main__':
    import sys    
    main(sys.argv) 
    
#==============================================================================
#EOF
    