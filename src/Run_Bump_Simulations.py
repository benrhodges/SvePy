#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 12:29:51 2018

@author: brh
"""

def main(argv):
    
    import SimulationRun as sr
    import numpy as np
    
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
    custom = {'user_case_name'     : 'flow_over_a_bump_0032_faceoutput', 
              'txtout_directory'   : 'outtxt_bump_0032_faceoutput',
              'binsave_directory'  : 'outbin_bump_0032_faceoutput',
              'flow_over_a_bump_NX': 32,
              'heightBC': 0.262,
              'time_total': 400,
              'plot_time_interval' : 2,
              'print_time_header_step_interval':50,
              'cfl_max':  0.7,
              'cfl_increase_dt': 0.6,
              'mannings_n_global_default': 0.0,
              'mannings_n_for_buffer': 0.03,
              'method_momentum' :  'T00',
              'method_hydjump_face': 'momentum_match',
              'method_Q_adjust_Vshape' : True,
              'method_Q_adjust_Vshape_coef' :1.0,
              'plot_element_limits': [0,16],
              'plot_types'         : ['free_surface','flowrate'],
              'swashes_bin'        : '/Users/brh/Box Sync/AA_Sync/CODE/SWASHES-1.03.00/bin/swashes'}

    #--------------------------------------------------------------------------
    # run the simulation
    [trun, setting] = sr.SimulationRun().simulation_toplevel(custom)
    
    del trun
    del setting
    
     
    custom = {'user_case_name'     : 'flow_over_a_bump_0064_faceoutput', 
              'txtout_directory'   : 'outtxt_bump_0064_faceoutput',
              'binsave_directory'  : 'outbin_bump_0064_faceoutput',
              'flow_over_a_bump_NX': 64,
              'heightBC': 0.262,
              'time_total': 400,
              'plot_time_interval' : 4,
              'print_time_header_step_interval':100,
              'cfl_max':  0.7,
              'cfl_increase_dt': 0.6,
              'mannings_n_global_default': 0.0,
              'mannings_n_for_buffer': 0.03,
              'method_momentum' :  'T00',
              'method_hydjump_face': 'momentum_match',
              'method_Q_adjust_Vshape' : True,
              'method_Q_adjust_Vshape_coef' :1.0,
              'plot_element_limits': [0,32],
               'plot_types'         : ['free_surface','flowrate'],
              'swashes_bin'        : '/Users/brh/Box Sync/AA_Sync/CODE/SWASHES-1.03.00/bin/swashes'}

    #--------------------------------------------------------------------------
    # run the simulation
    [trun, setting] = sr.SimulationRun().simulation_toplevel(custom)
    
    del trun
    del setting
    
    #--------------------------------------------------------------------------
    # new setup
    custom = {'user_case_name'     : 'flow_over_a_bump_0128_faceoutput', 
              'txtout_directory'   : 'outtxt_bump_0128_faceoutput',
              'binsave_directory'  : 'outbin_bump_0128_faceoutput',
              'flow_over_a_bump_NX': 128,
              'heightBC': 0.262,
              'time_total': 4300,
              'plot_time_interval' : 4,
              'print_time_header_step_interval':100,
              'cfl_max':  0.7,
              'cfl_increase_dt': 0.6,
              'mannings_n_global_default': 0.0,
              'mannings_n_for_buffer': 0.03,
              'method_momentum' :  'T00',
              'method_hydjump_face': 'momentum_match',
              'method_Q_adjust_Vshape' : True,
              'method_Q_adjust_Vshape_coef' :1.0,
              'plot_element_limits': [0,64],
              'plot_types'         : ['free_surface','flowrate'],
              'swashes_bin'        : '/Users/brh/Box Sync/AA_Sync/CODE/SWASHES-1.03.00/bin/swashes'}
   
    [trun, setting] = sr.SimulationRun().simulation_toplevel(custom)
    
    del trun
    del setting

    #--------------------------------------------------------------------------
    # new setup
    custom = {'user_case_name'     : 'flow_over_a_bump_0256_faceoutput', 
              'txtout_directory'   : 'outtxt_bump_0256_faceoutput',
              'binsave_directory'  : 'outbin_bump_0256_faceoutput',
              'flow_over_a_bump_NX': 256,
              'heightBC': 0.262,
              'time_total': 400,
              'plot_time_interval' : 4,
              'print_time_header_step_interval':100,
              'cfl_max':  0.7,
              'cfl_increase_dt': 0.6,
              'mannings_n_global_default': 0.0,
              'mannings_n_for_buffer': 0.03,
              'method_momentum' :  'T00',
              'method_hydjump_face': 'momentum_match',
              'method_Q_adjust_Vshape' : True,
              'method_Q_adjust_Vshape_coef' :1.0,
              'plot_element_limits': [0,128],
              'plot_types'         : ['free_surface','flowrate'],
              'swashes_bin'        : '/Users/brh/Box Sync/AA_Sync/CODE/SWASHES-1.03.00/bin/swashes'}
   
    [trun, setting] = sr.SimulationRun().simulation_toplevel(custom)
    
    del trun
    del setting

    #--------------------------------------------------------------------------
    # new setup
    custom = {'user_case_name'     : 'flow_over_a_bump_0512_faceoutput', 
              'txtout_directory'   : 'outtxt_bump_0512_faceoutput',
              'binsave_directory'  : 'outbin_bump_0512_faceoutput',
              'flow_over_a_bump_NX': 512,
              'heightBC': 0.262,
              'time_total': 400,
              'plot_time_interval' : 4,
              'print_time_header_step_interval':100,
              'cfl_max':  0.7,
              'cfl_increase_dt': 0.6,
              'mannings_n_global_default': 0.0,
              'mannings_n_for_buffer': 0.03,
              'method_momentum' :  'T00',
              'method_hydjump_face': 'momentum_match',
              'method_Q_adjust_Vshape' : True,
              'method_Q_adjust_Vshape_coef' :1.0,
              'plot_element_limits': [0,256],
              'plot_types'         : ['free_surface','flowrate'],
              'swashes_bin'        : '/Users/brh/Box Sync/AA_Sync/CODE/SWASHES-1.03.00/bin/swashes'}
   
    [trun, setting] = sr.SimulationRun().simulation_toplevel(custom)
    
    del trun
    del setting
    
#    #--------------------------------------------------------------------------
#    # new setup
#    custom = {'user_case_name'     : 'flow_over_a_bump_1024', 
#              'txtout_directory'   : 'outtxt_bump_1024',
#              'binsave_directory'  : 'outbin_bump_1024',
#              'flow_over_a_bump_NX': 1024,
#              'heightBC': 0.262,
#              'time_total': 300,
#              'plot_time_interval' : 4,
#              'print_time_header_step_interval':100,
#              'cfl_max':  0.7,
#              'cfl_increase_dt': 0.6,
#              'mannings_n_global_default': 0.0,
#              'mannings_n_for_buffer': 0.03,
#              'method_momentum' :  'nonuniform',
#              'method_hydjump_face': 'simple',
#              'method_Q_damp_oscillation_coef': 0.1,
#              'method_Q_damp_oscillation_type': 'linear_simple', 
#              'method_Q_adjust_Vshape' : True,
#              'method_Q_adjust_Vshape_coef' :1.0,
#              'plot_element_limits': [0,512],
#               'plot_types'         : ['free_surface','flowrate'],
#              'swashes_bin'        : '/Users/brh/Box Sync/AA_Sync/CODE/SWASHES-1.03.00/bin/swashes'}
#   
#    [trun, setting] = sr.SimulationRun().simulation_toplevel(custom)
#    
#    del trun
#    del setting
#
    

if __name__ == '__main__':
    import sys    
    main(sys.argv) 
    
#==============================================================================
#EOF
