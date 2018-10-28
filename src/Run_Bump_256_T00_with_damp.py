#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 11:48:49 2018

Simulations for flow over a bump
Using 64 grid cells (32 in domain)
Modifying CFL limits to look at convergence with decreasing dt

Using rectangular hydraulic jump

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
    
    #--------------------------------------------------------------------------
    # new setup
    

    custom = {'user_case_name'     : 'flow_over_a_bump_0256_damp01', 
              'txtout_directory'   : 'outtxt_bump_0256_damp01',
              'binsave_directory'  : 'outbin_bump_0256_damp01',
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
              'method_Q_damp_oscillation_coef': 0.1,
              'method_Q_damp_oscillation_type': 'linear_simple', 
              'method_Q_adjust_Vshape' : True,
              'method_Q_adjust_Vshape_coef' :1.0,
              'method_use_cosangle' : True,
              'plot_element_limits': [10,80],            
              'plot_types'         : ['free_surface','flowrate'],
              'swashes_bin'        : '/Users/brh/Box Sync/AA_Sync/CODE/SWASHES-1.03.00/bin/swashes'}

    #--------------------------------------------------------------------------
    # run the simulation
    [trun, setting] = sr.SimulationRun().simulation_toplevel(custom)
    
    del trun
    del setting

    
#    #--------------------------------------------------------------------------
#    # new setup
#    custom = {'user_case_name'     : 'flow_over_a_bump_0032bR', 
#              'txtout_directory'   : 'outtxt_bump_0032bR',
#              'binsave_directory'  : 'outbin_bump_0032bR',
#              'flow_over_a_bump_NX': 32,
#              'heightBC': 0.262,
#              'time_total': 400,
#              'plot_time_interval' : 4,
#              'print_time_header_step_interval':100,
#              'cfl_max':  0.35,
#              'cfl_increase_dt': 0.3,
#              'mannings_n_global_default': 0.0,
#              'mannings_n_for_buffer': 0.03,
#              'method_momentum' :  'T00',
#              'method_hydjump_face': 'momentum_match',
#              'method_Q_adjust_Vshape' : True,
#              'method_Q_adjust_Vshape_coef' :1.0,
#              'method_use_cosangle' : True,
#              'plot_element_limits': [12,20],            
#              'plot_types'         : ['free_surface','flowrate'],
#              'swashes_bin'        : '/Users/brh/Box Sync/AA_Sync/CODE/SWASHES-1.03.00/bin/swashes'}
#
#    #--------------------------------------------------------------------------
#    # run the simulation
#    [trun, setting] = sr.SimulationRun().simulation_toplevel(custom)
#    
#    del trun
#    del setting
#
#    #--------------------------------------------------------------------------
#    # new setup
#
#    custom = {'user_case_name'     : 'flow_over_a_bump_0032cR', 
#              'txtout_directory'   : 'outtxt_bump_0032cR',
#              'binsave_directory'  : 'outbin_bump_0032cR',
#              'flow_over_a_bump_NX': 32,
#              'heightBC': 0.262,
#              'time_total': 400,
#              'plot_time_interval' : 4,
#              'print_time_header_step_interval':500,
#              'cfl_max':  0.1,
#              'cfl_increase_dt': 0.085,
#              'mannings_n_global_default': 0.0,
#              'mannings_n_for_buffer': 0.03,
#              'method_momentum' :  'T10',
#              'method_hydjump_face': 'momentum_match',
#              'method_Q_adjust_Vshape' : True,
#              'method_Q_adjust_Vshape_coef' :1.0,
#              'method_use_cosangle' : True,
#              'plot_element_limits': [12,20],            
#              'plot_types'         : ['free_surface','flowrate'],
#              'swashes_bin'        : '/Users/brh/Box Sync/AA_Sync/CODE/SWASHES-1.03.00/bin/swashes'}
#
#    #--------------------------------------------------------------------------
#    # run the simulation
#    [trun, setting] = sr.SimulationRun().simulation_toplevel(custom)
#    
#    del trun
#    del setting
#
#    #--------------------------------------------------------------------------
#    # new setup
#
#    custom = {'user_case_name'     : 'flow_over_a_bump_0032dR', 
#              'txtout_directory'   : 'outtxt_bump_0032dR',
#              'binsave_directory'  : 'outbin_bump_0032dR',
#              'flow_over_a_bump_NX': 32,
#              'heightBC': 0.262,
#              'time_total': 400,
#              'plot_time_interval' : 4,
#              'print_time_header_step_interval':1000,
#              'cfl_max':  0.035,
#              'cfl_increase_dt': 0.03,
#              'mannings_n_global_default': 0.0,
#              'mannings_n_for_buffer': 0.03,
#              'method_momentum' :  'T20',
#              'method_hydjump_face': 'momentum_match',
#              'method_Q_adjust_Vshape' : True,
#              'method_Q_adjust_Vshape_coef' :1.0,
#              'method_use_cosangle' : True,
#              'plot_element_limits': [12,20],            
#              'plot_types'         : ['free_surface','flowrate'],
#              'swashes_bin'        : '/Users/brh/Box Sync/AA_Sync/CODE/SWASHES-1.03.00/bin/swashes'}
#
#    #--------------------------------------------------------------------------
#    # run the simulation
#    [trun, setting] = sr.SimulationRun().simulation_toplevel(custom)
#    
#    del trun
#    del setting


if __name__ == '__main__':
    import sys    
    main(sys.argv) 
    
#==============================================================================
#EOF
