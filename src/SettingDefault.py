#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 14:14:33 2018

@author: brh
"""
class SettingDefault:
    
    def __init__(self):
        '''
        Default definitions for settings in the model
        '''
        #----------------------------------------------------------------------                      
        # COMMON USER CONTROLS
        # set inoisy to False to suppress command line output 
        # inoisy = false will print a . for every time step (without newline)
        # to suppress the . print, set completely_silent to True
        # NOTE THAT IT IS USUALLY BETTER TO SET THESE IN custom_case
        self.inoisy = True
        self.completely_silent = False
        
        self.working_directory = None
        
        # geometry case (should be overwritten by custom case)
        self.geometry_case = 'flow_over_a_bump'
        self.flow_over_a_bump_NX = 128
        self.mannings_n_use_global_default = True
        self.inflowBC = 0.18
        self.heightBC = 0.33
        
        # baseline time interval
        self.dt = 1.0  
        # total simulation time
        self.time_total = 50000
        # units for simulation time
        self.time_units = 'seconds' 
        # maximum number of time steps in this simulation
        self.steps_max = 2
    
        # global mannings n (unless overwritten in SystemGeometry class)     
        #self.mannings_n_use_global_default = False
        self.mannings_n_global_default = 0.015
        # maximum number of width:depth pairs for general cross-section.
        # set this to 0 to save memory if general cross-section not used
        # set to 1 if you want the widthdepth pairs set by input files
        self.geometry_number_widthdepth_pairs = 0
         
        # T/F for writing binary output file at end of simulation
        # See BINARY OUTPUT CONTROLS
        self.binsave_writedata = True
        self.binsave = True
    
        # T/F for writing restart files of volume and flowrate
        # see RESTART OUTPUT CONTROLS
        self.txtout_writedata = True
        # header isn't required, but must of format '"......"\n'
        self.txtout_writedata_header = 'Restart File'      
    
        # Real-time print controls during simulation
        # time_header is various debug information as simulation progresses
        self.print_time_header_step_interval = 1
        # The debug_iterstart allows debugging prints to be silent until the
        #   time step is reached
        self.print_debug_iterstart = 1
    
        # controls for ploting at command line
        self.plot_time_interval = 0.01
        # iterstart allows plotting to be delayed so as to not waste time
        #   plotting during ramp-up time.
        self.plot_iterstart = 0
        self.plot_types = {}
    
        # ramp up of initial conditions.
        self.inflow_rampup = False
        self.inflow_rampup_time = 50
        self.inflow_rampup_time_units = 'seconds'
        self.flowrate_IC = 1.0# 0.01 # for rampup
    
        #----------------------------------------------------------------------                      
        # RESTART INPUT CONTROLS
    
        #self.IC_type = 'restart_file'
        #self.IC_type = 'elevation_file_interp'
        self.IC_type = 'custom' # standard use
        
        # Using volume and flowrate txt files (writedata) as restart files.
        #self.restart = False
        self.restart_time = 0
        self.restart_time_units = 'seconds'
        self.restart_volumefile   =   'volume_output_#.txt'
        self.restart_flowratefile = 'flowrate_output_#.txt'
        self.restart_directory = 'restart'
    

        # using a old restart file for interpolated initial conditions    
        #self.IC_interp_elevation = False
        self.IC_interp_elevation_time = 0
        self.IC_interp_elevation_folder = 'elevation folder name'
        self.IC_interp_elevation_filename = 'elevation_output_#.txt'
        self.IC_interp_xvalue_folder = 'xvalue folder name'
        self.IC_interp_xvalue_filename = 'xvalue_output_#.txt'
        
        #----------------------------------------------------------------------                      
        # RESTART OUTPUT CONTROLS
        
        # note that this is written every timeinterval after starttime
        self.txtout_writedata_starttime = 0
        self.txtout_volume_filename = 'volume_output.txt'
        self.txtout_flowrate_filename = 'flowrate_output.txt'
        self.txtout_flowrate_face_filename = 'flowrate_face_output.txt'
        self.txtout_elevation_filename = 'elevation_output.txt'
        self.txtout_elevation_faceM_filename = 'elevation_faceM_output.txt'
        self.txtout_elevation_faceP_filename = 'elevation_faceP_output.txt'
        self.txtout_hyddepth_filename = 'hyddepth_output.txt'
        self.txtout_area_filename = 'area_output.txt'
        self.txtout_timestep_filename = 'timestep_output.txt'
        
        self.txtout_writedata_timeinterval = 20.0
        self.txtout_writedata_time_units = 'seconds'
        self.txtout_directory = 'output_text'
        self.txtout_items_on_line_maximum = 100
        
        # outputs that generally only is provided at the start
        self.txtout_zbottom_filename = 'zbottom_output.txt'
        self.txtout_zbottom_face_filename = 'zbottom_face_output.txt'
        self.txtout_length_filename = 'length_output.txt'
        self.txtout_xvalue_filename = 'xvalue_output.txt'
        self.txtout_xvalue_face_filename = 'xvalue_face_output.txt'
        self.txtout_manningsn_filename = 'manningsn_output.txt'
        self.txtout_breadth_filename = 'breadth_output.txt'
        
        # T/F for control is in the COMMON USER CONTROLS
        #self.output_restart_writedata = True
        #----------------------------------------------------------------------                      
        # SIMULATION SETTING OUTPUT
        
        self.setting_filename = 'control_setting.txt'
        
        #----------------------------------------------------------------------                      
        # COMMON ERROR CHECKING AND DEBUGGING
        
        # error checking for volume at some interval of steps (not time)
        self.check_total_volume_stepinterval = 1
    
        #----------------------------------------------------------------------                      
        # BINARY OUTPUT CONTROLS       
        
        # Data saved to python binary file at end of simulation.
        # Note that all data are stored in a single file
        #  except for the time array, which gets a separate file.
        #HACK - need to develope an output functionality that writes at 
        #  intervals the present approach stores, but doesn't write until  the
        #  end, so all data is lost on an error.
        self.binsave_data_filename   = 'output_bindata.npy'
        self.binsave_time_filename   = 'output_bintime.npy'
       
        self.binsave_element = \
            ['volume','flowrate','eta','froude']
        self.binsave_face    =  \
            []
        self.binsave_timeinterval = 1000 
        self.binsave_time_units = 'seconds'
        self.binsave_directory = 'output_binary'
        # T/F for control is in the COMMON USER CONTROLS 
        #self.binsave_data = True
    
        
        #----------------------------------------------------------------------                      
        # AD HOC ERROR CHECKING AND LIMITERs
        # see also the time controls section
        
        # element for checking during debug
        self.check_element =  None#[701, 702, 703, 704, 705, 706]# [124, 125, 126, 127]# [78, 79, 80, 81] #[43, 44, 45, 46, 47]# None#[880, 881, 882,883] #None# [310, 311, 312] # [242, 243, 244, 245] #[231,232,233,234]# [315,316]# [242, 243, 244, 245]#[240, 241, 242] # # [269, 270, 271, 272] # [242, 243, 244, 245]# [266, 267, 268, 269] #None# [90, 91, 92, 93, 94]# [265, 266, 267, 268]# [242, 243, 244, 245]# [82, 83,84,85, 86]# [242, 243, 244, 245] #
        #[91,92, 93, 94, 95, 96, 97] # None #[265, 266, 267, 268] #[242, 243, 244, 245] #[91,92,93,94]# None #[115, 116, 117, 118] # [91,92,93]
                
        # Small value used as switch to small-depth momentum model
        self.depth_small_value = 0.01 #0.2
        # mannin's n used in small-depth model when 0 is provided
        self.depth_small_value_default_n = 0.01
        
        # Large depth for checking geometry
        self.depth_maximum_expected = 5.0
        
        self.face_area_min = 0.01
        # zero value is used in place of zero or negative values
        self.volume_zero_value = 1e-6
        self.area_zero_value = 1e-7
        self.depth_zero_value = 1e-4
        self.topwidth_zero_value = 1e-4
        self.flowrate_zero_value = 1e-4
        self.velocity_zero_value = 1e-4
        
        # limit the momentum solution for small volumes
        self.velocity_limit_small_volume = 0.2
        self.froude_limit_small_volume = 0.5        
        self.velocity_limit = 10# 3.0
               
        # maximum area used in width-depth geometry. Needs to be larger than
        # any expected area.
        self.area_maximum = 2.0
        
        # angle (in degrees) that is considered small for purposes of
        # preventing problems with vertical or flat surfaces. 
        # Do not use 0.0, as there are problems with quadratic formulas
        # etc. when perfectly vertical or flat channel edges encountered.
        self.angle_minimum = 0.1
        
        # time scale limits for tau-based interpolation
        self.time_scale_maximum = 1e6
        self.time_scale_minimum = 1e-6
    
        #----------------------------------------------------------------------                      
        # GEOMETRY ADJUSTMENT
        # allows automated fixes for small inconsistencies in width depth
        
        # top-level on/off for geometry consistency adjustment
        self.geometry_adjust = True
        # maximum allowable fractional change in a cell width or depth value 
        self.geometry_adjust_fraction_max = 0.05
        # fractional change imposed on that is inconsistent with pair below
        # and upper pair is also inconsistent
        self.geometry_adjust_fraction = 0.01
        # small value relative to channel width
        self.geometry_small_width = 1e-8
        # add a downstream buffer cell 
        # normally set this to false unless customizing for geometry that
        # has severe restrictions near outflow BC
        self.geometry_add_downstream_buffer = False
        
        self.geometry_downstream_minimum_length = 0.0
      
        #----------------------------------------------------------------------                      
        # TIME CONTROLS (rarely changed)
        
        # minimum time step
        self.dt_min = 1e-7 
        # maximum allowed CFL. dt is decreased above this
        self.cfl_max = 0.7
        # CFL below which dt can be increased.
        self.cfl_increase_dt = 0.5  
        # time step interval between dt changes:
        self.cfl_increase_stepinterval_from_decrease = 10  
    
        #----------------------------------------------------------------------                      
        # CONSTANTS
        self.gravity = 9.81
               
        #----------------------------------------------------------------------                      
        # STORAGE
        # these are variables stored in the setting for convenience
    
        # counter for increase_stepinterval in TIME CONTROLS 
        self.cfl_increase_stepcounter = 0    
        # Initial conditions storage during ramp-up time.        
        self.inflowBCsaved = 0 
        # storage for volume checking
        self.total_volume_element = 0;
        self.total_volume_face = 0;
        
        #----------------------------------------------------------------------                      
        # MISC THAT ARE NOT GENERALLY USER CONTROLLED
        # these are set within the code
        self.steady_solution_exists = False
        self.steady_solution_type = 'eta'
        
        # swashes bin is address of the binary file for SWASHES solution. 
        # Note that on OSX or Linux that you may need to do a chmod +x swashes 
        # file to make the file executable
        self.swashes_bin                                                      \
            = '/Users/brh/Box Sync/AA_Sync/CODE/SWASHES-1.03.00/bin/swashes'
        self.swashes_stype = 1 # Bump
        self.swashes_domain = 1 # L=25
        self.swashes_choice = 3 # transcritical with shock

        
        #----------------------------------------------------------------------                              
        # METHODS
        #self.method_area_interpolation = 'equalweight' 
        #self.method_area_interpolation = 'upwind'
        #self.method_area_interpolation = 'linear'
        self.method_area_interpolation = 'timescale'
        #self.method_area_interpolation = 'froude'
        
        # NOTE: linear is strongly recommended for eta.
        #self.method_eta_interpolation = 'equalweight' 
        #self.method_eta_interpolation = 'upwind'
        self.method_eta_interpolation = 'linear'
        #self.method_eta_interpolation = 'timescale'
    
        #self.method_topwidth_interpolation = 'equalweight' 
        #self.method_topwidth_interpolation = 'upwind'
        #self.method_topwidth_interpolation = 'linear'
        self.method_topwidth_interpolation = 'timescale'
     
    
        #self.method_perimeter_interpolation = 'equalweight' 
        #self.method_perimeter_interpolation = 'upwind'
        #self.method_perimeter_interpolation = 'linear'
        self.method_perimeter_interpolation = 'timescale'
    
       
        #self.method_flowrate_interpolation = 'upwind'
        #self.method_flowrate_interpolation = 'equalweight'
        #self.method_flowrate_interpolation = 'linear'
        self.method_flowrate_interpolation = 'timescale'
        
            
        #self.method_momentum = 'uniform'
        self.method_momentum = 'T00'
        #self.method_momentum = 'T10'
        #self.method_momentum = 'T20'
         
        self.method_rungekutta = 'rk4_classic'
        #self.method_rungekutta = 'rk4_3/8'
        #self.method_rungekutta = 'ssp_(3,3)'
        #self.method_rungekutta = 'ssp_(4,3)'
        #self.method_rungekutta = 'ssp_(5,3)'
        #self.method_rungekutta = 'ssp_(5,4)'
        #self.method_rungekutta = 'ssp_(6,3)'
        #self.method_rungekutta = 'ssp_(6,4)'
                
        # adjusts element flowrate and velocity when both faces are
        # of opposite sign (indicateor of incipient instability)
        self.method_Q_adjust_inconsistent = False
        
        # adjusts to damp V-shaped flowrates across element
        self.method_Q_adjust_Vshape = True
        self.method_Q_adjust_Vshape_coef = 0.1
        

        # Using Froude number for hydraulic jump faces
        self.method_froude_interpolation = True

        # interpolation/extrapolation accros a hydraulic jump
        self.method_hydjump_face = 'simple'
        #self.method_hydjump_face = 'momentum_match'
        #self.method_hydjump_face = 'extrapolate_surface'
        #self.method_hydjump_face = 'energy_limit'
        #self.method_hydjump_face = 'rectangular'
        
        
        self.method_hydjump_Froude_epsilon = 0.1
        self.method_hydjump_Momentum_epsilon = 0.05
              
        self.method_use_cosangle = True
                                
        return
        
#==============================================================================
        
    def custom_overwrite(self,setting, custom):
        
        import sys
        
        for thisitem, thisvalue in custom.items():
            setattr(setting,thisitem,thisvalue)
            if setting.inoisy:
                print('CUSTOM: setting.',thisitem,' = ',thisvalue)
            #endif    
        #endfor
        
        return setting
#==============================================================================

#EOF    