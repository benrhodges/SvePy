#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 08:13:32 2018

@author: brh
"""

class CaseFlowOverBump:
    
    def setting_CaseFlowOverBump(self, setting):
        
        
        setting.IC_type = 'custom'
        
        setting.cfl_max =  0.7# 1.7# 1.7
        
        # CFL below which dt can be increased.
        setting.cfl_increase_dt = 0.6 #1.5#  1.5 

        # FLOW OVER A BUMP
        # geometry case
        setting.geometry_case = 'flow_over_a_bump'
        setting.flow_over_a_bump_NX = 128
        setting.mannings_n_use_global_default = True
        setting.mannings_n_global_default = 0.0
        setting.inflowBC = 0.18
        setting.heightBC = 0.33
        
        setting.IC_type = 'custom'
        
        setting.depth_small_value = 0.001 
               
        setting.method_Q_damp_oscillation_coef = 0.1 
        setting.method_Q_damp_oscillation_type =  None 
        
        setting.method_Q_adjust_Vshape = True
        setting.method_Q_adjust_Vshape_coef = 1.0

        setting.dt = 0.01 
        # total simulation time
        setting.time_total = 300
        # units for simulation time
        setting.time_units = 'seconds' 
        # maximum number of time steps in this simulation
        setting.steps_max = 50000
    
        setting.print_debug_iterstart = 1

        # Real-time print controls during simulation
        # time_header is various debug information as simulation progresses
        setting.print_time_header_step_interval = 1# 50

        # controls for ploting at command line
        setting.plot_time_interval =  0.01# 1.0
        # iterstart allows plotting to be delayed so as to not waste time
        #   plotting during ramp-up time.
        setting.plot_iterstart = 0

        # ramp up of initial conditions.
        setting.inflow_rampup = True
        setting.inflow_rampup_time = 100
        setting.inflow_rampup_time_units = 'seconds'
        setting.flowrate_IC = 0.01 # for rampup         

        setting.txtout_writedata_timeinterval = 2.0

        setting.binsave_timeinterval = 5
 
        setting.method_hydjump_face = 'momentum_match'
        
        setting.method_eta_interpolation =  'linear'

        setting.method_area_interpolation = 'timescale' 
        setting.method_flowrate_interpolation = 'timescale'
        setting.method_perimeter_interpolation = 'timescale'
        setting.method_topwidth_interpolation = 'timescale'
        
        setting.method_rungekutta = 'rk4_classic'#'ssp_(6,4)'#'rk4_3/8'# 'rk4_3/8'
       
        return setting

    
    def get_number_of_cells(self,setting):
        '''
        define the number of cells in a reach
        '''
        NX = setting.flow_over_a_bump_NX

        return NX

    def get_number_of_widthdepth_pairs(self,setting):
        '''
        define the number of widthdepth pairs in a reach
        '''
        # HACK need to match to case geometry for widthdepth channel definition
        npair = 0
     
        return npair
   
    def define_geometry(self, geo, setting, NX):    
        '''
        flow over a bump as found in paper by Catella, Paris, and Solari (2008)
        '''
        import sys
        import numpy as np
               
        # note that only the rectangular channel has a defined solution
        # in the code
        use_rectangular_channel = True
        
        use_uniform_lengths = True
        #use_uniform_lengths = False

        #length of domain (m)
        setting.geometry_total_length = 50# 25 #50
        
        # rectangular channel
        if use_rectangular_channel == True:
            geo['etype'][:] = 'rectangular_channel'
            setting.geometry_channel_type = 'rectangular'
            geo['breadth'][:] = 1.0;   
        else:
            print('error')
            sys.exit()
            # test of a parabolic channel where 
            # z(y) = Ay^2  or y(z) = sqrt( z / A)
            # A = ??
            #geo['etype'][:] = 'parabolic_channel'
            #geo['parabolic_value'][:] = 0.587 # 0.587provides same area as rectangle at 0.33 depth
            
            # test of trapezoidal channel            
            geo['etype'][:] = 'trapezoidal_channel'
            setting.geometry_channel_type = 'trapezoidal'
            geo['breadth'][:] = 0.835 # bottom breadth
            geo['trapezoid_angle'][:] = 60.0  #angle, in degrees              
            geo['trapezoid_angle'][:] = np.deg2rad(geo['trapezoid_angle'][:]) # convert to radians
            
            # test of general channel
            # make sure that setting.geometry_number_widthdepth_pairs
            # is equal to the maximum size defined.
#            geo['etype'][:] = 'widthdepth_pair'
#            setting.geometry_channel_type = 'widthdepth_pair'
#            for ii in range(0,NX):
#                geo['widthdepth'][ii,0,0:2] = [0.7,  0.0]
#                geo['widthdepth'][ii,1,0:2] = [0.75,  0.1]
#                geo['widthdepth'][ii,2,0:2] = [0.75 , 0.2]
#                geo['widthdepth'][ii,3,0:2] = [1.0 , 0.3]
#                geo['widthdepth'][ii,4,0:2] = [1.1 , 0.4]
#                geo['widthdepth'][ii,5,0:2] = [1.2 , 1.0]
                     
        if use_uniform_lengths == True:
            # uniform element lengths over reach
            geo['length'][:] = setting.geometry_total_length / NX
        else:
            baselength = setting.geometry_total_length / NX
            dx = [+0.2, 0.0, +0.05, -0.15, 0.0, +0.18, 0.0, 0.0, -0.08, -0.20, 0.0 ]
            kk = 0
            for ii in range(0,NX):
                #print(ii,NX,kk, len(dx))
                geo['length'][ii] = baselength * (1.0 + dx[kk])
                kk = kk+1
                if kk >= len(dx):
                    kk = 0
                #endif---------------------------------------------------------    
                tlen = sum(geo['length'][0:ii])
                if tlen > setting.geometry_total_length:
                    if ii == NX-1:
                        tlen = sum(geo['length'][0:ii-1])
                        geo['length'][ii] = setting.geometry_total_length -tlen
                    else:
                        print(ii, NX, tlen)
                        print('problem in setup of nonuniform length')
                        sys.exit()
                    #endif-----------------------------------------------------    
                #endif---------------------------------------------------------    
            #endfor------------------------------------------------------------      
        #endif-----------------------------------------------------------------
        
        
        # x values measured from upstream face as x=0
        geo['xvalue'][0] = geo['length'][0]/2
        for ii in range(1,NX):
            geo['xvalue'][ii] = geo['xvalue'][ii-1]                             \
                + 0.5*(geo['length'][ii-1] + geo['length'][ii])
        
        # z values for Catella et al 2008 flow over bump   
        geo['zbottom'][:] = 0.2 - 0.05 * ((geo['xvalue'][:] - 10.0)**2.0)
    
        aa = geo['xvalue'][:] < 8.0
        geo['zbottom'][aa] = 0
        bb = geo['xvalue'][:] > 12.0
        geo['zbottom'][bb] = 0
        
    
        return [geo, setting]   
#==============================================================================
#EOF
        
        
        
        
        
        
        