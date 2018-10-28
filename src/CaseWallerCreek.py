#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 15:38:02 2018

@author: brh
"""
class CaseWallerCreek:
 
#==============================================================================
    def setting_CaseWallerCreek(self, setting):

        setting.IC_type = 'custom'
        
        # controls for printing and plotting at command line
        setting.print_time_header_step_interval = 1
        setting.plot_time_interval = 0.01

        setting.geometry_add_downstream_buffer = False
        setting.geometry_downstream_minimum_length = 100.0

        setting.check_element = None# [315, 316, 317]
        
        setting.cfl_max = 0.7 #0.6
        # CFL below which dt can be increased.
        setting.cfl_increase_dt =  0.6 #0.5

        setting.geometry_case = 'Waller_Creek'        
        setting.Waller_Creek_filename = 'WLR_WidthDepthList5.txt'
        setting.Waller_Creek_input_subdirectory = 'input'  
        setting.Waller_Creek_initial_depth = 0.5# 0.2 #0.2 #0.5 had problems     
        setting.Waller_Creek_cellsize_target =  None # 10 (use 0 or None for default)
        setting.mannings_n_use_global_default = False
        setting.inflowBC = 1.0 #0.18
        setting.heightBC = 132.0 #0.33
        
        # baseline time interval
        setting.dt = 1.0
        # total simulation time
        setting.time_total = 100000
        # units for simulation time
        setting.time_units = 'seconds' 
        # maximum number of time steps in this simulation
        setting.steps_max = 100000
 
        # viscosity used for damping grid-scale energy spikes
        setting.viscous_damping_KE  =  0.0 # 1.0e-5
        setting.viscous_damping_KE_type = 'CFL' 
        
        setting.viscous_damping_Q  = 0.0
        setting.viscous_damping_Q_type = None #'linear_normalized' 

        setting.method_flowrate_interpolation = 'timescale'
        setting.method_adjust_for_inconsistent_flowrate = True
       
        # Small value used as switch to small-depth momentum model
        setting.depth_small_value = 0.1
        
        # Large depth for checking geometry
        setting.depth_maximum_expected = 5.0
        
        setting.face_area_min = 0.01
        # zero value is used in place of zero or negative values
        setting.volume_zero_value = 1e-6
        setting.area_zero_value = 1e-7
        setting.depth_zero_value = 1e-4
        setting.topwidth_zero_value = 1e-4
        setting.flowrate_zero_value = 1e-4
        setting.velocity_zero_value = 1e-4
        
        # limit the momentum solution for small volumes
        setting.velocity_limit_small_volume = 10.0
        setting.froude_limit_small_volume = 0.5
        
        setting.velocity_limit = 10


        return setting        
   
#==============================================================================
    def get_number_of_cells(self,setting):
        '''
        obtains the number of cells in a reach from a file
        '''        
        import FileReadWrite2        
        frw  = FileReadWrite2.FileReadWrite()
        
        # get the filename for the width-depth pairs
        widthdepth_filename = self.set_input_file(setting)
        
        NX = frw.read_number_of_cells(setting, widthdepth_filename)
                   
        if setting.geometry_add_downstream_buffer:
            NX = NX+1
        
        return NX

#==============================================================================
    def get_number_of_widthdepth_pairs(self,setting):
        '''
        define the number of widthdepth pairs in a reach
        '''
        import FileReadWrite2        
        frw  = FileReadWrite2.FileReadWrite()

        # get the filename for the width-depth pairs
        widthdepth_filename = self.set_input_file(setting)

        npair = frw.read_max_number_pairs(setting, widthdepth_filename)
        
        return npair

#==============================================================================
    def define_geometry(self, geo, setting, NX):    
        '''
        reads the geometry file for Waller Creek and customizes
        '''
        import numpy as np
        import sys

        import matplotlib.pyplot as plt     
        
        import DataStorage2
        data = DataStorage2.DataStorage()
 
        import SystemGeometry2       
        sg = SystemGeometry2.SystemGeometry()  
        
        import FileReadWrite2
        frw  = FileReadWrite2.FileReadWrite()
        
        #----------------------------------------------------------------------                       
        # get the filename for the width-depth pairs
        widthdepth_filename = self.set_input_file(setting)
        
        #----------------------------------------------------------------------                       
        # read the data into a structure
        geo = frw.read_widthdepth_pairs(geo, setting, NX, widthdepth_filename)

        #----------------------------------------------------------------------                       
        # HACK stretching out the last cell (index 316)
        if setting.geometry_downstream_minimum_length > 0:
            oldX = geo['xvalue'][316]
            oldL = geo['length'][316]
            if setting.geometry_downstream_minimum_length > oldL:
                geo['length'][316] = setting.geometry_downstream_minimum_length
                geo['xvalue'][316] = oldX                                     \
                   + 0.5*(setting.geometry_downstream_minimum_length - oldL) 
       
        #----------------------------------------------------------------------                       
        # creating buffer cell at the end of the domain
        if setting.Waller_Creek_cellsize_target != None:
            if setting.Waller_Creek_cellsize_target > 0:
                buffer_length = setting.Waller_Creek_cellsize_target              
            else:
                buffer_length = 1000.0
            #endif
        else:
            buffer_length = 1000.0
        #endif
        
        #----------------------------------------------------------------------                       
        # set the buffer cell
        if setting.geometry_add_downstream_buffer:
            geo = sg.add_buffer_cell (geo, setting, 317, \
                             buffer_length, \
                             buffer_breadth = 11.0,   \
                             buffer_depth = 40.0)
       
        #----------------------------------------------------------------------                       
        # splitting domain into smaller cells
        if setting.Waller_Creek_cellsize_target != None:
            ncell = 0
            NXold = NX
            if setting.Waller_Creek_cellsize_target > 0:
                #temporary creation of zbottom and xvalue on face for 
                #establishing interpolated zbottom throughout the smaller cells
                zbottom = np.zeros(NXold+1, dtype=np.float64)
                xvalue  = np.zeros(NXold+1, dtype=np.float64)
                fa = {'zbottom':zbottom, 'xvalue':xvalue}   
                [fa, rkf] = sg.face_geometry_initialization(fa, None, geo, NX) 
        
                if setting.geometry_number_widthdepth_pairs > 0:
                    # check widthdepth pair geometry for consistency
                    geo = sg.widthdepth_pair_loop (geo, setting)                
                    # compute additional geometry data for widthdepth pairs
                    [geo, setting] = sg.widthdepth_pair_auxiliary(geo, setting)
                
                # cycle through to find the number of cells to add at
                # each cross-section
                for ii in range(0,NXold):
                    if geo['length'][ii] > 1.5                                \
                            * setting.Waller_Creek_cellsize_target:
                        nadd = int(geo['length'][ii]                          \
                                   // setting.Waller_Creek_cellsize_target)
                        ncell = ncell + nadd
                    else:
                        ncell = ncell + 1
                    #endif    
                #endfor

                NX = ncell
                setting.NX = NX
                kk = 0
                # define the newgeometry data space for the new NX
                [newgeo, setting] = data.geometry(NX, setting)
                # cycle through to get cell length of each new cell
                for ii in range(0,NXold):
                    # location of face
                    if ii < NXold:
                        xface = geo['xvalue'][ii] - 0.5 * geo['length'][ii]
                    else:
                        xface = geo['xvalue'][ii-1] + geo['length'][ii-1]
                    #endif
                    # get the number of cells to add and their length
                    if geo['length'][ii] > 1.5                                \
                            * setting.Waller_Creek_cellsize_target:
                        nadd = int(geo['length'][ii]                          \
                                   // setting.Waller_Creek_cellsize_target)
                        dx = geo['length'][ii] / float(nadd)
                    else:
                        nadd = 1
                        dx = geo['length'][ii]
                    #endif
                    for gg in range(0,nadd):
                        newgeo['length'][kk]  = dx
                        newgeo['xvalue'][kk] = xface + 0.5*dx
                        newgeo['manningsn'][kk] = geo['manningsn'][ii]
                        newgeo['breadth'][kk]   = geo['breadth'][ii]
                        newgeo['ID'][kk] = kk+1
                        newgeo['etype'][kk] = geo['etype'][ii]
                        # linear interpolation from faces for zbottom
                        zup = fa['zbottom'][ii]
                        zct = geo['zbottom'][ii]
                        zdn = fa['zbottom'][ii+1]
                        xup = fa['xvalue'][ii]
                        xct = geo['xvalue'][ii]
                        xdn = fa['xvalue'][ii+1]
                        xval = newgeo['xvalue'][kk]
                        
                        newgeo['widthdepth'][kk,:,:]                          \
                            = geo['widthdepth'][ii,:,:]
                        newgeo['auxWD1'][kk,:] = geo['auxWD1'][ii,:]
                        newgeo['auxWD2'][kk,:] = geo['auxWD2'][ii,:]
                        newgeo['npair'][kk] = geo['npair'][ii]
                        
                        if xval <= xct:
                            newgeo['zbottom'][kk] =                           \
                                (  (xval - xup) * zct                         \
                                 + (xct - xval) * zup 
                                ) / (xct - xup) 
                        else:
                            newgeo['zbottom'][kk] =                           \
                                (  (xval - xct) * zdn                         \
                                 + (xdn - xval) * zct 
                                ) / (xdn - xct) 
                            
                        xface = xface + dx
                        kk=kk+1
                #endfor
                
                # reset the geometry to the new data array
                geo = newgeo
            #endif
        #endif
                
        return [geo, setting]

#==============================================================================        
    def set_input_file(self, setting):
        '''
        input file for Waller Creek Case
        Requires setting.working_directory to be previously defined
        '''
        import os
        
        widthdepth_filename = setting.Waller_Creek_filename
        input_directory     = setting.Waller_Creek_input_subdirectory
        
        widthdepth_filename = os.path.join \
            (setting.working_directory, input_directory, widthdepth_filename)
        
        return widthdepth_filename