#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 09:22:10 2018

@author: brh
"""

class InitialConditions:
#==============================================================================
    def custom_initial_condition(self, el, geo, BC, setting, NX):
        
        import sys
        import ElementGeometry
        import BoundaryConditions2
        
        bc = BoundaryConditions2.BoundaryConditions()  
        eg = ElementGeometry.ElementGeometry()
        
        if setting.geometry_case == 'flow_over_a_bump':
            # set elevation to the output height
            el['eta'][:] = bc.height(setting)
            # get the volume for the elevation
            el = eg.get_volume_from_elevation(el, geo, setting, NX)
            # set the flowrate to the BC value
            # note do not use bc.flow, instead use BC.inflow which has
            #   been adjusted for rampup time
            el['flowrate'][:] = BC.inflowrate
            
        elif setting.geometry_case == 'Waller_Creek':            
            depth = setting.Waller_Creek_initial_depth
            el['eta'][:] = depth + geo['zbottom'][:]
            el['eta'][NX-1] = setting.heightBC
            
            # ensure starting condition is monotonic free surface
            for ii in range(NX-1,1, -1):
                #print(ii)
                if el['eta'][ii-1] < el['eta'][ii]:
                    el['eta'][ii-1] = el['eta'][ii]
                #endif
            #endfor ii    
            
            # get the volume for these elevation
            el = eg.get_volume_from_elevation(el, geo, setting, NX)
            
            # set flowrate to BC value (including ramp-up, if used)
            el['flowrate'][:] = BC.inflowrate
            
        else:
            print('error: unknown value for setting.geometry_case of ', \
                  setting.geometry_case)
            sys.exit()
        #endif setting.geometry_case   
        
        return el
    
#==============================================================================
    def initial_flow_and_volume(self, el, geo, BC, setting, NX):
        '''
        Set up initial conditions for flow and volume on elements
        '''
        import sys
        
        import BoundaryConditions2
        import ElementGeometry

        bc = BoundaryConditions2.BoundaryConditions()  
        eg = ElementGeometry.ElementGeometry()
                
        if setting.IC_type == 'custom':
            el = self.custom_initial_condition(el, geo, BC, setting, NX)     
            
        elif setting.IC_type == 'elevation_file_interp':
            # get locations of the elevation and xvalue files
            [el, setting] = self.interp_elevation_from_file                   \
                            (el, geo, setting, NX)
            el = eg.get_volume_from_elevation(el, geo, setting, NX)
              
            # Flow initial conditions consistent with BC 
            el['flowrate'][:] = BC.inflow
            
        else:
            print('error, unknown value for setting.IC_type of ',              \
                  setting.IC_type)
            sys.exit()
        #endif
        
        return el
    
#==============================================================================
    def interp_elevation_from_file(self, el, geo, setting, NX):    

        import sys
        import numpy as np
        
        # read the elevation file and xvalues
        [elevation, xvalue, setting] = self.elevation_interp_file_read(setting) 

        # check that the xvalues and elevation are consistent      
        if xvalue.size != elevation.size:
            print('error, mismatching size in X and Eta interpolation input')
            sys.exit()
        #endif            
            
        # cycle through to interpolate elevation values    
        for ii in range(0,NX):
            thisx = geo['xvalue'][ii]
            xdiff = abs(xvalue[:] - thisx)
            xx = np.argmin(xdiff)
            xlo = xvalue[xx-1]
            elo = elevation[xx-1]
            xmd = xvalue[xx]
            emd = elevation[xx]   
            if ii < NX-1:
                xhi = xvalue[xx+1]
                ehi = elevation[xx+1]
            else:
                xhi = xvalue[xx]
                ehi = elevation[xx]
            #endif    
            if (thisx >= xlo) & (thisx <= xmd):
                thisEta = elo + (emd - elo)  \
                    * (thisx - xlo) / (xmd - xlo)
            elif (thisx > xmd) & (thisx <= xhi):        
                thisEta = emd + (ehi - emd)  \
                    * (thisx - xmd) / (xhi - xmd)
            elif thisx > xhi:
                # retain old eta
                pass 
            else:
                print('error in logic')                        
                sys.exit()
            #endif
            el['eta'][ii] = thisEta
        #endfor

        return [el, setting]

#==============================================================================
    def elevation_interp_file_read(self, setting):
        
        import os
        import FileReadWrite2       
        frw = FileReadWrite2.FileReadWrite()
    
        # get locations of the elevation and xvalue files
        setting.IC_interp_elevation_folder = frw.absolute_path                \
            (setting.IC_interp_elevation_folder, setting.working_directory)  

        setting.IC_interp_xvalue_folder = frw.absolute_path                   \
            (setting.IC_interp_xvalue_folder, setting.working_directory)  

        setting.IC_interp_elevation_filename = os.path.join                   \
              (setting.IC_interp_elevation_folder,                            \
               setting.IC_interp_elevation_filename) 

        setting.IC_interp_xvalue_filename = os.path.join                      \
              (setting.IC_interp_xvalue_folder,                               \
               setting.IC_interp_xvalue_filename) 

        # open the files    
        file_elevation = open(setting.IC_interp_elevation_filename,'r')
        file_xvalue    = open(setting.IC_interp_xvalue_filename,'r')
        
        # read the line that is equal to or after the restart time
        [temp_time,elevation] = frw.get_oneline_data_matching_time \
            (setting.IC_interp_elevation_time, file_elevation, False)

        [temp_time,xvalue] = frw.get_oneline_data_matching_time \
            (0, file_xvalue, False)
         
        # close the files    
        file_elevation.close()
        file_xvalue.close()
        
        return [elevation, xvalue, setting]
    
#==============================================================================
#EOF
