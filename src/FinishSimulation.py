#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 17:14:17 2018

@author: brh
"""

class FinishSimulation:

    def final_housekeeping(self, setting, trun):
        
        # error checking before stopping                           
        if  trun.present_step >= setting.steps_max and  \
            trun.present_time < setting.time_total:
                
            if setting.inoisy == True:
                print()
                print('Warning, only completed ',trun.present_time, \
                      ' seconds of simulation')
                print('that as planned for ',setting.time_total,' seconds.' )
                print('Simulation stopped based on ', \
                      setting.steps_max,' steps as setting.steps.max')        
       
        # writing data that has been saved in binary file
        if setting.binsave_writedata == True: 
            NS = trun.index_next_save_element
            file_output_data = open(setting.binsave_data_filename,'wb')
            trun.element_data_saved[:,0:NS].tofile(file_output_data)
            file_output_data.close()
        
            file_output_time = open(setting.binsave_time_filename,'wb')
            trun.element_time_saved[:,0:NS].tofile(file_output_time)
            file_output_time.close()
            
        # close restart files    
        if setting.txtout_writedata == True:
            setting.txtout_file_flowrate.close()
            setting.txtout_file_volume.close()
            setting.txtout_file_hyddepth.close()
            setting.txtout_file_area.close()
            setting.txtout_file_elevation.close()
    
        if setting.inoisy == True:
            print()
            print('...completed.')
 
        return [setting]