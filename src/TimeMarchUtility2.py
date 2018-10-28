#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 13:04:54 2018

@author: brh
"""

class TimeMarchUtility:

#==============================================================================
    def plot_at_command_line (self, el, fa, geo, setting, trun, NX ):
        '''
        Provides debug plotting at command line during run
        '''
        import matplotlib.pyplot as plt     
        import sys
                        
        if trun.present_step >= setting.plot_iterstart:
            #note that these should be customized, depending on simulation case
            if (trun.present_time >                                           \
                       (trun.time_next_plotdata - 0.499*setting.dt))          \
                    or (trun.present_step == setting.steps_max):
                        
                
                if setting.plot_element_limits == []:
                    pstartf = 0
                    pendf = NX
                    pstarte = 0
                    pende = NX
                else:
                    pstartf = setting.plot_element_limits[0] 
                    pendf =   setting.plot_element_limits[1] 
                    pstarte = setting.plot_element_limits[0]
                    pende =   setting.plot_element_limits[1] 
                #endif
                
                trun.time_next_plotdata = trun.present_time                   \
                            + setting.plot_time_interval
   
                for plot_type in setting.plot_types:
                    if plot_type == 'free_surface':
                        plt.plot(geo['xvalue'][pstarte:pende],geo['zbottom'][pstarte:pende],'k-')
                        plt.plot(geo['xvalue'][pstarte:pende],el['eta'][pstarte:pende],'bo')
                        plt.plot( fa['xvalue'][pstartf:pendf],fa['etaP'][pstartf:pendf],'gx')
                        plt.plot( fa['xvalue'][pstartf:pendf],fa['etaM'][pstartf:pendf],'r+')
                        if setting.steady_solution_exists:
                            fps = pstarte * setting.steady_solution_resolution
                            fpe = pende * setting.steady_solution_resolution
                            #plt.plot(geo['xvalue'][pstarte:pende],geo['eta_solution'][pstarte:pende])
                            plt.plot(geo['xvalue_fine_solution'][2*fps:2*fpe], geo['eta_fine_solution'][2*fps:2*fpe],'g-')
                        #endif
                        plt.title('elevation')
                        plt.show()
                        
                    elif plot_type == 'flowrate':    
                        plt.plot(geo['xvalue'][pstarte:pende],el['flowrate'][pstarte:pende],'k-o')
                        plt.plot( fa['xvalue'][pstartf:pendf],fa['flowrate'][pstartf:pendf],'rx')
                        plt.title('flowrate')
                        plt.show()  
                        
                    elif plot_type == 'froude':    
                        plt.plot(geo['xvalue'][pstarte:pende],el['froude'][pstarte:pende],'ko')
                        plt.title('Froude No.')
                        plt.show()  

                    elif plot_type == 'issmallvolume':    
                        plt.plot(geo['xvalue'][pstarte:pende],el['issmallvolume'][pstarte:pende],'ro')
                        plt.title('is small volume')
                        plt.show()  
                       
                    #endif
                    
                #endfor     

                #print(el['flowrate'][56:60])
                   
            #endif
        #endif
        
        return trun

#==============================================================================
    def initialize_outputs(self, trun, el, fa, geo, setting, NX):
        '''
        setup for writing files at each time step
        '''
        import FileReadWrite2
        
        frw = FileReadWrite2.FileReadWrite()
        
        # SETUP THE BINARY SAVE FILES
        trun = self.first_binary_save(trun, el, setting)

        # WRITE TEXT FILES INITIALIZATION
        if setting.txtout_writedata == True:
            trun.time_next_txtout_writedata = trun.present_time \
                + setting.txtout_writedata_timeinterval
        #endif
        
        # PLOT DATA    
        trun.time_next_plotdata = trun.present_time +setting.plot_time_interval
        
        # VOLUME CHECKING
        if setting.check_total_volume_stepinterval == 0:
            setting.check_total_volume_stepinterval = setting.steps_max + 1
        #endif
        
        # OUTPUT A FILE OF ALL THE SIMULATION SETTINGS
        frw.write_control_setting(setting)
        
        # OUTPUT A FILE OF GEOMETRY
        frw.geometry_txout_write(geo, fa, setting, trun, NX)

        return [trun, setting]
    
#==============================================================================
    def first_binary_save (self, trun, el, setting):
        '''
        first save of binary data files
        '''      
        import DataStorage2
        data = DataStorage2.DataStorage()
        
        if setting.binsave == True:
            # save the starting data
            trun.element_data_saved = data.saveslice                          \
                (el, trun.element_data_saved, setting.binsave_element,        \
                 trun.index_next_save_element)
                
            # save the starting time
            trun.element_time_saved = data.savetimes                          \
                            (trun.element_time_saved, trun.present_step,      \
                              trun.present_time, trun.index_next_save_element)    
                            
            # increment for the next save time             
            trun.time_next_save_element                                       \
                = trun.present_time +setting.binsave_timeinterval
                
            # increment the index    
            trun.index_next_save_element = trun.index_next_save_element+1
        #endif
        
        return trun
    
#==============================================================================
    def ramp_up_inflow_BC (self, present_time, BC, setting):
        '''
        Adjusts the inflow BC if during a ramp-up period
        '''        
        if present_time < setting.inflow_rampup_time:            
            BC.inflowrate = setting.inflowBCsaved * present_time \
                            / setting.inflow_rampup_time
                            
            if  BC.inflowrate < setting.flowrate_IC:
                BC.inflowrate = setting.flowrate_IC   
            #endif
            
        else:
            BC.inflowrate = setting.inflowBCsaved 
        #endif
        
        return BC

#==============================================================================
    def dynamic_time_step_size(self, el, fa, geo, setting, trun, NX):
        '''
        Increases the model time step if the CFL everywhere is low, and 
        decreases the model time step if the CFL everywhere is high.
        To prevent this from oscillating between increase and decrease on
        successive time steps, there is a stepcounter that keeps track of the
        number of time steps since the last decrease. There must be a number
        of steps (cfl_increase_stepinterval_from_decrease) without a decrease
        before the step size can be increased.
        '''
        import numpy as np
        
        maxcfl1 = max(el['CFL_up'][:])
        maxcfl2 = max(el['CFL_dn'][:])       
        maxcfl = max(maxcfl1,maxcfl2)
        setting.cfl_increase_stepcounter = setting.cfl_increase_stepcounter + 1
        
        # volume flow out of a cell (+ is out)
        facevolbck = -fa['flowrate'][0:NX]   * setting.dt
        facevolfwd =  fa['flowrate'][1:NX+1] * setting.dt
        
        # max out of any cell
        maxfromcell = np.maximum(facevolbck, facevolfwd)
        
        # the volume ratio is equivalent to a CFL
        volratio = maxfromcell / el['volume'] 
        
        # eliminate small volume cells from setting the CFL
        aa = el['issmallvolume'][:] == True
        volratio[aa] = 0.0
        aa[:] = False
        
        maxratio = max(volratio)        
        maxcfl   = max(maxcfl,maxratio)
                 
        # Reset to smaller dt if cfl limit is exceeded
        if maxcfl >= setting.cfl_max:
            setting.dt = setting.dt * 0.75 * setting.cfl_max / maxcfl
            if setting.dt < setting.dt_min:
                setting.dt = setting.dt_min
            if setting.inoisy == True:  
                if trun.present_step%setting.print_time_header_step_interval  \
                    == 0:
                    aa = el['CFL_up'][:] >= setting.cfl_max 
                    print('REDUCED dt TO ',setting.dt,' at ID ',geo['ID'][aa])
                #endif
            #endif
            setting.cfl_increase_stepcounter = 0
        #endif
        
        # reset to larger dt if enough time has gone by since the last
        # time step reduction
        if maxcfl <  setting.cfl_increase_dt and \
            setting.cfl_increase_stepcounter >  \
                setting.cfl_increase_stepinterval_from_decrease: 
            setting.dt = setting.dt * 1.25
            if setting.inoisy == True:
                if trun.present_step%setting.print_time_header_step_interval  \
                    == 0:    
                    print('increased dt to ',setting.dt)
                #endif
            #endif
        #endif
            
        return setting
#==============================================================================
    def total_volume_check(self, el, fa, present_step, BC, setting, NX):
        '''
        diagnostic printout of volume change over a time step
        '''
        import numpy as np
        
        if present_step%setting.print_time_header_step_interval == 0:
            totalvolume2e = np.sum(el['volume'])
            deltavolumeE = totalvolume2e - setting.total_volume_element
            inflow = BC.inflowrate * setting.dt
            outflow = fa['flowrate'][NX] * setting.dt
            volumeerrorE = (deltavolumeE - (inflow-outflow)) / totalvolume2e

            if setting.inoisy == True:
                print('Vol (in,out,cons): ', inflow, outflow, volumeerrorE)
                
        return setting
#==============================================================================
    def end_of_step_checks_and_output                                         \
        (self, el, fa, geo, BC, setting, trun, NX):
        
        import numpy as np
        import sys
        
        import FileReadWrite2
        frw = FileReadWrite2.FileReadWrite()
    
        
        # error checking for NaN
        if any(np.isnan(el['flowrate'])):
            print()
            print('ERROR NaN found ======================================')
            print('VOLUME')
            print(el['volume'][:])
            print()
            print('FLOWRATE')
            print(el['flowrate'][:])
            sys.exit()
    
        # dynamic time step size for high CFL 
        setting = self.dynamic_time_step_size(el, fa, geo, setting, trun, NX)
        
        #error check for volume
        if trun.present_step%setting.check_total_volume_stepinterval == 0:
            setting = self.total_volume_check \
                        (el, fa, trun.present_step, BC, setting, NX)
    
        # WRITE DATA (RESTART FILES)
        trun = frw.time_loop_txtout_write(el, fa, setting, trun, NX) 

        # SAVE DATA (BINARY)    
        trun = frw.time_loop_datasave(el, setting, trun)
        
        # print a . after a time step for a non-noisy simulation
        if setting.inoisy == False and setting.completely_silent == False:
            print('.', end='', flush=True)

        return [trun, setting]     
#==============================================================================
      
#EOF
