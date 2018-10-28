#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 14:20:56 2018

@author: brh
"""

class SimulationRun:
  
#==============================================================================
    
    def simulation_toplevel(self, custom):
        ''' 
        Required initialization for a simulation run.
        '''
        import sys

        import FileReadWrite2    
        frw  = FileReadWrite2.FileReadWrite()
        
        import SystemGeometry2
        sg   = SystemGeometry2.SystemGeometry()
        
        import FinishSimulation
        fs   = FinishSimulation.FinishSimulation()

        import SettingDefault as sd

        #----------------------------------------------------------------------               
        # Initialize the default settings for simulation
        setting = sd.SettingDefault() 
        
        # make sure we have the correct working directory
        setting = frw.set_working_directory(setting)  

        # Adjust the setting values for a custom values.
        # At this point, we really only need the case_description and
        setting.user_case_name = custom['user_case_name']
        
       # Adjust setting for a pre-defined case
        setting = self.predefined_case_setup(setting)
        
        # Overwrite the pre-defined case values with custom values
        setting.custom_overwrite(setting, custom)
        
        #----------------------------------------------------------------------           
        # Get the number of cells in the reach
        # HACK data storage and computation for multiple reaches has not 
        #   been implemented. Need to build data structures with different
        #   values of NX in each reach and connectivity between reaches.
        setting.NX = sg.get_number_of_cells(setting)
        NX = setting.NX
        
        #----------------------------------------------------------------------  
        # Boundary Conditions         
        # HACK need to develop dictionary of BC
        def BC(self):
            pass
                           
        #----------------------------------------------------------------------           
        # Setup geometry, bc, intial conditions, etc.
        [geo, el, fa, rkc, rke, rkf, setting, trun, BC, NX]                   \
            = self.setup_simulation(setting, BC)
                        
        if setting.inoisy == True:    
            print('Finished with preliminaries, starting time-marching loop')
            print()

        [el, fa] = self.time_marching_loop                                    \
                    (el, fa, geo, rke, rkf, rkc, setting, trun, BC, NX)
                    
        setting = fs.final_housekeeping(setting, trun)
    
        return [trun, setting]

#==============================================================================
    def predefined_case_setup(self, setting):
        '''
        These are pre-defined cases
        '''
        import sys
        
        if setting.user_case_name[0:16] == 'flow_over_a_bump':
            import CaseFlowOverBump
            cfob = CaseFlowOverBump.CaseFlowOverBump()
            setting = cfob.setting_CaseFlowOverBump(setting)
            
        elif setting.user_case_name == 'Waller_Creek_baseline':
            import CaseWallerCreek
            cwc = CaseWallerCreek.CaseWallerCreek()
            setting = cwc.setting_CaseWallerCreek(setting)
            # ensure that time 0 values (no restart) are used
            setting.IC_type = 'custom'

        else:
            print('error, generalized customization not yet created')
            sys.exit()
        #endif  
                
        return setting        

#==============================================================================
    def setup_simulation(self,setting, BC):
        '''
        Procedures that initialize a simulation
        '''
        import BoundaryConditions2
        import DataStorage2
        import FileReadWrite2    
        import InitialConditions2
        import SystemGeometry2
        import UpdateValues
        import TimeMarchUtility2
        import sys
#        
#        import sys
#        
        bc   = BoundaryConditions2.BoundaryConditions()
        data = DataStorage2.DataStorage()
        frw  = FileReadWrite2.FileReadWrite()
        ic   = InitialConditions2.InitialConditions()
        sg   = SystemGeometry2.SystemGeometry()
        tmu  = TimeMarchUtility2.TimeMarchUtility()    
        uv   = UpdateValues.UpdateValues()
#    
        # HACK
        def trun(self):
            pass
         
        #----------------------------------------------------------------------                       
        # DEFINE SYSTEM GEOMETRY    
        [geo, el, fa, rkc, rke, rkf, setting] = sg.set_all_geometry(setting)
        
        # presently designed for a single reach with NX elements
        NX = setting.NX
        
        #----------------------------------------------------------------------                       
        # PRELIMINARIES - FILES AND STORAGE    
        setting = frw.initialize_all_writefile_names(setting)
    
        # OUTPUT (RESTART) FILE SETUP
        # Add time stamp to output filenames so that they cannot be overwritten
        if setting.txtout_writedata == True:
            setting = frw.open_txtout_file(setting, NX)
        #endif    
                                
        # INPUT OF RESTART FILE
        if setting.IC_type == 'restart_file':
            trun = frw.read_restart_files(setting, trun)        
        else:
            trun.restart_time = []
        #endif    
            
        # STORAGE FOR BINARY OUTPUT
        if setting.binsave == True:    
            # storage for simulation data (for binary output at end of run)
            trun.element_data_saved = data.savedstorage                                \
                (NX, setting.binsave_element, setting, trun.restart_time)
            # storage for time array                    
            trun.element_time_saved = data.timestorage(setting)                   
            trun.index_next_save_element = 0
        #endif    
            
        #----------------------------------------------------------------------                       
        # DEFINE SYSTEM BOUNDARY CONDITIONS    
        # HACK this needs revision with new Boundary Condition storage
        BC.inflowrate = bc.flow(setting)
        BC.outheight  = bc.height(setting)    
        setting.inflowBCsaved = BC.inflowrate
    
        #----------------------------------------------------------------------                       
        # DEFINE ANALYTICAL SOLUTION
        # requires pyswashes package and swashes compiled executable.
        if setting.steady_solution_exists == True:  
            [geo, setting] = self.swashes_solution(geo, setting, NX)
        #endif
        
        #----------------------------------------------------------------------                       
        # DEFINE SYSTEM INITIAL CONDITIONS
        if setting.IC_type == 'restart_file':
            el['flowrate'][:] = trun.flowrate_restart[:]
            el['volume'][:]   = trun.volume_restart[:]
            trun.present_time = trun.restart_time 
                   
        else:
            trun.present_time = 0.0
            if setting.inflow_rampup == True:              
                # HACK need to revise BC approach
                BC = tmu.ramp_up_inflow_BC (trun.present_time, BC, setting)
            #endif
            
            el = ic.initial_flow_and_volume(el, geo, BC, setting, NX)       
        #endif
        
        # always start the present run from a zero step    
        trun.present_step = 0
        
        #----------------------------------------------------------------------                       
        # COMPUTE ELEMENT AND FACE AUXILIARY DATA
        
        # set the element relationships
        [el, fa] = uv.update_auxiliary_relationships                          \
                    (el, fa, geo, BC, setting, NX)
                      
        #----------------------------------------------------------------------                       
        # SETTING UP FOR VARIOUS OUTPUTS         
        [trun, setting] = tmu.initialize_outputs(trun, el, fa, geo, setting, NX)
        
        return [geo, el, fa, rkc, rke, rkf, setting, trun, BC, NX] 
    
#==============================================================================
    def time_marching_loop                                                    \
        (self, el, fa, geo, rke, rkf, rkc, setting, trun, BC, NX):
        
        import numpy as np
        import sys
        import RungeKutta
        import TimeMarchUtility2
        import UpdateValues
        rk  = RungeKutta.RungeKutta()
        tmu = TimeMarchUtility2.TimeMarchUtility()
        uv  = UpdateValues.UpdateValues()
        
        while trun.present_time <  setting.time_total and                     \
              trun.present_step <= setting.steps_max:

            # PLOT DATA WITHIN TIME-MARCHING LOOP        
            if setting.inoisy == True:
                trun = tmu.plot_at_command_line(el, fa, geo, setting, trun, NX)
            #endif
               
            # increment the time and step counters      
            trun.present_time = trun.present_time + setting.dt
            trun.present_step = trun.present_step + 1       
    
            # store present volume for volume check
            if trun.present_step%setting.check_total_volume_stepinterval == 0:
                setting.total_volume_element = np.sum(el['volume'])
            #endif
            
            # print header/debug informations
            if setting.inoisy == True:
                if trun.present_step%setting.print_time_header_step_interval  \
                    == 0:
                    print('==================================================')
                    print(trun.present_step, '= step; ', \
                          trun.present_time, '= seconds;', \
                          setting.dt, '= dt;')
                    print('max cfl =',max(el['CFL_dn']), max(el['CFL_up']))                    
            #endif
            
            # adjust BC during ramp-up time
            if setting.inflow_rampup == True:
                # HACK need to revise BC
                BC = tmu.ramp_up_inflow_BC (trun.present_time, BC, setting)
            #endif
                
            # HACK zero everything to prevent accidental carry-over
            # after complete debugging, this could be removed.
            [el,fa] = uv.zero_all_except_element_flow_and_volume(el, fa)

            [el, fa] = uv.update_auxiliary_relationships                      \
                        (el, fa, geo, BC, setting, NX)
                        
            [el,fa] = rk.RKstep(el, fa, geo, rke, rkf, rkc, BC, setting, NX )
               
            [trun,setting] = tmu.end_of_step_checks_and_output                \
                                (el, fa, geo, BC, setting, trun, NX) 
                                
        return [el, fa]
#==============================================================================
    def swashes_solution(self, geo, setting, NX)    :
        import pyswashes
        import sys
        import numpy as np
        import matplotlib.pyplot as plt     
        
        stype = setting.swashes_stype # 1 # Bump
        domain = setting.swashes_domain # 1 # L=25
        choice = setting.swashes_choice # 3 # transcritical with shock
        ncell = NX #  number of cells
        # address of the binary file. Note that on OSX or Linux that
        # you may need to do a chmod +x swashes to make the file executable
        swashes_bin = '/Users/brh/Box Sync/AA_Sync/CODE/SWASHES-1.03.00/bin/swashes'
        swashout = pyswashes.OneDimensional(stype,domain, choice, ncell, swashes_bin)
        temp = swashout.dataframe()['head']
        
        setting.steady_solution_type = 'eta'       
        geo['eta_solution'] = np.zeros(NX,dtype=np.float64)        
        geo['eta_solution'][:] = temp.values
        
        #print(geo['eta_solution'][:])
        #sys.exit()
              
        xvalues = temp.index
        #xx = temp.index
        #head = temp.values
        #print(xx)
        #print(head)
        
        flowrate = setting.inflowBC 
        
        ps = 1
        pe = NX

        
#        # conjugate depth computation
#        depth = geo['eta_solution'][:] - geo['zbottom'][:]
#        area  = depth[:] * geo['breadth'][:]
#        velocity = flowrate / area[:]
#        froude = velocity[:] / np.sqrt(setting.gravity * depth[:])
#        Cdepth = depth[:] * 0.5 * (-1  + np.sqrt(1.0 + 8.0 *(froude**2.0)))
#        Ceta = Cdepth[:] + geo['zbottom'][:]
#        plt.plot(geo['xvalue'][ps:pe], Ceta[ps:pe],'r-o')
        
        plt.plot(xvalues[ps:pe], geo['eta_solution'][ps:pe],'b-o')
        plt.plot(geo['xvalue'][ps:pe], geo['eta_solution'][ps:pe],'rx')
        plt.plot(geo['xvalue'][ps:pe], geo['zbottom'][ps:pe],'k-+')
        
        setting.steady_solution_resolution = 20
        fps = setting.steady_solution_resolution*ps
        fpe = setting.steady_solution_resolution*pe
        ncell = setting.steady_solution_resolution*NX
        swashout = pyswashes.OneDimensional(stype,domain, choice, ncell, swashes_bin)
        temp2 = swashout.dataframe()['head']
        
        geo['eta_fine_solution']    = np.zeros(ncell,dtype=np.float64) 
        geo['xvalue_fine_solution'] = np.zeros(ncell,dtype=np.float64) 
        
        geo['eta_fine_solution'][:] = temp2.values
        geo['xvalue_fine_solution'][:] = temp2.index
                
        plt.plot(geo['xvalue_fine_solution'][fps:fpe], geo['eta_fine_solution'][fps:fpe],'g-')
        plt.show()
    
        #sys.exit()
        
        return [geo, setting]


#EOF   