#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 15:56:08 2018

@author: brh
"""
class RungeKutta:
#==============================================================================
    def RKstep (self, el, fa, geo, rke, rkf, rkc, BC, setting, NX ):
        '''
        RK method with continuity and momentum colocated at cell center
        with different solution options
        '''
        import sys
        import numpy as np
        import Adjustments2
        import UpdateValues
        adj = Adjustments2.Adjustments()
        uv  = UpdateValues.UpdateValues()
                        
        [rke, rkf, rkc, kV, kQ] = self.initialize_rk                          \
                                    (rke, rkf, rkc, el, fa, setting)
                
        # HACK - note that the ButcherTableC is not used in the following as 
        # time-varying BC are not included

#        print()
#        print('0 ',0 , 0, el['flowrate'][56:59])         
        
        # cycle through the ButcherTable rows
        for ii in range(0,setting.ButcherTableA.shape[0]+1):
            # RHS of equation (without friction)
            rkc = self.volume_step(kV[ii], rkc,      rkf,      NX, setting)    
            rkc = self.flow_step  (kQ[ii], rkc, rke, rkf, geo, NX, setting )

#            print('a ',ii ,0, rkc[kQ[ii]][56:59])         
              
            if ii < (setting.ButcherTableA.shape[0]): 
                # reset to time n values
                rke['volume'][:]   = el['volume'][:]
                rke['flowrate'][:] = el['flowrate'][:]

                for kk in range(0,ii+1):
                    
                    if setting.ButcherTableA[ii,kk] != 0.0:
                        # rk for volume
                        rke['volume'][:]   = rke['volume'][:]  \
                                + setting.ButcherTableA[ii,kk] * setting.dt   \
                                * rkc[kV[kk]][:]
                                
                        # small volume limiter   
                        # This loses exact volume conservation to ensure 
                        # volume positivity.
                        rke = adj.negative_volume_reset(rke, setting)     
                             
                        # rk for flowrate, neglecting friction
                        rke['flowrate'][:] = rke['flowrate'][:] \
                                + setting.ButcherTableA[ii,kk] * setting.dt   \
                                * rkc[kQ[kk]][:]
 
#                        print('b ',ii , kk, rke['flowrate'][56:59])         
                               
                        # effective friction term (with limiter)
                        if ii==kk:                            
                            # friction is flowrate limited at ii=kk iteration
                            # HACK this makes friction inconsistent at different
                            # RK levels, which may be a problem. Need to experiment
                            # with only setting friction value at very start
                            # (after rkc) and keeping it the same throughout.
                            rke = self.friction_effect                        \
                                (rke, geo,setting.ButcherTableA[ii,kk],setting)
                                
                            rkc[kQ[ii]][:] = rkc[kQ[ii]][:]                   \
                                - rke['friction'][:] / geo['length'][:] 
                            
                        else:
                            # this friction was possibly changed when ii<kk
                            # but then is held constant.
                            rke['flowrate'][:] = rke['flowrate'][:] \
                                - setting.ButcherTableA[ii,kk] * setting.dt   \
                                * rke['friction'][:] / geo['length'][:]  
                        #endif
                        
#                        print('c ',ii , kk, rke['flowrate'][56:59])   
                        
                    #endif
                #endfor kk
                
                # HACK may need a small volume velocity limiter that works
                # with the ButcherTable approach    
                
#                print('d ',ii , 0, rke['flowrate'][56:59])         
                
                # update the values
                [rke, rkf] = uv.update_auxiliary_relationships                \
                                (rke, rkf, geo, BC, setting, NX)  

#                print('e ',ii , 0, rke['flowrate'][56:59])         

                       
#                # Reset the rkc storage for ad hoc flow rates
#                aa = rke['isadhocflowrate'][:] == True
#                rkc[kQ[ii]][aa] = (rke['flowrate'][aa] -  el['flowrate'][aa]) \
#                    / (setting.ButcherTableA[ii,ii] * setting.dt)
#                
#                #print('ad hoc ', rkc[kQ[ii]][aa])
#                aa[:] = False    
                
            #endif ii < (setting.ButcherTableA.shape[0])
        #endfor ii   
        
#        print (' - - ')
        # RK step updates
        
        for ii in range(0,setting.ButcherTableA.shape[0]+1):

#            print('f ',ii , 0, rkc[kQ[ii]][56:59])         
                        
            el['volume'][:] = el['volume'][:] \
                    + setting.dt * setting.ButcherTableB[ii] * rkc[kV[ii]][:]

            el['flowrate'][:] = el['flowrate'][:] \
                    + setting.dt * setting.ButcherTableB[ii] * rkc[kQ[ii]][:]
        
#            print(el['flowrate'][56:59])
        #endfor

#        print('g ',0 , 0, el['flowrate'][56:59])         
       
        #print('====')
        #print(el['flowrate'][56:60])
        
        # update all aux values
        [el, fa] = uv.update_auxiliary_relationships                          \
                                (el, fa, geo, BC, setting, NX)  


#        print('h ',0 , 0, el['flowrate'][56:59])         
                                                                  
#        print('after update')
#        print(el['flowrate'][56:60])

        if setting.method_Q_adjust_inconsistent:
            el = adj.adjust_spatial_inconsistent_flowrates(el, fa, setting, NX)
        #endif
        
#        print('==')
        #print(el['flowrate'][56:60])
        
        if setting.method_Q_adjust_Vshape:
            el = adj.adjust_Vshaped_flows(el, fa, setting, NX)
        #endif    

#        print(el['flowrate'][56:60])
#        print()
        #print('---')
        #print('     ', el['volume'][57:59])
        #print(fa['flowrate'][57:60])
        
        
        if setting.method_Q_damp_oscillation_type != None:
            el = adj.flowrate_oscillation_damping(el, fa, geo, setting, NX)  
        #endif    

#        for mm in range(0,NX-1):
#            if fa['jumptype'][ii+1] == +1:
#                el['flowrate'][ii] = el['flowrate'][ii] * 1.01
                
#        print('i ',0 , 0, el['flowrate'][56:59])         
               
        return [el, fa]

#==============================================================================
    def initialize_rk(self, rke, rkf, rkc, el, fa, setting):
        '''
        Sets up rke and rkf arrays with data from el and fa so that we can
        write functions that do not depend on el and fa for the first step
        '''
        import sys
        import UpdateValues
        uv = UpdateValues.UpdateValues() 
        
        [rke, rkf] = uv.zero_all_except_element_flow_and_volume (rke, rkf)
        
        rkf['flowrate'][:] = fa['flowrate'][:]
        rkf['velocityP'][:] = fa['velocityP'][:]
        rkf['velocityM'][:] = fa['velocityM'][:]
        rkf['etaP'][:] = fa['etaP'][:]
        rkf['etaM'][:] = fa['etaM'][:]
        rkf['areaP'][:] = fa['areaP'][:]
        rkf['areaM'][:] = fa['areaM'][:]
        rkf['xvalue'][:] = fa['xvalue'][:]
        rkf['zbottom'][:]= fa['zbottom'][:]
        
        rke['eta'][:] = el['eta'][:]
        rke['friction'][:] = el['friction'][:]
        rke['area'][:] = el['area'][:]
        rke['volume'][:]   = el['volume'][:]
        rke['flowrate'][:] = el['flowrate'][:]
        rke['isadhocflowrate'][:] = el['isadhocflowrate'][:]
        rke['issmallvolume'][:] = el['issmallvolume'][:]
        rke['smallvolume_ratio'][:] = el['smallvolume_ratio'][:]

        if setting.rungekutta_levels == 2:
            kV = ['k1v','k2v']
            kQ = ['k1q','k2q']
        if setting.rungekutta_levels == 3:
            kV = ['k1v','k2v','k3v']
            kQ = ['k1q','k2q','k3q']
        elif setting.rungekutta_levels == 4:
            kV = ['k1v','k2v','k3v','k4v']
            kQ = ['k1q','k2q','k3q','k4q']
        elif setting.rungekutta_levels == 5:
            kV = ['k1v','k2v','k3v','k4v','k5v']
            kQ = ['k1q','k2q','k3q','k4q','k5q']
        elif setting.rungekutta_levels == 6:
            kV = ['k1v','k2v','k3v','k4v','k5v','k6v']
            kQ = ['k1q','k2q','k3q','k4q','k5q','k6q']
        elif setting.rungekutta_levels == 7:
            kV = ['k1v','k2v','k3v','k4v','k5v','k6v','k7v']
            kQ = ['k1q','k2q','k3q','k4q','k5q','k6q','k7q']
        elif setting.rungekutta_levels == 8:
            kV = ['k1v','k2v','k3v','k4v','k5v','k6v','k7v','k8v']
            kQ = ['k1q','k2q','k3q','k4q','k5q','k6q','k7q','k8q']
        elif setting.rungekutta_levels == 9:
            kV = ['k1v','k2v','k3v','k4v','k5v','k6v','k7v','k8v','k9v']
            kQ = ['k1q','k2q','k3q','k4q','k5q','k6q','k7q','k8q','k9q']
        else:
            print('error,unknown value of ',setting.rungekutta_levels, \
                  ' for setting.rungekutta_levels')
            sys.exit()
        #endif

        for ii in range(0,setting.ButcherTableA.shape[0]+1):
            rkc[kV[ii]][:] = 0.0
            rkc[kQ[ii]][:] = 0.0
        #endif
        
        return [rke, rkf, rkc, kV, kQ]

#==============================================================================
    def continuity_RHS_explicit(self, flowrate, ii , setting):
        '''
        RHS for continuity in finite volume is simply the flow rate difference 
        across the cell.
        '''
        return ( flowrate[ii] - flowrate[ii+1] ) 

#==============================================================================
    def flow_step (self, kNq, rkc, rke, rkf, geo, NX, setting):
        '''
        The RHS of momentum for an RK step looped over the domain 
        '''
        for ii in range(0,NX):      
                      
            rkc[kNq][ii] = self.momentum_RHS_no_friction  \
                            ( rke, rkf, geo, ii, setting ) 
           
        return rkc
    
#==============================================================================
    def friction_effect (self, rke, geo, tw, setting):
        '''
        separate friction term in RKE advance with limiter to prevent
        unphysical reversals
        '''
        import numpy as np
                
        # friction's effect on flow rate.     
        rke['friction'][:] = tw * setting.dt * rke['friction'][:]             \
                            / geo['length'][:]
        
        # friction in correct direction, but too large
        # we reduce friction so that it cannot reverse velocity
        aa =  (np.sign(rke['friction'][:]) * np.sign(rke['flowrate'][:]) > 0) \
            & (abs(rke['friction'][:]) > abs(rke['flowrate'][:]) )
            
        rke['friction'][aa] = rke['flowrate'][aa]  \
                - np.sign(rke['flowrate'][aa]) * setting.flowrate_zero_value 
        aa[:] = False
                    
        # friction is in incorrect direction (flow reversal)
        # we set friction to the correct direction and give it a value
        # based on either limiting the friction smaller than velocity,
        # or the magnitude that was previously computed
        # HACK Need to think about this further.
        aa = np.sign(rke['friction'][:]) * np.sign(rke['flowrate'][:]) < 0        

        rke['friction'][aa] = np.sign(rke['flowrate'][aa]) \
                * np.minimum(abs(rke['friction'][aa]),  \
                             abs(rke['flowrate'][aa]) )
        aa[:] = False        
        
        # Update the flow rate for the friction term
        rke['flowrate'][:] = rke['flowrate'][:] - rke['friction'][:]  
        
        # return friction to its standard form
        rke['friction'][:] = rke['friction'][:] * geo['length'][:]            \
                            /  (tw * setting.dt)
        
        return rke
        
#==============================================================================
    def momentum_RHS_no_friction(self, el, fa, geo, ii, setting):  
        '''
        RHS of momentum equation without friction term.
        The uniform option is the simpler version of the pressure term that
        is effectively an approximation of uniform flow.
        The nonuniform option should be used except for testing.
        '''
        import numpy as np
        import sys
        
        if setting.method_momentum == 'uniform':
            fout = (  fa['flowrate'][ii]   * fa['velocityP'][ii]              \
                    * np.sign(fa['flowrate'][ii])                             \
                    - fa['flowrate'][ii+1] * fa['velocityM'][ii+1]            \
                         * np.sign(fa['flowrate'][ii+1])                      \
                    + setting.gravity * el['area'][ii]                        \
                        * ( fa['etaP'][ii] - fa['etaM'][ii+1] )               \
                     ) / geo['length'][ii] 
            
        elif setting.method_momentum == 'T00': 
           
            fout = ( fa['flowrate'][ii] * fa['velocityP'][ii]                 \
                         * np.sign(fa['flowrate'][ii]) * fa['cosangle'][ii]   \
                     - fa['flowrate'][ii+1] * fa['velocityM'][ii+1]           \
                         * np.sign(fa['flowrate'][ii+1]) *fa['cosangle'][ii+1]\
                     + setting.gravity * fa['areaP'][ii  ] * fa['etaP'][ii]   \
                         * fa['cosangle'][ii]                                 \
                     - setting.gravity * fa['areaM'][ii+1] * fa['etaM'][ii+1] \
                         * fa['cosangle'][ii+1]                               \
                     + setting.gravity * el['eta'][ii]                        \
                     * ( fa['areaM'][ii+1] -  fa['areaP'][ii])                \
                    ) / geo['length'][ii] 

        elif setting.method_momentum == 'T10': 
           
            fout = ( fa['flowrate'][ii] * fa['velocityP'][ii]                 \
                         * np.sign(fa['flowrate'][ii])* fa['cosangle'][ii]    \
                     - fa['flowrate'][ii+1] * fa['velocityM'][ii+1]           \
                         * np.sign(fa['flowrate'][ii+1])*fa['cosangle'][ii+1] \
                     + setting.gravity * (fa['cosangle'][ii] - 0.5)           \
                         * fa['areaP'][ii  ] * fa['etaP'][ii]                 \
                     - setting.gravity * (fa['cosangle'][ii+1] - 0.5)         \
                         * fa['areaM'][ii+1] * fa['etaM'][ii+1]               \
                     + setting.gravity * 0.5                                  \
                     * (  fa['areaM'][ii+1] * fa['etaP'][ii]                  \
                        - fa['areaP'][ii]   * fa['etaM'][ii+1])               \
                    ) / geo['length'][ii] 
                    
        elif setting.method_momentum == 'T20': 
           
            fout = ( fa['flowrate'][ii] * fa['velocityP'][ii]                 \
                         * np.sign(fa['flowrate'][ii]) * fa['cosangle'][ii]   \
                     - fa['flowrate'][ii+1] * fa['velocityM'][ii+1]           \
                         * np.sign(fa['flowrate'][ii+1])*fa['cosangle'][ii+1] \
                     + setting.gravity * (fa['cosangle'][ii] - 1.0/6.0)       \
                         * fa['areaP'][ii  ] * fa['etaP'][ii]                 \
                     - setting.gravity * (fa['cosangle'][ii+1] - 1.0/6.0)     \
                         * fa['areaM'][ii+1] * fa['etaM'][ii+1]               \
                     + setting.gravity * (1.0/6.0)                            \
                     * (  fa['areaM'][ii+1]                                   \
                            * (fa['etaP'][ii]   + 4.0 * el['eta'][ii] )       \
                        - fa['areaP'][ii]                                     \
                            * (fa['etaM'][ii+1] + 4.0 * el['eta'][ii] ) )     \
                    ) / geo['length'][ii] 
        else:
            print('error, unknonw setting.method_momentum ', \
                  setting.method_momentum)
            sys.exit()
        #endif
             
        return fout
#==============================================================================
    def volume_step(self, kNv, rkc, rkf, NX, setting):   
        '''
        The RHS of continuity for an RK step looped over the domain
        '''
        for ii in range(0,NX):
            rkc[kNv][ii] = self.continuity_RHS_explicit \
                            ( rkf['flowrate'], ii, setting ) 
        return rkc
        
#==============================================================================
#==============================================================================
#EOF    
