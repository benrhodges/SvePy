#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 10:37:01 2018

@author: brh
"""

class ElementDynamics:
    
#==============================================================================
    def element_dynamics(self, el, fa, geo, setting):
        '''
        update the element values that are dynamic
        Note that any use of fa[:] values here are lagged based on prior update
        because the fa[:] are updated after the element.
        '''
        import numpy as np
        import sys

        NX = setting.NX

        el['isadhocflowrate'][:] = False
        el['velocity'][:] = 0.0
                        
        #----------------------------------------------------------------------   
        # baseline velocity update                    
        el['velocity'][:] = el['flowrate'][:] / el['area'][:]
        
        #----------------------------------------------------------------------   
        # ad hoc adjustment for absolute velocity limit
        if setting.velocity_limit > 0:
            aa = abs(el['velocity'][:]) > setting.velocity_limit
            el['velocity'][aa] = np.sign(el['velocity'][aa])                  \
                            * 0.99 * setting.velocity_limit
            el['flowrate'][aa] = el['velocity'][aa] * el['area'][aa]
            el['isadhocflowrate'][aa] = True
            aa[:] = False
        #endif
        
        #----------------------------------------------------------------------   
        #HACK small volumes a simple surface gradeint argument
        # uses linear combination of previously computed and Chezy-Manning
        el = self.blended_smallvolume_velocity(el, fa, geo, setting, NX)

        #----------------------------------------------------------------------   
        # ad hoc adjustment for near zero volumes
        if setting.volume_zero_value > 0.0:  
            aa = el['volume'] <= setting.volume_zero_value 
            el['flowrate'][aa] = np.sign(el['flowrate'][aa])                  \
                            * setting.flowrate_zero_value
            el['velocity'][aa] = np.sign(el['flowrate'][aa])                  \
                            * setting.velocity_zero_value                           
            el['isadhocflowrate'][aa] = True
            aa[:] = False
        #endif
                        
        #----------------------------------------------------------------------  
        # dynamics based on time scales
        el = self.time_scale_computations(el, geo, setting)
        
        return el
    
#==============================================================================
    def blended_smallvolume_velocity(self, el, fa, geo,setting, NX):
        '''
        Blend the computed velocity with a Chezy-Manning approximation for
        small volumes
        '''
        import numpy as np
        
        # HACK - this needs to be turned into array-processed rather than loop
        for ii in range(0,NX):
            if el['issmallvolume'][ii] == True:
                # Slope of water surface, which is an approximation of the
                # energy grade line in a small volume where the velocity is
                # expected to be small.
                # Note that this is a time-lagged surface gradient, which should
                # be addequate for ad hoc small volume adjustment
                eslope = (fa['etaP'][ii] - fa['etaM'][ii+1]) /geo['length'][ii]
                
                # velocity by Chezy-Manning
                if geo['manningsn'][ii] > 0.0:
                    velocity = np.sign(eslope) * (1.0 /geo['manningsn'][ii])  \
                            * el['hydradius'][ii]**(2.0/3.0)                  \
                            * np.sqrt(abs(eslope))
                else:
                    # for n=0, set a default value
                    velocity = np.sign(eslope)                                \
                            * (1.0 /setting.depth_small_value_default_n)      \
                            * el['hydradius'][ii]**(2.0/3.0)                  \
                            * np.sqrt(abs(eslope))
                #endif
                
                # velocity blending by smallvolume_ratio
                velocity = el['smallvolume_ratio'][ii] * el['velocity'][ii]   \
                    + (1.0 - el['smallvolume_ratio'][ii]) * velocity
                
                # if raw and adjusted velocities have the same sign, then
                # use the smaller of the two. Otherwise, use the blended value.
                if np.sign(velocity) * np.sign(el['velocity'][ii]) > 0:
                    el['velocity'][ii] = np.sign(velocity)                    \
                        * min(abs(velocity),abs(el['velocity'][ii]))
                else:
                    el['velocity'][ii] = velocity
                #endif
    
                # adjust flowrate to match the adhoc velocity
                el['flowrate'][ii] = el['velocity'][ii] * el['area'][ii]
                el['isadhocflowrate'][ii] = True
                
            #endif            
        #endfor   
 
        return el
#==============================================================================
    def time_scale_computations(self, el, geo, setting):
        '''
        element dynamics based on time scales
        '''
        import numpy as np
        
        wavespeed = np.sqrt(setting.gravity * el['hyddepth'][:] )
        
        el['froude'][:]  = np.abs(el['velocity'][:]) / wavespeed[:] 
                            
        el['CFL_dn'][:] = ( setting.dt / geo['length'][:] )                  \
                        * ( el['velocity'][:] + wavespeed[:] )
       
        el['CFL_up'][:] = -( setting.dt / geo['length'][:] )                  \
                         * ( el['velocity'][:] - wavespeed[:] )
                         
        el['tscale_up'][:] = -geo['length'][:] \
                                / ( el['velocity'][:] - wavespeed[:] )                  

        el['tscale_dn'][:] = geo['length'][:] \
                                / ( el['velocity'][:] + wavespeed[:] )  
                                
        aa = el['tscale_dn'][:] < 0.0  
        el['tscale_dn'][aa] = setting.time_scale_maximum
        aa[:] = False

        aa = el['tscale_up'][:] < 0.0  
        el['tscale_up'][aa] = setting.time_scale_maximum
        aa[:] = False
        
        aa = el['tscale_dn'][:] < setting.time_scale_minimum  
        el['tscale_dn'][aa] = setting.time_scale_minimum 
        aa[:] = False

        aa = el['tscale_up'][:] < setting.time_scale_minimum  
        el['tscale_up'][aa] = setting.time_scale_minimum 
        aa[:] = False

        aa = el['tscale_dn'][:] > setting.time_scale_maximum  
        el['tscale_dn'][aa] = setting.time_scale_maximum
        aa[:] = False

        aa = el['tscale_up'][:] > setting.time_scale_maximum  
        el['tscale_up'][aa] = setting.time_scale_maximum
        aa[:] = False
        
        return el
    
#==============================================================================
#==============================================================================
    