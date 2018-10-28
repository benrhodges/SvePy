#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 16:27:13 2018

@author: brh
"""
class Adjustments:
#==============================================================================
    def adjust_spatial_inconsistent_flowrates(self, el, fa, setting, NX):
        '''
        Handle inconsistent flow rate that can result from non-smooth
        bathymetry. Ad hoc adjustment when both the faces have a different
        velocity sign than the element.
        '''
        import numpy as np
        
        # the cells where the flowrate sign is opposite from the face flow
        # rates
        aa = (np.sign(el['flowrate'][:]) *np.sign(fa['flowrate'][0:NX]) < 0 ) \
           & (np.sign(el['flowrate'][:]) *np.sign(fa['flowrate'][1:NX+1]) < 0)\
           & (abs(fa['flowrate'][0:NX])   > setting.flowrate_zero_value)      \
           & (abs(fa['flowrate'][1:NX+1]) > setting.flowrate_zero_value)

        tmp1 = fa['flowrate'][0:NX]
        tmp2 = fa['flowrate'][1:NX+1]   
                
        el['isadhocflowrate'][aa] = True

        # Set flowrate on element to a simple weighted value of faces
        el['flowrate'][aa] = 0.5 * (tmp1[aa] + tmp2[aa])
        el['velocity'][aa] = el['flowrate'][aa] / el['area'][aa]
       
        return el
    
#==============================================================================
    def adjust_Vshaped_flows(self, el, fa, setting, NX): 
        '''
        Damps a V-shaped flow rate across an element
        '''
        import numpy as np
        
        # detect V-shaped flow rate at a grid cell              
        isVshape = np.sign(fa['flowrate'][0:NX]   - el['flowrate'][:]) \
                 * np.sign(fa['flowrate'][1:NX+1] - el['flowrate'][:]) > 0
                 
#        # only adjust cells that are not part of a hydraulic jump        
#        isVshape = isVshape & (fa['jumptype'][0:NX] == 0) \
#                            & (fa['jumptype'][1:NX+1] == 0)        
        
        
        # TEST replace with simple average              
        #el['temp1'][:] = 0.5 * (fa['flowrate'][0:NX] + fa['flowrate'][1:NX+1])
        
        # replace with a time-scaled average
        el['temp1'][:] = (                                                    \
                    el['tscale_up'][:]  * fa['flowrate'][1:NX+1]              \
                  + el['tscale_dn'][:]  * fa['flowrate'][0:NX] )              \
                / ( el['tscale_up'][:]  + el['tscale_dn'][:])

        el['temp1'][0]    = el['flowrate'][0]
        el['temp1'][NX-1] = el['flowrate'][NX-1]
        
        #el['flowrate'][isVshape] = el['temp1'][isVshape]
        
        # Blending the V-adjusted velocity with the original velocity
        # if setting.method_Q_adjust_Vshape_coef = 1.0, then we get
        # just the V-adjusted velocity.
        el['flowrate'][isVshape]  \
            =  setting.method_Q_adjust_Vshape_coef * el['temp1'][isVshape]     \
            +  (1.0 - setting.method_Q_adjust_Vshape_coef)                     \
                * el['flowrate'][isVshape] 

        return el
    
#==============================================================================    
    def negative_volume_reset(self, el, setting):
        
        aa = el['volume'][:] <= setting.volume_zero_value 
        
        el['volume'][aa] = setting.volume_zero_value
        
        return el
#==============================================================================
    def flowrate_oscillation_damping(self, el, fa, geo, setting, NX):
        
        import sys
        
        # treating the damping coeffienct as a simple velocity fraction
        if setting.method_Q_damp_oscillation_type == 'linear_simple':
            el['temp1'][1:NX-1] = setting.method_Q_damp_oscillation_coef      \
                            * (       fa['flowrate'][1:NX-1]                  \
                               -2.0 * el['flowrate'][1:NX-1]                  \
                               +      fa['flowrate'][2:NX]  ) 
                            
            el['flowrate'][1:NX-1] = el['flowrate'][1:NX-1]                   \
                                      + el['temp1'][1:NX-1]  
          
        # treating the damping coefficient as a viscosity                    
        elif setting.method_Q_damp_oscillation_type == 'linear_viscous':
            el['temp1'][1:NX-1] = setting .method_Q_damp_oscillation_coef     \
                            * ( setting.dt / (geo['length'][1:NX-1]**2.0) )   \
                            * (       fa['flowrate'][1:NX-1]                  \
                               -2.0 * el['flowrate'][1:NX-1]                  \
                               +      fa['flowrate'][2:NX]  ) 
                            
            el['flowrate'][1:NX-1] = el['flowrate'][1:NX-1]                   \
                                      + el['temp1'][1:NX-1]  
        else:
            print(setting.setting.method_Q_damp_oscillation_type)
            print('error unknown value for setting.method_Q_damp_oscillation_type')
            sys.exit()
        #endif

        el['velocity'][:] = el['flowrate'][:] / el['area'][:]    

        return el        
#==============================================================================
    
