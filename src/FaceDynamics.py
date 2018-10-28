#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 11:07:50 2018

@author: brh
"""
class FaceDynamics:
#==============================================================================
       
    def face_dynamics(self, el, fa, geo,  BC, setting, NX):
        '''
        compute the dynamic values on the face and the element friction
        (which is affected by face velocities)
        '''
        
        fa = self.face_flowrate             (el, fa, geo, BC, setting, NX )
        fa = self.face_velocities           (el, fa, setting, NX)    
        
        if setting.method_froude_interpolation == True:            
            fa = self.face_hydraulic_jump(el, fa, geo, setting, NX)
        #endif
        
        fa = self.face_smallvolume_adjust   (el, fa, setting, NX)
        el = self.friction_with_face_effects(el, fa, geo, setting, NX)
       
        return [fa, el]    
    
#==============================================================================
    def face_flowrate(self, el, fa, geo, BC, setting, NX ):
        '''
        update face flowrate from elements
        '''
        import sys
        
        # simple equal weight interpolation -----------------------------------
        if setting.method_flowrate_interpolation == 'equalweight':
            fa['flowrate'][1:NX] \
                = 0.5 * (el['flowrate'][0:NX-1] + el['flowrate'][1:NX]) 

        # upwind interpolation ------------------------------------------------
        elif setting.method_flowrate_interpolation == 'upwind':
            # HACK assumes downstream flow only
            fa['flowrate'][1:NX] \
                = el['flowrate'][0:NX-1]
            print('error: code needs to be changed for ',                     \
                  'method_flowrate_interpolation = upwind')  
            sys.exit()
        
        # linear interpolation ------------------------------------------------
        elif setting.method_flowrate_interpolation == 'linear':
            fa['flowrate'][1:NX] = (                                          \
                  geo['length'][1:NX]   * el['flowrate'][0:NX-1]              \
                + geo['length'][0:NX-1] * el['flowrate'][1:NX] )              \
              / ( geo['length'][1:NX]   + geo['length'][0:NX-1] )
 
        # timescale interpolation ---------------------------------------------
        elif setting.method_flowrate_interpolation == 'timescale':
            fa['flowrate'][1:NX] = (                                          \
                    el['tscale_up'][1:NX]    * el['flowrate'][0:NX-1]         \
                  + el['tscale_dn'][0:NX-1]  * el['flowrate'][1:NX] )         \
                / ( el['tscale_dn'][0:NX-1]  + el['tscale_up'][1:NX])
        else:
            print('error: unknown value setting.method_face_interpolation ',  \
                  ' of ',setting.method_face_interpolation)
            print()
            sys.exit()
        #endif ----------------------------------------------------------------      
         
        #----------------------------------------------------------------------   
        # boundary conditions          
        fa['flowrate'][0]  = BC.inflowrate 
       
        # HACK simple extrapolation - may not work  at supercritical boundary
        fa['flowrate'][NX] = el['flowrate'][NX-1]
        
        return fa
        
#==============================================================================
    def face_velocities(self, el, fa, setting, NX):
        '''
        face velocities consistent with flows and areas
        '''
        
        fa['velocityM'][:] = fa['flowrate'][:] / fa['areaM'][:]
        fa['velocityP'][:] = fa['flowrate'][:] / fa['areaP'][:]
 
#        # HYDRAULIC JUMP        
#        if setting.method_froude_interpolation == True:
#            # HACK need to make this array processed.
#            for ii in range(0,NX-1):
#                if fa['jumptype'][ii+1] == +1:                        
#                    fa['velocityM'][ii+1] = el['velocity'][ii]                    
#                                    
#                elif fa['jumptype'][ii] == -1:
#                    fa['velocityP'][ii] = el['velocity'][ii]                   
#                #endif
#            #endfor
#        #endif    
       
        return fa
#==============================================================================
    def face_smallvolume_adjust(self, el, fa, setting, NX):
        '''
        ad hoc adjust faces for for small volumes on elements
        '''
        # HACK need to make this array processed.
        if any(el['issmallvolume']):
            for ii in range(0,NX-1):
                if el['issmallvolume'][ii]:                    
                    # if the volume is zero, then zero out the outflows
                    if el['volume'][ii] <= setting.volume_zero_value:
                        if fa['flowrate'][ii+1] >= 0.0:
                            fa['flowrate'][ii+1] = setting.flowrate_zero_value
                            fa['velocityM'][ii+1] = setting.velocity_zero_value
                            fa['velocityP'][ii+1] = setting.velocity_zero_value
                        elif fa['flowrate'][ii] <= 0.0:
                            fa['flowrate'][ii] = -setting.flowrate_zero_value
                            fa['velocityM'][ii] = -setting.velocity_zero_value
                            fa['velocityP'][ii] = -setting.velocity_zero_value                     
                        #endif
                    #endif
                #endif
            #endfor
        #endif
        
        return fa
    
#==============================================================================
    def friction_with_face_effects(self, el, fa, geo, setting, NX):  
        '''
        Friction on an element using the face velocities
        '''
        import numpy as np
        
        # velocities in the upper and lower sections of an element
        # by simple interpolation
        upperVelocity = 0.5*(el['velocity'][:] + fa['velocityP'][0:NX])
        lowerVelocity = 0.5*(el['velocity'][:] + fa['velocityM'][1:NX+1])
        
        # hydraulic radiums in the upper and lower sections of an element
        # by simple interpolation
        upperHydRadius = 0.5*(el['area'][:] + fa['areaP'][0:NX] )             \
                       / (0.5*(el['perimeter'][:] + fa['perimeter'][0:NX]))
        lowerHydRadius = 0.5*(el['area'][:] + fa['areaM'][1:NX+1] )           \
                       / (0.5*(el['perimeter'][:] + fa['perimeter'][1:NX+1]))
                       
        # volume in the upper and lower sections of an element 
        # by simple interpolation               
        upperVolume = 0.5*geo['length'][:]                                    \
                    * 0.5 *(el['area'][:] + fa['areaP'][0:NX])               
        lowerVolume = 0.5*geo['length'][:]                                    \
                    * 0.5 *(el['area'][:] + fa['areaM'][1:NX+1])               
                
        # element friction term that sums the two sections        
        el['friction'][:] = np.sign(el['velocity'][:]) * setting.gravity      \
            *( geo['manningsn'][:]**2.0 )                                     \
            *(  upperVolume * (upperVelocity**2.0) / upperHydRadius**(4.0/3.0)\
              + lowerVolume * (lowerVelocity**2.0) / lowerHydRadius**(4.0/3.0)) 
        
        return el
    
#==============================================================================
    def face_hydraulic_jump(self, el, fa, geo, setting, NX):
        '''
        compute hydraulic jumps on faces
        '''
        import numpy as np
        import sys
        
        fa['jumptype'][:] = 0
        
        # HACK this needs to be rewritten for array processing
        
        for ii in range(0,NX-1):

            #------------------------------------------------------------------      
            # jump with downstream flow 
            # note tha face flowrate is not updated yet, so this includes
            # a time-lagged discrimator the effectively prevents the
            # system from oscillating with jumps.
            if (el['froude'][ii  ] > 1.0                                      \
                    + setting.method_hydjump_Froude_epsilon) and                                 \
               (el['froude'][ii+1] < 1.0                                      \
                    - setting.method_hydjump_Froude_epsilon) and                                 \
               (fa['flowrate'][ii+1] > 0.0)                                   \
               and                                                            \
               (el['eta'][ii] < el['eta'][ii+1]):

               fa['jumptype'][ii+1] = +1
                   
               # baseline simple extrapolation    
               fa['etaM'][ii+1]  = el['eta'][ii] 
               fa['areaM'][ii+1] = el['area'][ii]
               fa['velocityM'][ii+1] = fa['flowrate'][ii+1] / fa['areaM'][ii+1]
               
               fa['etaP'][ii+1]  = el['eta'][ii+1]
               fa['areaP'][ii+1] = el['area'][ii+1]
               fa['velocityP'][ii+1] = fa['flowrate'][ii+1] / fa['areaP'][ii+1]
               
               # using depths on elements
               depthU  = el['hyddepth'][ii]
               froudeU = el['froude'][ii]
               
               depthD  = el['hyddepth'][ii+1]
               froudeD = el['froude'][ii+1]


               # Note that el['flowrate'] < Qnorm implies that we are above
               # the normal depth (i.e. decreasing depth is needed to get
               # a normal Q that meets Chezy-Manning). Hence for Q < Qnorm
               # we are on an S2 curve and expect an immediate jump. In
               # contrast, el['flowrate'] > Qnorm implies that increasing
               # the depth is needed to match C-M at the observed flow, so
               # we are on an S3, M3, C3, H3 or A3 curve. With a #3 curve
               # we use the downstream depth as the known jump depth and 
               # compute the upstream depth. With an S2 curve we use the
               # upstream depth as the known depth and compute the 
               # downstream depth.

               # use conjugate depth relationships for rectangular channel
               slope = (fa['zbottom'][ii]- fa['zbottom'][ii+1])           \
                       / geo['length'][ii]                  
               # normal flow rate
               if slope > 0.0:
                   if (geo['manningsn'][ii] == 0.0):
                       Qnorm = 1e12
                   else:
                       Qnorm = (1.0 / geo['manningsn'][ii])               \
                               * el['area'][ii]                           \
                               * el['hydradius'][ii]**(2.0/3.0)           \
                               * np.sqrt(slope)
                   #endif
               else:
                   Qnorm = 0
               #endif

               # approximation of the discontinuity at the jump    
               if setting.method_hydjump_face == 'simple':
                   pass
               
               elif setting.method_hydjump_face == 'extrapolate_surface':
                   # extrapolate the surface gradients on either side
                   # of the jump to estimate the surface and areas
                   # just upstream and just downstream of jump
                   # note that this option is not no-neighbor compliant
                   
                   fa['etaM'][ii+1] = el['eta'][ii]                           \
                       - (fa['etaP'][ii] - el['eta'][ii])
                       
                   fa['areaM'][ii+1] = el['area'][ii]                         \
                       - (fa['etaP'][ii] - el['eta'][ii])*fa['topwidth'][ii+1]    
                       
                   fa['velocityM'][ii+1] = fa['flowrate'][ii+1]               \
                                           / fa['areaM'][ii+1]
                        
                   fa['etaP'][ii+1] = el['eta'][ii+1]                         \
                       - (fa['etaP'][ii+2] - el['eta'][ii+1])
                       
                   fa['areaP'][ii+1] = el['area'][ii+1]                       \
                       - (fa['etaP'][ii+2] - el['eta'][ii+1])                 \
                       * fa['topwidth'][ii+1]

                   fa['velocityP'][ii+1] = fa['flowrate'][ii+1]               \
                                           / fa['areaP'][ii+1]

               elif setting.method_hydjump_face=='energy_limit':                   
                    # check for energy increase across jump and
                    # set limit on downstream behavior

                    EnergyUp = fa['areaM'][ii+1] \
                            * ( 0.5*(fa['velocityM'][ii+1]**2.0)  \
                                   +setting.gravity * fa['etaM'][ii+1])
                    EnergyDn = fa['areaP'][ii+1] \
                            * ( 0.5*(fa['velocityP'][ii+1]**2.0)  \
                                   +setting.gravity * fa['etaP'][ii+1])
                    
                    if EnergyUp < EnergyDn:
                        Hd = fa['areaP'][ii+1] / fa['topwidth'][ii+1]
                        
                        p0 = 1.0
                        
                        p1 = - (fa['etaP'][ii+1] + 2.0 * Hd)
                        
                        p2 = 2.0 * Hd * fa['etaP'][ii+1] + Hd**2.0 \
                            + EnergyUp / (setting.gravity * fa['topwidth'][ii+1])
                            
                        p3 = EnergyUp * Hd / (setting.gravity * fa['topwidth'][ii+1]) \
                            - (Hd**2.0) * fa['etaP'][ii+1] \
                            - 0.5 * (fa['flowrate'][ii+1]**2.0) \
                               / (setting.gravity * (fa['topwidth'][ii+1]**2.0))
                        
                        tdelta = np.roots([p0,p1,p2,p3])                        
                        tdelta = tdelta.real[abs(tdelta.imag) < 1.e-5]                        
                        etaDelta = fa['etaP'][ii+1] - fa['etaM'][ii+1]                      
                        tdelta = tdelta[(tdelta > 0) and (tdelta < etaDelta)]                       
                        tdelta = min(tdelta)
                        if tdelta > 0.0:
                            fa['etaP'][ii+1] = fa['etaP'][ii+1] - tdelta
                            fa['areaP'][ii+1] = fa['areaP'][ii+1]             \
                                - tdelta *fa['topwidth'][ii+1] 
                            fa['velocityP'][ii+1] = fa['flowrate'][ii+1]      \
                                / fa['areaP'][ii+1]
                    #endif      

               elif setting.method_hydjump_face=='momentum_match':                   
                    # match the downstream momentum with the upstream
                    # momentum.
                    
                    #upstream depth is above normal depth    
                    centroidDepthUp = el['eta'][ii] - geo['zbottom'][ii] \
                        - 0.5 * el['hyddepth'][ii]
                    
                    centroidDepthDn = el['eta'][ii+1] - geo['zbottom'][ii+1] \
                        - 0.5 * el['hyddepth'][ii+1]
                        
                    MomUp = fa['flowrate'][ii+1] * fa['velocityM'][ii+1] \
                        + setting.gravity * fa['areaM'][ii+1] \
                        * centroidDepthUp
                         
                    MomDn = fa['flowrate'][ii+1] * fa['velocityP'][ii+1] \
                        + setting.gravity * fa['areaP'][ii+1] \
                        * centroidDepthDn
                         
                    MomDiff = MomUp - MomDn     
                    if MomDiff > MomUp                                   \
                                    * setting.method_hydjump_Momentum_epsilon:
                                                               
                        #print(ii, ' fixing jump')
                        msign = np.sign(MomDiff)                
                                        
                        Hd = fa['areaP'][ii+1] / fa['topwidth'][ii+1]
                        
                        p0 = 1.0
                    
                        #p1 = msign * 3.0 * Hd
                        p1 = 2.0 * msign * (Hd + centroidDepthDn)
                        
                        #p2 = 3.0 * (Hd**2.0) - 2.0 * MomUp \
                        #    / (setting.gravity * fa['topwidth'][ii+1])
                            
                        p2 = (Hd**2.0)  + 4.0 * Hd * centroidDepthDn          \
                            - 2.0 * MomUp                                     \
                            / (setting.gravity * fa['topwidth'][ii+1])
                            
                        #p3 = msign *(Hd**3.0)                                 \
                        #    + 2.0 * msign                                     \
                        #         *  ( (  (fa['flowrate'][ii+1]**2.0)          \
                        #                - fa['areaP'][ii+1] * MomUp) )        \
                        #       / (setting.gravity *(fa['topwidth'][ii+1]**2.0)) 
    
                        p3 = 2.0 * msign                                      \
                            * (  (fa['flowrate'][ii+1]**2.0)                  \
                                - fa['areaP'][ii+1] * MomUp                   \
                                + setting.gravity * (fa['areaP'][ii+1]**2.0)  \
                                  * centroidDepthDn                           \
                                )                                             \
                            / (setting.gravity *(fa['topwidth'][ii+1]**2.0))
    
                        tdelta = np.roots([p0,p1,p2,p3])                        
                        tdelta = tdelta.real[abs(tdelta.imag) < 1.e-5]                        
                        tdelta = tdelta[(tdelta > 0)]                       
                        tdelta = min(tdelta)
                        if tdelta > 0.0:
                            fa['etaP'][ii+1] = fa['etaP'][ii+1] + msign*tdelta
                            fa['areaP'][ii+1] = fa['areaP'][ii+1]             \
                                + msign * tdelta *fa['topwidth'][ii+1] 
                            fa['velocityP'][ii+1] = fa['flowrate'][ii+1]      \
                                / fa['areaP'][ii+1]
                    #endif                              
              
               elif setting.method_hydjump_face == 'rectangular':
                   
                   
                   if el['flowrate'][ii] <= Qnorm:
                       #upstream depth is above normal depth    
                       depthD = depthU * 0.5 * (-1.0                          \
                          + np.sqrt(1.0 + 8.0 *(froudeU**2.0)))
                       fa['areaP'][ii+1] =  depthD * fa['topwidth'][ii+1]                                         
                       fa['etaP'][ii+1] = fa['zbottom'][ii+1] + depthD
                       fa['velocityP'][ii+1] = fa['flowrate'][ii+1]           \
                                               / fa['areaP'][ii+1]
                   else:
                       depthU = depthD * 0.5 * (-1.0                          \
                          + np.sqrt(1.0 + 8.0 *(froudeD**2.0)))
                       fa['areaM'][ii+1] =  depthU * fa['topwidth'][ii+1]                                         
                       fa['etaM'][ii+1]  = fa['zbottom'][ii+1] + depthU
                       fa['velocityM'][ii+1] = fa['flowrate'][ii+1]           \
                                               / fa['areaM'][ii+1]
                   #endif
                   
                                           

               else:
                   print('error, ',                                           \
                         'unknown value for setting.method_hydjump_face of ', \
                         setting.method_hydjump_face)
                   sys.exit()
               #endif
                             
            #endif
            
            #------------------------------------------------------------------      
            # jump with flow in the nominal upstream direction
            if  (el['froude'][ii-1] < 1.0                                     \
                     + setting.method_hydjump_Froude_epsilon) and                                \
                (el['froude'][ii] > 1.1                                       \
                     - setting.method_hydjump_Froude_epsilon) and                                  \
                (fa['flowrate'][ii] < 0.0)                                    \
                 and                                                          \
                (el['eta'][ii-1] > el['eta'][ii]):

               fa['jumptype'][ii] = -1                    

               # baseline simple extrapolation    
               fa['etaM'][ii]  = el['eta'][ii-1] 
               fa['areaM'][ii] = el['area'][ii-1]
               fa['velocityM'][ii] = fa['flowrate'][ii] / fa['areaM'][ii]
               
               fa['etaP'][ii]  = el['eta'][ii]
               fa['areaP'][ii] = el['area'][ii]
               fa['velocityP'][ii] = fa['flowrate'][ii] / fa['areaP'][ii]

               # using depths on elements with U and D as jump upstream
               # and downstream
               depthU  = el['hyddepth'][ii]
               froudeU = el['froude'][ii]
               
               depthD  = el['hyddepth'][ii-1]
               froudeD = el['froude'][ii-1]

               # use conjugate depth relationships for rectangular channel
               # note this slope is positive if bottom decreasing upstream
               slope = (fa['zbottom'][ii+1]- fa['zbottom'][ii])           \
                       / geo['length'][ii]                  
               # normal flow rate
               if slope > 0.0:
                   # give this Q a negative because el is < 0
                   if (geo['manningsn'][ii] == 0.0):
                       Qnorm = -1e12
                   else:
                       Qnorm = -(1.0 / geo['manningsn'][ii])               \
                               * el['area'][ii]                           \
                               * el['hydradius'][ii]**(2.0/3.0)           \
                               * np.sqrt(slope)
                   #endif
               else:
                   Qnorm = 0
               #endif

               # approximation of the discontinuity at the jump    
               if setting.method_hydjump_face == 'simple':
                   # simple extrapolation of element values
                   # as jump values    
                   pass
               
               elif setting.method_hydjump_face == 'extrapolate_surface':
                   
                   fa['etaM'][ii] = el['eta'][ii-1]                           \
                       - (fa['etaP'][ii-1] - el['eta'][ii-1])
                       
                   fa['areaM'][ii] = el['area'][ii-1]                         \
                       - (fa['etaP'][ii-1] - el['eta'][ii-1])                 \
                       *fa['topwidth'][ii]    
                       
                   fa['velocityM'][ii] = fa['flowrate'][ii]                   \
                                           / fa['areaM'][ii]
                       
  
                   fa['etaP'][ii] = el['eta'][ii]                             \
                       - (fa['etaP'][ii+1] - el['eta'][ii])
                       
                   fa['areaP'][ii] = el['area'][ii]                           \
                       - (fa['etaP'][ii+1] - el['eta'][ii])                   \
                       * fa['topwidth'][ii]

                   fa['velocityP'][ii] = fa['flowrate'][ii]                   \
                                           / fa['areaP'][ii]
                   
               elif setting.method_hydjump_face == 'rectangular':
                   
                   # Note that el['flowrate'] < Qnorm implies that we are above
                   # the normal depth (i.e. decreasing depth is needed to get
                   # a normal Q that meets Chezy-Manning). Hence for Q < Qnorm
                   # we are on an S2 curve and expect an immediate jump. In
                   # contrast, el['flowrate'] > Qnorm implies that increasing
                   # the depth is needed to match C-M at the observed flow, so
                   # we are on an S3, M3, C3, H3 or A3 curve. With a #3 curve
                   # we use the downstream depth as the known jump depth and 
                   # compute the upstream depth. With an S2 curve we use the
                   # upstream depth as the known depth and compute the 
                   # downstream depth.
                   
                   if el['flowrate'][ii] >= Qnorm: # note that el < 0 here
                       #upstream depth is above normal depth    
                       depthD = depthU * 0.5 * (-1                            \
                          + np.sqrt(1.0 + 8.0 *(froudeU**2.0)))
                       fa['areaM'][ii] =  depthD * fa['topwidth'][ii]                                         
                       fa['etaM'][ii]   = fa['zbottom'][ii] + depthD
                       fa['velocityM'][ii] = fa['flowrate'][ii]               \
                                               / fa['areaM'][ii]
                   else:
                       depthU = depthD * 0.5 * (-1                            \
                          + np.sqrt(1.0 + 8.0 *(froudeD**2.0)))
                       fa['areaP'][ii] =  depthU * fa['topwidth'][ii]                                         
                       fa['etaP'][ii]  = fa['zbottom'][ii] + depthU
                       fa['velocityP'][ii] = fa['flowrate'][ii]               \
                                               / fa['areaP'][ii]
                   #endif
                   
                                           
               else:
                   print('error, ',                                           \
                         'unknown value for setting.method_hydjump_face of ', \
                         setting.method_hydjump_face)
                   sys.exit()
               #endif                                                                        
            #endif
            
        #endfor
          
        return fa
#==============================================================================
