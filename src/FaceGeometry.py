#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 11:07:50 2018

@author: brh
"""

class FaceGeometry:
#==============================================================================
    def face_geometry(self, el, fa, geo, BC, setting, NX): 
        '''
        compute all the face geometry from the element geometry 
        '''                
        # HACK - the following functions could be simplified and combined
        # into a single face interpolation function
        fa = self.face_topwidth         (el, fa, geo, setting, NX)
        fa = self.face_wetted_perimeter (el, fa, geo, setting, NX)  
        fa = self.face_free_surface     (el, fa, geo, BC, setting, NX)
        fa = self.face_area             (el, fa, geo, BC, setting, NX)
        
#        if setting.method_froude_interpolation == True:            
#            fa = self.face_hydraulic_jump(el, fa, geo, setting, NX)
#        #endif
        
        return fa
    
#==============================================================================
    def face_topwidth(self, el, fa, geo, setting, NX):
        '''
        compute face topwidth by interpolation
        '''       
        import sys
        
        # equal weight interpolation ------------------------------------------
        if setting.method_topwidth_interpolation == 'equalweight':
            fa['topwidth'][1:NX] \
                = 0.5 * (el['topwidth'][0:NX-1] + el['topwidth'][1:NX]) 
                
        # upwind interpolation ------------------------------------------------
        elif setting.method_topwidth_interpolation == 'upwind':
            # HACK assumes flow direction
            fa['topwidth'][1:NX] \
                = el['topwidth'][0:NX-1] 
            print('error: code needs to be changed for ',                     \
                  'method_topwidth_interpolatio n= upwind')  
            sys.exit()
            
        # linear interpolation ------------------------------------------------
        elif setting.method_topwidth_interpolation == 'linear':
            fa['topwidth'][1:NX] = (                                          \
                  geo['length'][1:NX]   * el['topwidth'][0:NX-1]              \
                + geo['length'][0:NX-1] * el['topwidth'][1:NX] )              \
              / ( geo['length'][1:NX]   + geo['length'][0:NX-1] )
                           
        # timescale interpolation ---------------------------------------------
        elif setting.method_topwidth_interpolation == 'timescale':
            fa['topwidth'][1:NX] = (                                          \
                    el['tscale_up'][1:NX]    * el['topwidth'][0:NX-1]         \
                  + el['tscale_dn'][0:NX-1]  * el['topwidth'][1:NX] )         \
                / ( el['tscale_dn'][0:NX-1]  + el['tscale_up'][1:NX])
                
        else:
            print('error: unknown value ',                                    \
                  'setting.method_topwidth_interpolation ',                   \
                  ' of ',setting.method_topwdith_interpolation)
            print()
            sys.exit()            
        #endif  
        
        #----------------------------------------------------------------------   
        # boundary conditions 
        fa['topwidth'][0]  = el['topwidth'][0]  
        fa['topwidth'][NX] = el['topwidth'][NX-1]  
 
        return fa
    
#==============================================================================
    def face_wetted_perimeter(self, el, fa, geo, setting, NX):      
        '''
        compute face wetted perimeter by interpolation
        '''
        import sys
        
        # equal weight interpolation ------------------------------------------
        if setting.method_perimeter_interpolation == 'equalweight':
            fa['perimeter'][1:NX] \
                = 0.5 * (el['perimeter'][0:NX-1] + el['perimeter'][1:NX]) 
                
        # upwind interpolation ------------------------------------------------
        elif setting.method_perimeter_interpolation == 'upwind':
            # HACK assumes flow direction
            fa['perimeter'][1:NX] \
                = el['perimeter'][0:NX-1] 
            print('error: code needs to be changed for ',                     \
                  'method_perimeter_interpolation = upwind')  
            sys.exit()
                
        # linear interpolation ------------------------------------------------
        elif setting.method_perimeter_interpolation == 'linear':
            fa['perimeter'][1:NX] = (                                         \
                  geo['length'][1:NX]   * el['perimeter'][0:NX-1]             \
                + geo['length'][0:NX-1] * el['perimeter'][1:NX] )             \
              / ( geo['length'][1:NX]   + geo['length'][0:NX-1] )
                           
        # timescale interpolation ---------------------------------------------
        elif setting.method_perimeter_interpolation == 'timescale':
            fa['perimeter'][1:NX] = (                                         \
                    el['tscale_up'][1:NX]    * el['perimeter'][0:NX-1]        \
                  + el['tscale_dn'][0:NX-1]  * el['perimeter'][1:NX] )        \
                / ( el['tscale_dn'][0:NX-1]  + el['tscale_up'][1:NX])
        else:
            print('error: unknown value setting.method_perimeter_',           \
                  'interpolation of ',setting.method_perimeter_interpolation)
            print()
            sys.exit()
        #endif  

        #----------------------------------------------------------------------   
        # boundary conditions 
        fa['perimeter'][0]  = el['perimeter'][0]  
        fa['perimeter'][NX] = el['perimeter'][NX-1]  
    
        return fa
    
#==============================================================================
    def face_free_surface(self, el, fa, geo, BC, setting, NX):
        '''
        compute face free surface elevation by interpolation
        '''
        import sys
        
        # equal weight interpolation ------------------------------------------
        if setting.method_eta_interpolation == 'equalweight':
            fa['etaM'][1:NX] \
                = 0.5 * (el['eta'][0:NX-1] + el['eta'][1:NX])
                
        # upwind interpolation ------------------------------------------------
        elif setting.method_eta_interpolation == 'upwind':
            # HACK assumes flow direction
            fa['etaM'][1:NX] = el['eta'][0:NX-1]
            print('error: code needs to be changed for ',                     \
                  'method_eta_interpolation = upwind')  
            sys.exit()

        # linear interpolation ------------------------------------------------
        elif setting.method_eta_interpolation == 'linear': 
            fa['etaM'][1:NX] =  (                                             \
                  geo['length'][1:NX]   * el['eta'][0:NX-1]                   \
                + geo['length'][0:NX-1] * el['eta'][1:NX] )                   \
              / ( geo['length'][1:NX]   + geo['length'][0:NX-1] )

        elif setting.method_eta_interpolation == 'timescale':
            fa['etaM'][1:NX] = (                                              \
                    el['tscale_up'][1:NX]    * el['eta'][0:NX-1]              \
                  + el['tscale_dn'][0:NX-1]  * el['eta'][1:NX] )              \
                / ( el['tscale_dn'][0:NX-1]  + el['tscale_up'][1:NX])
                
        else:
            print('error: unknown value setting.method_eta_interpolation ',   \
                  ' of ',setting.method_eta_interpolation)
            print()
            sys.exit()
        #endif   
        
        #----------------------------------------------------------------------   
        # boundary conditions 
        fa['etaM'][0] = el['eta'][0] 
                        
        fa['etaM'][NX] = BC.outheight 
              
        #----------------------------------------------------------------------   
        # set a single value of free surface (changed elsewhere)
        fa['etaP'][:] = fa['etaM'][:]
        
        return fa
    
#==============================================================================
    def face_area(self, el, fa, geo, BC, setting, NX):
        '''
        compute face area by interpolation
        '''
        import sys
        
        # equal weight interpolation ------------------------------------------
        if setting.method_area_interpolation == 'equalweight':
            fa['areaM'][1:NX] = 0.5 * (el['area'][0:NX-1] + el['area'][1:NX]) 
                
        # upwind interpolation ------------------------------------------------
        elif setting.method_area_interpolation == 'upwind':
            # HACK assumes flow direction
            fa['areaM'][1:NX] = el['area'][0:NX-1]
            print('error: code needs to be changed for ',                     \
                  'method_area_interpolation = upwind')  
            sys.exit()

        # linear interpolation ------------------------------------------------
        elif setting.method_area_interpolation == 'linear':
            fa['areaM'][1:NX] = (                                             \
                  geo['length'][1:NX]   * el['area'][0:NX-1]                  \
                + geo['length'][0:NX-1] * el['area'][1:NX] )                  \
              / ( geo['length'][1:NX]   + geo['length'][0:NX-1] )
              
        # timescale interpolation ---------------------------------------------
        elif setting.method_area_interpolation == 'timescale':
            fa['areaM'][1:NX] = (                                             \
                    el['tscale_up'][1:NX]    * el['area'][0:NX-1]             \
                  + el['tscale_dn'][0:NX-1]  * el['area'][1:NX] )             \
                / ( el['tscale_dn'][0:NX-1]  + el['tscale_up'][1:NX])
                
        else:
            print('error: unknown value setting.method_area_interpolation ',  \
                  ' of ',setting.method_area_interpolation)
            print()
            sys.exit()
        #endif-----------------------------------------------------------------    
        
        #----------------------------------------------------------------------   
        # Boundary conditions
        fa = self.face_area_BC(el, fa, geo, BC, setting, NX)

        #----------------------------------------------------------------------   
        # require a single face area at each face (overwritten elsewhere)
        fa['areaP'][:]  = fa['areaM'][:]

        return fa
    
#==============================================================================
    def face_area_BC(self, el, fa, geo, BC, setting, NX):
        '''
        BC on face areas
        '''
        import numpy as np
        import sys
        
        #----------------------------------------------------------------------   
        # boundary conditions 
        # upstream boundary by extrapolation
        fa['areaM'][0]  = el['area'][0]
        
        # downstream boundary depth
        depthB = BC.outheight - fa['zbottom'][NX]    
        if depthB < setting.depth_zero_value:
            print()
            print(depthB, BC.outheight, fa['zbottom'][NX],NX)
            print('error: downstream boundary condition height makes depth ', \
                  'less than setting.depth_min')
            sys.exit()
        #endif

        #----------------------------------------------------------------------   
        # downstream BC depends on channel type        
        if geo['etype'][NX-1] == 'rectangular_channel':
            fa['areaM'][NX] = geo['breadth'][NX-1] * depthB   
            
        elif geo['etype'][NX-1] == 'parabolic_channel':
            fa['areaM'][NX] = (4.0/3.0)*np.sqrt((depthB**3.0)                 \
                  /geo['parabola_value'][NX-1])  
            
        elif geo['etype'][NX-1] == 'trapezoidal_channel':
            fa['areaM'][NX] = geo['breadth'][NX-1] * depthB \
                        + (depthB**2.0) / np.tan(geo['trapezoid_angle'][NX-1])   
                        
        elif geo['etype'][NX-1] == 'widthdepth_pair':    
            # aliases for indexes
            widthAtLayerTop                                                   \
                = setting.geometry_widthdepth_values['widthAtLayerTop']
            depthAtLayerTop                                                   \
                = setting.geometry_widthdepth_values['depthAtLayerTop']
            areaTotalBelowThisLayer                                           \
                = setting.geometry_widthdepth_values['areaTotalBelowThisLayer']
            Dwidth = setting.geometry_widthdepth_values['Dwidth']
            Ddepth = setting.geometry_widthdepth_values['Ddepth']
            angle =  setting.geometry_widthdepth_values['angle']
            
            # depth sets bracketing each level
            depthLow = geo['widthdepth'][NX-1,:,depthAtLayerTop] \
                        - geo['widthdepth'][NX-1,:,Ddepth]
            
            depthAbove =  geo['widthdepth'][NX-1,:,depthAtLayerTop]
            
            # mask for level containing depth
            aa = (depthB > depthLow) & (depthB <= depthAbove)
            
            # the difference between local depth and the lower widthdepth level
            deltaD = depthB - depthLow[aa]

            # the area below + the trapezoid at the level
            fa['areaM'][NX]                                                   \
                = geo['widthdepth'][NX-1,aa,areaTotalBelowThisLayer]          \
                + (geo['widthdepth'][NX-1,aa,widthAtLayerTop]                 \
                   - geo['widthdepth'][NX-1,aa,Dwidth]                        \
                   )* deltaD                                                  \
                + (deltaD**2.0) / np.tan(geo['widthdepth'][NX-1,aa,angle])
        
            aa[:] = False
                
        else:       
            print('error: unknown value for etype of ',geo['etype'][NX-1])
            sys.exit()
        #endif-----------------------------------------------------------------         
                
        return fa   
     
##==============================================================================
#    def face_hydraulic_jump(self, el, fa, geo, setting, NX):
#        '''
#        compute hydraulic jumps on faces
#        '''
#        import numpy as np
#        import sys
#        
#        fa['jumptype'][:] = 0
#        
#        # HACK this needs to be rewritten for array processing
#        
#        for ii in range(0,NX-1):
#
#            #------------------------------------------------------------------      
#            # jump with downstream flow 
#            # note tha face flowrate is not updated yet, so this includes
#            # a time-lagged discrimator the effectively prevents the
#            # system from oscillating with jumps.
#            if (el['froude'][ii  ] > 1.1) and                                 \
#               (el['froude'][ii+1] < 0.9) and                                 \
#               (fa['flowrate'][ii+1] > 0.0)                                   \
#               and                                                            \
#               (el['eta'][ii] < el['eta'][ii+1]):
#
#               fa['jumptype'][ii+1] = +1
#                   
#               # baseline simple extrapolation    
#               fa['etaM'][ii+1]  = el['eta'][ii] 
#               fa['areaM'][ii+1] = el['area'][ii]
#               fa['velocityM'][ii+1] = fa['flowrate'][ii+1] / fa['areaM'][ii+1]
#               
#               fa['etaP'][ii+1]  = el['eta'][ii+1]
#               fa['areaP'][ii+1] = el['area'][ii+1]
#               fa['velocityP'][ii+1] = fa['flowrate'][ii+1] / fa['areaP'][ii+1]
#               
#               # using depths on elements
#               depthU  = el['hyddepth'][ii]
#               froudeU = el['froude'][ii]
#               
#               depthD  = el['hyddepth'][ii+1]
#               froudeD = el['froude'][ii+1]
#
#               # approximation of the discontinuity at the jump    
#               if setting.method_hydjump_face == 'simple':
#                   # simple extrapolation of element values
#                   # as jump values       
##                   print()
##                   print(ii, \
##                         fa['flowrate'][ii+1]*fa['velocityM'][ii+1] \
##                             + setting.gravity * fa['areaM'][ii+1] * fa['etaM'][ii+1], \
##                         fa['flowrate'][ii+1]*fa['velocityP'][ii+1] \
##                             + setting.gravity * fa['areaP'][ii+1] * fa['etaP'][ii+1], \
##                         )
#                   
##                   fa['velocityP'][ii+1]  \
##                       = (( fa['flowrate'][ii+1]*fa['velocityM'][ii+1] \
##                          + setting.gravity * fa['areaM'][ii+1] * fa['etaM'][ii+1]) \
##                          - setting.gravity * fa['areaP'][ii+1] * fa['etaP'][ii+1]) \
##                          / fa['flowrate'][ii+1]
#
##                   print(ii, \
##                         fa['flowrate'][ii+1]*fa['velocityM'][ii+1] \
##                             + setting.gravity * fa['areaM'][ii+1] * fa['etaM'][ii+1], \
##                         fa['flowrate'][ii+1]*fa['velocityP'][ii+1] \
##                             + setting.gravity * fa['areaP'][ii+1] * fa['etaP'][ii+1], \
##                         )
#                             
#                   #print(ii, 'area: ', fa['areaP'][ii+1], fa['areaM'][ii+1] )
#                   #print('flow,vP,vM: ', fa['flowrate'][ii+1], fa['velocityP'][ii+1], fa['velocityM'][ii+1] )
#                   #print('momentum: ',fa['flowrate'][ii+1]*fa['velocityP'][ii+1], fa['flowrate'][ii+1]*fa['velocityM'][ii+1] )
#                  
#                   pass
#               
#               elif setting.method_hydjump_face == 'rectangular':
#                   # use conjugate depth relationships for rectangular channel
#                   slope = (fa['zbottom'][ii]- fa['zbottom'][ii+1])           \
#                           / geo['length'][ii]                  
#                   # normal flow rate
#                   if slope > 0.0:
#                       if (geo['manningsn'][ii] == 0.0):
#                           Qnorm = 1e12
#                       else:
#                           Qnorm = (1.0 / geo['manningsn'][ii])               \
#                                   * el['area'][ii]                           \
#                                   * el['hydradius'][ii]**(2.0/3.0)           \
#                                   * np.sqrt(slope)
#                       #endif
#                   else:
#                       Qnorm = 0
#                   #endif
#                   
#                   # Note that el['flowrate'] < Qnorm implies that we are above
#                   # the normal depth (i.e. decreasing depth is needed to get
#                   # a normal Q that meets Chezy-Manning). Hence for Q < Qnorm
#                   # we are on an S2 curve and expect an immediate jump. In
#                   # contrast, el['flowrate'] > Qnorm implies that increasing
#                   # the depth is needed to match C-M at the observed flow, so
#                   # we are on an S3, M3, C3, H3 or A3 curve. With a #3 curve
#                   # we use the downstream depth as the known jump depth and 
#                   # compute the upstream depth. With an S2 curve we use the
#                   # upstream depth as the known depth and compute the 
#                   # downstream depth.
#                   
#                   if el['flowrate'][ii] <= Qnorm:
#                       #upstream depth is above normal depth    
#                       depthD = depthU * 0.5 * (-1.0                            \
#                          + np.sqrt(1.0 + 8.0 *(froudeU**2.0)))
#                       fa['areaP'][ii+1] =  depthD * fa['topwidth'][ii+1]                                         
#                       fa['etaP'][ii+1] = fa['zbottom'][ii+1] + depthD
#                   else:
#                       depthU = depthD * 0.5 * (-1.0                            \
#                          + np.sqrt(1.0 + 8.0 *(froudeD**2.0)))
#                       fa['areaM'][ii+1] =  depthU * fa['topwidth'][ii+1]                                         
#                       fa['etaM'][ii+1]  = fa['zbottom'][ii+1] + depthU
#                   #endif
#               else:
#                   print('error, ',                                           \
#                         'unknown value for setting.method_hydjump_face of ', \
#                         setting.method_hydjump_face)
#                   sys.exit()
#               #endif
#                             
#            #endif
#            
#            #------------------------------------------------------------------      
#            # jump with flow in the nominal upstream direction
#            if  (el['froude'][ii-1] < 0.9) and                                \
#                (el['froude'][ii] > 1.1) and                                  \
#                (fa['flowrate'][ii] < 0.0)                                    \
#                 and                                                          \
#                (el['eta'][ii-1] > el['eta'][ii]):
#
#               fa['jumptype'][ii] = -1                    
#
#               # baseline simple extrapolation    
#               fa['etaM'][ii]  = el['eta'][ii-1] 
#               fa['areaM'][ii] = el['area'][ii-1]
#               fa['velocityM'][ii] = fa['flowrate'][ii] / fa['areaM'][ii]
#               
#               fa['etaP'][ii]  = el['eta'][ii]
#               fa['areaP'][ii] = el['area'][ii]
#               fa['velocityP'][ii] = fa['flowrate'][ii] / fa['areaP'][ii]
#
#               # using depths on elements with U and D as jump upstream
#               # and downstream
#               depthU  = el['hyddepth'][ii]
#               froudeU = el['froude'][ii]
#               
#               depthD  = el['hyddepth'][ii-1]
#               froudeD = el['froude'][ii-1]
#
#               # approximation of the discontinuity at the jump    
#               if setting.method_hydjump_face == 'simple':
#                   # simple extrapolation of element values
#                   # as jump values    
#                   pass
#               
#               elif setting.method_hydjump_face == 'rectangular':
#                   # use conjugate depth relationships for rectangular channel
#                   # note this slope is positive if bottom decreasing upstream
#                   slope = (fa['zbottom'][ii+1]- fa['zbottom'][ii])           \
#                           / geo['length'][ii]                  
#                   # normal flow rate
#                   if slope > 0.0:
#                       # give this Q a negative because el is < 0
#                       if (geo['manningsn'][ii] == 0.0):
#                           Qnorm = -1e12
#                       else:
#                           Qnorm = -(1.0 / geo['manningsn'][ii])               \
#                                   * el['area'][ii]                           \
#                                   * el['hydradius'][ii]**(2.0/3.0)           \
#                                   * np.sqrt(slope)
#                       #endif
#                   else:
#                       Qnorm = 0
#                   #endif
#                   
#                   # Note that el['flowrate'] < Qnorm implies that we are above
#                   # the normal depth (i.e. decreasing depth is needed to get
#                   # a normal Q that meets Chezy-Manning). Hence for Q < Qnorm
#                   # we are on an S2 curve and expect an immediate jump. In
#                   # contrast, el['flowrate'] > Qnorm implies that increasing
#                   # the depth is needed to match C-M at the observed flow, so
#                   # we are on an S3, M3, C3, H3 or A3 curve. With a #3 curve
#                   # we use the downstream depth as the known jump depth and 
#                   # compute the upstream depth. With an S2 curve we use the
#                   # upstream depth as the known depth and compute the 
#                   # downstream depth.
#                   
#                   if el['flowrate'][ii] >= Qnorm: # note that el < 0 here
#                       #upstream depth is above normal depth    
#                       depthD = depthU * 0.5 * (-1                                            \
#                          + np.sqrt(1.0 + 8.0 *(froudeU**2.0)))
#                       fa['areaM'][ii] =  depthD * fa['topwidth'][ii]                                         
#                       fa['etaM'][ii]   = fa['zbottom'][ii] + depthD
#                   else:
#                       depthU = depthD * 0.5 * (-1                                            \
#                          + np.sqrt(1.0 + 8.0 *(froudeD**2.0)))
#                       fa['areaP'][ii] =  depthU * fa['topwidth'][ii]                                         
#                       fa['etaP'][ii]  = fa['zbottom'][ii] + depthU
#                   #endif
#               else:
#                   print('error, ',                                           \
#                         'unknown value for setting.method_hydjump_face of ', \
#                         setting.method_hydjump_face)
#                   sys.exit()
#               #endif                                                                        
#            #endif
#            
#        #endfor
#          
#        return fa
##==============================================================================
#EOF