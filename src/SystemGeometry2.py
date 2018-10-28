#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 14:54:12 2018

@author: brh
"""

class SystemGeometry:
    
#==============================================================================    
    def geometry_case(self, geo, NX, setting):
        '''
        customize for different cases
        output should be entire geo structure
        '''
        import sys
        import numpy as np
               
        #----------------------------------------------------------------------                       
        if setting.geometry_case == 'flow_over_a_bump':
            import CaseFlowOverBump        
            cBump = CaseFlowOverBump.CaseFlowOverBump()
            [geo, setting] = cBump.define_geometry(geo, setting, NX) 
 
            # setup the roughness for geometry
            if setting.mannings_n_use_global_default:
                geo = self.roughness_initialization(geo, setting)   
                
            # use a manning's n to damp flow in buffer domain
            if geo['xvalue'][NX-1] > 25.0:
                aa = geo['xvalue'][:] > 27.0
                geo['manningsn'][aa] = setting.mannings_n_for_buffer
            #endif
            
            
            setting.steady_solution_exists = True
        
        #----------------------------------------------------------------------                       
        elif setting.geometry_case == 'Waller_Creek':
            import CaseWallerCreek
            wCreek = CaseWallerCreek.CaseWallerCreek()           
            [geo, setting] = wCreek.define_geometry(geo, setting, NX)
                       
            setting.steady_solution_exists = False
            
        #----------------------------------------------------------------------                       
        else:
            print('error: unknown setting.geometry_case of ', \
                  setting.geometry_case )
            sys.exit()
            
        #endif 
        
        return [geo, setting]
    

#==============================================================================
    def set_all_geometry(self, setting):
        '''
        Calls all the functions to set up the geometry of the system for
        a single reach with NX cells
        '''
        import numpy as np
        import sys

        import DataStorage2
        ds = DataStorage2.DataStorage()
        
        #----------------------------------------------------------------------                       
        # the desired number of cells (NX) of this reach based on user input
        thisNX = setting.NX

        #----------------------------------------------------------------------                       
        # get the maximum number of width-depth pairs
        setting.geometry_number_widthdepth_pairs                              \
            = self.get_number_of_widthdepth_pairs(setting)

        #----------------------------------------------------------------------                       
        # convert minimum angle to radians
        setting.angle_minimum = np.deg2rad(setting.angle_minimum)

        #----------------------------------------------------------------------                       
        # define the solution storage space
        [geo, el, fa, rkc, rke, rkf, setting]                                 \
            = ds.set_all_storage(thisNX, setting)

        #----------------------------------------------------------------------                       
        # get the geometry for the particular case (which may change NX!)
        [geo, setting] = self.geometry_case(geo, thisNX, setting)
        
        #----------------------------------------------------------------------                       
        #HACK reset NX if changed in case (adding buffer etc.)
        if thisNX != setting.NX:
            #HACK the NX was changed by the particular case, so reset storage
            thisNX = setting.NX
            [dummy, el, fa, rkc, rke, rkf, setting]                           \
                = ds.set_all_storage(thisNX, setting)
            del dummy
        #endif
        NX = thisNX
        
        #----------------------------------------------------------------------                       
        # compute face geometry
        [fa, rkf] = self.face_geometry_initialization                         \
                    (fa, rkf, geo, NX, setting) 
       
        #----------------------------------------------------------------------  
        # perform width-depth geometry computations                     
        if setting.geometry_number_widthdepth_pairs > 0:
            geo = self.widthdepth_pair_loop (geo, setting)
            
            # compute additional geometry data for widthdepth pairs
            [geo, setting] = self.widthdepth_pair_auxiliary(geo, setting)
        #endif
        
        #----------------------------------------------------------------------  
        # initialize the smallvolumes consistent with small depths
        geo = self.smallvolume_initialization(geo, setting)
                                    
        return [geo, el, fa, rkc, rke, rkf, setting]

#==============================================================================        
    def face_geometry_initialization(self, fa, rkf, geo, NX, setting):   
        '''
        provides initial values for face arrays from geometry array
        '''
        import sys
        import numpy as np
        
        fa = self.face_xvalue (fa, geo, NX)
                    
        fa = self.face_zbottom (fa, geo, NX)
        
        fa = self.face_cosangle (fa, geo, NX, setting)
        
        # for rkf
        if rkf != None:
            rkf['zbottom'][:] = fa['zbottom'][:]
            rkf['cosangle'][:] = fa['cosangle'][:]
        #endf
        
        return [fa, rkf]

#==============================================================================        
    def face_xvalue (self, fa, geo, NX):
        '''
        compute and store xvalues on faces as interpolation from geo lengths
        stored on centers
        '''
        import sys
        
        fa['xvalue'][0:NX] = geo['xvalue'][:] - 0.5 * geo['length'][:]
                
        fa['xvalue'][NX] = fa['xvalue'][NX-1] + geo['length'][NX-1]
        
        #error checking
        xdiff = fa['xvalue'][1:NX+1] - fa['xvalue'][0:NX]
        aa = xdiff <= 0.0
        
        if any(aa):
            print()
            print(xdiff[aa])
            print('error, some of the xvalues in the geometry are not monotonic')
            sys.exit()

        return fa

#==============================================================================        
    def face_zbottom (self, fa, geo, NX):
        '''
        Compute z bottom at faces as interpolation from zbottom from geo array
        '''
        fa['zbottom'][1:NX] \
            = (   geo['zbottom'][0:NX-1] * geo['length'][1:NX  ]              \
                + geo['zbottom'][1:NX  ] * geo['length'][0:NX-1] )            \
              / ( geo['length'] [1:NX  ] + geo['length'][0:NX-1])

                
        fa['zbottom'][0]  = geo['zbottom'][0]
        fa['zbottom'][NX] = geo['zbottom'][NX-1]
        
        return fa

#==============================================================================        
    def face_cosangle (self, fa, geo, NX, setting):
        '''
        Computes the cosine of slope angle across the face
        Defined with a downstream slope as positive
        '''
        import numpy as np
        
        if setting.method_use_cosangle:
        
            zdist = geo['zbottom'][1:NX] - geo['zbottom'][0:NX-1]
            xdist = geo['xvalue'][1:NX]  - geo['xvalue'][0:NX-1]
            sdist = np.sqrt(zdist**2.0 + xdist**2.0)
            fa['cosangle'][1:NX] = xdist / sdist
            
            fa['cosangle'][0] = 1.0
            fa['cosangle'][NX] = 1.0
        else:
            fa['cosangle'][:] = 1.0
        #endif
               
        return fa
#==============================================================================
    def get_number_of_cells(self,setting):
        '''
        Determine the target number of cells in a reach, which depends on case.
        The goal is allow for subdivision of input reaches.
        '''      
        import sys
        
        if setting.geometry_case == 'flow_over_a_bump':
            import CaseFlowOverBump        
            cBump = CaseFlowOverBump.CaseFlowOverBump()            
            targetNX = cBump.get_number_of_cells(setting)
            
        elif setting.geometry_case == 'Waller_Creek':
            import CaseWallerCreek
            wCreek = CaseWallerCreek.CaseWallerCreek()
            targetNX = wCreek.get_number_of_cells(setting)
            
        else:
            print('error: unknown setting.geometry_case of ', \
                  setting.geometry_case )
            sys.exit()
        #endif setting.geometry_case == 'flow_over_a_bump'
        
        return targetNX

#==============================================================================    
    def get_number_of_widthdepth_pairs(self, setting):
        '''
        determine the maximum number of widthdepth pairs in a reach
        which is needed to set geometry array
        '''
        import sys
        
        if setting.geometry_case == 'flow_over_a_bump':
            import CaseFlowOverBump        
            cBump = CaseFlowOverBump.CaseFlowOverBump()            
            npair = cBump.get_number_of_widthdepth_pairs(setting)
            
        elif setting.geometry_case == 'Waller_Creek':
            import CaseWallerCreek
            wCreek = CaseWallerCreek.CaseWallerCreek()
            npair = wCreek.get_number_of_widthdepth_pairs(setting)
            
        else:
            print('error: unknown setting.geometry_case of ', \
                  setting.geometry_case )
            sys.exit()
            
        #endif 
        
        return npair

#==============================================================================
    def roughness_initialization(self, geo, setting):   
        '''
        stores the global mannings n across the entire reach
        '''
        # manning's n
        geo['manningsn'][:] = setting.mannings_n_global_default
        
        return geo

#==============================================================================
    def smallvolume_initialization(self, geo, setting):
        '''
        Defines the small volume consistent with the setting.depth_small_value
        value for each element.
        '''
        import numpy as np
        import sys
        
        #----------------------------------------------------------------------  
        # aliases for handling width-depth pairs
        if setting.geometry_number_widthdepth_pairs > 0:
            # aliases for indexes
            widthAtLayerTop                                                   \
                = setting.geometry_widthdepth_values['widthAtLayerTop']
            depthAtLayerTop                                                   \
                = setting.geometry_widthdepth_values['depthAtLayerTop']
            #areaThisLayer = setting.geometry_widthdepth_values['areaThisLayer'] 
            
            # note this does not include the layer defined by depth and width
            areaTotalBelowThisLayer                                           \
                = setting.geometry_widthdepth_values['areaTotalBelowThisLayer']
            Dwidth = setting.geometry_widthdepth_values['Dwidth']
            Ddepth = setting.geometry_widthdepth_values['Ddepth']
            #angle =  setting.geometry_widthdepth_values['angle']
            #perimeterBelowThisLayer                                           \ 
            #    = setting.geometry_widthdepth_values['perimeterBelowThisLayer']
        #endif
        
        #----------------------------------------------------------------------  
        # Masks for different types of geometry
        rr = geo['etype'][:] == 'rectangular_channel'
        pp = geo['etype'][:] == 'parabolic_channel'
        tt = geo['etype'][:] == 'trapezoidal_channel'
        wd = geo['etype'][:] == 'widthdepth_pair'

        #----------------------------------------------------------------------  
        # initialization
        geo['smallvolume'][:] = 0.0
        
        #----------------------------------------------------------------------  
        # Small Volumes for various channel types
        # for rectangular channel
        geo['smallvolume'][rr] = setting.depth_small_value                    \
                             * geo['breadth'][rr] * geo['length'][rr]
        #----------------------------------------------------------------------  
        #for parabolic channel
        geo['smallvolume'][pp] = geo['length'][pp] * (4.0/3.0)                \
                * np.sqrt( setting.depth_small_value**3.0                     \
                         / geo['parabola_value'][pp] ) 
                
        #----------------------------------------------------------------------  
        # for trapezoidal channel
        geo['smallvolume'][tt] = geo['length'][tt] * (                        \
                            setting.depth_small_value * geo['breadth'][tt]    \
                          +(setting.depth_small_value**2.0 )                  \
                          / np.tan(geo['trapezoid_angle'][tt]) )

        #----------------------------------------------------------------------  
        # for widthpair channel based on the two lowest widthdepth pairs
        # HACK needs to be cleaned up and compressed
        if setting.geometry_number_widthdepth_pairs > 0:
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            for ii in range(0,setting.NX):
                #print(ii, setting.depth_small_value, geo['npair'][ii],geo['etype'][ii])
                # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                if geo['etype'][ii] == 'widthdepth_pair':
                    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
                    # if the small depth is in the lowest widthdepth pair
                    if setting.depth_small_value                              \
                        <= geo['widthdepth'][ii,0,depthAtLayerTop]:
                            
                        # compute a small cross-sectional area    
                        geo['smallvolume'][ii]                                \
                            = 0.5 * (setting.depth_small_value**2)            \
                            * geo['widthdepth'][ii,0,widthAtLayerTop]         \
                            / geo['widthdepth'][ii,0,depthAtLayerTop]
                    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
                    else:
                        #cycle through the width-depth apirs
                        for kk in range(1,geo['npair'][ii]):
                            #print('       ',kk, geo['widthdepth'][ii,kk,depthAtLayerTop] )
                            if    (setting.depth_small_value <= geo['widthdepth'][ii,kk  ,depthAtLayerTop]) \
                                & (setting.depth_small_value >  geo['widthdepth'][ii,kk-1,depthAtLayerTop]): 
                                    
                                depthinc = setting.depth_small_value - geo['widthdepth'][ii,kk-1,depthAtLayerTop]   

                                # compute a small cross-sectional area
                                geo['smallvolume'][ii]                                         \
                                    = geo['widthdepth'][ii,kk,areaTotalBelowThisLayer]         \
                                    + depthinc * geo['widthdepth'][ii,kk-1,widthAtLayerTop]    \
                                    +(depthinc**2.0 )                                          \
                                    * (geo['widthdepth'][ii,kk,Dwidth])                        \
                                    / (geo['widthdepth'][ii,kk,Ddepth]) 
                                #print(geo['smallvolume'][ii], geo['widthdepth'][ii,kk,areaBelowThisLayer], geo['widthdepth'][ii,kk+1,areaBelowThisLayer])    
                                #/print('here ',kk,  geo['smallvolume'][ii])    
                                #print( geo['widthdepth'][ii,kk-1,areaTotalBelow])    
                                break #kk
                            else:
                                pass
                            #endif
                        #endfor kk    
                    #endif
                    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
                    # the "volume" above is actually an area!
                    geo['smallvolume'][ii] = geo['smallvolume'][ii] \
                                           * geo['length'][ii]
                        
                    if geo['smallvolume'][ii] == 0.0:
                        print('error no small volume found in ',ii)
                        sys.exit()
                    #endif
                #endif
                # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            #endfor
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        #endif
        
        return geo

#==============================================================================
    def widthdepth_pair_loop (self, geo, setting):
        '''
        loops through the width-depth pairs to checkc consistency
        '''
        import sys
        
        # check widthdepth pair geometry for consistency
        nfix = 1
        mm = 0
        while nfix > 0:
            [geo, nfix] = self.widthdepth_pair_consistency(geo, setting)
            mm=mm+1
            #print()
            #print(mm,' time', nfix)      
            #print()

            if (nfix > 0) and (mm > 1000):
                print('widthdepth pair inconsistency appears stuck in loop')
                sys.exit()
            #endif
        #endwhile
        
        return geo
   
#==============================================================================
    def widthdepth_pair_consistency(self, geo, setting):
        '''
        Check to see that widthdepth pairs have width and depth both increasing
        with the vertical iteration variable (2nd index).
        Calls for widthdepth_pair_fix if there is a problem
        '''
        import sys
        import numpy as np
        
        NX = setting.NX
                
        aa = geo['etype'][:] == 'widthdepth_pair'
        
        width = setting.geometry_widthdepth_values['widthAtLayerTop']
        depth = setting.geometry_widthdepth_values['depthAtLayerTop']

        #----------------------------------------------------------------------                               
        if any(aa):    
            # Note that we must loop over all cells (slow) because the 
            # number of valid widthdepth pairs may be different in each cell.
            # Thus, the invalid pair storage space has zeros and will trip
            # a test for consistency
            nfix = 0
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            for ii in range(0,NX):
                # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                if geo['etype'][ii] == 'widthdepth_pair':
                    # the number of pairs at this element
                    npair = geo['npair'][ii]
                    # compute the difference in width across each level
                    dWidth = geo['widthdepth'][ii,1:npair  ,width] \
                           - geo['widthdepth'][ii,0:npair-1,width]
                    # compute difference in depth across eacg level
                    dDepth = geo['widthdepth'][ii,1:npair  ,depth] \
                           - geo['widthdepth'][ii,0:npair-1,depth]
                    
                    # negative values indicate non-monotonic behavior that
                    # can be fixed.
                    aa = (dWidth < 0.0) | (dDepth < 0.0)
                    if any(aa):
                        nfix = nfix + np.count_nonzero(aa)
                        if setting.geometry_adjust:
                            geo = self.widthdepth_pair_fix(setting, geo, ii)
                        #endif setting.geometry_adjust
                    #endif any(aa) 
                
                    # check that the width-depth pairs cover enough depth
                    # and fix with vertical walls
                    if max(geo['widthdepth'][ii,:,depth]) < setting.depth_maximum_expected:
                        if setting.inoisy:
                            print('insufficient depth in widthdepth pair at cell ID',geo['ID'][ii])
                            print('increasing depth available by vertical walls')
                            #print(geo['widthdepth'][ii,:,depth])
                            #print(geo['npair'][ii])
                            #print(geo['widthdepth'][ii,geo['npair'][ii-1],depth])
                        #endif
            
                        geo['widthdepth'][ii,geo['npair'][ii],depth]  \
                            = 2.0 * setting.depth_maximum_expected
                
                        geo['widthdepth'][ii,geo['npair'][ii],width] \
                            = geo['widthdepth'][ii,geo['npair'][ii-1],width]
                
                        # an additional pair has been added at this element
                        geo['npair'][ii] = geo['npair'][ii] + 1           
                    #endif  
                    
                #endif geo['etype'][ii] == 'widthdepth_pair'
                # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            #endfor
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        #----------------------------------------------------------------------                       
        else:
            print('error: no widthdepth channel sections found, but ',
                  'setting.geometry_number_widthdepth_pairs > 0. ',
                  'Set to 0 to reduce memory use.')
            sys.exit()
        #endif any(aa)
        #----------------------------------------------------------------------                       
       
        return [geo, nfix]
    
#==============================================================================
    def widthdepth_pair_fix(self, setting, geo, ii):
        '''
        Fix inconsistency in widthdepth pairs for ii cell;
        i.e. when upper width or depth is smaller than width or depth below
        Note that inconsistency is upper width < lower width, or
        upper depth <= lower depth. That is, the depth must always increase
        (cannot stay the same), but the width may either stay the same or 
        increase.
        '''
        import sys
        
        # max number of pairs for this cells
        npair = geo['npair'][ii]      
        # alias for geo index
        width = setting.geometry_widthdepth_values['widthAtLayerTop']
        depth = setting.geometry_widthdepth_values['depthAtLayerTop']
        
        #print('fixing geometry here for cell ',ii)
        
        changed_this_depth = False
        changed_this_width = False
        
        #----------------------------------------------------------------------                       
        # cycle through the pairs in the ii cell
        for kk in range(0,npair-1):
            # width 2 layers above
            up2W = geo['widthdepth'][ii,kk+2,width]
            # width 1 layer above
            up1W = geo['widthdepth'][ii,kk+1,width]
            # width this layer (below)
            lowW = geo['widthdepth'][ii,kk  ,width]
            
            # depths named similar to above
            up2D = geo['widthdepth'][ii,kk+2,depth]
            up1D = geo['widthdepth'][ii,kk+1,depth]
            lowD = geo['widthdepth'][ii,kk  ,depth]
            
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Fix inconsistent depth above
            if up1D <= lowD:
                changed_this_depth = True
                if up2D <= lowD:
                    # multi-level inconsistency - simple expansion
                    geo['widthdepth'][ii,kk+1,depth] = lowD                   \
                        * (1.0 + setting.geometry_adjust_fraction)
                else:
                    # single-level inconsistency
                    if (up1W > lowW) and (up2W > up1W):
                        # width is consistent - use linear interpolation
                        geo['widthdepth'][ii,kk+1,depth] = lowD               \
                            + (up2D-lowD) * (up1W-lowW) / (up2W-lowW)
                    else:
                        # both depth and width are inconsistent - simple
                        # expansion
                        geo['widthdepth'][ii,kk+1,depth] = lowD               \
                            * (1.0 + setting.geometry_adjust_fraction)
                    #endif (up1W > lowW) and (up2W > up1W)
                #endif up2D <= lowD     
            #endif up1D <= lowD    
            
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Fix inconsistent widht above
            if up1W < lowW:
                changed_this_width = True
                if up2W < lowW:
                    # multi-level inconsistency - simple expansion
                    geo['widthdepth'][ii,kk+1,width] = lowW                   \
                        * (1.0 + setting.geometry_adjust_fraction)
                else:
                    # single-level inconsistency
                    if (up1D > lowD) and (up2D > up1D):
                        # depth is consistent - use linear interpolation
                        geo['widthdepth'][ii,kk+1,width] = lowW               \
                            + (up2W-lowW) * (up1D-lowD) / (up2D-lowD)
                    else:
                        # both depth and width are inconsistent - simple
                        # expansion
                        geo['widthdepth'][ii,kk+1,width] = lowW               \
                            * (1.0 + setting.geometry_adjust_fraction)
                    #endif (up1D > lowD) and (up2D > up1D)
                #endif up2W < lowW     
            #endif up1W < lowW   
            
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # printouts
            if setting.inoisy:
                if changed_this_width:                   
                    print('NOTE: widthdepth geometry inconsistent in cell ',  \
                          ii,' for pair ',kk)
                    print('      adjusted width from ',up1W,' to ',           \
                          geo['widthdepth'][ii,kk+1,width])  
                    changed_this_width = False
                #endif changed_this_width    
                
                if changed_this_depth:                            
                    print('NOTE: widthdepth geometry inconsistent in cell ',  \
                          ii,' for pair ',kk)
                    print('      adjusted depth from ',up1D,' to ',           \
                          geo['widthdepth'][ii,kk+1,depth])  
                    changed_this_depth = False
                #endif changed_this_depth
                
            #endif setting.inoisy
        #endfor kk
            
        return geo
    
#==============================================================================
    def widthdepth_pair_auxiliary(self, geo, setting):
        '''
        Defines the auxiliary values between two width-depth pair sets.  
        Note that the point depth = 0, width = 0 is an implicit value
        point. Thus, the depth at level ii is the depth to top of that
        level and the width is to the top of that level.
        '''
        import sys
        import numpy as np
        
        #----------------------------------------------------------------------                       
        # aliases for indexes
        widthAtLayerTop = setting.geometry_widthdepth_values['widthAtLayerTop']
        depthAtLayerTop = setting.geometry_widthdepth_values['depthAtLayerTop']
        areaThisLayer = setting.geometry_widthdepth_values['areaThisLayer']         
        areaTotalBelowThisLayer                                               \
            = setting.geometry_widthdepth_values['areaTotalBelowThisLayer']
        Dwidth = setting.geometry_widthdepth_values['Dwidth']
        Ddepth = setting.geometry_widthdepth_values['Ddepth']
        angle =  setting.geometry_widthdepth_values['angle']
        perimeterBelowThisLayer \
            = setting.geometry_widthdepth_values['perimeterBelowThisLayer']

        eIn1 = geo['widthdepth'].shape[1]

        wp = geo['etype'][:] == 'widthdepth_pair'
        
        #----------------------------------------------------------------------                       
        # ensure aux data are all zero other than width-depth values
        geo['widthdepth'][:,:,2:8] = 0.0    
        
        #----------------------------------------------------------------------                       
        # lowest layer is triangular
        geo['widthdepth'][wp,0,areaThisLayer] = 0.5 * \
                geo['widthdepth'][wp,0,widthAtLayerTop] \
             *  geo['widthdepth'][wp,0,depthAtLayerTop]
               
        # the area in this layer   
        geo['widthdepth'][wp,1:eIn1,areaThisLayer] = \
              + 0.5 *(  geo['widthdepth'][wp,1:eIn1,  widthAtLayerTop]     \
                      + geo['widthdepth'][wp,0:eIn1-1,widthAtLayerTop])    \
                    *(  geo['widthdepth'][wp,1:eIn1,  depthAtLayerTop]     \
                      - geo['widthdepth'][wp,0:eIn1-1,depthAtLayerTop])

        #----------------------------------------------------------------------                       
        # set areas to zero above the uppermost pair
        for ii in range(0,geo['etype'].shape[0]):
            if geo['etype'][ii] == 'widthdepth_pair':
                geo['widthdepth'][ii,geo['npair'][ii]:eIn1,areaThisLayer] = 0.0
            #endif
        #endfor                 
        
        geo['widthdepth'][wp,0,Dwidth] =geo['widthdepth'][wp,0,widthAtLayerTop]          
        geo['widthdepth'][wp,0,Ddepth] =geo['widthdepth'][wp,0,depthAtLayerTop]          

        #----------------------------------------------------------------------       
        # cycle over the valid layers and store width and depth differences                
        for ii in range(1,geo['widthdepth'].shape[1]-1):
            # delta width between top and bottom of this layer
            geo['widthdepth'][wp,ii,Dwidth]                                   \
                = geo['widthdepth'][wp,ii  ,widthAtLayerTop]                  \
                - geo['widthdepth'][wp,ii-1,widthAtLayerTop]

            # delta depth between top and bottom of this layer
            geo['widthdepth'][wp,ii,Ddepth]                                   \
                = geo['widthdepth'][wp,ii  ,depthAtLayerTop]                  \
                - geo['widthdepth'][wp,ii-1,depthAtLayerTop] 
        #endfor

        #----------------------------------------------------------------------  
        # Trapezoidal area computations        
          
        # angle of trapezoid - first for all Dwidth > 0
        # mask for the width-depth pair cells
        wpWD      = geo['auxWD2'][:,:]
        wpWD[:,:] = 0.0
        wpWD[wp,:] = 1.0
        # pairs that are not 90 degree angles
        aa = (geo['widthdepth'][:,:,Dwidth] > setting.geometry_small_width)  \
            & (wpWD[:,:] > 0.0)
        geo['widthdepth'][aa,angle] = np.arctan(                              \
           2.0 * geo['widthdepth'][aa,Ddepth] / geo['widthdepth'][aa,Dwidth] )

        # A near-90 degree angle will have an infinite tangent. 
        # We handle this case by setting all these angles to pi/2 - small value           
        aa = (geo['widthdepth'][:,:,Dwidth] <= setting.geometry_small_width) \
            & (wpWD[:,:] > 0.0)
        geo['widthdepth'][aa,angle] = np.pi/2.0 - setting.angle_minimum 
 
        # accumulated area of the trapezoids to the ii width-depth level
        geo['widthdepth'][wp,0,areaTotalBelowThisLayer] = 0.0
        for kk in range(1,eIn1):
            geo['widthdepth'][wp,kk,areaTotalBelowThisLayer] = \
                  geo['widthdepth'][wp,kk-1,areaTotalBelowThisLayer]   \
                + geo['widthdepth'][wp,kk-1,areaThisLayer]
        #endfor
                           
        #----------------------------------------------------------------------       
        #check that the setting.area_maximum value is greater than any 
        # accumulated area at the uppermost level
        thismax = np.amax(geo['widthdepth'][wp,eIn1-1,areaTotalBelowThisLayer])
        if setting.area_maximum < thismax:
            setting.area_maximum = 2.0*thismax
        #endif
        
        #----------------------------------------------------------------------       
        # perimeter below this layer
        geo['widthdepth'][wp,0, perimeterBelowThisLayer] = 0.0 
        for ii in range(1,eIn1):
            geo['widthdepth'][wp,ii, perimeterBelowThisLayer] =               \
                geo['widthdepth'][wp,ii-1, perimeterBelowThisLayer]           \
                + 2.0 * np.sqrt(      geo['widthdepth'][wp,ii-1,Ddepth]**2.0  \
                                +(0.5*geo['widthdepth'][wp,ii-1,Dwidth])**2.0 )
        #endfor
        
        return [geo, setting]
