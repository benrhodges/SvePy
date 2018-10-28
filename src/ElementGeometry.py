#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 08:11:05 2018

@author: brh
"""

class ElementGeometry:
    
#==============================================================================
    def element_geometry(self, el, geo, setting):
        '''
        Compute auxiliary values on element center
        Assumes that volume, length, zbottom, breadth are defined
        '''             
        import sys
        
        NX = setting.NX

        # updated values are all zeroed to make sure we don't miss anything        
        el['eta'][:] = 0.0
        el['topwidth'][:] = 0.0
        el['hyddepth'][:] = 0.0
        el['perimeter'][:] = 0.0
        el['hydradius'][:] = 0.0
        el['issmallvolume'][:] = False
        
        # Boolean Masks for different types of geometry
        rr = geo['etype'][:] == 'rectangular_channel'
        pp = geo['etype'][:] == 'parabolic_channel'
        tt = geo['etype'][:] == 'trapezoidal_channel'
        wd = geo['etype'][:] == 'widthdepth_pair'
                   
        #----------------------------------------------------------------------                       
        # handle small volumes
        el = self.smallvolume_treatment(el, geo, setting)      
        
        #----------------------------------------------------------------------                       
        # area
        el = self.area_computation(el, geo, setting)    
        
        #----------------------------------------------------------------------                       
        # other geometry
        if sum(rr) > 0:
            el = self.rectangular_channel_geometry(el, geo, setting, rr)
        #endif
        if sum(pp) > 0:
            el = self.parabolic_channel_geometry(el, geo, setting, pp)
        #endif
        if sum(tt) > 0:
            el = self.trapezoidal_channel_geometry(el, geo, setting, tt)
        #endif    
        if setting.geometry_number_widthdepth_pairs > 0:
            el = self.widthdepth_pair_channel_geometry                        \
                (el, geo, setting, wd, NX)    
        #endif
        
        #----------------------------------------------------------------------                       
        # reset for small or negative values 
        aa = el['eta'][:] - geo['zbottom'][:] <= setting.depth_zero_value 
        el['eta'][aa] = geo['zbottom'][aa] + setting.depth_zero_value
        aa[:] = False

        # reset for small or negative values
        aa = el['topwidth'][:] <= setting.topwidth_zero_value 
        el['topwidth'][aa] = setting.topwidth_zero_value
        aa[:] = False

        # reset for small or negative values 
        aa = el['perimeter'][:] <= setting.topwidth_zero_value 
        el['perimeter'][aa] = setting.topwidth_zero_value
        aa[:] = False
                  
        #----------------------------------------------------------------------                       
        # hydraulic depth       
        el['hyddepth'][:]   = el['area'][:] / el['topwidth'][:]

        # reset for small or negative values 
        aa = el['hyddepth'][:] <= setting.depth_zero_value 
        el['hyddepth'][aa] = setting.depth_zero_value
        aa[:] = False

        #----------------------------------------------------------------------                       
        # hydraulic radius 
        el['hydradius'][:] = el['area'][:] / el['perimeter'][:]
 
       # reset for small or negative values 
        aa = el['hydradius'][:] <= setting.topwidth_zero_value 
        el['hydradius'][aa] = setting.topwidth_zero_value
        aa[:] = False

        # HACK handle near-zero volumes 
        if any(el['issmallvolume']):
            aa = el['volume'][:] <= setting.volume_zero_value
            el['area'][aa] = setting.area_zero_value
            el['eta'][aa] = geo['zbottom'][aa] + setting.depth_zero_value
            el['topwidth'][aa] = setting.topwidth_zero_value
            el['hyddepth'][aa] = setting.depth_zero_value
            el['perimeter'][aa] = setting.topwidth_zero_value
            el['hydradius'][aa] = setting.topwidth_zero_value
            aa[:] = False
        #endif
        
        return el
         
#==============================================================================
    def smallvolume_treatment(self, el, geo, setting):
        
        # zero out negative volumes
        # Note this is a source of volume loss
        aa = el['volume'][:] < 0.0
        el['volume'][aa] = 0.0
        aa[:] = False
        
        # Compute smallvolume_ratio where the local volume is small
        # The ratio fraction used as a blending function in small volumes
        # and is 1.0 where the element volume is exactly the small volume
        aa = el['volume'][:] < geo['smallvolume'][:]
        el['issmallvolume'][aa] = True
        el['smallvolume_ratio'][aa] = el['volume'][aa] / geo['smallvolume'][aa]
        aa[:] = False
              
        # set the ratio to zero for volumes below the user zero value
        aa = el['volume'][:] <= setting.volume_zero_value
        el['smallvolume_ratio'][aa] = 0.0
        aa[:] = False

        return el
#==============================================================================
    def area_computation(self, el, geo, setting):        
            # AREA
        el['area'][:]   = el['volume'][:] / geo['length'][:]
        
        # reset for small or negative values ---------------------------------
        aa = el['area'][:] <= setting.area_zero_value 
        el['area'][aa] = setting.area_zero_value
        aa[:] = False

        return el
    
#==============================================================================
    def rectangular_channel_geometry(self, el, geo, setting, rr):

        # free surface
        el['eta'][rr] = geo['zbottom'][rr] \
                        + el['area'][rr] / geo['breadth'][rr]

        # top width
        el['topwidth'][rr] = geo['breadth'][rr] 

        # wetted perimeter
        el['perimeter'][rr] = geo['breadth'][rr]                              \
            + 2.0 * (el['eta'][rr] - geo['zbottom'][rr])
        
        return el
#==============================================================================
    def parabolic_channel_geometry(self, el, geo, setting, pp):
        '''
        geometry for a parabolic channel cross-section
        '''
        import numpy as np
        
        # free surface
        el['eta'][pp] = geo['zbottom'][pp] \
            + geo['parabola_value'][pp]**(1/3) * (0.75 * el['area'][pp])**(2/3)
            
        # top width
        el['topwidth'][pp] = 2.0 * np.sqrt(                                   \
                       (el['eta'][pp] - geo['zbottom'][pp]) / geo['gA'][pp] )
       
        # wetted perimeter
        el['perimeter'][pp] = (2.0 / (3.0 * geo['parabola_value'][pp]) )      \
              *(    (1 + geo['parabola_value'][pp]                            \
                        * (el['eta'][pp] - geo['zbottom'][pp])                \
                     )**(1.5)  - 1.0 )
              
        return el
    
#==============================================================================
    def trapezoidal_channel_geometry(self, el, geo, setting, tt):
        '''
        Quadratic solution for geometry from cross-sectional area in a
        trapezoidal channel
        '''
        import numpy as np
        
        trapz_tanTheta =  np.tan(geo['trapezoid_angle'][tt])
        
        CC = - el['area'][tt] * trapz_tanTheta  
        
        BB = geo['breadth'][tt] * trapz_tanTheta
        
        trapz_depth = -0.5 * (BB - np.sqrt(BB**2.0 - 4.0 * CC))
        
        # free surface
        el['eta'][tt] = trapz_depth + geo['zbottom'][tt]
        
        # top width
        el['topwidth'][tt] = geo['breadth'][tt]                               \
                            + 2.0 * trapz_depth  / trapz_tanTheta
       
        # wetted perimeter
        el['perimeter'][tt] = geo['breadth'][tt]                              \
            + 2.0 * trapz_depth / np.sin(geo['trapezoid_angle'][tt])
                            
        return el

#==============================================================================
    def widthdepth_pair_channel_geometry(self, el, geo, setting, wd, NX):
        '''
        Computes all geometry from cross-section area for a channel
        cross-section defined by width-depth pairs.
        '''
        import numpy as np
                
        # aliases for indexes
        widthAtLayerTop                                                       \
            = setting.geometry_widthdepth_values['widthAtLayerTop']
        depthAtLayerTop                                                       \
            = setting.geometry_widthdepth_values['depthAtLayerTop']
        areaThisLayer = setting.geometry_widthdepth_values['areaThisLayer'] 
        areaTotalBelowThisLayer                                               \
            = setting.geometry_widthdepth_values['areaTotalBelowThisLayer']
        Dwidth = setting.geometry_widthdepth_values['Dwidth']
        Ddepth = setting.geometry_widthdepth_values['Ddepth']
        angle =  setting.geometry_widthdepth_values['angle']
        perimeterBelowThisLayer                                               \
            = setting.geometry_widthdepth_values['perimeterBelowThisLayer']

        # HACK - this needs to be reworked so that it is array-processed
        # rather than iterated 
               
        for ii in range(0,NX):
            if wd[ii]:
                # difference between the cumulative area at each level
                # and the present area
                area_difference = el['area'][ii]                              \
                    - geo['widthdepth'][ii,:,areaTotalBelowThisLayer]
                    
                # difference btwn the remaining area and the area of each layer
                local_difference = area_difference                            \
                    - geo['widthdepth'][ii,:,areaThisLayer]
                    
                # the index is where the area difference changes sign
                ind = np.sign(area_difference) * np.sign(local_difference) < 0
                    
                # Coefficients for a quadratic solution
                AA = 1.0 / np.tan(geo['widthdepth'][ii,ind,angle])   
                BB = geo['widthdepth'][ii,ind,widthAtLayerTop]                \
                    - geo['widthdepth'][ii,ind,Dwidth]           
                CC = - area_difference[ind]      
                # quadratic solution
                DD = (-BB + np.sqrt(BB**2.0 - 4.0*AA*CC)) / (2.0 * AA)
                
                # total elevation from the depth                   
                el['eta'][ii] = geo['zbottom'][ii] + DD \
                    + geo['widthdepth'][ii,ind,depthAtLayerTop] \
                    - geo['widthdepth'][ii,ind,Ddepth]

                el['topwidth'][ii] = geo['widthdepth'][ii,ind,widthAtLayerTop]\
                    - (geo['widthdepth'][ii,ind,Ddepth] - DD) \
                    *  geo['widthdepth'][ii,ind,Dwidth]  \
                    /  geo['widthdepth'][ii,ind,Ddepth]
                
                el['perimeter'][ii]                                           \
                    = geo['widthdepth'][ii,ind,perimeterBelowThisLayer] \
                   + 2.0 * DD / np.sin(geo['widthdepth'][ii,ind,angle])

            #endif
        #endfor
        
# THE FOLLOWING NEEDS TO BE CHECKED TO BE SURE IT WORKS AS AN ARRAY BROADCAST
# REPLACEMENT FOR THE LOOP ABOVE            
#            geo['auxWD1'][:,:] = 0.0
#            #eIn1 = geo['widthdepth'].shape[1]
#            
#            # store the duplicates of the present area in the same shape
#            # as the geo[widthdepth] array for array processing
#            for ii in range(0,geo['widthdepth'].shape[1]):
#                geo['auxWD1'][wd,ii] = el['area'][wd]
#            # difference between the cumulative area at each level
#            # and the present volume
#            cumulative_area_difference =  \
#                      geo['auxWD1'][wd,:] - geo['widthdepth'][wd,:,areaTotalBelow]
#                      
#            # negative values are empty trapezoids.
#            aa = cumulative_area_difference <= 0.0
#            cumulative_area_difference[aa] = setting.area_maximum
#                        
#            # ind is index that provides the cells in 2D width-depth array where
#            # the free surface resides
#            ind = cumulative_area_difference - geo['widthdepth'][wd,:,areaUp] <= 0.0
#            
#            # quadratic to find the depth in this layer  
#            AA = 1.0 / np.tan(geo['widthdepth'][ind,angle])   
#            BB = geo['widthdepth'][ind,widthAtLayerTop]           
#            CC = - cumulative_area_difference[ind]        
#            DD = (-BB + np.sqrt(BB**2.0 - 4.0*AA*CC)) / (2.0 * AA)
#            # total elevation from the depth
#            el['eta'][wd] = geo['zbottom'][wd] + geo['widthdepth'][ind,depthAtLayerTop] + DD
#
#            geo['auxWD1'][:,:] = 0.0
        
        return el        

#==============================================================================
    def get_volume_from_elevation(self, el, geo, setting, NX):
        '''
        Compute an element volume from elevation.
        Typically used as part of initial condition computation
        '''        
        import numpy as np
        
        # Masks for different types of geometry
        rr = geo['etype'][:] == 'rectangular_channel'
        pp = geo['etype'][:] == 'parabolic_channel'
        tt = geo['etype'][:] == 'trapezoidal_channel'
               
        # rectangular channel
        el['volume'][rr] = (el['eta'][rr] - geo['zbottom'][rr])               \
                                * geo['length'][rr]
    
        # parabolic   channel 
        el['volume'][pp] = (4.0 / 3.0) * geo['length'][pp]                    \
            / np.sqrt( (el['eta'][pp] - geo['zbottom'][pp])**3.0              \
                      / geo['parabola_value'][pp])
        # trapezoidal channel
        el['volume'][tt] = (geo['length'][tt] / geo['trapezoid_angle'][tt] )  \
                            * (el['eta'][tt] - geo['zbottom'][tt])**2.0       \
                         + geo['breadth'][tt] * geo['length'][tt]             \
                            * (el['eta'][tt] - geo['zbottom'][tt])
        
        # width-depth indexed channel
        if setting.geometry_number_widthdepth_pairs > 0:
            el = self.widthdepth_volume_from_elevation(el, geo, setting, NX)  
        #endif
        
        return el
    
#==============================================================================
    def widthdepth_volume_from_elevation(self, el, geo, setting, NX):
        '''
        computes the volume from a given elevation for a width-depth pair
        cross-section. Typically used in initial conditions.
        '''
        import sys
        
        wd = geo['etype'][:] == 'widthdepth_pair'
        
        #----------------------------------------------------------------------                       
        # aliases for indexes
        width  = setting.geometry_widthdepth_values['widthAtLayerTop']
        depth  = setting.geometry_widthdepth_values['depthAtLayerTop']
        areaTotalBelow                                                        \
            = setting.geometry_widthdepth_values['areaTotalBelowThisLayer']
        Dwidth = setting.geometry_widthdepth_values['Dwidth']
        Ddepth = setting.geometry_widthdepth_values['Ddepth']
                    
        #----------------------------------------------------------------------  
        # define the depth
        thisdepth = el['eta'][wd] - geo['zbottom'][wd]
        
        #----------------------------------------------------------------------  
        # error checking
        aa = thisdepth <= setting.depth_zero_value
        if any(aa):
            print()
            print(thisdepth)
            print()
            print('error: initial conditions has depth ',                     \
                  '< setting.depth_zero_value ')
            sys.exit()
        #endif any(aa)
        
        #----------------------------------------------------------------------  
        # HACK if this is going to be used other than initial conditions
        # it needs to be rewritten as array-processed            
        for ii in range(0,NX):
            if wd[ii] == True:
                thisarea = 0.0
                tdepth = el['eta'][ii] - geo['zbottom'][ii]
                # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                # check if the elevation is in the lowest pair
                if tdepth <= geo['widthdepth'][ii,0,depth]:
                    
                    thisarea = 0.5 * geo['widthdepth'][ii,0,width]            \
                        * (tdepth / geo['widthdepth'][ii,0,depth])            \
                        * tdepth

                # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                else:
                    # Cycle through upper pairs to find the free surface.
                    # Note that we don't have to compute areas of each of the
                    # width-depth pairs because they've been stored already.
                    for kk in range(1,geo['npair'][ii]):
                        
                        if     (geo['widthdepth'][ii,kk-1,depth] < tdepth)    \
                             & (geo['widthdepth'][ii,kk  ,depth  ] >= tdepth):
                                 
                            # depth increment above the lower depth
                            dinc = tdepth - geo['widthdepth'][ii,kk-1,depth]
                            # width at the water depth
                            thiswidth = geo['widthdepth'][ii,kk,Dwidth]       \
                                * dinc /  geo['widthdepth'][ii,kk,Ddepth]     
                            # add area increment to area below
                            thisarea = geo['widthdepth'][ii,kk,areaTotalBelow]\
                                + 0.5 * (geo['widthdepth'][ii,kk-1,width]     \
                                       + thiswidth) * dinc
                            # if you're here, you've found the top!             
                            break   
                        
                        #endif 
                    #endfor    
                    
                    el['volume'][ii] = thisarea * geo['length'][ii]
                #endif depth
                # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                
            #endif wd
        #endfor ii
        
        return el
#==============================================================================
