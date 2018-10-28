#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 15:06:32 2018

@author: brh
"""

class DataStorage:

#==============================================================================
    
    def set_all_storage(self, NX, setting):
                  
        [geo, setting] = self.geometry(NX,setting)
        el  = self.element(NX)
        fa  = self.face(NX)
        [rkc, setting] = self.rk4storage(NX, setting)
        rke = self.element(NX)
        rkf = self.face(NX)
        
        setting.geometry_types =     \
            ['rectangular_channel',  \
             'parabolic_channel',    \
             'trapezoidal_channel',   \
             'widthdepth_pair']
            
        return [geo, el, fa, rkc, rke, rkf, setting]

#==============================================================================
    def geometry(self, NX, setting): 
        '''
        geometry for each element
        '''
        import numpy as np      
        import sys
        
        #----------------------------------------------------------------------                       
        # fixed geometry of elements for solution
        # element length
        length     = np.zeros(NX,dtype=np.float64 )
        # element bottom elevation
        zbottom    = np.zeros(NX,dtype=np.float64 ) 
        # element Manning's n
        manningsn  = np.zeros(NX,dtype=np.float64 ) 
        # element X distance
        xvalue     = np.zeros(NX,dtype=np.float64 ) 
        # small volume for this element
        smallvolume= np.zeros(NX,dtype=np.float64 ) 
        
        # HACK - we can later combine the unique geometry into one or
        # two common storage locations. Keeping them separate for debugging
        
        #  element bottom breadth (square channel)
        breadth    = np.zeros(NX,dtype=np.float64 )
        # angle (radians) for a simple trapezoidal section
        trapezoid_angle = np.zeros(NX,dtype=np.float64 )  
        # value to characterize a parabolic section
        parabola_value = np.zeros(NX,dtype=np.float64 )  
#        gB         = np.zeros(NX,dtype=np.float64 )  
        #common01   = np.zeros(NX,dtype=np.float64 ) 
        ID         = np.zeros(NX,dtype=int)
        etype      = ['-----------------------------']*NX    
        etype      = np.array(etype)    
        geo = {'length':length, 'zbottom':zbottom, 
                'manningsn':manningsn, 
                'breadth':breadth,
                'trapezoid_angle':trapezoid_angle,
                'parabola_value':parabola_value,
                'smallvolume':smallvolume,               
                'xvalue':xvalue, 'etype':etype, 'ID':ID}

#                'gA': gA, 'gB':gB, 'common01':common01,

        #----------------------------------------------------------------------               
        # setting up additiona data storage for general geometry
        # defined by width-depth pairs.
        if setting.geometry_number_widthdepth_pairs > 0:
            # note that this is the maximum number
            # this data storage may be reset based on input files
            NY = setting.geometry_number_widthdepth_pairs 
            widthdepth = np.zeros([NX,NY,9], dtype=np.float64)
            auxWD1 = np.zeros([NX,NY], dtype=np.float64)
            auxWD2 = np.zeros([NX,NY], dtype=np.float64)
            npair = np.zeros(NX, dtype=int)
            #auxWD2 = np.zeros([NX,NY], dtype=np.float64)
            geo['widthdepth'] = widthdepth
            geo['auxWD1'] = auxWD1
            geo['auxWD2'] = auxWD2
            geo['npair'] = npair
                        
            setting.geometry_widthdepth_values = \
                {'widthAtLayerTop':0, \
                 'depthAtLayerTop':1, \
                 'areaThisLayer':2,  \
                 'areaTotalBelowThisLayer':3, \
                 'Dwidth':4, \
                 'Ddepth':5, \
                 'angle':6, \
                 'perimeterBelowThisLayer':7}
        
        return [geo, setting]
    
#==============================================================================
    def element(self,NX):
        '''
        element data. 
        _up are values in the nominal upstream direction
        _dn are values oin the nominal downstream direction
        '''
        import numpy as np
        # dynamic variables on elements      
        area      = np.zeros(NX, dtype=np.float64) # Cross-sectional area
        CFL_up    = np.zeros(NX, dtype=np.float64) # time-scale CFL
        CFL_dn    = np.zeros(NX, dtype=np.float64) # time-scale CFL
        eta       = np.zeros(NX, dtype=np.float64) # Free surface elevation
        friction  = np.zeros(NX, dtype=np.float64) # Friction on cell        
        flowrate  = np.zeros(NX, dtype=np.float64) # flow rate
        froude    = np.zeros(NX, dtype=np.float64) # Froude number
        hyddepth  = np.zeros(NX, dtype=np.float64) # hydraulic depth
        hydradius = np.zeros(NX, dtype=np.float64) # hydraulic radius
        perimeter = np.zeros(NX, dtype=np.float64) # wetted perimeter
        topwidth  = np.zeros(NX, dtype=np.float64) # Top Width
        tscale_up = np.zeros(NX, dtype=np.float64) # time scale to upstream face
        tscale_dn = np.zeros(NX, dtype=np.float64) # time scale to downstream face
        velocity  = np.zeros(NX, dtype=np.float64) # velocity
        volume    = np.zeros(NX, dtype=np.float64) # volume
        smallvolume_ratio  = np.zeros(NX, dtype=np.float64)
        issmallvolume = np.zeros(NX, dtype=bool)
        isadhocflowrate = np.zeros(NX, dtype=bool)
        aux01      = np.zeros(NX, dtype=np.float64)
        temp1      = np.zeros(NX, dtype=np.float64)

        return {'area':area, 
                'CFL_up':CFL_up, 'CFL_dn':CFL_dn, 
                'eta':eta, 
                'friction':friction, 'flowrate':flowrate,
                'froude':froude, 'hyddepth':hyddepth,
                'hydradius':hydradius, 'perimeter':perimeter, 
                'topwidth':topwidth, 'tscale_up':tscale_up, 
                'tscale_dn':tscale_dn, 'velocity':velocity, 'volume':volume,
                'smallvolume_ratio':smallvolume_ratio,
                'issmallvolume':issmallvolume, 
                'isadhocflowrate':isadhocflowrate,
                'aux01':aux01, 'temp1':temp1}

#==============================================================================
    def face(self,NX):
        '''
        face values
        *P are values to the positive (downstream) side of face ii
        *M are values to the negative (upstream) side of face ii
        '''
        import numpy as np
        # face values, P is downstream, M is upstream
        areaP     = np.zeros(NX+1, dtype=np.float64)
        areaM     = np.zeros(NX+1, dtype=np.float64)
        etaM      = np.zeros(NX+1, dtype=np.float64)
        etaP      = np.zeros(NX+1, dtype=np.float64)
        flowrate  = np.zeros(NX+1, dtype=np.float64)
        topwidth  = np.zeros(NX+1, dtype=np.float64)
        perimeter = np.zeros(NX+1, dtype=np.float64)
        velocityP = np.zeros(NX+1, dtype=np.float64)
        velocityM = np.zeros(NX+1, dtype=np.float64)  
        zbottom   = np.zeros(NX+1, dtype=np.float64)
        cosangle  = np.zeros(NX+1, dtype=np.float64)
        xvalue    = np.zeros(NX+1, dtype=np.float64)
        aux1      = np.zeros(NX+1, dtype=np.float64)
        temp1     = np.zeros(NX+1, dtype=np.float64)
        temp2     = np.zeros(NX+1, dtype=np.float64)
        jumptype  = np.zeros(NX+1, dtype=int)
        
        return { 
                 'areaP':areaP, 'areaM':areaM, 
                 'etaP':etaP, 'etaM':etaM,
                 'flowrate':flowrate, 'topwidth':topwidth, 'perimeter':perimeter,
                 'velocityP':velocityP, 'velocityM':velocityM,
                 'zbottom':zbottom, 'cosangle':cosangle, 'xvalue':xvalue, 
                 'aux1':aux1,'temp1':temp1, 'temp2':temp2, 
                 'jumptype':jumptype}

#==============================================================================
    def rk4storage(self,NX, setting): 
        '''
        storage for rk4 time march
        The RK methods are defined from
        Colin Barr Macdonald, "Constructing High-Order Runge-Kutta Methods with
            Embedded Strong-Stability-Preserving Pairs"
            MS Thesis, Department of Mathematics, Simon Frasier University, 2003; 90 pgs.
            https://www.math.ubc.ca/~cbm/mscthesis/cbm-mscthesis.pdf
            Downloaded 20180217
        '''
        import numpy as np
        import sys
        
        if setting.method_rungekutta == 'rk4_classic':
            setting.rungekutta_levels = 4
            setting.ButcherTableA = np.array( \
                                    [[0.5    , 0.0    , 0.0   ], \
                                     [0.0    , 0.5    , 0.0   ], \
                                     [0.0    , 0.0    , 1.0   ]] )    
            setting.ButcherTableB = np.array( \
                                        [1.0/6.0, 2.0/6.0, 2.0/6.0, 1.0/6.0] )
            
            setting.ButcherTableC = np.array([0, 0.5, 0.5, 1.0])
                                             
        elif setting.method_rungekutta == 'rk4_3/8':
            setting.rungekutta_levels = 4
            setting.ButcherTableA = np.array( \
                                    [[ 1.0/3.0,  0.0    , 0.0 ], \
                                     [-1.0/3.0,  1.0    , 0.0 ], \
                                     [ 1.0    , -1.0    , 1.0 ]])
            setting.ButcherTableB = np.array( \
                                     [ 1.0/8.0,  3.0/8.0, 3.0/8.0, 1.0/8.0] ) 

            setting.ButcherTableC = np.array([1.0/3.0, 2.0/3.0, 1.0])
                
        elif setting.method_rungekutta == 'ssp_(3,3)':
            setting.rungekutta_levels = 3
            setting.ButcherTableA = np.array( \
                                    [[ 1.0 ,    0.0    ], \
                                     [ 0.25,    0.25   ]])
            setting.ButcherTableB = np.array( \
                                     [ 1.0/6.0, 1.0/6.0, 2.0/6.0] )   
            setting.ButcherTableC = np.array([0.5, 1.0])
                
        elif setting.method_rungekutta == 'ssp_(4,3)':
            setting.rungekutta_levels = 4
            setting.ButcherTableA = np.array( \
                                    [[ 0.5,     0.0,     0.0    ], \
                                     [ 0.5,     0.5,     0.0    ], \
                                     [ 1.0/6.0, 1.0/6.0, 1.0/6.0]])
            setting.ButcherTableB = np.array( \
                                     [ 1.0/6.0, 1.0/6.0, 1.0/6.0, 3.0/6.0] )
            setting.ButcherTableC = np.array([0.5, 1.0, 0.5])
            
        elif setting.method_rungekutta == 'ssp_(5,3)':
            setting.rungekutta_levels = 5
            setting.ButcherTableA = np.array( \
                                    [[ 0.37727, 0.0,     0.0    , 0.0    ], \
                                     [ 0.37727, 0.37727, 0.0    , 0.0    ], \
                                     [ 0.24300, 0.24300, 0.24300, 0.0    ], \
                                     [ 0.15359, 0.15359, 0.15359, 0.23846]])
            setting.ButcherTableB = np.array( \
                                     [ 0.20673, 0.20673, 0.11710, 0.18180, 0.28763] )
            setting.ButcherTableC = np.array([0.37727, 0.75454, 0.72899, 0.69923])
                
        elif setting.method_rungekutta == 'ssp_(5,4)':
            setting.rungekutta_levels = 5
            setting.ButcherTableA = np.array( \
                                    [[ 0.39175, 0.0,     0.0    , 0.0    ], \
                                     [ 0.21767, 0.36841, 0.0    , 0.0    ], \
                                     [ 0.082692,0.13996, 0.25189, 0.0    ], \
                                     [ 0.067966,0.11503, 0.20703, 0.54497]])
            setting.ButcherTableB = np.array( \
                                     [ 0.14681, 0.24848, 0.10426, 0.27444, 0.22601])
            setting.ButcherTableC = np.array([0.39175, 0.58608, 0.47454,  0.93501])
            
        elif setting.method_rungekutta == 'ssp_(6,3)':
            setting.rungekutta_levels = 6
            setting.ButcherTableA = np.array( \
                                    [[ 0.28422, 0.0,     0.0    , 0.0    , 0.0,   ], \
                                     [ 0.28422, 0.28422, 0.0    , 0.0    , 0.0    ], \
                                     [ 0.23071, 0.23071, 0.23071, 0.0    , 0.0    ], \
                                     [ 0.13416, 0.13416, 0.13416, 0.16528, 0.0    ], \
                                     [ 0.13416, 0.13416, 0.13416, 0.16528, 0.28422]])
            setting.ButcherTableB = np.array( \
                                     [ 0.17016, 0.17016, 0.10198, 0.12563, 0.21604, 0.21604 ] )
            setting.ButcherTableC = np.array([0.28422, 0.56844, 0.69213, 0.56776, 0.85198])

        elif setting.method_rungekutta == 'ssp_(6,4)':
            setting.rungekutta_levels = 6
            setting.ButcherTableA = np.array( \
                                    [[ 0.35530, 0.0,     0.0    , 0.0    , 0.0,   ], \
                                     [ 0.27049, 0.33179, 0.0    , 0.0    , 0.0    ], \
                                     [ 0.1224 , 0.15014, 0.19721, 0.0    , 0.0    ], \
                                     [ 0.076343,0.093643,0.123  , 0.27182, 0.0    ], \
                                     [ 0.076343,0.093643,0.123  , 0.27182, 0.43582]])
            setting.ButcherTableB = np.array( \
                                     [ 0.15225, 0.18675, 0.15554, 0.13485, 0.2162 , 0.15442 ] )
            setting.ButcherTableC = np.array([0.35530, 0.60227, 0.46975, 0.56481, 1.00063])
        else:
            print('error, unknown value of ',setting.method_rungekutta, \
                  ' for setting.method_rungekutta')
            sys.exit()
        #endif
          
        k1v       = np.zeros(NX, dtype=np.float64 ) 
        k2v       = np.zeros(NX, dtype=np.float64 )
        k1q       = np.zeros(NX, dtype=np.float64 ) 
        k2q       = np.zeros(NX, dtype=np.float64 ) 
        dout = {'k1v':k1v, 'k2v':k2v, 'k1q':k1q, 'k2q':k2q }
        if setting.rungekutta_levels > 2:
            k3v       = np.zeros(NX, dtype=np.float64 ) 
            k3q       = np.zeros(NX, dtype=np.float64 ) 
            dout['k3v'] = k3v
            dout['k3q'] = k3q
        #endif
        if setting.rungekutta_levels > 3:
            k4v       = np.zeros(NX, dtype=np.float64 ) 
            k4q       = np.zeros(NX, dtype=np.float64 )      
            dout['k4v'] = k4v
            dout['k4q'] = k4q
        #endif
        if setting.rungekutta_levels > 4:
            k5v       = np.zeros(NX, dtype=np.float64 ) 
            k5q       = np.zeros(NX, dtype=np.float64 )      
            dout['k5v'] = k5v
            dout['k5q'] = k5q
        #endif
        if setting.rungekutta_levels > 5:
            k6v       = np.zeros(NX, dtype=np.float64 ) 
            k6q       = np.zeros(NX, dtype=np.float64 )      
            dout['k6v'] = k6v
            dout['k6q'] = k6q
        #endif
        if setting.rungekutta_levels > 6:
            k7v       = np.zeros(NX, dtype=np.float64 ) 
            k7q       = np.zeros(NX, dtype=np.float64 )      
            dout['k7v'] = k7v
            dout['k7q'] = k7q
        #endif
        if setting.rungekutta_levels >7:
            k8v       = np.zeros(NX, dtype=np.float64 ) 
            k8q       = np.zeros(NX, dtype=np.float64 )      
            dout['k8v'] = k8v
            dout['k8q'] = k8q
        #endif
        if setting.rungekutta_levels >8:
            k9v       = np.zeros(NX, dtype=np.float64 ) 
            k9q       = np.zeros(NX, dtype=np.float64 )      
            dout['k9v'] = k9v
            dout['k9q'] = k9q
        #endif
        if setting.rungekutta_levels >9:
            print('error, unknown value of ',setting.rungekutta_levels, \
                  ' for setting.rungekutta_levels')
            sys.exit()
        #endif
                         
        return [dout, setting]
   
#==============================================================================
    def savedstorage(self,NX, data_saved_type, setting, present_time):
        '''
        creates storage for all data expected to be saved during simulation
        '''
        import math
        import numpy as np
        import sys
        
        if setting.IC_type == 'restart_file':
            simulation_time = setting.time_total - present_time
        else:
            simulation_time = setting.time_total
        #endif
        
        #print(simulation_time, setting.time_total, present_time)
        #sys.exit()
        
        storeslices = math.ceil( 1.1 *  simulation_time                       \
                                     / setting.binsave_timeinterval )
        data_saved = np.zeros                                                 \
           ((NX,len(data_saved_type),storeslices), dtype=np.float64)
           
        data_saved[:] = np.nan
        
        
        return data_saved 

#==============================================================================
    def timestorage(self,setting):
        '''
        creates storage for time step index and time clock for stored data
        '''
        import math
        import numpy as np
        
        storeslices = math.ceil( 1.1 * setting.time_total                     \
                                     / setting.binsave_timeinterval )
        
        data_saved = np.zeros((2,storeslices), dtype=np.float64)
    
        data_saved[:] = np.nan
        
        return data_saved
    
#==============================================================================
    def saveslice(self, data, data_saved, data_saved_type, thisindex):
        '''
        From the data array, extracts the types in data_saved_type at
        present time and stores in the data_saved array at slice thisindex
        '''
        for ii in range(0,len(data_saved_type)):
            data_saved[:,ii,thisindex] = data[data_saved_type[ii]][:]
                         
        return data_saved
#==============================================================================
    def savetimes(self, data_saved, thisstep, thistime, thisindex):
        '''
        Saves the time step and time clock associated with the saveslice data
        '''
        data_saved[0,thisindex] = thisstep
        data_saved[1,thisindex] = thistime
                
        return data_saved
#==============================================================================
#EOF
