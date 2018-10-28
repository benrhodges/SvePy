#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 10:30:19 2018

@author: brh
"""

class UpdateValues:

#==============================================================================
    def update_auxiliary_relationships(self, el, fa, geo, BC, setting, NX):
        
        import sys
        
        import ElementGeometry        
        eg = ElementGeometry.ElementGeometry()
        
        import ElementDynamics
        ed = ElementDynamics.ElementDynamics()
        
        import FaceGeometry
        fg = FaceGeometry.FaceGeometry()
        
        import FaceDynamics
        fd = FaceDynamics.FaceDynamics()
 
        #----------------------------------------------------------------------   
        # update the element geometry
        el = eg.element_geometry(el, geo, setting)

        #----------------------------------------------------------------------   
        # update element dynamics
        # note that using fa in element is a time lag
        el = ed.element_dynamics(el, fa, geo, setting)

        #----------------------------------------------------------------------   
        # update the face geometry
        fa = fg.face_geometry(el, fa, geo, BC, setting, NX)

        #----------------------------------------------------------------------   
        # update the face dynamics
        [fa,el] = fd.face_dynamics(el, fa, geo, BC, setting, NX)

        #----------------------------------------------------------------------   
        # debugging printouts
        self.check_values(el, fa, setting, NX)
        
        return [el, fa]
    
#==============================================================================
    def check_values(self, el, fa, setting, NX):
        
        if setting.check_element != None:
            te = setting.check_element 
            print()
            print('ID = ',geo['ID'][te])
            print('element ',te)
            print('issmallV = ',el['issmallvolume'][te])
            print('Vratio   = ',el['smallvolume_ratio'][te])
            print('xvalue   = ',geo['xvalue'][te])
            print('zbottom  = ',geo['zbottom'][te])
            print('depth    = ',el['eta'][te] - geo['zbottom'][te])
            print('hyddepth = ',el['hyddepth'][te])
            print('CFL_up   = ',el['CFL_up'][te])
            print('CFL_dn   = ',el['CFL_dn'][te])
            print('tscale_up= ',el['tscale_up'][te])
            print('tscale_dn= ',el['tscale_dn'][te])
            print('volume   = ',el['volume'][te])
            print('area     = ',el['area'][te])
            print('flowrate = ',el['flowrate'][te])
            print('friction = ',el['friction'][te])
            print('eta      = ',el['eta'][te])
            print('velocity = ',el['velocity'][te])
            print('energy   = ',el['energy_dia'][te])
            print('froude   = ',el['froude'][te])
            print()

        if setting.check_element != None:
            te = setting.check_element
            print()
            print('ID = ',geo['ID'][te])
            #te = [te[0],te[1],317]
            print('fa ',te)            
            print('xvalue    = ',fa['xvalue'][te])
            print('zbottom   = ',fa['zbottom'][te])
            print('areaM     = ',fa['areaM'][te])
            print('areaP     = ',fa['areaP'][te])
            print('flowrate  = ',fa['flowrate'][te])
            print('topwidth  = ',fa['topwidth'][te])
            print('etaM      = ',fa['etaM'][te])
            print('etaP      = ',fa['etaP'][te])
            print('velocityM = ',fa['velocityM'][te])
            print('velocityP = ',fa['velocityP'][te])
            print('jumptype  = ',fa['jumptype'][te])
            print()
            
        return
    
#==============================================================================
    def zero_all_except_element_flow_and_volume(self, el, fa):
        '''
        Set el and fa storage to zero (except fixe geometry in face)
        This is useful for testing, but should be removed in production code.
        '''
        
        el['area'][:]= 0.0
        el['CFL_up'][:]= 0.0 
        el['CFL_dn'][:]= 0.0
        el['eta'][:]= 0.0 
        el['friction'][:]= 0.0
        el['froude'][:]= 0.0
        el['hyddepth'][:]= 0.0
        el['hydradius'][:]= 0.0
        el['perimeter'][:]= 0.0
        el['topwidth'][:]= 0.0
        el['tscale_up'][:]= 0.0 
        el['tscale_dn'][:]= 0.0 
        el['velocity'][:]= 0.0
        el['smallvolume_ratio'][:]= 0.0
        el['issmallvolume'][:]= False
        el['aux01'][:]= 0.0
        el['temp1'][:]= 0.0

        fa['areaP'][:]= 0.0
        fa['areaM'][:]= 0.0
        fa['etaP'][:]= 0.0
        fa['etaM'][:]= 0.0
        fa['flowrate'][:]= 0.0
        fa['topwidth'][:]= 0.0
        fa['velocityP'][:]= 0.0
        fa['velocityM'][:]= 0.0
        fa['aux1'][:]= 0.0
        fa['temp1'][:]= 0.0
        fa['temp2'][:]= 0.0
        fa['jumptype'][:]= 0              
        
        return [el, fa]
#==============================================================================

#EOF
    