#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 09:13:52 2018

@author: brh
"""

class BoundaryConditions:
    
    def flow(self, setting):
        
        flow = setting.inflowBC
        
        return flow
    
    def height(self, setting):
        
        height = setting.heightBC
        
        return height