#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""

import numpy as np
import astra

class AstraOper:
    """A simple 2D parallel beam projection/backprojection class based on ASTRA toolbox"""
    def __init__(self, data, DetectorsDim, AnglesVec, ObjSize):        
        self.data = data
        self.DetectorsDim = DetectorsDim
        self.AnglesVec = AnglesVec
        self.ObjSize = ObjSize
        self.proj_geom = astra.create_proj_geom('parallel', 1.0, DetectorsDim, AnglesVecP)
		self.vol_geom = astra.create_vol_geom(ObjSize, ObjSize)		
    def forwproj(self):
        """Applying forward projection"""
		#proj_id = astra.create_projector('cuda', self.proj_geom, self.vol_geom) # for GPU
		proj_id = astra.create_projector('strip', self.proj_geom, self.vol_geom) # for CPU
		sinogram_id, sinogram = astra.create_sino(self.data, proj_id)
		astra.data2d.delete(sinogram_id)
		astra.data2d.delete(proj_id)
        return sinogram
    def backproj(self):
        """Applying backprojection"""
        #proj_id = astra.create_projector('cuda', self.proj_geom, self.vol_geom) # for GPU
		proj_id = astra.create_projector('strip', self.proj_geom, self.vol_geom) # for CPU
        rec_id, X = astra.create_backprojection(self.data, proj_id)
		astra.data2d.delete(proj_id)
		astra.data2d.delete(rec_id)		
		return X
