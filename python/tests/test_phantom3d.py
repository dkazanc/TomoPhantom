import unittest
import numpy as np
import tomophantom.phantom3d

class TestTomophantom3D(unittest.TestCase):
    def test_create_phantom3d(self):
        data = phantom3d.build_phantom_3d('models/Phantom3DLibrary.dat',1,256)
        self.assetEqual(data.shape, (256,256))