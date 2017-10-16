import unittest
import numpy as np
import tomophantom
import tomophantom.phantom3d
import os
class TestTomophantom3D(unittest.TestCase):
    def test_create_phantom3d(self):
        [tpath, filename] = os.path.split(os.path.abspath(tomophantom.__file__))
        libpath = os.path.join(tpath,'models/Phantom3DLibrary.dat')

        data = tomophantom.phantom3d.build_volume_phantom_3d(libpath,1,256)
        self.assertEqual(data.shape, (256,256,256))
        self.assertNotEqual(data[128,128,128], 0.0)
        self.assertEqual(data[0,0,0],0.0)
        
    def test_create_singoram_phantom3d(self):
        [tpath, filename] = os.path.split(os.path.abspath(tomophantom.__file__))
        libpath = os.path.join(tpath,'models/Phantom3DLibrary.dat')

        angles = np.linspace(0,180, 64, dtype='float32')
        centering = 1 #astra
        data = tomophantom.phantom3d.build_sinogram_phantom_3d(libpath,1,256, 256, angles, centering)
        self.assertEqual(data.shape, (64, 256,256))
        
if __name__ == "__main__":
    unittest.main()