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
        
    def test_create_phantom3d_single(self):
        [tpath, filename] = os.path.split(os.path.abspath(tomophantom.__file__))
        libpath = os.path.join(tpath,'models/Phantom3DLibrary.dat')

        data = tomophantom.phantom3d.build_volume_phantom_3d(libpath,1,256)
        params = np.array([(1, 1.00, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.0),], dtype=[('Obj', np.int), ('C0', np.float32), ('x0', np.float32), ('y0',np.float32), ('z0', np.float32),('a',np.float32), ('b', np.float32), ('c', np.float32), ('phi', np.float32)])
        data_single = tomophantom.phantom3d.build_volume_phantom_3d_params(256, params)
        self.assertEqual(np.allclose(data, data_single), True)
        params = np.array([(2, 1.00, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.0),], dtype=[('Obj', np.int), ('C0', np.float32), ('x0', np.float32), ('y0',np.float32), ('z0', np.float32),('a',np.float32), ('b', np.float32), ('c', np.float32), ('phi', np.float32)])
        data_single_not = tomophantom.phantom3d.build_volume_phantom_3d_params(256, params)
        self.assertEqual(np.allclose(data, data_single_not), False)        
        
    def test_create_singoram_phantom3d(self):
        [tpath, filename] = os.path.split(os.path.abspath(tomophantom.__file__))
        libpath = os.path.join(tpath,'models/Phantom3DLibrary.dat')

        angles = np.linspace(0,180, 64, dtype='float32')
        centering = 1 #astra
        data = tomophantom.phantom3d.build_sinogram_phantom_3d(libpath,1,256, 256, angles, centering)
        self.assertEqual(data.shape, (64, 256,256))
        
    def test_create_sinogram_phantom3d_single(self):
        [tpath, filename] = os.path.split(os.path.abspath(tomophantom.__file__))
        libpath = os.path.join(tpath,'models/Phantom3DLibrary.dat')

        angles = np.linspace(0,180, 64, dtype='float32')
        centering = 1 #astra
        data = tomophantom.phantom3d.build_sinogram_phantom_3d(libpath,1,256, 256, angles, centering)
        params = np.array([(1, 1.00, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.0),], dtype=[('Obj', np.int), ('C0', np.float32), ('x0', np.float32), ('y0',np.float32), ('z0', np.float32),('a',np.float32), ('b', np.float32), ('c', np.float32), ('phi', np.float32)])
        data_single = tomophantom.phantom3d.build_sinogram_phantom_3d_params(256, 256, angles, centering, params)
        self.assertEqual(np.allclose(data, data_single), True)  
        params = np.array([(2, 1.00, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.0),], dtype=[('Obj', np.int), ('C0', np.float32), ('x0', np.float32), ('y0',np.float32), ('z0', np.float32),('a',np.float32), ('b', np.float32), ('c', np.float32), ('phi', np.float32)])
        data_single_not = tomophantom.phantom3d.build_sinogram_phantom_3d_params(256, 256, angles, centering, params)
        self.assertEqual(np.allclose(data, data_single_not), False)        
        
        
if __name__ == "__main__":
    unittest.main()