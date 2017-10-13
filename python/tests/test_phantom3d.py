import unittest
import numpy as np
import tomophantom
import tomophantom.phantom3d
import os
class TestTomophantom3D(unittest.TestCase):
    def test_create_phantom3d(self):
        [tpath, filename] = os.path.split(os.path.abspath(tomophantom.__file__))
        libpath = os.path.join(tpath,'models/Phantom3DLibrary.dat')
        print(libpath)
        data = tomophantom.phantom3d.build_phantom_3d(libpath,1,256)
        self.assertEqual(data.shape, (256,256,256))
        self.assertNotEqual(data[128,128,128], 0.0)
        self.assertEqual(data[0,0,0],0.0)
        
        
if __name__ == "__main__":
    unittest.main()