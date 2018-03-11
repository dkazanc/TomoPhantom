## Prerequisites
This software depends on cython and numpy.

## Install   
Installation using conda (preferred):
```
conda build conda-recipe --numpy 1.12 --python 3.5
conda install tomophantom --use-local --force
```
Stand-alone
```
python setup.py build_ext --inplace
python setup.py install --user
```

## Usage:

```python

from tomophantom import TomoP2D
from tomophantom import TomoP3D
model = 1
sizeN = 256
#This will generate 256x256 phantom from model no.1
phantom = TomoP2D.Model(model, sizeN,'path/to/the/Phantom2DLibrary.dat')

angles_num = int(0.5*np.pi*N_size); # angles number
angles = np.linspace(0,180,angles_num,dtype='float32')
angles_rad = angles*(np.pi/180)
P = int(np.sqrt(2)*N_size) #detectors
#This will generate a sinogram of model no.1
sino_an = TomoP2D.ModelSino(model, N_size, P, angles, 'path/to/the/Phantom2DLibrary.dat')

# Lets build 3D phantom without using the library
N3D = 256 # size of the phantom
params = np.array([('gaussian', 1.00, -0.25, 0.1, 0.0, 0.2, 0.35, 0.7, 30.0, 60.0, -25.0),], dtype=[('Obj',  '|S22'), ('C0', np.float32), ('x0', np.float32), 
('y0',np.float32), ('z0',np.float32), ('a',np.float32), ('b', np.float32), ('c', np.float32), ('psi1', np.float32),('psi2', np.float32),('psi3', np.float32)])
Object3D = TomoP3D.Object(N3D, params)
```
## see Demos for more information
DemoModel.py - shows how to create 2D/3D models and their tomographic projections

DemoModel_temporal.py - demonstrate temporal capabilities of the software creating 3D/4D models 

DemoObject.py - shows how to generate use software without library (*dat) files
