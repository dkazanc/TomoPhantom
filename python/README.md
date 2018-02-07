## Prerequisites
This software depends on cython and numpy.

## Install   
Installation using conda (preferred):
```
conda build conda-recipe --numpy 1.12 --python 3.5
conda install tomophantom --use-local --force
```
or stand-alone
```
python setup.py build_ext --inplace
python setup.py install --user
```

## Generating 2D/3D phantoms
```python
tomophantom.phantom2d.buildPhantom2D( <model id:int>, <phantom size:int>, '<model parameters files>')
tomophantom.phantom3d.buildPhantom3D( <model id:int>, <phantom size:int>, '<model parameters files>')
```

```python
# see DemoTomoPhantom.py in Demos for more information
from tomophantom import phantom3d
#This will generate 256x256x256 phantom
model = 1
sizeN = 256
data = phantom3d.buildPhantom3D(model, sizeN,'models/Phantom3DLibrary.dat')
```

## Usage
```python
### under development
tomophantom.phantom3d.buildSinogram('<model parameters files>', <model id:int>, <volume size:int>, <detector size:int>, <angles: numpy array float>, <centering: int 0-radon,1-astra>)
```

```python
### under development
from tomophantom import phantom3d
#This will generate 64x256x256 phantom with astra centering
data = phantom3d.build_sinogram_phantom_3d('models/Phantom3DLibrary.dat', 1, 256, 256, numpy.linspace(0,180,64,dtype='float32'), 1)
```
