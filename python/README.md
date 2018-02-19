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
# see DemoTomoPhantom.py in Demos for more information
from tomophantom import phantom3d
#This will generate 256x256x256 phantom
model = 1
sizeN = 256
data = phantom3d.buildPhantom3D(model, sizeN,'models/Phantom3DLibrary.dat')
```
