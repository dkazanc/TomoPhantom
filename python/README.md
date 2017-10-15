## Prerequisites
This software depends on cython and numpy.

## Install
   
```
python setup.py build_ext 
python setup.py install
```

To test the code

```
python setup.py test
```

## Usage
```python
tomophantom.phantom3d.build_phantom_3d('<model parameters files>', <model id:int>, <phantom size:int>)
```

```python
from tomophantom import phanton3d
#This will generate 256x256x256 phantom
data = phantom3d.build_phantom_3d('models/Phantom3DLibrary.dat', 1, 256)
```
