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
tomophantom.phantom3d.build_volume_phantom_3d('<model parameters files>', <model id:int>, <phantom size:int>)
```

```python
from tomophantom import phantom3d
#This will generate 256x256x256 phantom
data = phantom3d.build_volume_phantom_3d(1, 256,'models/Phantom3DLibrary.dat')
```


## Usage
```python
tomophantom.phantom3d.build_sinogram_phantom_3d('<model parameters files>', <model id:int>, <volume size:int>, <detector size:int>, <angles: numpy array float>, <centering: int 0-radon,1-astra>)
```

```python
from tomophantom import phantom3d
#This will generate 64x256x256 phantom with astra centering
data = phantom3d.build_sinogram_phantom_3d('models/Phantom3DLibrary.dat', 1, 256, 256, numpy.linspace(0,180,64,dtype='float32'), 1)
```


