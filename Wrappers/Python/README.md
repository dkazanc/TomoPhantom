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

## see Demos for more information
