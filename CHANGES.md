# ChangeLog

## v3.0.1 (2024.12)

* Documentation has been updated in January 2025 with a lot of additional information, API links, tutorials, etc.
* Docstrings of the modules were updated with the references to documentation. 


## v3.0 (2023.12)
* Project reorganised into two parts: a library that is build as a shared object using Cmake and Ctypes bindings and pure Python part that can be 
installed separately. 
* In general, ver.3.* should be fully compatabile to ver2.*, however, please note that paths for some functions have been changed. For instance,
```python 
from tomophantom.artefacts import artefacts_mix
```
* Documentation is now available, please check all API for TomoPhantom.
* Demos were updated
