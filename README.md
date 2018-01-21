<div align="center">
  <img src="docs/img/TomoPhantomLogo.jpg" height="350"><br>
</div>

****************
**TomoPhantom** is a toolbox to generate customisable analytical 2D and 3D phantoms for various image processing tasks (reconstruction, denoising, deblurring, etc.).
****************

### Detailed description:

**TomoPhantom** is recommended for various image processing tasks that require extensive numerical testing: image reconstruction, denoising, deblurring, etc. 
The software is mainly suited for tomographic image reconstruction (TIR) cases. For TIR algorithms testing, the popular [Shepp-Logan phantom](https://en.wikipedia.org/wiki/Shepp%E2%80%93Logan_phantom) is not always a 
good choice due to the piecewise-constant nature. This toolbox provides a simple modular approach to build customisable 2D and 3D phantoms consisting of 
piecewise-constant and also smooth analytical objects. The objects include: gaussians, parabolas, ellipses, cones, rectangulars, etc. The exact Radon
Transform (currently parallel beam geometry only available) can be obtained, therefore producing an analytical sinogram. The sinograms can be used for TIR testing purposes
without so-called the ['Inverse Crime'](http://www.sciencedirect.com/science/article/pii/S0377042705007296). TomoPhantom is also compatable with 
['ASTRA-toolbox'](http://www.astra-toolbox.com/) and therefore the generated data can be directly reconstructed using ASTRA-toolbox (see provided examples). 

### Package contents:

**TomoPhantom** is available for MATLAB and Python (core modules are written in C-OMP)
- **Phantom2DGeneratorDemo.m** and **Phantom3DGeneratorDemo.m** are MATLAB demo scripts
- **SpectralPhantomDemo.m** MATLAB script to generate spectral phantom with 4 dedicated materials
- **DemoTomoPhantom.py** Python demo script showing 2D and 3D cases
- **Phantom2DLibrary.dat** and **Phantom3DLibrary.dat** are editable text files with parametrised models

### Installation:
- For MATLAB, run **compile_mex.m** to compile MEX-ed C functions
- For Python, see ReadMe in python 'directory'

### License:
- The project uses Apache License v.2, but some demos where ['ASTRA-toolbox'](http://www.astra-toolbox.com/) is used are of GPLv3 license


### If software is used, please cite the paper:

[1] [D. Kazantsev, V. Pickalov "New iterative reconstruction methods for fan-beam tomography", IPSE, 2017](https://ccpforge.cse.rl.ac.uk/gf/download/frsrelease/582/8704/GP_IPSE.pdf)

For any questions, please e-mail daniil.kazantsev@manchester.ac.uk 

