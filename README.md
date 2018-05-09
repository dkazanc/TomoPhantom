<div align="center">
  <img src="docs/img/TomoPhantomLogo.png" height="350"><br>
  <img src="docs/img/models2Dtime/2DtModel14.gif" height="175"><img src="docs/img/models4D/model11_4D.gif "height="175" width="200"><br>
</div>

****************
**TomoPhantom is a toolbox to generate customisable 2D/3D/4D phantoms (with a temporal capability) and their 
analytical projection datae for various image processing tasks (reconstruction, denoising, deblurring, etc.).**

<a href="https://zenodo.org/badge/latestdoi/95991001"><img src="https://zenodo.org/badge/95991001.svg" alt="DOI"></a>
****************    
   
 <div class="post-content">
        <h3 class="post-title">About TomoPhantom </h3>
        <p> **TomoPhantom** is recommended for various image processing tasks that require extensive numerical testing: image reconstruction, denoising, deblurring, etc. 
Specifically, the software is best-suited for testing various tomographic image reconstruction (TIR) methods. For TIR algorithms testing, the popular <a href="https://en.wikipedia.org/wiki/Shepp%E2%80%93Logan_phantom">Shepp-Logan phantom</a> 
is not always a good choice due to its piecewise-constant nature. This toolbox provides a simple modular approach to efficently build customisable 2D-4D phantoms consisting of 
piecewise-constant and also smooth analytical objects. The objects include: Gaussians, parabolas, ellipses, cones, rectangulars. The exact tomographic projections (sinograms) as a result of applying Radon
Transform to analytical objects can be obtained. The sinograms can be used for TIR benchmarking purposes
without so-called the <a href="http://www.sciencedirect.com/science/article/pii/S0377042705007296">'Inverse Crime'</a>. TomoPhantom is also compatable with 
<a href="http://www.astra-toolbox.com/">ASTRA-toolbox</a> and <a href="http://tomopy.readthedocs.io/en">TomoPy</a> packages. Generated data can be directly reconstructed using 
toolboxes (see examples). Additionally, **TomoPhantom** provides simple temporal extension, therefore a capability of creating 2D+time and 3D+time objects.   
        </p>
 </div>

## **Tomophantom** prerequisites: 

 * [MATLAB](www.mathworks.com/products/matlab/) OR
 * Python (tested ver. 3.5); Cython
 * C compilers (GCC/MinGW/Visual Studio)

### Other dependencies (reconstruction):
 * [ASTRA-toolbox](http://www.astra-toolbox.com/)
 * [TomoPy](http://tomopy.readthedocs.io/en)

## Installation:

### Python (conda-build preferrable)
```
	conda build conda-recipe --numpy 1.12 --python 3.5
	conda install tomophantom --use-local --force
```
### Matlab
```
	run compile_mex.m % to compile CPU modules
```

### Package modules:
- **Phantom2DLibrary.dat** and **Phantom3DLibrary.dat** are editable text files with parametrised models
- See MATLAB and Python demos

### License:
- The project uses Apache License v.2, but some demos where ['ASTRA-toolbox'](http://www.astra-toolbox.com/) is used are of GPLv3 license

### Related software projects on GitHub:
- [xdesign](https://github.com/tomography/xdesign) XDesign is an open-source Python package for generating configurable simulation phantoms for benchmarking tomographic image reconstruction.
- [syris](https://github.com/ufo-kit/syris) Syris (synchrotron radiation imaging simulation) is a framework for simulations of X-ray absorption and phase contrast dynamic imaging experiments, like time-resolved radiography, tomography or laminography.

### References:

[1] D. Kazantsev et al. 2018, *TomoPhantom, a software package to generate 2D-4D analytical phantoms for CT image reconstruction algorithm benchmarks*, Software X (accepted)

[2] [D. Kazantsev, V. Pickalov "New iterative reconstruction methods for fan-beam tomography", IPSE, 2017](https://doi.org/10.1080/17415977.2017.1340946)

Software related questions/comments please e-mail to Daniil Kazantsev at dkazanc@hotmail.com
