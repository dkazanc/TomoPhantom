<div align="center">
  <img src="docs/img/TomoPhantomLogo.png" height="350"><br>
  <img src="docs/img/models2Dtime/2DtModel14.gif" height="175"><img src="docs/img/models4D/model11_4D.gif "height="175" width="200"><br>
</div>

****************
**TomoPhantom [1] is a toolbox to generate customisable 2D-4D phantoms (with a temporal capability) and their 
analytical projection data for various image processing tasks (reconstruction, denoising, deblurring, etc.).**

<a href="https://zenodo.org/badge/latestdoi/95991001"><img src="https://zenodo.org/badge/95991001.svg" alt="DOI"></a>
****************    
   
 <div class="post-content">
        <h3 class="post-title">About TomoPhantom </h3>
        <p> TomoPhantom is recommended for various image processing tasks that require extensive numerical testing: image reconstruction, denoising, deblurring, etc. Specifically,TomoPhantom is best-suited for testing various tomographic image reconstruction (TIR) methods. For TIR algorithms testing, the popular <a href="https://en.wikipedia.org/wiki/Shepp%E2%80%93Logan_phantom">Shepp-Logan phantom</a> 
is not always a good choice due to its piecewise-constant nature. This toolbox provides a simple modular approach to efficently build customisable 2D-4D phantoms consisting of 
piecewise-constant, piecewise-smooth, and smooth analytical objects. The objects include: Gaussians, parabolas, ellipses, cones, rectangulars. The exact tomographic projections (analytical sinograms) can be also obtained. The sinograms can be used for TIR benchmarking purposes, which enables reconstruction without so-called <a href="http://www.sciencedirect.com/science/article/pii/S0377042705007296">'Inverse Crime'</a>. TomoPhantom is also compatable with <a href="http://www.astra-toolbox.com/">ASTRA-toolbox</a> and <a href="http://tomopy.readthedocs.io/en">TomoPy</a> packages. Generated data can be directly reconstructed using the toolboxes (see examples). Additionally, TomoPhantom provides a simple temporal extension, i.e. a capability of creating 2D+time and 3D+time models.   
        </p>
 </div>

### **TomoPhantom** prerequisites: 

 * [MATLAB](www.mathworks.com/products/matlab/) OR
 * Python (tested ver. 3.5); Cython
 * C compilers: GCC/MinGW/[TDM-GCC](http://tdm-gcc.tdragon.net/)/Visual Studio

### Other dependencies (reconstruction):
 * [ASTRA-toolbox](http://www.astra-toolbox.com/)
 * [TomoPy](http://tomopy.readthedocs.io)

## Installation:

### Python (conda-build preferrable)
```
	conda build conda-recipe --numpy 1.12 --python 3.5
	conda install tomophantom --use-local --force
```
### Matlab
```
	run compile_mex_linux.m % to compile CPU modules on linux
	run compile_mex_windows.m % to compile CPU modules on Windows
```

### Package modules:
- **Phantom2DLibrary.dat** and **Phantom3DLibrary.dat** are editable text files with parametrised models (2D/3D versions of Shepp-Logan, Defrise, and QRM phantoms are included). The generation of new phantoms is highly encouraged, please submit them through pull requests or via e-mail bellow. 
- See MATLAB and Python demos

### License:
- TomoPhantom is released under [Apache License v.2](http://www.apache.org/licenses/LICENSE-2.0). Note that some demos where ['ASTRA-toolbox'](http://www.astra-toolbox.com/) is used are of GPLv3 license and also BSD-3 license for [TomoPy](http://tomopy.readthedocs.io/en) package.

### Related software projects on GitHub:
- [xdesign](https://github.com/tomography/xdesign) XDesign is an open-source Python package for generating configurable simulation phantoms for benchmarking tomographic image reconstruction.
- [syris](https://github.com/ufo-kit/syris) Syris (synchrotron radiation imaging simulation) is a framework for simulations of X-ray absorption and phase contrast dynamic imaging experiments, like time-resolved radiography, tomography or laminography.

### References:

[1] [D. Kazantsev et al. 2018, *TomoPhantom, a software package to generate 2D-4D analytical phantoms for CT image reconstruction algorithm benchmarks*, Software X, Volume 7, January–June 2018, Pages 150–155](https://doi.org/10.1016/j.softx.2018.05.003)

[2] [D. Kazantsev, V. Pickalov "New iterative reconstruction methods for fan-beam tomography", IPSE, 2017](https://doi.org/10.1080/17415977.2017.1340946)

### Applications: 
* [Regularised FISTA-type iterative reconstruction algorithm for X-ray tomographic reconstruction with highly inaccurate measurements](https://github.com/dkazanc/FISTA-tomo)
* [Joint image reconstruction method with correlative multi-channel prior for X-ray spectral computed tomography](https://github.com/dkazanc/multi-channel-X-ray-CT)

Software related questions/comments please e-mail to Daniil Kazantsev at dkazanc@hotmail.com
