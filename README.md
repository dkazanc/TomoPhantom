# TomoPhantom (ver 1.1)
****************
Software to generate simple analytical phantoms for image processing (reconstruction, denoising, deblurring...)
****************
This software is recommended for various image processing tasks that require numerical testing: image reconstruction, denoising, deblurring, etc. The popular Shepp-Logan phantom is not always fit for testing due to its piecewise-constant nature. This package provides a modular approach to build customisable phantoms consisting of piecewise constant and smooth analytical objects, such as: gaussians,  parabolas, ellipses, cones, rectangulars, etc. Additionally, due to analytical approach of calculating phantoms there is no loss of the spatial resolution due to the absence of the interpolation.

Package contents: PhantomGeneratorDemo.m : demo to run; PhantomLibrary.dat : text file with some models already created, please
feel free to design own models or edit existing ones; buildPhantom.m: function to build phantoms; DeformObject.m and supp files: function to deform objects according to nonlinear geometrical transformation proposed in [1]

If used please refer to:
[1] [D. Kazantsev, V. Pickalov "New iterative reconstruction methods for fan-beam tomography", IPSE, 2017](https://ccpforge.cse.rl.ac.uk/gf/download/frsrelease/582/8704/GP_IPSE.pdf)

Future plans: python code; 3D version (required?)
If any questions, please e-mail daniil.kazantsev@manchester.ac.uk 
