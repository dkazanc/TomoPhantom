# TomoPhantom (ver 0.1)
****************
Software to generate simple analytical phantoms for image processing (reconstruction, denoising, deblurring...)
****************
This software is ideally suitable for various image processing tasks that require numerical testing: image reconstruction, denoising, deblurring, etc. The popular Shepp-Logan phantom is not always suitable for testing due to its piecewise-constant nature. This package provides a modular approach to build customisable phantoms consisting of piecewise constant and smooth analytical objects, such as: gaussians,  parabolas, ellipses, cones, rectangulars, etc. Additionally, due to analytical approach of calculating phantoms, any discrete grid can be used without interpolation and therefore the loss of the spatial resolution.

Package contents:
PhantomGeneratorDemo.m : demo to run
PhantomLibrary.dat : text file with some models already created, please
feel free to design own models or edit existing ones.
buildPhantom:function to build phantoms
DeformObject and supp files : function to deform objects according to
nonlinear geometrical transformation proposed in [1]

Future plans: 
*python code
*Analytical sinograms (Radon transforms)
*3D version (required?)

Refrences:
[1] D. Kazantsev, V. Pickalov "New iterative reconstruction methods for 
fan-beam tomography", IPSE, 2017

If any questions, please e-mail daniil.kazantsev@manchester.ac.uk 
