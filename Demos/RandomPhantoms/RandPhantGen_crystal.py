#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
i23 reandom phantoms generator
"""
import numpy as np
import matplotlib.pyplot as plt

from tomophantom import TomoP2D
from tomophantom.TomoP2D import Objects2D
from scipy.ndimage import gaussian_filter
from tomobar.methodsDIR import RecToolsDIR
from tomophantom.artefacts import artefacts_mix
import random


# create a 2D object explicitly without using parameters file
N_size = 512  # define the grid

# A PHANTOM WITHOUT ARTEFACTS
a_el1_min = 0.7
a_el1_max = 0.9
a_el1 = random.uniform(a_el1_min, a_el1_max)
b_el1_min = 0.6
b_el1_max = 0.75
b_el1 = random.uniform(b_el1_min, b_el1_max)

el1 = {
    "Obj": Objects2D.ELLIPSE,
    "C0": 0.7,
    "x0": 0.0,
    "y0": 0.0,
    "a": a_el1,
    "b": b_el1,
    "phi": 0.0,
}

a_el2_min = 0.6
a_el2_max = a_el1
a_el2 = random.uniform(a_el2_min, a_el2_max)
b_el2_min = 0.6
b_el2_max = b_el1
b_el2 = random.uniform(b_el2_min, b_el2_max)

el2 = {
    "Obj": Objects2D.ELLIPSE,
    "C0": -0.4,
    "x0": 0.0,
    "y0": 0.0,
    "a": a_el2,
    "b": b_el2,
    "phi": 0.0,
}

C0_min = 0.01
C0_max = 0.2
C_0 = random.uniform(C0_min, C0_max)
a_el3_min = 0.1
a_el3_max = 0.7
a_el3 = random.uniform(a_el3_min, a_el3_max)
b_el3_min = 0.1
b_el3_max = 0.7
b_el3 = random.uniform(b_el3_min, b_el3_max)
phi_min = 0.0
phi_max = 90.0
phi = random.uniform(phi_min, phi_max)

el3 = {
    "Obj": Objects2D.RECTANGLE,
    "C0": C_0,
    "x0": 0.0,
    "y0": 0.0,
    "a": a_el3,
    "b": b_el3,
    "phi": phi,
}

GROUND_TRUTH = TomoP2D.Object(N_size, [el1, el2, el3])

plt.close("all")
plt.figure(1)
plt.rcParams.update({"font.size": 21})
plt.imshow(GROUND_TRUTH, vmin=0, vmax=1, cmap="BuPu")
plt.colorbar(ticks=[0, 0.5, 1], orientation="vertical")
plt.title("{}".format("Ground Truth Object"))
# %%

# define all objects bellow:
el1 = {
    "Obj": Objects2D.ELLIPSE,
    "C0": 0.7,
    "x0": 0.0,
    "y0": 0.0,
    "a": a_el1,
    "b": b_el1,
    "phi": 0.0,
}

el2 = {
    "Obj": Objects2D.ELLIPSE,
    "C0": 0.7,
    "x0": 0.0,
    "y0": 0.0,
    "a": a_el1 - 0.015,
    "b": b_el1 - 0.015,
    "phi": 0.0,
}

Object1 = TomoP2D.Object(N_size, [el1])
Object2 = TomoP2D.Object(N_size, [el2])

Object_loop = Object1 - Object2
Object_loop = gaussian_filter(Object_loop, sigma=3)
Object_loop = gaussian_filter(Object1 + Object_loop, sigma=3)


el3 = {
    "Obj": Objects2D.ELLIPSE,
    "C0": -0.4,
    "x0": 0.0,
    "y0": 0.0,
    "a": a_el2,
    "b": b_el2,
    "phi": 0.0,
}

Object3 = TomoP2D.Object(N_size, [el3])

Object_liquor = gaussian_filter(Object3, sigma=3)

Object = Object_loop + Object_liquor

el6 = {
    "Obj": Objects2D.RECTANGLE,
    "C0": 2.0 * C_0,
    "x0": 0.0,
    "y0": 0.0,
    "a": a_el3 + 0.02,
    "b": b_el3 + 0.02,
    "phi": phi,
}

Object6 = TomoP2D.Object(N_size, [el6])
Object_crystal_edge = gaussian_filter(Object6, sigma=3)
Object += Object_crystal_edge

el5 = {
    "Obj": Objects2D.RECTANGLE,
    "C0": C_0,
    "x0": 0.0,
    "y0": 0.0,
    "a": a_el3,
    "b": b_el3,
    "phi": phi,
}

Object5 = TomoP2D.Object(N_size, [el5])
Object_crystal = gaussian_filter(Object5, sigma=3)

Object -= Object_crystal

# forming dictionaries with artifact types
_noise_ = {"type": "Gaussian", "sigma": 0.1, "seed": None}  # noise amplitude

# adding zingers and stripes
_zingers_ = {"percentage": 0.5, "modulus": 50}

_stripes_ = {
    "percentage": 1.2,
    "maxthickness": 3.0,
    "intensity": 0.25,
    "type": "partial",
    "variability": 0.005,
}

Object = artefacts_mix(Object, _noise_, _zingers_, _stripes_, _sinoshifts_={})

Object = gaussian_filter(Object, sigma=3)


plt.figure(2)
plt.rcParams.update({"font.size": 21})
plt.imshow(Object, vmin=0, vmax=1, cmap="BuPu")
plt.colorbar(ticks=[0, 0.5, 1], orientation="vertical")
plt.title("{}".format("Distorted Phantom"))

# %%
# Generate projection data of distorted phantom

angles_num = int(np.pi * N_size)
# angles number
angles = np.linspace(0.0, 179.9, angles_num, dtype="float32")
angles_rad = angles * (np.pi / 180.0)
P = N_size  # detectors

Rectools = RecToolsDIR(
    DetectorsDimH=P,  # Horizontal detector dimension
    DetectorsDimV=None,  # Vertical detector dimension (3D case)
    CenterRotOffset=0.0,  # Center of Rotation scalar
    AnglesVec=angles_rad,  # A vector of projection angles in radians
    ObjSize=N_size,  # Reconstructed object dimensions (scalar)
    device_projector="gpu",
)

sino_num = Rectools.FORWPROJ(Object)

_noise_ = {}
_zingers_ = {}
_sinoshifts_ = {}

_stripes_ = {
    "percentage": 0.75,
    "maxthickness": 2.0,
    "intensity": 0.15,
    "type": "mix",
    "variability": 0.005,
}

sino_num_artifacts = artefacts_mix(
    sino_num, _noise_, _zingers_, _stripes_, _sinoshifts_
)

plt.figure(3)
plt.rcParams.update({"font.size": 21})
plt.imshow(sino_num_artifacts, cmap="BuPu")
plt.title("{}".format("Distorted Phantom"))
# %%
FBPrec = Rectools.FBP(sino_num_artifacts)  # perform FBP reconstruction

plt.figure()
plt.rcParams.update({"font.size": 20})
plt.imshow(FBPrec, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation="vertical")
plt.title("FBP reconstruction")
# %%
from tomobar.methodsIR import RecToolsIR

RectoolsIR = RecToolsIR(
    DetectorsDimH=P,  # DetectorsDimH # detector dimension (horizontal)
    DetectorsDimV=None,  # DetectorsDimV # detector dimension (vertical) for 3D case only
    CenterRotOffset=None,  # Center of Rotation (CoR) scalar (for 3D case only)
    AnglesVec=angles_rad,  # array of angles in radians
    ObjSize=N_size,  # a scalar to define reconstructed object dimensions
    datafidelity="LS",  # data fidelity, choose LS, PWLS (wip), GH (wip), Student (wip)
    device_projector="gpu",
)

_data_ = {"projection_norm_data": sino_num_artifacts}  # data dictionary
lc = RectoolsIR.powermethod(
    _data_
)  # calculate Lipschitz constant (run once to initialise)
_algorithm_ = {"iterations": 350, "lipschitz_const": lc}

# adding regularisation using the CCPi regularisation toolkit
_regularisation_ = {
    "method": "PD_TV",
    "regul_param": 0.001,
    "iterations": 150,
    "device_regulariser": "gpu",
}

# Run FISTA reconstrucion algorithm with regularisation
RecFISTA_reg = RectoolsIR.FISTA(_data_, _algorithm_, _regularisation_)

plt.figure()
plt.imshow(RecFISTA_reg, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation="vertical")
plt.title("TV-Regularised FISTA reconstruction")
plt.show()

# %%
# Segment the result of FBPrec while using the GROUND_TRUTH as ideally segmented phantom

# %%
##############################################################################
####################3D phantom generator######################################
##############################################################################
import numpy as np
import matplotlib.pyplot as plt

from tomophantom import TomoP3D
from tomophantom.TomoP3D import Objects3D
from scipy.ndimage import gaussian_filter
from tomobar.methodsDIR import RecToolsDIR
from tomophantom.artefacts import artefacts_mix
import random


N_size = 256
# A PHANTOM WITHOUT ARTEFACTS
a_el1_min = 0.7
a_el1_max = 0.95
a_el1 = random.uniform(a_el1_min, a_el1_max)
b_el1_min = 0.6
b_el1_max = 0.75
b_el1 = random.uniform(b_el1_min, b_el1_max)
c_el1_min = 0.6
c_el1_max = 0.85
c_el1 = random.uniform(c_el1_min, c_el1_max)

el1 = {
    "Obj": Objects3D.ELLIPSOID,
    "C0": 0.7,
    "x0": 0.0,
    "y0": 0.0,
    "z0": 0.0,
    "a": a_el1,
    "b": b_el1,
    "c": c_el1,
    "phi1": 0.0,
}


a_el2_min = 0.6
a_el2_max = a_el1
a_el2 = random.uniform(a_el2_min, a_el2_max)
b_el2_min = 0.6
b_el2_max = b_el1
b_el2 = random.uniform(b_el2_min, b_el2_max)
c_el2_min = 0.6
c_el2_max = c_el1
c_el2 = random.uniform(c_el2_min, c_el2_max)

el2 = {
    "Obj": Objects3D.ELLIPSOID,
    "C0": -0.4,
    "x0": 0.0,
    "y0": 0.0,
    "z0": 0.0,
    "a": a_el2,
    "b": b_el2,
    "c": c_el2,
    "phi1": 0.0,
}

C0_min = 0.01
C0_max = 0.2
C_0 = random.uniform(C0_min, C0_max)
a_el3_min = 0.1
a_el3_max = 0.7
a_el3 = random.uniform(a_el3_min, a_el3_max)
b_el3_min = 0.1
b_el3_max = 0.7
b_el3 = random.uniform(b_el3_min, b_el3_max)
c_el3_min = 0.1
c_el3_max = 0.7
c_el3 = random.uniform(c_el3_min, c_el3_max)
x0_rand = random.uniform(-0.15, 0.15)
y0_rand = random.uniform(-0.15, 0.15)
z0_rand = random.uniform(-0.15, 0.15)
phi_min = 0.0
phi_max = 180.0
phi1 = random.uniform(phi_min, phi_max)

el3 = {
    "Obj": Objects3D.CUBOID,
    "C0": C_0,
    "x0": x0_rand,
    "y0": y0_rand,
    "z0": z0_rand,
    "a": a_el3,
    "b": b_el3,
    "c": c_el3,
    "phi1": phi1,
}

GROUND_TRUTH = TomoP3D.Object(N_size, [el1, el2, el3])

GROUND_TRUTH[GROUND_TRUTH > 0.7] = 0.43336788
GROUND_TRUTH[(GROUND_TRUTH > 0.0) & (GROUND_TRUTH < 0.29999998)] = 0.43336788

sliceSel = (int)(N_size / 2)
plt.figure()
plt.subplot(131)
plt.imshow(GROUND_TRUTH[sliceSel, :, :])
plt.title("Ideal Phantom1")
plt.subplot(132)
plt.imshow(GROUND_TRUTH[:, sliceSel, :])
plt.title("Ideal Phantom2")
plt.subplot(133)
plt.imshow(GROUND_TRUTH[:, :, sliceSel])
plt.title("Ideal Phantom3")
plt.show()
# %%
"""
print("Generating artificial phantom")
el1_add = {'Obj': Objects3D.ELLIPSOID,
      'C0' : 0.7,
      'x0' : 0.0,
      'y0' : 0.0,
      'z0' : 0.0,
      'a'  : a_el1-0.015,
      'b'  : b_el1-0.015,
      'c'  : c_el1-0.015,
      'phi1': 0.0}

Object1 = TomoP3D.Object(N_size, [el1])
Object2 = TomoP3D.Object(N_size, [el1_add])

Object_loop = Object1 - Object2
Object_loop = gaussian_filter(Object_loop, sigma=2)
Object_loop = gaussian_filter(Object1+Object_loop, sigma=2)

Object3 = TomoP3D.Object(N_size, [el2])
Object_liquor = gaussian_filter(Object3, sigma=3)
Object = Object_loop + Object_liquor

el2_add = {'Obj': Objects3D.CUBOID,
      'C0' : 2.0*C_0,
      'x0' : x0_rand,
      'y0' : y0_rand,
      'z0' : z0_rand,
      'a'  : a_el3+0.02,
      'b'  :  b_el3+0.02,
      'c'  :  c_el3+0.02,
      'phi1': phi1}

Object6 = TomoP3D.Object(N_size, [el2_add])
Object_crystal_edge = gaussian_filter(Object6, sigma=2)
Object += Object_crystal_edge

Object5 = TomoP3D.Object(N_size, [el3])
Object_crystal = gaussian_filter(Object5, sigma=2)
Object -= Object_crystal
"""
# forming dictionaries with artifact types
_noise_ = {
    "noise_type": "Gaussian",
    "noise_amplitude": 0.02,  # noise amplitude
    "noise_seed": None,
}

# adding zingers and stripes
_zingers_ = {"zingers_percentage": 1.5, "zingers_modulus": 50}

_stripes_ = {
    "stripes_percentage": 1.2,
    "stripes_maxthickness": 3.0,
    "stripes_intensity": 0.25,
    "stripes_type": "partial",
    "stripes_variability": 0.005,
}

Object = artefacts_mix(GROUND_TRUTH, **_noise_, **_zingers_, **_stripes_)


sliceSel = (int)(N_size / 2)
plt.figure()
plt.subplot(131)
plt.imshow(Object[:, sliceSel, :])
plt.title("Distorted Phantom1")
plt.subplot(132)
plt.imshow(Object[sliceSel, :, :])
plt.title("Distorted Phantom2")
plt.subplot(133)
plt.imshow(Object[:, :, sliceSel])
plt.title("Distorted Phantom3")
plt.show()

# %%
print(
    "Simulate synthetic flat fields, add flat field background to the projections and add noise"
)
from tomophantom.flatsgen import synth_flats
from tomophantom.artefacts import artefacts_mix

I0 = 75000
# Source intensity
flatsnum = 20  # the number of the flat fields required

angles_num = int(np.pi * N_size)
# angles number
angles = np.linspace(0.0, 179.9, angles_num, dtype="float32")
angles_rad = angles * (np.pi / 180.0)
P = N_size  # detectors

Rectools = RecToolsDIR(
    DetectorsDimH=P,  # Horizontal detector dimension
    DetectorsDimV=N_size,  # Vertical detector dimension (3D case)
    CenterRotOffset=0.0,  # Center of Rotation scalar
    AnglesVec=angles_rad,  # A vector of projection angles in radians
    ObjSize=N_size,  # Reconstructed object dimensions (scalar)
    device_projector="gpu",
)

projection_data3D = Rectools.FORWPROJ(Object)
intens_max_clean = np.max(projection_data3D)

_fresnel_propagator_ = {
    "fresnel_dist_observation": 40,
    "fresnel_scale_factor": 10,
    "fresnel_wavelenght": 0.007,
}
projection_data3D_fresnel = artefacts_mix(projection_data3D, **_fresnel_propagator_)
# %%
[projData3D_noisy, flatsSIM] = synth_flats(
    projection_data3D_fresnel,
    source_intensity=I0,
    source_variation=0.01,
    arguments_Bessel=(1, 10, 10, 12),
    specklesize=15,
    kbar=0.3,
    jitter=0.1,
    sigmasmooth=3,
    flatsnum=flatsnum,
)
# del projData3D_analyt

print("Normalise projections using ToMoBAR software")
from tomobar.supp.suppTools import normaliser

# normalise the data, the required format is [detectorsX, Projections, detectorsY]
projData3D_norm = normaliser(
    projData3D_noisy, flatsSIM, darks=None, log="true", method="mean"
)

# del projData3D_noisy
intens_max = 0.3 * np.max(projData3D_norm)
sliceSel = 150
plt.figure()
plt.subplot(131)
plt.imshow(projData3D_norm[:, sliceSel, :], vmin=0, vmax=intens_max)
plt.title("Normalised 2D Projection (erroneous)")
plt.subplot(132)
plt.imshow(projData3D_norm[sliceSel, :, :], vmin=0, vmax=intens_max)
plt.title("Sinogram view")
plt.subplot(133)
plt.imshow(projData3D_norm[:, :, sliceSel], vmin=0, vmax=intens_max)
plt.title("Tangentogram view")
plt.show()

# %%
from tomobar.methodsDIR import RecToolsDIR

RectoolsDIR = RecToolsDIR(
    DetectorsDimH=P,  # DetectorsDimH # detector dimension (horizontal)
    DetectorsDimV=N_size,  # DetectorsDimV # detector dimension (vertical) for 3D case only
    CenterRotOffset=None,  # Center of Rotation (CoR) scalar (for 3D case only)
    AnglesVec=angles_rad,  # array of angles in radians
    ObjSize=N_size,  # a scalar to define reconstructed object dimensions
    device_projector="gpu",
)

print("Reconstruction using FBP from tomobar")
recNumerical = RectoolsDIR.FBP(projData3D_norm)  # FBP reconstruction
recNumerical *= intens_max_clean

sliceSel = int(0.5 * N_size)
max_val = 1
# plt.gray()
plt.figure()
plt.subplot(131)
plt.imshow(recNumerical[sliceSel, :, :], vmin=0, vmax=max_val)
plt.title("3D Reconstruction, axial view")

plt.subplot(132)
plt.imshow(recNumerical[:, sliceSel, :], vmin=0, vmax=max_val)
plt.title("3D Reconstruction, coronal view")

plt.subplot(133)
plt.imshow(recNumerical[:, :, sliceSel], vmin=0, vmax=max_val)
plt.title("3D Reconstruction, sagittal view")
plt.show()

# %%
