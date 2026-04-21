#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPLv3 license (ASTRA toolbox)
Note that the TomoPhantom package is released under Apache License, Version 2.0

* Script to generate 3D analytical phantoms and their projection data using TomoPhantom
* Projection data is also generated numerically and reconstructed using ToMoBAR

>>>>> Dependencies (reconstruction): <<<<<
1. ASTRA toolbox
2. ToMoBAR
3. CuPy

@author: Daniil Kazantsev
"""
import timeit
import os
import matplotlib.pyplot as plt
import numpy as np
import tomophantom
from tomophantom import TomoP3D
from tomophantom.qualitymetrics import QualityTools

print("Building 3D phantom using TomoPhantom software")
tic = timeit.default_timer()
model = 13  # select a model number from the library
N_size = 256  # Define phantom dimensions using a scalar value (cubic phantom)
path = os.path.dirname(tomophantom.__file__)
path_library3D = os.path.join(path, "phantomlib", "Phantom3DLibrary.dat")
# This will generate a N_size x N_size x N_size phantom (3D)
phantom_tm = TomoP3D.Model(model, N_size, path_library3D)
toc = timeit.default_timer()
Run_time = toc - tic
print("Phantom has been built in {} seconds".format(Run_time))

sliceSel = int(0.5 * N_size)
# plt.gray()
plt.figure()
plt.subplot(131)
plt.imshow(phantom_tm[sliceSel, :, :], vmin=0, vmax=1)
plt.title("3D Phantom, axial view")

plt.subplot(132)
plt.imshow(phantom_tm[:, sliceSel, :], vmin=0, vmax=1)
plt.title("3D Phantom, coronal view")

plt.subplot(133)
plt.imshow(phantom_tm[:, :, sliceSel], vmin=0, vmax=1)
plt.title("3D Phantom, sagittal view")
plt.show()

# Projection geometry related parameters:
Horiz_det = int(2 * N_size)  # detector column count (horizontal)
Vert_det = N_size  # detector row count (vertical) (no reason for it to be > N)
angles_num = int(0.5 * np.pi * N_size)
# angles number
angles = np.linspace(0.0, 179.9, angles_num, dtype="float32")  # in degrees
angles_rad = angles * (np.pi / 180.0)
# %%
print("Building 3D analytical projection data with TomoPhantom")
projData3D_analyt = TomoP3D.ModelSino(
    model, N_size, Horiz_det, Vert_det, angles, path_library3D
)

intens_max = 70
sliceSel = int(0.5 * N_size)
plt.figure()
plt.subplot(131)
plt.imshow(projData3D_analyt[:, sliceSel, :], vmin=0, vmax=intens_max)
plt.title("2D Projection (analytical)")
plt.subplot(132)
plt.imshow(projData3D_analyt[sliceSel, :, :], vmin=0, vmax=intens_max)
plt.title("Sinogram view")
plt.subplot(133)
plt.imshow(projData3D_analyt[:, :, sliceSel], vmin=0, vmax=intens_max)
plt.title("Tangentogram view")
plt.show()
# %%
print("Adding noise to projection data")
from tomophantom.artefacts import artefacts_mix

# forming dictionaries with artifact types
_noise_ = {
    "noise_type": "Poisson",
    "noise_amplitude": 10000,  # noise amplitude
    "noise_seed": 0,
}

_stripes_ = {
    "stripes_percentage": 1.2,
    "stripes_maxthickness": 3,
    "stripes_intensity": 0.25,
    "stripes_type": "mix",
    "stripes_variability": 0.005,
}


projData3D_analyt_noisy = artefacts_mix(projData3D_analyt, **_noise_, **_stripes_)

intens_max = np.max(projData3D_analyt_noisy)
sliceSel = int(0.5 * N_size)
plt.figure()
plt.subplot(131)
plt.imshow(projData3D_analyt_noisy[:, sliceSel, :], vmin=0, vmax=intens_max)
plt.title("2D noisy Projection (analytical)")
plt.subplot(132)
plt.imshow(projData3D_analyt_noisy[sliceSel, :, :], vmin=0, vmax=intens_max)
plt.title("Noisy sinogram view")
plt.subplot(133)
plt.imshow(projData3D_analyt_noisy[:, :, sliceSel], vmin=0, vmax=intens_max)
plt.title("Noisy tangentogram view")
plt.show()
# %%
print("Reconstruction using FBP from tomobar")
# initialise tomobar DIRECT reconstruction class ONCE
from tomobar.methodsDIR import RecToolsDIR

RectoolsDIR = RecToolsDIR(
    DetectorsDimH=Horiz_det,  # DetectorsDimH # detector dimension (horizontal)
    DetectorsDimH_pad=0,  # Padding size of horizontal detector
    DetectorsDimV=Vert_det,  # DetectorsDimV # detector dimension (vertical) for 3D case only
    CenterRotOffset=None,  # The Center of Rotation (CoR) scalar
    AnglesVec=angles_rad,  # array of angles in radians
    ObjSize=N_size,  # a scalar to define reconstructed object dimensions
    device_projector="gpu",
)

recNumerical = RectoolsDIR.FBP(projData3D_analyt_noisy)  # FBP reconstruction

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

# calculate errors
Qtools = QualityTools(phantom_tm, recNumerical)
RMSE = Qtools.rmse()
print("Root Mean Square Error is {}".format(RMSE))
# %%
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print("Reconstructing with FISTA-OS-TV method using tomobar")
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
# ! you will need to have CuPy installed to be able to run iterative methods
import cupy as cp 
from tomobar.methodsIR_CuPy import RecToolsIRCuPy
input_data_labels = ["detY", "angles", "detX"]


Rectools = RecToolsIRCuPy(
    DetectorsDimH=Horiz_det,  # DetectorsDimH # detector dimension (horizontal)
    DetectorsDimH_pad=0,  # Padding size of horizontal detector
    DetectorsDimV=Vert_det,  # DetectorsDimV # detector dimension (vertical) for 3D case only
    CenterRotOffset=0.0,  # Center of Rotation (CoR) scalar (for 3D case only)
    AnglesVec=angles_rad,  # array of angles in radians
    ObjSize=N_size,  # a scalar to define reconstructed object dimensions
    device_projector=0,
    OS_number=10,  # The number of ordered subsets
)
# prepare dictionaries with parameters:
_data_ = {
    "data_fidelity": "LS",
    "projection_data":  cp.asarray(projData3D_analyt_noisy),  # Normalised projection data
    "data_axes_labels_order": input_data_labels,
    }
lc = Rectools.powermethod(
    _data_
)  # calculate Lipschitz constant (run once to initialise)

# Run FISTA-OS reconstrucion algorithm without regularisation
_algorithm_ = {"iterations": 15, "lipschitz_const": lc}

_regularisation_ = {
    "method": "PD_TV",  # Selected regularisation method
    "regul_param": 0.000075,  # Regularisation parameter
    "iterations": 40,  # The number of regularisation iterations
    "half_precision": True,  # enabling half-precision calculation
}


RecFISTA_os_reg = Rectools.FISTA(_data_, _algorithm_, _regularisation_)

RecFISTA_os_reg = cp.asnumpy(RecFISTA_os_reg)

sliceSel = int(0.5 * N_size)
max_val = 1
plt.figure()
plt.subplot(131)
plt.imshow(RecFISTA_os_reg[sliceSel, :, :], vmin=0, vmax=max_val)
plt.title("3D FISTA Reconstruction, axial view")

plt.subplot(132)
plt.imshow(RecFISTA_os_reg[:, sliceSel, :], vmin=0, vmax=max_val)
plt.title("3D FISTA Reconstruction, coronal view")

plt.subplot(133)
plt.imshow(RecFISTA_os_reg[:, :, sliceSel], vmin=0, vmax=max_val)
plt.title("3D FISTA Reconstruction, sagittal view")
plt.show()
# %%
