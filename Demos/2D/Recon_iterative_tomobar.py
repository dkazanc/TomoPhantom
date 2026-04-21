#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPLv3 license (ASTRA toolbox)
Note that the TomoPhantom package is released under Apache License, Version 2.0

Script to generate 2D analytical phantoms and their sinograms with added noise and artifacts 
then reconstructed using ToMoBAR

>>>>> Dependencies (reconstruction): <<<<<
1. ASTRA toolbox
2. ToMoBAR
3. CuPy

@author: Daniil Kazantsev
"""
import numpy as np
import matplotlib.pyplot as plt
from tomophantom import TomoP2D
import os
import tomophantom
from tomophantom.qualitymetrics import QualityTools

model = 15  # select a model
N_size = 256  # set dimension of the phantom
path = os.path.dirname(tomophantom.__file__)
path_library2D = os.path.join(path, "phantomlib", "Phantom2DLibrary.dat")
phantom_2D = TomoP2D.Model(model, N_size, path_library2D)

plt.close("all")
plt.figure(1)
plt.rcParams.update({"font.size": 21})
plt.imshow(phantom_2D, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation="vertical")
plt.title("{}" "{}".format("2D Phantom using model no.", model))

# create sinogram analytically
angles_num = int(0.5 * np.pi * N_size)
# angles number
angles = np.linspace(0.0, 179.9, angles_num, dtype="float32")
angles_rad = angles * (np.pi / 180.0)
P = N_size  # detectors

sino_an = TomoP2D.ModelSino(model, N_size, P, angles, path_library2D)

plt.figure(2)
plt.rcParams.update({"font.size": 21})
plt.imshow(sino_an, cmap="gray")
plt.colorbar(ticks=[0, 150, 250], orientation="vertical")
plt.title("{}" "{}".format("Analytical sinogram of model no.", model))

# %%
# Adding artefacts and noise
from tomophantom.artefacts import artefacts_mix

plt.close("all")
# forming dictionaries with artifact types
_noise_ = {"noise_type": "Poisson", "noise_amplitude": 10000}
# misalignment dictionary
_sinoshifts_ = {"datashifts_maxamplitude_pixel": 10}
[noisy_sino_misalign, shifts] = artefacts_mix(sino_an, **_noise_, **_sinoshifts_)

# adding zingers and stripes
_zingers_ = {"zingers_percentage": 2, "zingers_modulus": 10}

_stripes_ = {
    "stripes_percentage": 0.8,
    "stripes_maxthickness": 2,
    "stripes_intensity": 0.25,
    "stripes_type": "full",
    "stripes_variability": 0.002,
}

noisy_zing_stripe = artefacts_mix(sino_an, **_noise_, **_zingers_, **_stripes_)

plt.figure()
plt.rcParams.update({"font.size": 21})
plt.imshow(noisy_zing_stripe, cmap="gray")
plt.colorbar(ticks=[0, 150, 250], orientation="vertical")
plt.title("{}" "{}".format("Analytical noisy sinogram with artefacts.", model))
# %%
# initialise tomobar DIRECT reconstruction class ONCE
from tomobar.methodsDIR import RecToolsDIR

RectoolsDIR = RecToolsDIR(
    DetectorsDimH=P,  # DetectorsDimH # detector dimension (horizontal)
    DetectorsDimH_pad=0,  # Padding size of horizontal detector
    DetectorsDimV=None,  # DetectorsDimV # detector dimension (vertical) for 3D case only
    CenterRotOffset=None,  # Center of Rotation (CoR) scalar (for 3D case only)
    AnglesVec=angles_rad,  # array of angles in radians
    ObjSize=N_size,  # a scalar to define reconstructed object dimensions
    device_projector="cpu",
)
# %%
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print("Reconstructing analytical sinogram using Fourier Slice method")
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
RecFourier = RectoolsDIR.FOURIER(sino_an, method="linear")
plt.figure()
plt.imshow(RecFourier, vmin=0, vmax=1, cmap="BuPu")
plt.colorbar(ticks=[0, 0.5, 1], orientation="vertical")
plt.title("Fourier slice reconstruction")
# %%
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print("Reconstructing analytical sinogram using FBP (tomobar)...")
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
FBPrec_ideal = RectoolsDIR.FBP(sino_an)  # ideal reconstruction
FBPrec_error = RectoolsDIR.FBP(noisy_zing_stripe)  # reconstruction with artifacts
FBPrec_misalign = RectoolsDIR.FBP(
    noisy_sino_misalign
)  # reconstruction with misalignment

plt.figure()
plt.subplot(131)
plt.imshow(FBPrec_ideal, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation="vertical")
plt.title("Ideal FBP reconstruction")
plt.subplot(132)
plt.imshow(FBPrec_error, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation="vertical")
plt.title("Erroneous data FBP Reconstruction")
plt.subplot(133)
plt.imshow(FBPrec_misalign, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation="vertical")
plt.title("Misaligned noisy FBP Reconstruction")
plt.show()

plt.figure()
plt.imshow(abs(FBPrec_ideal - FBPrec_error), vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation="vertical")
plt.title("FBP reconsrtuction differences")
# %%
# ! you will need to have CuPy installed to be able to run iterative methods
import cupy as cp 
from tomobar.methodsIR_CuPy import RecToolsIRCuPy
input_data_labels = ["angles", "detX"]

Rectools = RecToolsIRCuPy(
    DetectorsDimH=P,  # DetectorsDimH # detector dimension (horizontal)
    DetectorsDimH_pad=0,  # Padding size of horizontal detector
    DetectorsDimV=None,  # DetectorsDimV # detector dimension (vertical) for 3D case only
    CenterRotOffset=0.0,  # Center of Rotation (CoR) scalar (for 3D case only)
    AnglesVec=angles_rad,  # array of angles in radians
    ObjSize=N_size,  # a scalar to define reconstructed object dimensions
    device_projector=0,
    OS_number=1,  # The number of ordered subsets
)

# %%
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print("Reconstructing analytical sinogram using SIRT (tomobar)...")
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
# prepare dictionaries with parameters:
_data_ = {"projection_data": cp.asarray(sino_an),
          "data_axes_labels_order": input_data_labels}  # data dictionary
_algorithm_ = {"iterations": 250}

SIRTrec_ideal = cp.asnumpy(Rectools.SIRT(_data_, _algorithm_))  # ideal reconstruction


_data_ = {"projection_data": cp.asarray(noisy_zing_stripe),
          "data_axes_labels_order": input_data_labels}  # data dictionary

SIRTrec_error = cp.asnumpy(Rectools.SIRT(_data_, _algorithm_))  # error reconstruction

plt.figure()
plt.subplot(121)
plt.imshow(SIRTrec_ideal[0,:,:], vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation="vertical")
plt.title("Ideal SIRT reconstruction (ASTRA)")
plt.subplot(122)
plt.imshow(SIRTrec_error[0,:,:], vmin=0, vmax=3, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation="vertical")
plt.title("Erroneous data SIRT Reconstruction (ASTRA)")
plt.show()

plt.figure()
plt.imshow(abs(SIRTrec_ideal[0,:,:] - SIRTrec_error[0,:,:]), vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation="vertical")
plt.title("SIRT reconsrtuction differences")
# %%
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print("Reconstructing using FISTA-OS-TV method (tomobar)")
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
Rectools = RecToolsIRCuPy(
    DetectorsDimH=P,  # DetectorsDimH # detector dimension (horizontal)
    DetectorsDimH_pad=0,  # Padding size of horizontal detector
    DetectorsDimV=None,  # DetectorsDimV # detector dimension (vertical) for 3D case only
    CenterRotOffset=0.0,  # Center of Rotation (CoR) scalar (for 3D case only)
    AnglesVec=angles_rad,  # array of angles in radians
    ObjSize=N_size,  # a scalar to define reconstructed object dimensions
    device_projector=0,
    OS_number=10,  # The number of ordered subsets
)

# prepare dictionaries with parameters:
_data_ = {
    "projection_data": cp.asarray(noisy_zing_stripe),
    "data_axes_labels_order": input_data_labels
}  # data dictionary

lc = Rectools.powermethod(
    _data_
)  # calculate Lipschitz constant (run once to initialise)
_algorithm_ = {"iterations": 15, "lipschitz_const": lc}

_regularisation_ = {
    "method": "PD_TV",  # Selected regularisation method
    "regul_param": 0.0005,  # Regularisation parameter
    "iterations": 40,  # The number of regularisation iterations
    "half_precision": True,  # enabling half-precision calculation
}


# Run FISTA reconstrucion algorithm with regularisation
RecFISTA_reg = Rectools.FISTA(_data_, _algorithm_, _regularisation_)


RecFISTA_reg = cp.asnumpy(RecFISTA_reg[0, :, :])
plt.figure()
plt.imshow(RecFISTA_reg, vmin=0, vmax=2, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 2], orientation="vertical")
plt.title("TV-Regularised FISTA reconstruction")
plt.show()

# calculate errors
Qtools = QualityTools(phantom_2D, RecFISTA_reg)
RMSE_FISTA_reg = Qtools.rmse()
print("RMSE for regularised FISTA is {}".format(RMSE_FISTA_reg))
# %%
from tomophantom.artefacts import artefacts_mix

# forming dictionaries with artefact types
_noise_ = {
    "noise_type": "Poisson",
    "noise_sigma": 200000,  # noise amplitude
    "noise_seed": 0,
}

# partial volume effect dictionary
_pve_ = {"pve_strength": 1}
_fresnel_propagator_ = {
    "fresnel_dist_observation": 10,
    "fresnel_scale_factor": 10,
    "fresnel_wavelenght": 0.003,
}

noisy_sino_pve = artefacts_mix(sino_an, **_noise_, **_pve_)

FBPrec_pve = RectoolsDIR.FBP(noisy_sino_pve)

plt.figure()
plt.imshow(FBPrec_pve, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation="vertical")
plt.title("FBP reconstruction from PVE sinogram")
plt.show()
# %%
from tomophantom.artefacts import artefacts_mix

# forming dictionaries with artefact types
_noise_ = {
    "noise_type": "Poisson",
    "noise_amplitude": 200000,  # noise amplitude
    "noise_seed": 0,
}
_fresnel_propagator_ = {
    "fresnel_dist_observation": 20,
    "fresnel_scale_factor": 10,
    "fresnel_wavelenght": 0.003,
}

noisy_sino_fresnel = artefacts_mix(sino_an, **_noise_, **_fresnel_propagator_)

FBPrec_fresnel = RectoolsDIR.FBP(noisy_sino_fresnel)

plt.figure()
plt.imshow(FBPrec_fresnel, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation="vertical")
plt.title("FBP reconstruction from sinogram with Fresnel propagator")
plt.show()
# %%
