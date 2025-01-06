#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to generate 2D analytical phantoms and their sinograms
If one needs to modify/add phantoms, please edit Phantom2DLibrary.dat

@author: Daniil Kazantsev
"""
import numpy as np
import timeit
import matplotlib.pyplot as plt
import os
import tomophantom
from tomophantom import TomoP2D
from tomophantom.qualitymetrics import QualityTools

model = 4  # select a model number from the library file (Phantom2DLibrary)
N_size = 512  # set the desired dimension of the phantom
path = os.path.dirname(tomophantom.__file__)
path_library2D = os.path.join(path, "phantomlib", "Phantom2DLibrary.dat")

# Generate a N_size x N_size phantom (2D)
phantom_2D = TomoP2D.Model(model, N_size, path_library2D)

plt.close("all")
plt.figure()
plt.rcParams.update({"font.size": 21})
plt.imshow(phantom_2D, vmin=0, vmax=1, cmap="BuPu")
plt.colorbar(ticks=[0, 0.5, 1], orientation="vertical")
plt.title("{}" "{}".format("2D Phantom using model no.", model))
plt.show()

# %%
# Parameters to generate a sinogram
angles_num = int(0.5 * np.pi * N_size)
# angles number
angles = np.linspace(0.0, 179.9, angles_num, dtype="float32")
P = int(np.sqrt(2) * N_size)  # detectors


tic = timeit.default_timer()

# create sinogram analytically
sino_an = TomoP2D.ModelSino(model, N_size, P, angles, path_library2D)
toc = timeit.default_timer()
Run_time = toc - tic
print("Analytical sinogram has been generated in {} seconds".format(Run_time))

plt.figure()
plt.rcParams.update({"font.size": 21})
plt.imshow(sino_an, vmin=0, vmax=300, cmap="BuPu")
plt.colorbar(ticks=[0, 150, 300], orientation="vertical")
plt.title("{}" "{}".format("Analytical sinogram of model no.", model))
plt.show()
# # %%
# # generate numerical sinogram
# angles_num = int(0.5 * np.pi * N_size)
# # angles number
# angles = np.linspace(0.0, 179.9, angles_num, dtype="float32")
# P = int(np.sqrt(2) * N_size)  # detectors

# tic = timeit.default_timer()
# sino_num = TomoP2D.SinoNum(phantom_2D, P, angles)
# toc = timeit.default_timer()
# Run_time = toc - tic
# print("Numerical sinogram has been generated in {} seconds".format(Run_time))

# plt.figure()
# plt.rcParams.update({"font.size": 21})
# plt.imshow(sino_num, vmin=0, vmax=300, cmap="BuPu")
# plt.colorbar(ticks=[0, 150, 300], orientation="vertical")
# plt.title("{}" "{}".format("Numerical sinogram of model no.", model))

# # plt.figure()
# # plt.imshow(abs(sino_an-sino_num), vmin=0, vmax=0.001, cmap="BuPu")
# # plt.colorbar(ticks=[0, 0.02, 0.05], orientation='vertical')
# # plt.title('Analytical vs Numerical singograms')
# # %%


# ###################################################################
# from tomobar.methodsDIR import RecToolsDIR

# angles_rad = angles * (np.pi / 180.0)

# # get numerical sinogram (ASTRA-toolbox)
# Rectools = RecToolsDIR(
#     DetectorsDimH=P,  # DetectorsDimH # detector dimension (horizontal)
#     DetectorsDimV=None,  # DetectorsDimV # detector dimension (vertical) for 3D case only
#     CenterRotOffset=0.0,  # Center of Rotation (CoR) scalar (for 3D case only)
#     AnglesVec=angles_rad,  # array of angles in radians
#     ObjSize=N_size,  # a scalar to define reconstructed object dimensions
#     device_projector="cpu",
# )

# tic = timeit.default_timer()
# sino_num_ASTRA = Rectools.FORWPROJ(phantom_2D)  # generate numerical sino (Ax)
# toc = timeit.default_timer()
# Run_time = toc - tic
# print("Numerical (ASTRA) sinogram has been generated in {} seconds".format(Run_time))

# plt.figure()
# plt.rcParams.update({"font.size": 21})
# plt.imshow(sino_num_ASTRA, vmin=0, vmax=150, cmap="BuPu")
# plt.colorbar(ticks=[0, 150, 250], orientation="vertical")
# plt.title("{}" "{}".format("Numerical sinogram (ASTRA) of model no.", model))
# # %%
# print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
# print("Reconstructing analytical sinogram using Fourier Slice method")
# print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

# RecFourier = Rectools.FOURIER(sino_an, "linear")

# plt.figure()
# plt.imshow(RecFourier, vmin=0, vmax=1, cmap="BuPu")
# plt.colorbar(ticks=[0, 0.5, 1], orientation="vertical")
# plt.title("Fourier slice reconstruction")
# # %%
# print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
# print("Reconstructing analytical sinogram using FBP (tomobar)...")
# print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
# # x = Atools.backproj(sino_an) # generate backprojection (A'b)

# plt.figure()
# plt.subplot(121)
# plt.imshow(sino_an, cmap="BuPu")
# plt.title("Analytical sinogram")
# plt.subplot(122)
# plt.imshow(sino_num_ASTRA, cmap="BuPu")
# plt.title("Numerical sinogram")
# plt.show()
# # calculate norm
# # rmse1 = np.linalg.norm(sino_an - sino_num_ASTRA)/np.linalg.norm(sino_num_ASTRA)

# print("Reconstructing analytical sinogram using FBP (astra TB)...")
# FBPrec1 = Rectools.FBP(sino_an)
# plt.figure()
# plt.imshow(FBPrec1, vmin=0, vmax=1, cmap="BuPu")
# plt.colorbar(ticks=[0, 0.5, 1], orientation="vertical")
# plt.title("FBP Reconstructed Phantom (analyt)")

# print("Reconstructing numerical sinogram using FBP (astra TB)...")
# FBPrec2 = Rectools.FBP(sino_num_ASTRA)

# plt.figure()
# plt.imshow(FBPrec2, vmin=0, vmax=1, cmap="BuPu")
# plt.colorbar(ticks=[0, 0.5, 1], orientation="vertical")
# plt.title("FBP Reconstructed Phantom (numeric)")

# plt.figure()
# plt.imshow(abs(FBPrec1 - FBPrec2), vmin=0, vmax=0.05, cmap="BuPu")
# plt.colorbar(ticks=[0, 0.02, 0.05], orientation="vertical")
# plt.title("FBP rec differences")
# # rmse2 = np.linalg.norm(FBPrec1 - FBPrec2)/np.linalg.norm(FBPrec2)

# Qtools = QualityTools(phantom_2D, FBPrec1)
# RMSE_FBP1 = Qtools.rmse()
# Qtools = QualityTools(phantom_2D, FBPrec2)
# RMSE_FBP2 = Qtools.rmse()
# print("RMSE for FBP (analyt) {}".format(RMSE_FBP1))
# print("RMSE for FBP (numeric) {}".format(RMSE_FBP2))
