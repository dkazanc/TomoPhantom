#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPLv3 license (ASTRA toolbox)
Note that the TomoPhantom package is released under Apache License, Version 2.0

Script to generate 2D analytical phantoms and their sinograms with added noise and artifacts
Sinograms then reconstructed using ASTRA TOOLBOX 

>>>>> Dependencies (reconstruction): <<<<<
1. ASTRA toolbox: conda install -c astra-toolbox astra-toolbox
2. tomobar: conda install -c dkazanc tomobar
or install from https://github.com/dkazanc/ToMoBAR

This demo demonstrates frequent inaccuracies which are accosiated with X-ray imaging:
zingers, rings and noise

@author: Daniil Kazantsev
"""
import numpy as np
import matplotlib.pyplot as plt
from tomophantom import TomoP2D
import os
import tomophantom
from tomophantom.supp.qualitymetrics import QualityTools

model = 15 # select a model
N_size = 512 # set dimension of the phantom
# one can specify an exact path to the parameters file
# path_library2D = '../../../PhantomLibrary/models/Phantom2DLibrary.dat'
path = os.path.dirname(tomophantom.__file__)
path_library2D = os.path.join(path, "Phantom2DLibrary.dat")
phantom_2D = TomoP2D.Model(model, N_size, path_library2D)

plt.close('all')
plt.figure(1)
plt.rcParams.update({'font.size': 21})
plt.imshow(phantom_2D, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('{}''{}'.format('2D Phantom using model no.',model))

# create sinogram analytically
angles_num = int(0.5*np.pi*N_size); # angles number
angles = np.linspace(0.0,179.9,angles_num,dtype='float32')
angles_rad = angles*(np.pi/180.0)
P = N_size #detectors

sino_an = TomoP2D.ModelSino(model, N_size, P, angles, path_library2D)

plt.figure(2)
plt.rcParams.update({'font.size': 21})
plt.imshow(sino_an,  cmap="gray")
plt.colorbar(ticks=[0, 150, 250], orientation='vertical')
plt.title('{}''{}'.format('Analytical sinogram of model no.',model))

#%%
# Adding artifacts and noise
from tomophantom.supp.artifacts import _Artifacts_

# forming dictionaries with artifact types
_noise_ =  {'type' : 'Poisson',
            'sigma' : 10000, # noise amplitude
            'seed' : 0}
# misalignment dictionary
_sinoshifts_ = {'maxamplitude' : 10}
[noisy_sino_misalign,shifts] = _Artifacts_(sino_an, _noise_, {}, {}, _sinoshifts_)

# adding zingers and stripes
_zingers_ = {'percentage' : 0.25,
             'modulus' : 10}

_stripes_ = {'percentage' : 0.8,
             'maxthickness' : 2.0,
             'intensity' : 0.25,
             'type' : 'full',
             'variability' : 0.002}

noisy_zing_stripe = _Artifacts_(sino_an, _noise_, _zingers_, _stripes_, _sinoshifts_= {})

plt.figure()
plt.rcParams.update({'font.size': 21})
plt.imshow(noisy_zing_stripe,cmap="gray")
plt.colorbar(ticks=[0, 150, 250], orientation='vertical')
plt.title('{}''{}'.format('Analytical noisy sinogram with artifacts.',model))
#%%
# initialise tomobar DIRECT reconstruction class ONCE
from tomobar.methodsDIR import RecToolsDIR
RectoolsDIR = RecToolsDIR(DetectorsDimH = P,  # DetectorsDimH # detector dimension (horizontal)
                    DetectorsDimV = None,  # DetectorsDimV # detector dimension (vertical) for 3D case only
                    CenterRotOffset = None, # Center of Rotation (CoR) scalar (for 3D case only)
                    AnglesVec = angles_rad, # array of angles in radians
                    ObjSize = N_size, # a scalar to define reconstructed object dimensions
                    device_projector='cpu')
#%%
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print ("Reconstructing analytical sinogram using Fourier Slice method")
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
RecFourier = RectoolsDIR.FOURIER(sino_an,'linear') 
plt.figure() 
plt.imshow(RecFourier, vmin=0, vmax=1, cmap="BuPu")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('Fourier slice reconstruction')
#%%
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print ("Reconstructing analytical sinogram using FBP (tomobar)...")
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
FBPrec_ideal = RectoolsDIR.FBP(sino_an)  # ideal reconstruction
FBPrec_error = RectoolsDIR.FBP(noisy_zing_stripe) # reconstruction with artifacts
FBPrec_misalign = RectoolsDIR.FBP(noisy_sino_misalign) # reconstruction with misalignment

plt.figure()
plt.subplot(131)
plt.imshow(FBPrec_ideal, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('Ideal FBP reconstruction')
plt.subplot(132)
plt.imshow(FBPrec_error, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('Erroneous data FBP Reconstruction')
plt.subplot(133)
plt.imshow(FBPrec_misalign, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('Misaligned noisy FBP Reconstruction')
plt.show()

plt.figure() 
plt.imshow(abs(FBPrec_ideal-FBPrec_error), vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('FBP reconsrtuction differences')
#%%
# initialise tomobar ITERATIVE reconstruction class ONCE
from tomobar.methodsIR import RecToolsIR
RectoolsIR = RecToolsIR(DetectorsDimH = P,  # DetectorsDimH # detector dimension (horizontal)
                    DetectorsDimV = None,  # DetectorsDimV # detector dimension (vertical) for 3D case only
                    CenterRotOffset = None, # Center of Rotation (CoR) scalar (for 3D case only)
                    AnglesVec = angles_rad, # array of angles in radians
                    ObjSize = N_size, # a scalar to define reconstructed object dimensions
                    datafidelity='LS',# data fidelity, choose LS, PWLS (wip), GH (wip), Student (wip)
                    device_projector='gpu')
#%%
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print ("Reconstructing analytical sinogram using SIRT (tomobar)...")
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
# prepare dictionaries with parameters:
_data_ = {'projection_norm_data' : sino_an} # data dictionary
_algorithm_ = {'iterations' : 250}
SIRTrec_ideal = RectoolsIR.SIRT(_data_,_algorithm_) # ideal reconstruction
_data_ = {'projection_norm_data' : noisy_zing_stripe} # data dictionary
SIRTrec_error = RectoolsIR.SIRT(_data_,_algorithm_) # error reconstruction

plt.figure()
plt.subplot(121)
plt.imshow(SIRTrec_ideal, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('Ideal SIRT reconstruction (ASTRA)')
plt.subplot(122)
plt.imshow(SIRTrec_error, vmin=0, vmax=3, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('Erroneous data SIRT Reconstruction (ASTRA)')
plt.show()

plt.figure() 
plt.imshow(abs(SIRTrec_ideal-SIRTrec_error), vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('SIRT reconsrtuction differences')
#%%
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print ("Reconstructing using ADMM method (tomobar)")
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
# Run ADMM reconstrucion algorithm with regularisation 
_data_ = {'projection_norm_data' : noisy_zing_stripe} # data dictionary
_algorithm_ = {'iterations' : 15,
               'ADMM_rho_const' : 5000.0}

# adding regularisation using the CCPi regularisation toolkit
_regularisation_ = {'method' : 'PD_TV',
                    'regul_param' : 0.1,
                    'iterations' : 150,
                    'device_regulariser': 'gpu'}

RecADMM_reg = RectoolsIR.ADMM(_data_, _algorithm_, _regularisation_)

plt.figure()
plt.imshow(RecADMM_reg, vmin=0, vmax=2, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 3], orientation='vertical')
plt.title('ADMM reconstruction')
plt.show()

# calculate errors 
Qtools = QualityTools(phantom_2D, RecADMM_reg)
RMSE_ADMM_reg = Qtools.rmse()
print("RMSE for regularised ADMM is {}".format(RMSE_ADMM_reg))
#%%
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print ("Reconstructing using FISTA method (tomobar)")
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
# prepare dictionaries with parameters:
_data_ = {'projection_norm_data' : noisy_zing_stripe,
          'huber_threshold' : 4.5,
          'ring_weights_threshold' : 10.0,
          'ring_tuple_halfsizes': (9,5,0)
          } # data dictionary

lc = RectoolsIR.powermethod(_data_) # calculate Lipschitz constant (run once to initialise)
_algorithm_ = {'iterations' : 350,
               'lipschitz_const' : lc}

# adding regularisation using the CCPi regularisation toolkit
_regularisation_ = {'method' : 'PD_TV',
                    'regul_param' : 0.001,
                    'iterations' : 150,
                    'device_regulariser': 'gpu'}

# Run FISTA reconstrucion algorithm with regularisation 
RecFISTA_reg = RectoolsIR.FISTA(_data_, _algorithm_, _regularisation_)

plt.figure()
plt.imshow(RecFISTA_reg, vmin=0, vmax=2, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 2], orientation='vertical')
plt.title('TV-Regularised FISTA reconstruction')
plt.show()

# calculate errors 
Qtools = QualityTools(phantom_2D, RecFISTA_reg)
RMSE_FISTA_reg = Qtools.rmse()
print("RMSE for regularised FISTA is {}".format(RMSE_FISTA_reg))
#%%