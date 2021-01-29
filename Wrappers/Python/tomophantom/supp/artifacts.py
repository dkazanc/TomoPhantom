#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#    Note that the TomoPhantom package is released under Apache License, Version 2.0
#    @author: Daniil Kazantsev
#    latest update: 26.01.2021

import numpy as np
import random

def _Artifacts_(data, **artefacts_dict):
    """
    The following keys in dictionaries can be provided to simulate noise and artefacts for 2D and 3D-projection data:
    - ['noise_type']: 'Poisson' or 'Gaussian'
    - ['noise_amplitude']: photon flux (Poisson) or variance for Gaussian
    - ['noise_seed']: seeds for noise (integers), None means random generation
    - ['noise_prelog']: set to True if  prelog (raw) data required
    - ['zingers_percentage']: the amount of zingers to be added to the data (streaks in reconstruction)
    - ['zingers_modulus']: modulus to control the amount of 4/6 pixel clusters to be added
    - ['stripes_percentage']: the amount of stripes in the data (rings in reconstruction)
    - ['stripes_maxthickness']: maxthickness defines the maximal thickness of a stripe
    - ['stripes_intensity']: controls the intensity levels of stripes
    - ['stripe_type']: stripe types can be 'partial' or 'full'
    - ['stripes_variability']: variability multiplier to incorporate the change of intensity in stripe
    - ['sinoshifts_maxamplitude']: misalignment (in int pixels) defines the maximal amplitude of each angular deviation
    - ['pve_strength']: the strength of partial volume effect, responsible for the limited resolution
    """
    ####### NOISE DICTIONARY########
    _noise_ = {}
    if ('noise_type' not in artefacts_dict):
        _noise_['noise_type'] = None
    else:
        _noise_['noise_type'] = artefacts_dict['noise_type']
    if ('noise_amplitude' not in artefacts_dict):
        _noise_['noise_amplitude'] = 10000
    else:
         _noise_['noise_amplitude'] = artefacts_dict['noise_amplitude']
    if ('noise_seed' not in artefacts_dict):
        _noise_['noise_seed'] = None
    else:
        _noise_['noise_seed'] = artefacts_dict['noise_seed']
    if ('noise_prelog' not in artefacts_dict):
        _noise_['noise_prelog'] = None
    else:
        _noise_['noise_prelog'] = artefacts_dict['noise_prelog']
    ####### ZINGERS ########
    _zingers_ = {}
    if ('zingers_percentage' not in artefacts_dict):
        _zingers_['zingers_percentage'] = None
    else:
        _zingers_['zingers_percentage'] = artefacts_dict['zingers_percentage']
    if ('zingers_modulus' not in artefacts_dict):
        _zingers_['zingers_modulus'] = 10
    else:
        _zingers_['zingers_modulus'] = artefacts_dict['zingers_modulus']
    ####### STRIPES ########
    _stripes_ = {}
    if ('stripes_percentage' not in artefacts_dict):
        _stripes_['stripes_percentage'] = None
    else:
        _stripes_['stripes_percentage'] = artefacts_dict['stripes_percentage']
    if ('stripes_maxthickness' not in artefacts_dict):
        _stripes_['stripes_maxthickness'] = 1.0
    else:
        _stripes_['stripes_maxthickness'] = artefacts_dict['stripes_maxthickness']
    if ('stripes_intensity' not in artefacts_dict):
        _stripes_['stripes_intensity'] = 0.1
    else:
        _stripes_['stripes_intensity'] = artefacts_dict['stripes_intensity']
    if ('stripes_type' not in artefacts_dict):
        _stripes_['stripes_type'] = 'full'
    else:
        _stripes_['stripes_type'] = artefacts_dict['stripes_type']
    if ('stripes_variability' not in artefacts_dict):
        _stripes_['stripes_variability'] = 0.0
    else:
        _stripes_['stripes_variability'] = artefacts_dict['stripes_variability']
    ####### SINOSHIFTS ########
    _sinoshifts_ = {}
    if ('sinoshifts_maxamplitude' not in artefacts_dict):
        _sinoshifts_['sinoshifts_maxamplitude'] = None
    else:
        _sinoshifts_['sinoshifts_maxamplitude'] = artefacts_dict['sinoshifts_maxamplitude']
    ####### PVE ########
    _pve_ = {}
    if ('pve_strength' not in artefacts_dict):
        _pve_['pve_strength'] = None
    else:
        _pve_['pve_strength'] = artefacts_dict['pve_strength']
    _fresnel_propagator_ = {}
    if ('fresnel_dist_observation' not in artefacts_dict):
        _fresnel_propagator_['fresnel_dist_observation'] = None
    else:
        _fresnel_propagator_['fresnel_dist_observation'] = artefacts_dict['fresnel_dist_observation']
    if ('fresnel_scale_factor' not in artefacts_dict):
        _fresnel_propagator_['fresnel_scale_factor'] = 10
    else:
        _fresnel_propagator_['fresnel_scale_factor'] = artefacts_dict['fresnel_scale_factor']
    if ('fresnel_wavelenght' not in artefacts_dict):
        _fresnel_propagator_['fresnel_wavelenght'] = 0.0001
    else:
        _fresnel_propagator_['fresnel_wavelenght'] = artefacts_dict['fresnel_wavelenght']
    ###########################################################################
    ################Applying artefacts and noise to the data###################
    ###########################################################################
    # PARTIAL VOLUME EFFECT
    if _pve_['pve_strength'] is not None:
        sino_artifacts = pve(data=data, pve_strength=_pve_['pve_strength'])
        print("Partial volume effect (PVE) have been simulated.")
    else:
        sino_artifacts = data
    # FRESNEL PROPAGATOR
    if _fresnel_propagator_['fresnel_dist_observation'] is not None:
        sino_artifacts = fresnel_propagator(data=sino_artifacts, dist_observation=_fresnel_propagator_['fresnel_dist_observation'], scale_factor = _fresnel_propagator_['fresnel_scale_factor'], wavelenght = _fresnel_propagator_['fresnel_wavelenght'])
        print("Fresnel propagator have been simulated.")
    # ZINGERS
    if _zingers_['zingers_percentage'] is not None:
        sino_artifacts = zingers(data=sino_artifacts, percentage=_zingers_['zingers_percentage'], modulus = _zingers_['zingers_modulus'])
        print("Zingers have been added to the data.")
    # STRIPES
    if _stripes_['stripes_percentage'] is not None:
        sino_artifacts = stripes(data=sino_artifacts,\
                                 percentage=_stripes_['stripes_percentage'],\
                                 maxthickness = _stripes_['stripes_maxthickness'],\
                                 intensity_thresh = _stripes_['stripes_intensity'],\
                                 stripe_type = _stripes_['stripes_type'],\
                                 variability = _stripes_['stripes_variability'])
        print("Stripes have been simulated.")
    # SINOSHIFTS
    if _sinoshifts_['sinoshifts_maxamplitude'] is not None:
        [sino_artifacts, shifts] = sinoshifts(data=sino_artifacts, maxamplitude = _sinoshifts_['sinoshifts_maxamplitude'])
        print("Sinogram shifts have been simulated.")
    # NOISE
    if _noise_['noise_type']  is not None:
        sino_artifacts = noise(data=sino_artifacts, sigma=_noise_['noise_amplitude'], noisetype = _noise_['noise_type'], seed = _noise_['noise_seed'], prelog=_noise_['noise_prelog'])
        print("{} noise have been added to the data.".format(_noise_['noise_type']))

    if _sinoshifts_['sinoshifts_maxamplitude'] is not None:
        return [sino_artifacts,shifts]
    else:
        return sino_artifacts

def stripes(data, percentage, maxthickness, intensity_thresh, stripe_type, variability):
    # Function to add stripes (constant offsets) to sinogram which results in rings in the
    # reconstructed image
    # - percentage defines the amount of stripes in the data
    # - maxthickness defines the maximal thickness of a stripe
    # - stripe_type can be 'partial' or 'full'
    # - variability multiplier to incorporate change of intensity in the stripe
    if (data.ndim == 2):
        (anglesDim, DetectorsDimH) = np.shape(data)
    else:
        (DetectorsDimV, anglesDim, DetectorsDimH) = np.shape(data)
    if 0 < percentage <= 100:
        pass
    else:
        raise ("percentage must be larger than zero but smaller than 100")
    if 0 <= maxthickness <= 10:
        pass
    else:
        raise ("maximum thickness must be in [0,10] range")
    if ((stripe_type != 'partial')):
        stripe_type = 'full'
    sino_stripes = data.copy()
    max_intensity = np.max(sino_stripes)
    range_detect = int((np.float32(DetectorsDimH))*(np.float32(percentage)/100.0))
    if (data.ndim == 2):
        for x in range(range_detect):
            for mm in range(0,20):
                randind = random.randint(0,DetectorsDimH-1) # generate random index
                if (sino_stripes[0,randind] != 0.0):
                    break
            if (stripe_type == 'partial'):
                randind_ang1 =  random.randint(0,anglesDim)
                randind_ang2 =  random.randint(0,anglesDim)
            else:
                randind_ang1 = 0
                randind_ang2 = anglesDim
            randthickness = random.randint(0,maxthickness) #generate random thickness
            randintens = random.uniform(-1.0, 0.5) # generate random multiplier
            intensity = max_intensity*randintens*intensity_thresh
            if ((randind > 0+randthickness) & (randind < DetectorsDimH-randthickness)):
                for x1 in range(-randthickness,randthickness+1):
                    if (variability != 0.0):
                        intensity_off = variability*max_intensity
                    else:
                        intensity_off = 0.0
                    for ll in range(randind_ang1,randind_ang2):
                        sino_stripes[ll,randind+x1] += intensity + intensity_off
                        intensity_off += [-1,1][random.randrange(2)]*variability*max_intensity
                        #sino_stripes[randind_ang1:randind_ang2,randind+x1] += intensity
    else:
        for j in range(DetectorsDimV):
            for x in range(range_detect):
                for mm in range(0,20):
                    randind = random.randint(0,DetectorsDimH-1) # generate random index
                    if (sino_stripes[j,0,randind] != 0.0):
                        break
                if (stripe_type == 'partial'):
                    randind_ang1 =  random.randint(0,anglesDim)
                    randind_ang2 =  random.randint(0,anglesDim)
                else:
                    randind_ang1 = 0
                    randind_ang2 = anglesDim
                randthickness = random.randint(0,maxthickness) #generate random thickness
                randintens = random.uniform(-1, 0.5) # generate random multiplier
                intensity = max_intensity*randintens*intensity_thresh
                if ((randind > 0+randthickness) & (randind < DetectorsDimH-randthickness)):
                    for x1 in range(-randthickness,randthickness+1):
                        if (variability != 0.0):
                            intensity_off = variability*max_intensity
                        else:
                            intensity_off = 0.0
                        for ll in range(randind_ang1,randind_ang2):
                            sino_stripes[j,ll,randind+x1] += intensity + intensity_off
                            intensity_off += [-1,1][random.randrange(2)]*variability*max_intensity
    return sino_stripes

def zingers(data, percentage, modulus):
    # adding zingers (zero single pixels or small 4 pixels clusters) to data or 6 voxels to 3D projection data
    # percentage = 0.5, #- percentage - the amount of zingers to be added to the data
    # modulus = 10 # modulus to control the amount of 4/6 pixel clusters to be added
    if (data.ndim == 2):
        (anglesDim, DetectorsDimH) = np.shape(data)
    else:
        (DetectorsDimV, anglesDim, DetectorsDimH) = np.shape(data)
    if 0.0 < percentage <= 100.0:
        pass
    else:
        raise ("percentage must be larger than zero but smaller than 100")
    if (modulus > 0):
        pass
    else:
        raise ("Modulus integer must be positive")
    sino_zingers = data.copy()
    length_sino = np.size(sino_zingers)
    num_values = int((length_sino)*(np.float32(percentage)/100.0))
    sino_zingers_fl = sino_zingers.flatten()
    for x in range(num_values):
        randind = random.randint(0,length_sino-1) # generate random index
        sino_zingers_fl[randind] = 0
        if ((x % int(modulus)) == 0):
            if (data.ndim == 2):
                if ((randind > DetectorsDimH) & (randind < length_sino-DetectorsDimH)):
                    sino_zingers_fl[randind+1] = 0
                    sino_zingers_fl[randind-1] = 0
                    sino_zingers_fl[randind+DetectorsDimH] = 0
                    sino_zingers_fl[randind-DetectorsDimH] = 0
            else:
                if ((randind > DetectorsDimH*DetectorsDimV) & (randind < length_sino-DetectorsDimH*DetectorsDimV)):
                    sino_zingers_fl[randind+1] = 0
                    sino_zingers_fl[randind-1] = 0
                    sino_zingers_fl[randind+DetectorsDimH] = 0
                    sino_zingers_fl[randind-DetectorsDimH] = 0
                    sino_zingers_fl[randind+DetectorsDimH*DetectorsDimV] = 0
                    sino_zingers_fl[randind-DetectorsDimH*DetectorsDimV] = 0
    sino_zingers[:] = sino_zingers_fl.reshape(sino_zingers.shape)
    return sino_zingers

def noise(data, sigma, noisetype, seed, prelog):
    # Adding random noise to data
    # noisetype = None, # 'Gaussian', 'Poisson' or None
    # sigma = 10000, # photon intensity (Poisson) or variance for Gaussian
    # seed = 0, # seeds for noise
    # prelog: None or True (get the raw pre-log data)
    sino_noisy = data.copy()
    if noisetype == 'Gaussian':
        # add normal Gaussian noise
        if seed is not None:
            np.random.seed(int(seed))
        sino_noisy += np.random.normal(loc = 0.0, scale = sigma, size = np.shape(sino_noisy))
        sino_noisy[sino_noisy<0] = 0
    elif noisetype == 'Poisson':
        # add Poisson noise
        maxSino = np.max(data)
        if maxSino > 0:
            sino_noisy = data/maxSino
            dataExp = sigma*np.exp(-sino_noisy)  # noiseless raw data
            sino_raw = np.random.poisson(dataExp) #adding Poisson noise
            div_res = np.float32(sino_raw)/np.max(sino_raw)
            sino_noisy = -np.log(div_res)*maxSino # log corrected data
            sino_noisy[sino_noisy<0] = 0
    else:
        print ("Select 'Gaussian' or 'Poisson' for noise type")
    if prelog is True:
        return [sino_noisy,sino_raw]
    else:
        return sino_noisy

def sinoshifts(data, maxamplitude):
    # A function to add random shifts to sinogram rows (an offset for each angular position)
    # maxamplitude (integer pixels) defines the maximal amplitude of each angular deviation
    if (data.ndim == 2):
        (anglesDim, DetectorsDimH) = np.shape(data)
        shifts = np.zeros(anglesDim, dtype='int8') # the vector of shifts
    else:
        (DetectorsDimV, anglesDim, DetectorsDimH) = np.shape(data)
        shifts = np.zeros([anglesDim,2], dtype='int8') # the 2D vector of shifts

    sino_shifts = np.zeros(np.shape(data),dtype='float32')
    non = lambda s: s if s<0 else None
    mom = lambda s: max(0,s)
    for x in range(anglesDim):
        rand_shift = random.randint(-maxamplitude,maxamplitude)  #generate random shift
        if (data.ndim == 2):
            shifts[x] = rand_shift
            projection = data[x,:] # extract 1D projection
            projection_shift = np.zeros(np.shape(projection),dtype='float32')
            projection_shift[mom(rand_shift):non(rand_shift)] = projection[mom(-rand_shift):non(-rand_shift)]
            sino_shifts[x,:] = projection_shift
        else:
            rand_shift2 = random.randint(-maxamplitude,maxamplitude)  #generate random shift
            shifts[x,0] = rand_shift
            shifts[x,1] = rand_shift2
            projection2D = data[:,x,:] # extract 2D projection
            projection2D_shift = np.zeros(np.shape(projection2D),dtype='float32')
            projection2D_shift[mom(rand_shift):non(rand_shift), mom(rand_shift2):non(rand_shift2)] = projection2D[mom(-rand_shift):non(-rand_shift),mom(-rand_shift2):non(-rand_shift2)]
            sino_shifts[:,x,:] = projection2D_shift
    return [sino_shifts,shifts]


def pve(data, pve_strength):
    from scipy.ndimage import gaussian_filter
    data_pve = data.copy()
    if (pve_strength > 0):
        pass
    else:
        raise ("Smoothing kernel must be positive")
    if (data.ndim == 2):
        (anglesDim, DetectorsDimH) = np.shape(data)
        for x in range(anglesDim):
            data_pve[x,:] = gaussian_filter(data_pve[x,:], pve_strength)
    else:
        (DetectorsDimV, anglesDim, DetectorsDimH) = np.shape(data)
        for x in range(anglesDim):
            data_pve[:,x,:] = gaussian_filter(data_pve[:,x,:], pve_strength)
    return data_pve


def fresnel_propagator(data, dist_observation, scale_factor, wavelenght):
    # adapted from the script by Adrián Carbajal-Domínguez, adrian.carbajal@ujat.mx
    data_fresnel = data.copy()
    if (data.ndim == 2):
        (anglesDim, DetectorsDimH) = np.shape(data)
        n1 = DetectorsDimH*0.5
        #Define the angular spectrum coordinates
        u = np.arange(-n1, n1, 1)
        #Define the propagation matrix
        propagator=np.exp(2*np.pi*1j*(dist_observation/scale_factor)*np.sqrt((1/wavelenght)**2-(u/10)**2))
        #### Compute the Fast Fourier Transform of each 1D projection
        for x in range(anglesDim):
            f=np.fft.fft(data_fresnel[x,:])
            #Correct the low and high frequencies
            fshift=np.fft.fftshift(f)
            #multiply both matrices: Fourier transform and the propagator matrices.
            field=fshift*propagator
            #Calculate the inverse Fourier transform
            field2=np.fft.ifft(field)
            data_fresnel[x,:] = np.abs(field2)
    else:
        (DetectorsDimV, anglesDim, DetectorsDimH) = np.shape(data)
        ####Define the size of the propagation function p(u,v). It has to be of the same size of the image.
        n1 = DetectorsDimV*0.5
        n2 = DetectorsDimH*0.5
        #Define the angular spectrum coordinates
        u = np.arange(-n1, n1, 1)
        v = np.arange(-n2, n2, 1)
        U,V = np.meshgrid(u,v)
        #Define the propagation matrix
        propagator=np.exp(2*np.pi*1j*(dist_observation/scale_factor)*np.sqrt((1/wavelenght)**2-(U/scale_factor)**2-(V/scale_factor)**2))
        #### Compute the Fast Fourier Transform of each 1D projection
        for x in range(anglesDim):
            f=np.fft.fft2(data_fresnel[:,x,:])
            #Correct the low and high frequencies
            fshift=np.fft.fftshift(f)
            #multiply both matrices: Fourier transform and the propagator matrices.
            field=fshift*propagator
            #Calculate the inverse Fourier transform
            field2=np.fft.ifft2(field)
            data_fresnel[:,x,:] = np.abs(field2)
    return data_fresnel
