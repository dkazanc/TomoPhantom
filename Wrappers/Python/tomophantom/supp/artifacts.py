#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Various methods to simulate noise and artifacts for sinograms (2D) and for 3D-projection data
 
What can be simulated: 
-- noise (Poisson or Gaussian)
-- zingers (result in streaks in the reconstruction)
-- stripes (result in rings in the reconstruction)
---- variable intensity
---- partial
-- shifts - data misalignment (result in blur in the reconstruction)

Note that the TomoPhantom package is released under Apache License, Version 2.0
@author: Daniil Kazantsev
"""

def _Artifacts_(sinogram, 
                _noise_, # noise related dictionary
                _zingers_, # zingers dictionary
                _stripes_, # stripes dictionary
                _sinoshifts_ #sinoshifts
                ):
    ####### NOISE ########
    if ('type' not in _noise_):
        _noise_['type'] = 'Poisson'
    if ((_noise_['type'] != 'Gaussian') and (_noise_['type'] != 'Poisson')):
        print("Unknown noise type, select Gaussian or Poisson, now set to Poisson")
        _noise_['type'] = 'Poisson'
    if ('sigma' not in _noise_):
        # photon intensity (Poisson) or variance for Gaussian 
        _noise_['sigma'] = 10000
    if ('seed' not in _noise_):
        # seeds for noise (integers), None means random generation
        _noise_['seed'] = None
    if ('prelog' not in _noise_):
        # set to True in order to get the prelog (raw) data 
        _noise_['prelog'] = None
    ####### ZINGERS ########
    if ('percentage' not in _zingers_):
        # # percentage the amount of zingers to be added to the data
        _zingers_['percentage'] = None
    if ('modulus' not in _zingers_):
        # modulus to control the amount of 4/6 pixel clusters to be added
        _zingers_['modulus'] = 10
    ####### STRIPES ########
    if ('percentage' not in _stripes_):
        # percentage defines the amount of stripes in the data
        _stripes_['percentage'] = None
    if ('maxthickness' not in _stripes_):
        # maxthickness defines the maximal thickness of a stripe
        _stripes_['maxthickness'] = 1.0
    if ('intensity' not in _stripes_):
        # controls the intensity levels of stripes
        _stripes_['intensity'] = 0.1
    if ('type' not in _stripes_):
        # stripe types can be 'partial' or 'full'
        _stripes_['type'] = 'full'
    # variability multiplier to incorporate change of intensity in the stripe
    if ('variability' not in _stripes_):
        _stripes_['variability'] = 0.0
    ####### SINOSHIFTS ########
    if ('maxamplitude' not in _sinoshifts_):
        # maxamplitude (in int pixels) defines the maximal amplitude of each angular deviation
        _sinoshifts_['maxamplitude'] = None
    'Applying artifacts and noise to the data'
    # ZINGERS
    if _zingers_['percentage'] is not None:
        sino_artifacts = zingers(sinogram=sinogram, percentage=_zingers_['percentage'], modulus = _zingers_['modulus'])
        print("Zingers have been added to the data.")
    else:
        sino_artifacts = sinogram
    # STRIPES
    if _stripes_['percentage'] is not None:
        sino_artifacts = stripes(sinogram=sino_artifacts,\
                                 percentage=_stripes_['percentage'],\
                                 maxthickness = _stripes_['maxthickness'],\
                                 intensity_thresh = _stripes_['intensity'],\
                                 stripe_type = _stripes_['type'],\
                                 variability = _stripes_['variability'])
        print("Stripes have been added to the data.")
    # SINOSHIFTS
    if _sinoshifts_['maxamplitude'] is not None:
        [sino_artifacts, shifts] = sinoshifts(sinogram=sino_artifacts, maxamplitude = _sinoshifts_['maxamplitude'])
        print("Sinogram shifts have been added to the data.")
    # NOISE
    sino_artifacts = noise(sinogram=sino_artifacts, sigma=_noise_['sigma'], noisetype = _noise_['type'], seed = _noise_['seed'], prelog=_noise_['prelog'])
    print("{} noise have been added to the data.".format(_noise_['type']))
    
    if _sinoshifts_['maxamplitude'] is not None:
        return [sino_artifacts,shifts]
    else:
        return sino_artifacts

def stripes(sinogram, percentage, maxthickness, intensity_thresh, stripe_type, variability):
    # Function to add stripes (constant offsets) to sinogram which results in rings in the 
    # reconstructed image
    # - percentage defines the amount of stripes in the data
    # - maxthickness defines the maximal thickness of a stripe
    # - stripe_type can be 'partial' or 'full'
    # - variability multiplier to incorporate change of intensity in the stripe
    import numpy as np
    import random
    if (sinogram.ndim == 2):
        (anglesDim, DetectorsDimH) = np.shape(sinogram)
    else:
        (DetectorsDimV, anglesDim, DetectorsDimH) = np.shape(sinogram)
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
    sino_stripes = sinogram.copy()
    max_intensity = np.max(sino_stripes)
    range_detect = int((np.float32(DetectorsDimH))*(np.float32(percentage)/100.0))
    if (sinogram.ndim == 2):
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

def zingers(sinogram, percentage, modulus):
    # adding zingers (zero single pixels or small 4 pixels clusters) to data or 6 voxels to 3D projection data
    # percentage = 0.5, #- percentage - the amount of zingers to be added to the data
    # modulus = 10 # modulus to control the amount of 4/6 pixel clusters to be added
    import numpy as np
    import random
    if (sinogram.ndim == 2):
        (anglesDim, DetectorsDimH) = np.shape(sinogram)
    else:
        (DetectorsDimV, anglesDim, DetectorsDimH) = np.shape(sinogram)
    if 0.0 < percentage <= 100.0:
        pass
    else:
        raise ("percentage must be larger than zero but smaller than 100")
    if (modulus > 0):
        pass
    else:
        raise ("Modulus integer must be positive")
    sino_zingers = sinogram.copy()
    length_sino = np.size(sino_zingers)
    num_values = int((length_sino)*(np.float32(percentage)/100.0))
    sino_zingers_fl = sino_zingers.flatten()
    for x in range(num_values):
        randind = random.randint(0,length_sino) # generate random index 
        sino_zingers_fl[randind] = 0
        if ((x % int(modulus)) == 0):
            if (sinogram.ndim == 2):
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

def noise(sinogram, sigma, noisetype, seed, prelog):
    # Adding random noise to data
    # noisetype = None, # 'Gaussian', 'Poisson' or None
    # sigma = 10000, # photon intensity (Poisson) or variance for Gaussian 
    # seed = 0, # seeds for noise
    # prelog: None or True (get the raw pre-log data) 
    import numpy as np
    sino_noisy = sinogram.copy()
    if noisetype == 'Gaussian':
        # add normal Gaussian noise
        if seed is not None:
            np.random.seed(int(seed))
        sino_noisy += np.random.normal(loc = 0.0, scale = sigma, size = np.shape(sino_noisy))
        sino_noisy[sino_noisy<0] = 0
    elif noisetype == 'Poisson':
        # add Poisson noise
        maxSino = np.max(sinogram)
        if maxSino > 0:
            sino_noisy = sinogram/maxSino
            dataExp = sigma*np.exp(-sino_noisy)  # noiseless raw data
            sino_raw = np.random.poisson(dataExp) #adding Poisson noise
            div_res = np.float32(sino_raw)/np.max(sino_raw)
            sino_noisy = -np.log(div_res)*maxSino # log corrected data -> sinogram
            sino_noisy[sino_noisy<0] = 0
    else:
        print ("Select 'Gaussian' or 'Poisson' for noise type")
    if prelog is True:
        return [sino_noisy,sino_raw]
    else:
        return sino_noisy

def sinoshifts(sinogram, maxamplitude):
    # A function to add random shifts to sinogram rows (an offset for each angular position)
    # maxamplitude (integer pixels) defines the maximal amplitude of each angular deviation
    import numpy as np
    import random
    if (sinogram.ndim == 2):
        (anglesDim, DetectorsDimH) = np.shape(sinogram)
        shifts = np.zeros(anglesDim, dtype='int8') # the vector of shifts
    else:
        (DetectorsDimV, anglesDim, DetectorsDimH) = np.shape(sinogram)
        shifts = np.zeros([anglesDim,2], dtype='int8') # the 2D vector of shifts

    sino_shifts = np.zeros(np.shape(sinogram),dtype='float32')
    non = lambda s: s if s<0 else None
    mom = lambda s: max(0,s)
    for x in range(anglesDim):
        rand_shift = random.randint(-maxamplitude,maxamplitude)  #generate random shift
        if (sinogram.ndim == 2):
            shifts[x] = rand_shift
            projection = sinogram[x,:] # extract 1D projection
            projection_shift = np.zeros(np.shape(projection),dtype='float32')
            projection_shift[mom(rand_shift):non(rand_shift)] = projection[mom(-rand_shift):non(-rand_shift)]
            sino_shifts[x,:] = projection_shift
        else:
            rand_shift2 = random.randint(-maxamplitude,maxamplitude)  #generate random shift
            shifts[x,0] = rand_shift
            shifts[x,1] = rand_shift2
            projection2D = sinogram[:,x,:] # extract 2D projection
            projection2D_shift = np.zeros(np.shape(projection2D),dtype='float32')
            projection2D_shift[mom(rand_shift):non(rand_shift), mom(rand_shift2):non(rand_shift2)] = projection2D[mom(-rand_shift):non(-rand_shift),mom(-rand_shift2):non(-rand_shift2)]
            sino_shifts[:,x,:] = projection2D_shift
    return [sino_shifts,shifts]
