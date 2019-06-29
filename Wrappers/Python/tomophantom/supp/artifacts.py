#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Various methods to simulate noise and artifacts for sinograms (2D) and for 3D-projection data
 
What can be simulated: 
-- noise (Poisson or Gaussian)
-- zingers (result in streaks in the reconstruction)
-- stripes (result in rings in the reconstruction)
-- shifts - data misalignment (result in blur in the reconstruction)

Note that the TomoPhantom package is released under Apache License, Version 2.0
@author: Daniil Kazantsev
"""

def _Artifacts_(sinogram, 
                noise_type = None, # 'Gaussian', 'Poisson' or None
                noise_sigma = 10000, # photon intensity (Poisson) or variance for Gaussian 
                noise_seed = None, # seeds for noise (integers), None means random generation
                zingers_percentage = None, # percentage the amount of zingers to be added to the data
                zingers_modulus = 10, # modulus to control the amount of 4/6 pixel clusters to be added
                stripes_percentage = None, # percentage defines the amount of stripes in the data
                stripes_maxthickness = 1.0, # maxthickness defines the maximal thickness of a stripe
                sinoshifts_maxamplitude = None # maxamplitude (in int pixels) defines the maximal amplitude of each angular deviation
                ):
    # ZINGERS
    if zingers_percentage is not None:
        sino_artifacts = zingers(sinogram=sinogram, percentage=zingers_percentage, modulus = zingers_modulus)
        print("Zingers have been added to the data.")
    else:
        sino_artifacts = sinogram
    # STRIPES
    if stripes_percentage is not None:
        sino_artifacts = stripes(sinogram=sino_artifacts, percentage=stripes_percentage, maxthickness = stripes_maxthickness)
        print("Stripes have been added to the data.")
    # SINOSHIFTS
    if sinoshifts_maxamplitude is not None:
        sino_artifacts = sinoshifts(sinogram=sino_artifacts, maxamplitude = sinoshifts_maxamplitude)
        print("Sinogram shifts have been added to the data.")
    # NOISE
    if noise_type is not None:
        sino_artifacts = noise(sinogram=sino_artifacts, sigma=noise_sigma, noisetype = noise_type, seed = noise_seed)
        print("{} noise have been added to the data.".format(noise_type))
    
    return sino_artifacts

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

def noise(sinogram, sigma, noisetype, seed):
    # Adding random noise to data
    # noisetype = None, # 'Gaussian', 'Poisson' or None
    # sigma = 10000, # photon intensity (Poisson) or variance for Gaussian 
    # seed = 0, # seeds for noise
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
            sino_noisy = np.random.poisson(dataExp) #adding Poisson noise
            div_res = np.float32(sino_noisy)/np.max(sino_noisy)
            sino_noisy = -np.log(div_res)*maxSino # log corrected data -> sinogram
            sino_noisy[sino_noisy<0] = 0
    else:
        print ("Select 'Gaussian' or 'Poisson' for noise type")
    return sino_noisy

def stripes(sinogram, percentage, maxthickness):
    # Function to add stripes (constant offsets) to sinogram which results in rings in the 
    # reconstructed image
    #- percentage defines the amount of stripes in the data
    #- maxthickness defines the maximal thickness of a stripe
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
    sino_stripes = sinogram.copy()
    max_intensity = np.max(sino_stripes)
    range_detect = int((np.float32(DetectorsDimH))*(np.float32(percentage)/100.0))
    if (sinogram.ndim == 2):
        for x in range(range_detect):
            randind = random.randint(0,DetectorsDimH) # generate random index
            randthickness = random.randint(0,maxthickness) #generate random thickness
            randintens = random.uniform(-1.0, 0.5) # generate random multiplier
            intensity = max_intensity*randintens
            if ((randind > 0+randthickness) & (randind < DetectorsDimH-randthickness)):
                for x1 in range(-randthickness,randthickness+1):
                    sino_stripes[:,randind+x1] += intensity
    else:
        for j in range(DetectorsDimV):
            for x in range(range_detect):
                randind = random.randint(0,DetectorsDimH) # generate random index
                randthickness = random.randint(0,maxthickness) #generate random thickness
                randintens = random.uniform(-1, 0.5) # generate random multiplier
                intensity = max_intensity*randintens
                if ((randind > 0+randthickness) & (randind < DetectorsDimH-randthickness)):
                    for x1 in range(-randthickness,randthickness+1):
                        sino_stripes[j,:,randind+x1] += intensity
    return sino_stripes

def sinoshifts(sinogram, maxamplitude):
    #A  function to add random shifts to sinogram rows (an offset for each angular position)
    # maxamplitude (in int pixels) defines the maximal amplitude of each angular deviation
    import numpy as np
    import random
    if (sinogram.ndim == 2):
        (anglesDim, DetectorsDimH) = np.shape(sinogram)
    else:
        (DetectorsDimV, anglesDim, DetectorsDimH) = np.shape(sinogram)

    sino_shifts = np.zeros(np.shape(sinogram),dtype='float32')
    non = lambda s: s if s<0 else None
    mom = lambda s: max(0,s)
    for x in range(anglesDim):
        rand_shift = random.randint(-maxamplitude,maxamplitude)  #generate random shift
        if (sinogram.ndim == 2):
            projection = sinogram[x,:] # extract 1D projection
            projection_shift = np.zeros(np.shape(projection),dtype='float32')
            projection_shift[mom(rand_shift):non(rand_shift)] = projection[mom(-rand_shift):non(-rand_shift)]
            sino_shifts[x,:] = projection_shift
        else:
            rand_shift2 = random.randint(-maxamplitude,maxamplitude)  #generate random shift
            projection2D = sinogram[:,x,:] # extract 2D projection
            projection2D_shift = np.zeros(np.shape(projection2D),dtype='float32')
            projection2D_shift[mom(rand_shift):non(rand_shift), mom(rand_shift2):non(rand_shift2)] = projection2D[mom(-rand_shift):non(-rand_shift),mom(-rand_shift2):non(-rand_shift2)]
            sino_shifts[:,x,:] = projection2D_shift
    return sino_shifts