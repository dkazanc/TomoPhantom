# this is a collection of speckle generating routines taken from:
# https://github.com/fperakis/Speckles-simulations

import numpy as np
from scipy.ndimage import gaussian_filter

def sample2ddist(asicshape,kbar,dist):
    """
    samples photons from a 2D distribution
    """
    #Definitions
    nasic = asicshape[0]*asicshape[1]
    nphotons = np.round(kbar*nasic).astype(int) 
    photons = np.zeros(asicshape)   
    
    #sample photon positions from given distribution
    sampl = np.random.choice(nasic,nphotons,p=dist)

    #convert to 2d and accumulate photons
    for i in range (nphotons):
        ix = int(sampl[i]/asicshape[1])
        iy = int(sampl[i]%asicshape[1])
        photons[ix,iy] += 1 
    
    return photons

def model_speckles(asicshape,specklesize):
    """
    models a speckle pattern with a given speckle size (integer) based on random
    phasors (see JC Goodman - Speckle Phenomena in Optics - Appendix G)
    """
    #definitions speckles
    randphasors = np.zeros([asicshape[0],asicshape[1]],dtype = complex)
    nphasorsx = np.round(asicshape[0]/specklesize).astype(int)
    nphasorsy = np.round(asicshape[1]/specklesize).astype(int)

    #random phasors
    for i in range(nphasorsx):
        for j in range(nphasorsy):
            randphase = np.random.random_sample()*2.j*np.pi
            randphasors[i,j] = np.exp(randphase)

    #speckles
    specklefield = np.fft.fft2(randphasors)
    speckleint = np.absolute(specklefield)*np.absolute(specklefield)
    speckleint/=np.sum(speckleint) #normalise
    
    return speckleint   

def model_speckles_modes(asicshape,specklesize,modes):
    """
    models a speckle pattern with a given speckle size (integer) and
    contrast (1/modes) 
    """
   
    #definitions
    nx,ny = asicshape[0],asicshape[1]

    # add speckle patterns as many times as the number of modes
    dist = model_speckles([nx,ny],specklesize)
    for m in range(modes-1):
        dist += model_speckles([nx,ny],specklesize)
    dist/=np.sum(dist) # normalise

    return  dist

def simulate_speckles_with_shot_noise(asicshape,modes,specklesize,kbar):
    """
    returns a 2D speckle pattern with given speckle size, contrast
    (1/modes) and photon density kbar (photons/pixel). 
    """
    #definitions
    nx,ny = asicshape[0],asicshape[1]
    dist = model_speckles_modes(asicshape,specklesize,modes).flatten()
    # sample photons from the distribution
    sig =sample2ddist([nx,ny],kbar,dist)

    return sig

def simulate_shot_noise(asicshape,kbar):
    """
    returns a 2d shot noise pattern with given dimensions and
    photon density kbar (photons/pixel)
    """
    #definitions
    nx,ny = asicshape[0],asicshape[1]

    # flat 2d distribution
    dist = np.ones((nx,ny)).flatten()
    dist/=np.sum(dist)

    # sample photons from the distribution
    sig =sample2ddist([nx,ny],kbar,dist)

    return np.array(sig,dtype = float)

def simulate_charge_sharing(asicshape,modes,specklesize,kbar,sigma=1):
    """
    returns a 2D speckle pattern with given speckle size, contrast
    (1/modes) and photon density kbar (photons/pixel)
    including charge sharing simulated with a gaussian at each photon position with with sigma
    """
    #Definitions
    nasic = asicshape[0]*asicshape[1]
    nphotons = np.round(kbar*nasic).astype(int) 
    photons = np.zeros(asicshape)
    photons_blur = np.zeros(asicshape)
    dist = model_speckles_modes(asicshape,specklesize,modes).flatten()

    #sample photon positions from given distribution
    sampl = np.random.choice(nasic,nphotons,p=dist)

    #convert to 2d and accumulate photons
    for i in range (nphotons):
        ix = int(sampl[i]/asicshape[1])
        iy = int(sampl[i]%asicshape[1])
        photons_blur = np.zeros(asicshape)
        photons_blur[ix,iy] = 1
        photons_blur = gaussian_filter(photons_blur,sigma=sigma)
        photons += photons_blur

    return photons