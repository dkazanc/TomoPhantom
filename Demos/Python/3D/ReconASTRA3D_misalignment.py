import timeit
import os
import matplotlib.pyplot as plt
import numpy as np
import tomophantom
from tomophantom import TomoP3D


print ("Building 3D phantom using TomoPhantom software")
tic=timeit.default_timer()
model = 13 # select a model number from the library
N_size = 128 # Define phantom dimensions using a scalar value (cubic phantom)
path = os.path.dirname(tomophantom.__file__)
path_library3D = os.path.join(path, "Phantom3DLibrary.dat")
#This will generate a N_size x N_size x N_size phantom (3D)
phantom_tm = TomoP3D.Model(model, N_size, path_library3D)
toc=timeit.default_timer()
Run_time = toc - tic
print("Phantom has been built in {} seconds".format(Run_time))

sliceSel = int(0.5*N_size)
#plt.gray()
plt.figure() 
plt.subplot(131)
plt.imshow(phantom_tm[sliceSel,:,:],vmin=0, vmax=1)
plt.title('3D Phantom, axial view')

plt.subplot(132)
plt.imshow(phantom_tm[:,sliceSel,:],vmin=0, vmax=1)
plt.title('3D Phantom, coronal view')

plt.subplot(133)
plt.imshow(phantom_tm[:,:,sliceSel],vmin=0, vmax=1)
plt.title('3D Phantom, sagittal view')
plt.show()

# Projection geometry related parameters:
Horiz_det = int(2*N_size) # detector column count (horizontal)
Vert_det = N_size # detector row count (vertical) (no reason for it to be > N)
angles_num = int(0.5*np.pi*N_size); # angles number
angles = np.linspace(0.0,179.9,angles_num,dtype='float32') # in degrees
angles_rad = angles*(np.pi/180.0)


#%%
print ("Building 3D analytical projection data with TomoPhantom")
projData3D_analyt= TomoP3D.ModelSino(model, N_size, Horiz_det, Vert_det, angles, path_library3D)

intens_max = 70
sliceSel = int(0.5*N_size)
plt.figure() 
plt.subplot(131)
plt.imshow(projData3D_analyt[:,sliceSel,:],vmin=0, vmax=intens_max)
plt.title('2D Projection (analytical)')
plt.subplot(132)
plt.imshow(projData3D_analyt[sliceSel,:,:],vmin=0, vmax=intens_max)
plt.title('Sinogram view')
plt.subplot(133)
plt.imshow(projData3D_analyt[:,:,sliceSel],vmin=0, vmax=intens_max)
plt.title('Tangentogram view')
plt.show()
#%%
print ("Adding noise to projection data")
from tomophantom.supp.artifacts import _Artifacts_

# forming dictionaries with artifact types
_noise_ =  {'noise_type' : 'Poisson',
            'noise_amplitude' : 10000, # noise amplitude
            'noise_seed' : 0}

_datashifts_ = {'datashifts_maxamplitude_subpixel' : 5} # subpixel misalignment

[prj, shifts_exact] = _Artifacts_(projData3D_analyt, **_noise_, **_datashifts_)

intens_max = 70
sliceSel = int(0.5*N_size)
plt.figure() 
plt.subplot(131)
plt.imshow(prj[:,sliceSel,:],vmin=0, vmax=intens_max)
plt.title('2D Projection shifted (analytical)')
plt.subplot(132)
plt.imshow(prj[sliceSel,:,:],vmin=0, vmax=intens_max)
plt.title('Sinogram shifted view')
plt.subplot(133)
plt.imshow(prj[:,:,sliceSel],vmin=0, vmax=intens_max)
plt.title('Tangentogram shifted view')
plt.show()
#%%
# perform reconstruction
from tomobar.methodsDIR import RecToolsDIR
RectoolsDIR = RecToolsDIR(DetectorsDimH = Horiz_det, # Horizontal detector dimension
                    DetectorsDimV = Vert_det,        # Vertical detector dimension (3D case)
                    CenterRotOffset = -shifts_exact, # Center of Rotation scalar
                    AnglesVec = angles_rad,          # A vector of projection angles in radians
                    ObjSize = N_size,                # Reconstructed object dimensions (scalar)
                    device_projector='gpu')

FBPrec = RectoolsDIR.FBP(prj)

sliceSel = int(0.5*N_size)
#plt.gray()
plt.figure() 
plt.subplot(131)
plt.imshow(FBPrec[sliceSel,:,:])
plt.title('3D Recon, axial view')

plt.subplot(132)
plt.imshow(FBPrec[:,sliceSel,:])
plt.title('3D Recon, coronal view')

plt.subplot(133)
plt.imshow(FBPrec[:,:,sliceSel])
plt.title('3D Recon, sagittal view')
plt.show()
#%%