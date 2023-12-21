import timeit
import os
import matplotlib.pyplot as plt
import numpy as np
import tomophantom
from tomophantom import TomoP3D


print("Building 3D phantom using TomoPhantom software")
tic = timeit.default_timer()
model = 13  # select a model number from the library
N_size = 128  # Define phantom dimensions using a scalar value (cubic phantom)
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
angles_num = int(0.5 * np.pi * N_size)  # angles number
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

_datashifts_ = {"datashifts_maxamplitude_subpixel": 10}  # subpixel misalignment

[prj, shifts_exact] = artefacts_mix(projData3D_analyt, **_noise_, **_datashifts_)

intens_max = 70
sliceSel = int(0.5 * N_size)
plt.figure()
plt.subplot(131)
plt.imshow(prj[:, sliceSel, :], vmin=0, vmax=intens_max)
plt.title("2D Projection shifted (analytical)")
plt.subplot(132)
plt.imshow(prj[sliceSel, :, :], vmin=0, vmax=intens_max)
plt.title("Sinogram shifted view")
plt.subplot(133)
plt.imshow(prj[:, :, sliceSel], vmin=0, vmax=intens_max)
plt.title("Tangentogram shifted view")
plt.show()
# %%
# scale the data first
scl = max(abs(prj.max()), abs(prj.min()))
prj /= scl
# %%
from tomobar.methodsDIR import RecToolsDIR

RectoolsDIR = RecToolsDIR(
    DetectorsDimH=Horiz_det,  # Horizontal detector dimension
    DetectorsDimV=Vert_det,  # Vertical detector dimension (3D case)
    CenterRotOffset=0.0,  # Center of Rotation scalar
    AnglesVec=angles_rad,  # A vector of projection angles in radians
    ObjSize=N_size,  # Reconstructed object dimensions (scalar)
    device_projector="gpu",
)

recon_nocorr = RectoolsDIR.FBP(prj)

sliceSel = int(0.5 * N_size)
plt.figure()
plt.subplot(131)
plt.imshow(recon_nocorr[sliceSel, :, :])
plt.title("3D Recon no correction, axial view")

plt.subplot(132)
plt.imshow(recon_nocorr[:, sliceSel, :])
plt.title("3D Recon no correction, coronal view")

plt.subplot(133)
plt.imshow(recon_nocorr[:, :, sliceSel])
plt.title("3D Recon no correction, sagittal view")
plt.show()

# %%
######################################################
############# Identify projection shifts #############
######################################################
# Based on alignment algorithm by D. Gursoy (see TomoPy package)
# use phase cross correlation to identify shifts
from skimage.registration import phase_cross_correlation

iter_number = 20  # the number of algorithm iterations
shifts_estimated = np.zeros([angles_rad.size, 2], dtype="float32")
recon = recon_nocorr.copy()
for iteration in range(0, iter_number):
    # re-project the reconstructed image
    reproj = RectoolsDIR.FORWPROJ(recon)

    for m in range(0, angles_rad.size):
        # Register current projection in sub-pixel precision
        shift, error, diffphase = phase_cross_correlation(
            prj[:, m, :], reproj[:, m, :], upsample_factor=20
        )
        shifts_estimated[m, 0] += shift[1]
        shifts_estimated[m, 1] += shift[0]

    # reconstruct with the estimated shifts
    RectoolsDIR = RecToolsDIR(
        DetectorsDimH=Horiz_det,  # Horizontal detector dimension
        DetectorsDimV=Vert_det,  # Vertical detector dimension (3D case)
        CenterRotOffset=shifts_estimated,  # Center of Rotation scalar + shifts passed
        AnglesVec=angles_rad,  # A vector of projection angles in radians
        ObjSize=N_size,  # Reconstructed object dimensions (scalar)
        device_projector="gpu",
    )

    recon = RectoolsDIR.FBP(prj)

sliceSel = int(0.5 * N_size)
plt.figure()
plt.subplot(131)
plt.imshow(recon[sliceSel, :, :])
plt.title("3D estimated Recon, axial view")

plt.subplot(132)
plt.imshow(recon[:, sliceSel, :])
plt.title("3D estimated Recon, coronal view")

plt.subplot(133)
plt.imshow(recon[:, :, sliceSel])
plt.title("3D estimated Recon, sagittal view")
plt.show()
# %%
# scale back the data
prj *= scl
# %%
# perform reconstruction with EXACT shifts
from tomobar.methodsDIR import RecToolsDIR

RectoolsDIR = RecToolsDIR(
    DetectorsDimH=Horiz_det,  # Horizontal detector dimension
    DetectorsDimV=Vert_det,  # Vertical detector dimension (3D case)
    CenterRotOffset=-shifts_exact,  # Center of Rotation scalar
    AnglesVec=angles_rad,  # A vector of projection angles in radians
    ObjSize=N_size,  # Reconstructed object dimensions (scalar)
    device_projector="gpu",
)

FBPrec_exact = RectoolsDIR.FBP(prj)

sliceSel = int(0.5 * N_size)
plt.figure()
plt.subplot(131)
plt.imshow(FBPrec_exact[sliceSel, :, :])
plt.title("3D exact Recon, axial view")

plt.subplot(132)
plt.imshow(FBPrec_exact[:, sliceSel, :])
plt.title("3D exact Recon, coronal view")

plt.subplot(133)
plt.imshow(FBPrec_exact[:, :, sliceSel])
plt.title("3D exact Recon, sagittal view")
plt.show()
# %%
