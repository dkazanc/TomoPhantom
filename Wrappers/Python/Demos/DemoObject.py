"""
GPLv3 license (ASTRA toolbox)
Note that the TomoPhantom package is released under Apache License, Version 2.0

Script to generate 2D/3D analytical objects and their sinograms
Recursively adding objects and sinos one can build a required model

>>>> Prerequisites: ASTRA toolbox, if one needs to do reconstruction <<<<<

Run demo from the folder "Demos"

@author: Daniil Kazantsev
"""
import numpy as np
import matplotlib.pyplot as plt
from tomophantom import TomoP2D 
from tomophantom.TomoP2D import Objects2D

# create a 2D object explicitly without using parameters file
N_size = 512 # define the grid

# define all objects bellow:
pp = {'Obj': Objects2D.GAUSSIAN,
      'C0' : 1.00,
      'x0' : 0.25,
      'y0' : -0.3,
      'a'  : 0.15,
      'b'  :  0.3,
      'phi': -30.0}

pp1 = {'Obj': Objects2D.RECTANGLE,
      'C0' : 1.00,
      'x0' : -0.2,
      'y0' : 0.2,
      'a'  : 0.25,
      'b'  :  0.4,
      'phi': 60.0}

myObjects = [pp, pp1] # dictionary of objects
Object1 = TomoP2D.Object(N_size, myObjects)

plt.close('all')
plt.figure(1)
plt.rcParams.update({'font.size': 21})
plt.imshow(Object1, vmin=0, vmax=1, cmap="BuPu")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('{}'.format('2D Object'))

# create sinogram analytically without using the parameters file
angles_num = int(0.5*np.pi*N_size); # angles number
angles = np.linspace(0,180,angles_num,dtype='float32')
angles_rad = angles*(np.pi/180)
P = int(np.sqrt(2)*N_size) #detectors

sino_an = TomoP2D.ObjectSino(N_size, P, angles, myObjects)

plt.figure(2)
plt.rcParams.update({'font.size': 21})
plt.imshow(sino_an, cmap="BuPu")
plt.colorbar(ticks=[0, 150, 250], orientation='vertical')
plt.title('{}'.format('Analytical sinogram of an object'))
#%%
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print ("Reconstructing analytical sinogram using FBP (ASTRA-TOOLBOX)...")
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
from tomophantom.supp.astraOP import AstraTools
Atools = AstraTools(P, angles_rad, N_size, 'cpu') # initiate a class object
FBPrec = Atools.fbp2D(sino_an)

plt.figure(3) 
plt.imshow(FBPrec, vmin=0, vmax=1, cmap="BuPu")
plt.title('FBP Reconstructed Model')
#%%
"""
# similarly one can create 3D objects
from tomophantom import TomoP3D
from tomophantom.TomoP3D import Objects3D
import matplotlib.pyplot as plt

N3D_size = 256

obj3D = {'Obj': Objects3D.GAUSSIAN, 
      'C0' : 1.00,
      'x0' :-0.25,
      'y0' : 0.1,
      'z0' : 0.0,
      'a'  : 0.2,
      'b'  :  0.35,
      'c'  :  0.7,
      'phi1'  :  30.0,
      'phi2'  :  60.0,
      'phi3': -25.0}

Object3D = TomoP3D.Object(N3D_size, [obj3D])

sliceSel = int(0.5*N3D_size)
#plt.gray()
plt.figure(7) 
plt.subplot(121)
plt.imshow(Object3D[sliceSel,:,:],vmin=0, vmax=1)
plt.title('3D Phantom, axial view')

plt.subplot(122)
plt.imshow(Object3D[:,sliceSel,:],vmin=0, vmax=1)
plt.title('3D Phantom, coronal view')
plt.show()
"""
#%%

