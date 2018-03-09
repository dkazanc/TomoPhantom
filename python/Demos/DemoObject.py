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
from tomophantom import TomoP3D
from astraOP import AstraTools
#%%
# create an object explicitly without using parameters file 
N_size = 512
params = np.array([('gaussian', 1.00, -0.25, 0.3, 0.15, 0.3, 30.0),], 
                  dtype=[('Obj',  '|S16'), ('C0', np.float32), ('x0', np.float32), ('y0',np.float32),('a',np.float32), ('b', np.float32),  ('phi', np.float32)])
pp = {'Obj': 'gaussian', 
      'C0' : 1.00, 
      'x0' : -0.25,
      'y0' : 0.3,
      'a'  : 0.15,
      'b'  :  0.3,
      'phi': 30.0}
pp1 = {'Obj': 'circle', 
      'C0' : 1.00, 
      'x0' : -0.25,
      'y0' : 0.3,
      'a'  : 0.15,
      'b'  :  0.3,
      'phi': 30.0}
myObjects = [pp, pp1]

Object1 = TomoP2D.Object(N_size, params)
#Object2 = TomoP2D.Object2(N_size, [pp])


plt.figure(1)
plt.rcParams.update({'font.size': 21})
plt.imshow(Object1, vmin=0, vmax=1)
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('{}'.format('2D Object.'))
#%%
# create sinogram analytically without using the parameters file
angles_num = int(0.5*np.pi*N_size); # angles number
angles = np.linspace(0,180,angles_num,dtype='float32')
angles_rad = angles*(np.pi/180)
P = int(np.sqrt(2)*N_size) #detectors

sino_an1 = TomoP2D.ObjectSino(N_size, P, angles, params)

plt.figure(2)
plt.rcParams.update({'font.size': 21})
plt.imshow(sino_an1)
plt.colorbar(ticks=[0, 150, 250], orientation='vertical')
plt.title('{}'.format('Analytical sinogram of an object'))
#%%
# Create another object
params = np.array([('rectangle', 1.00, 0.2, -0.2, 0.25, 0.4, 60.0),], dtype=[('Obj',  '|S16'), ('C0', np.float32), ('x0', np.float32), ('y0',np.float32),('a',np.float32), ('b', np.float32),  ('phi', np.float32)])
Object2 = TomoP2D.Object(N_size, params)


plt.figure(3)
plt.rcParams.update({'font.size': 21})
plt.imshow(Object2, vmin=0, vmax=1)
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('{}'.format('2D Object.'))
#%%
# create sinogram analytically without using the parameters file
angles_num = int(0.5*np.pi*N_size); # angles number
angles = np.linspace(0,180,angles_num,dtype='float32')
angles_rad = angles*(np.pi/180)
P = int(np.sqrt(2)*N_size) #detectors

sino_an2 = TomoP2D.ObjectSino(N_size, P, angles, params)

plt.figure(4)
plt.rcParams.update({'font.size': 21})
plt.imshow(sino_an2)
plt.colorbar(ticks=[0, 150, 250], orientation='vertical')
plt.title('{}'.format('Analytical sinogram of an object'))
#%%
# Now add objects into a model and reconstruct added sinos
Model = Object1 + Object2 # composite model
SinoModel = sino_an1 + sino_an2 # composite sino

plt.figure(5) 
plt.subplot(121)
plt.imshow(Model)
plt.title('New model')
plt.subplot(122)
plt.imshow(SinoModel)
plt.title('Sinogram of the model')
plt.show()
#%%
# lets reconstruct using ASTRA
Atools = AstraTools(P, angles_rad + 0.5*np.pi, N_size, 'cpu') # initiate a class object
FBPrec = Atools.fbp2D(SinoModel)

plt.figure(6) 
plt.imshow(FBPrec, vmin=0, vmax=1)
plt.title('FBP Reconstructed Model')

#%%
# similarly one can create 3D objects explicitly calling to object functions
N3D = 256
params = np.array([('gaussian', 1.00, -0.25, 0.1, 0.0, 0.2, 0.35, 0.7, 30.0, 60.0, -25.0),], dtype=[('Obj',  '|S22'), ('C0', np.float32), ('x0', np.float32), ('y0',np.float32), ('z0',np.float32), ('a',np.float32), ('b', np.float32), ('c', np.float32), ('psi1', np.float32),('psi2', np.float32),('psi3', np.float32)])
Object3D = TomoP3D.Object(N3D, params)

sliceSel = int(0.5*N3D)
#plt.gray()
plt.figure(7) 
plt.subplot(121)
plt.imshow(Object3D[sliceSel,:,:],vmin=0, vmax=1)
plt.title('3D Phantom, axial view')

plt.subplot(122)
plt.imshow(Object3D[:,sliceSel,:],vmin=0, vmax=1)
plt.title('3D Phantom, coronal view')
plt.show()
#%%

