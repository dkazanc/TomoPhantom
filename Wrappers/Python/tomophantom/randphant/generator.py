#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This function generates random phantoms 

@author: Daniil Kazantsev
"""

import random
import numpy as np
import math
import matplotlib.pyplot as plt

from tomophantom import TomoP2D 
from tomophantom.TomoP2D import Objects2D

N_size = 512 # define the grid size
tot_objects = 100 # the total number of objects to generate
attemptsNo = 1000

x0min = -0.9
x0max = 0.9
y0min = -0.9
y0max = 0.9
c0min = 0.01
c0max = 1.0
ab_min = 0.01
ab_max = 0.25

X0 = np.float32(np.zeros(tot_objects))
Y0 = np.float32(np.zeros(tot_objects))
AB = np.float32(np.zeros(tot_objects))
C0_var = np.float32(np.zeros(tot_objects))

def rand_init2D(x0min, x0max, y0min, y0max, c0min, c0max, ab_min, ab_max):
    x0 = np.random.uniform(low = x0min, high=x0max)
    y0 = np.random.uniform(low = y0min, high=y0max)
    c0 = np.random.uniform(low = c0min, high=c0max)
    ab = np.random.uniform(low = ab_min, high=ab_max)
    return (x0,y0,c0,ab)

for i in range(0,tot_objects):
    (x0,y0,c0,ab) = rand_init2D(x0min, x0max, y0min, y0max, c0min, c0max, ab_min, ab_max)
    if (i > 0):
        breakj = False
        for j in range(0,attemptsNo):
            if breakj:
                (x0,y0,c0,ab) = rand_init2D(x0min, x0max, y0min, y0max, c0min, c0max, ab_min, ab_max)
                breakj = False
            else:
                for l in range(0,i): # checks consistency with previously created objects
                    dist = math.sqrt((X0[l]-x0)**2 + (Y0[l]-y0)**2)
                    if (dist < (ab + AB[l])) or ((abs(x0) + ab)**2 + (abs(y0) + ab)**2 > 1.0):
                        breakj = True
                        break
                if (breakj == False): # re-initialise if doesn't fit the criteria
                    X0[i] = x0
                    Y0[i] = y0
                    AB[i] = ab
                    C0_var[i] = c0
                    break
    else:
        X0[i] = x0
        Y0[i] = y0
        AB[i] = ab
        C0_var[i] = c0

myObjects = [] # dictionary of objects

for obj in range(0,len(X0)):
    #if ((abs(X0[obj]) + AB[obj])**2 + (abs(Y0[obj]) + AB[obj])**2 < 1.0):
    curr_obj = {'Obj': Objects2D.ELLIPSE,
                'C0' : C0_var[obj],
                'x0' : X0[obj],
                'y0' : Y0[obj],
                'a'  : AB[obj],
                'b'  : AB[obj],
                'phi': 0.0}
    myObjects.append(curr_obj)

Object = TomoP2D.Object(N_size, myObjects)


plt.figure()
plt.imshow(Object,  vmin=0, vmax=0.1, cmap="gray")
#%%
indeces_del = np.uint8(np.zeros(tot_objects))
for i in range(0,tot_objects):
    for j in range(0,tot_objects):
        if (i != j):
            # check if the distance between the current and previous object
            dist = math.hypot(X0[i] - X0[j], Y0[i]-Y0[j])
            if (dist < (AB[i] + AB[j])):
                # circles intersect -> remove i-th cirle
                indeces_del[i] = j
                break


#indeces_del2 = indeces_del[np.nonzero(indeces_del)]
indeces_del2 = indeces_del[indeces_del>0] 
X0 = np.delete(X0,indeces_del2)
Y0 = np.delete(Y0,indeces_del2)
AB = np.delete(AB,indeces_del2)
C0_var = np.delete(C0_var,indeces_del2)
#%%
