#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The set of functions to generate random 2D/3D phantoms

The TomoPhantom package is released under Apache License, Version 2.0
@author: Daniil Kazantsev
"""


def rand_init2D(x0min, x0max, y0min, y0max, c0min, c0max, ab_min, ab_max):
    import numpy as np

    x0 = np.random.uniform(low=x0min, high=x0max)
    y0 = np.random.uniform(low=y0min, high=y0max)
    c0 = np.random.uniform(low=c0min, high=c0max)
    ab = np.random.uniform(low=ab_min, high=ab_max)
    return (x0, y0, c0, ab)


def rand_init3D(x0min, x0max, y0min, y0max, z0min, z0max, c0min, c0max, ab_min, ab_max):
    import numpy as np

    x0 = np.random.uniform(low=x0min, high=x0max)
    y0 = np.random.uniform(low=y0min, high=y0max)
    z0 = np.random.uniform(low=z0min, high=z0max)
    c0 = np.random.uniform(low=c0min, high=c0max)
    ab = np.random.uniform(low=ab_min, high=ab_max)
    return (x0, y0, z0, c0, ab)


# Function to generate 2D foam-like structures using randomly located circles
def foam2D(
    x0min,
    x0max,
    y0min,
    y0max,
    c0min,
    c0max,
    ab_min,
    ab_max,
    N_size,
    tot_objects,
    object_type,
):
    import numpy as np
    import math
    import random

    # 2D functions
    from tomophantom import TomoP2D
    from tomophantom.TomoP2D import Objects2D

    attemptsNo = 2000  # the number of attempts to fit the object
    # objects accepted: 'ellipse', 'parabola', 'gaussian', 'mix'
    mix_objects = False
    if object_type == "ellipse":
        object_type = Objects2D.ELLIPSE
    elif object_type == "parabola":
        object_type = Objects2D.PARABOLA
    elif object_type == "gaussian":
        object_type = Objects2D.GAUSSIAN
    elif object_type == "mix":
        mix_objects = True
    else:
        raise TypeError("object_type can be only ellipse, parabola, gaussian or mix")
    X0 = np.float32(np.zeros(tot_objects))
    Y0 = np.float32(np.zeros(tot_objects))
    AB = np.float32(np.zeros(tot_objects))
    C0_var = np.float32(np.zeros(tot_objects))
    for i in range(0, tot_objects):
        (x0, y0, c0, ab) = rand_init2D(
            x0min, x0max, y0min, y0max, c0min, c0max, ab_min, ab_max
        )
        if i > 0:
            breakj = False
            for j in range(0, attemptsNo):
                if breakj == True:
                    (x0, y0, c0, ab) = rand_init2D(
                        x0min, x0max, y0min, y0max, c0min, c0max, ab_min, ab_max
                    )
                    breakj = False
                else:
                    for l in range(
                        0, i
                    ):  # checks consistency with previously created objects
                        dist = math.sqrt((X0[l] - x0) ** 2 + (Y0[l] - y0) ** 2)
                        if (dist < (ab + AB[l])) or (
                            (abs(x0) + ab) ** 2 + (abs(y0) + ab) ** 2 > 1.0
                        ):
                            breakj = True
                            break
                    if breakj == False:  # re-initialise if doesn't fit the criteria
                        X0[i] = x0
                        Y0[i] = y0
                        AB[i] = ab
                        C0_var[i] = c0
                        break
        if AB[i] == 0.0:
            X0[i] = x0
            Y0[i] = y0
            AB[i] = 0.0001
            C0_var[i] = c0

    myObjects = []  # dictionary of objects
    for obj in range(0, len(X0)):
        if mix_objects == True:
            rand_obj = random.randint(0, 2)
            if rand_obj == 0:
                object_type = Objects2D.ELLIPSE
            if rand_obj == 1:
                object_type = Objects2D.PARABOLA
            if rand_obj == 2:
                object_type = Objects2D.GAUSSIAN
        curr_obj = {
            "Obj": object_type,
            "C0": C0_var[obj],
            "x0": X0[obj],
            "y0": Y0[obj],
            "a": AB[obj],
            "b": AB[obj],
            "phi": 0.0,
        }
        myObjects.append(curr_obj)
    Object = TomoP2D.Object(N_size, myObjects)
    return (Object, myObjects)


# Function to generate 3D foam-like structures using randomly located spheres
def foam3D(
    x0min,
    x0max,
    y0min,
    y0max,
    z0min,
    z0max,
    c0min,
    c0max,
    ab_min,
    ab_max,
    N_size,
    tot_objects,
    object_type,
):
    import numpy as np
    import math
    import random

    # 3D functions
    from tomophantom import TomoP3D
    from tomophantom.TomoP3D import Objects3D

    attemptsNo = 2000
    # objects accepted: 'ellipsoid', 'paraboloid', 'gaussian', 'mix'
    mix_objects = False
    if object_type == "ellipsoid":
        object_type = Objects3D.ELLIPSOID
    elif object_type == "paraboloid":
        object_type = Objects3D.PARABOLOID
    elif object_type == "gaussian":
        object_type = Objects3D.GAUSSIAN
    elif object_type == "mix":
        mix_objects = True
    else:
        raise TypeError("object_type can be only ellipse, parabola, gaussian or mix")
    X0 = np.float32(np.zeros(tot_objects))
    Y0 = np.float32(np.zeros(tot_objects))
    Z0 = np.float32(np.zeros(tot_objects))
    AB = np.float32(np.zeros(tot_objects))
    C0_var = np.float32(np.zeros(tot_objects))
    for i in range(0, tot_objects):
        (x0, y0, z0, c0, ab) = rand_init3D(
            x0min, x0max, y0min, y0max, z0min, z0max, c0min, c0max, ab_min, ab_max
        )
        if i > 0:
            breakj = False
            for j in range(0, attemptsNo):
                if breakj:
                    (x0, y0, z0, c0, ab) = rand_init3D(
                        x0min,
                        x0max,
                        y0min,
                        y0max,
                        z0min,
                        z0max,
                        c0min,
                        c0max,
                        ab_min,
                        ab_max,
                    )
                    breakj = False
                else:
                    for l in range(
                        0, i
                    ):  # checks consistency with previously created objects
                        dist = math.sqrt(
                            (X0[l] - x0) ** 2 + (Y0[l] - y0) ** 2 + (Z0[l] - z0) ** 2
                        )
                        if (dist < (ab + AB[l])) or (
                            (abs(x0) + ab) ** 2
                            + (abs(y0) + ab) ** 2
                            + (abs(z0) + ab) ** 2
                            > 1.0
                        ):
                            breakj = True
                            break
                    if breakj == False:  # re-initialise if doesn't fit the criteria
                        X0[i] = x0
                        Y0[i] = y0
                        Z0[i] = z0
                        AB[i] = ab
                        C0_var[i] = c0
                        break
        if AB[i] == 0.0:
            X0[i] = x0
            Y0[i] = y0
            AB[i] = 0.0001
            C0_var[i] = c0

    myObjects = []  # dictionary of objects
    for obj in range(0, len(X0)):
        if mix_objects == True:
            rand_obj = random.randint(0, 2)
            if rand_obj == 0:
                object_type = Objects3D.ELLIPSOID
            if rand_obj == 1:
                object_type = Objects3D.PARABOLOID
            if rand_obj == 2:
                object_type = Objects3D.GAUSSIAN
        curr_obj = {
            "Obj": object_type,
            "C0": C0_var[obj],
            "x0": X0[obj],
            "y0": Y0[obj],
            "z0": Z0[obj],
            "a": AB[obj],
            "b": AB[obj],
            "c": AB[obj],
            "phi1": 0.0,
        }
        myObjects.append(curr_obj)
    Object3D = TomoP3D.Object(N_size, myObjects)
    return (Object3D, myObjects)
