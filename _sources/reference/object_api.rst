.. _ref_object_api:

Object API
==========

In TomoPhantom, we have a dedicated functionality that allows the user to build bespoke phantoms by stacking multiple objects together.
One should look into the :func:`Object` function of :func:`tomophantom.TomoP2D` or :func:`tomophantom.TomoP3D` modules if one wishes to build a 2D or 3D phantom, respectively.
If one needs also projection data to be generated, then please see the :func:`ObjectSino` function in the same modules. 

One can notice that one of the input parameters for :func:`Object` and :func:`ObjectSino` functions is :data:`obj_params` (*a dictionary or a list of dictionaries with the parameters for object(s)*).
The :ref:`ref_object_2d` section gives the details about the meaning of parameters for 2D objects and :ref:`ref_object_api_python2d` explains how to
use it in Python. For 3D phantoms, see also :ref:`ref_object_3d` and :ref:`ref_object_api_python3d` sections, respectively.

.. _ref_object_2d:

Parameters of 2D objects explained
----------------------------------

There are 7 parameters in total to describe a 2D object. Please refer to :numref:`fig_objects2d_notations` when exploring parameters bellow.

.. dropdown:: 1. object's name, type: string

	For 2D phantoms the accepted objects names are :code:`gaussian`, :code:`parabola`, :code:`parabola1`, :code:`ellipse`, :code:`cone` and :code:`rectangle`. See more in the TomoPhantom paper [SX2018]_ for analytical formulae. Please note that when working with :func:`Object` functions, the object names need to be defined in capital letters, see :ref:`ref_object_api_python2d`.

.. dropdown:: 2. :math:`C_{0}`: intensity level, type: float

	Parameter :math:`C_{0}` defines the grayscale intensity level, given as a floating point number. :math:`C_{0}` can be either negative or positive. Objects in the model are concatenated by summation, so one can do a subtraction of objects by defining negative intensities.  

.. dropdown:: 3. :math:`x_{0}` : object's center along the :math:`X` axis, type: float

	Parameter :math:`x_{0}`  defines the location of the object's center along the :math:`X` axis. The range of the parameter is normally within [-0.5, 0.5]. 

.. dropdown:: 4. :math:`y_{0}` : object's center along the :math:`Y` axis, type: float

	Parameter :math:`y_{0}` defines the location of the object's center along the :math:`Y` axis. The range of the parameter is normally within [-0.5, 0.5]. 

.. dropdown:: 5. :math:`a` : object's half-width, type: float

	Parameter :math:`a`  defines the half-width of the object. The range of the parameter is normally within (0, 2). 

.. dropdown:: 6. :math:`b` : object's half-length, type: float

	Parameter :math:`b`  defines the half-length of the object. The range of the parameter is normally within (0, 2). 

.. dropdown:: 7. :math:`\phi` : object's rotation angle in degrees, type: float

	Parameter :math:`\phi` defines the object's rotation angle in degrees. Positive values lead to clockwise rotation. 

.. _fig_objects2d_notations:
.. figure::  ../_static/api_objects/objects2d_notations.png
    :scale: 50 %
    :alt: Notation for a 2D object

    The coordinate system in which the objects are defined and an ellipsoid object as an example. Note the position of :math:`X` and :math:`Y` axes and the ranges for the visible field of view. Parameters describing a 2D object are: object's name, :math:`C_{0}`, :math:`x_{0}`, :math:`y_{0}`, :math:`a`, :math:`b` and the rotation angle :math:`\phi`, see detailed description above.

.. _ref_object_api2d:

Parameters of 2D objects in library file
----------------------------------------

More general information about what the library file can contain is given in :ref:`ref_library_files_api`. Here we look at the part where the object is described: 

.. code-block:: text

    Object : ellipse 1.0 0.0 0.0 0.69 0.92 0.0;

TomoPhantom converts this string for a 2D model into the following behind the scenes:

.. code-block:: python

    objects_name="ellipse" # object's name
    C0=1.0 # object's intensity level
    x0=0.0 # object's center along the X-axis
    y0=0.0 # object's center along the Y-axis
    a=0.69 # object's half-width
    b=0.92 # object's half-length
    phi=0.0 # object's rotation angle


.. _ref_object_api_python2d:

Defining 2D objects in Python
-----------------------------

In order to create a 2D object in Python one should use :func:`Object` function of :func:`tomophantom.TomoP2D` to create a phantom-object and 
:func:`ObjectSino` function to generate sinogram data for that object. 

You will need to create a dictionary for one object or a list of dictionaries for multiple objects. For example, to create a phantom with one object one can do:

.. code-block:: python

    from tomophantom import TomoP2D
    from tomophantom.TomoP2D import Objects2D
    phantom_size = 256
    object1 = {
        "Obj": Objects2D.GAUSSIAN,
        "C0": 1.00,
        "x0": 0.25,
        "y0": -0.3,
        "a": 0.15,
        "b": 0.3,
        "phi": -30.0,
    }
    phantom = TomoP2D.Object(phantom_size, object1) # which will generate 256^2 phantom

.. note:: To define object's name in the enum class :code:`Objects2D` one need to use capital letters for objects: :code:`GAUSSIAN`, :code:`PARABOLA`, :code:`PARABOLA1`, :code:`ELLIPSE`, :code:`CONE`, :code:`RECTANGLE`.

And for multiple objects stacked together, one can do 

.. code-block:: python

    from tomophantom import TomoP2D
    from tomophantom.TomoP2D import Objects2D
    phantom_size = 256
    object1 = {
        "Obj": Objects2D.GAUSSIAN,
        "C0": 1.00,
        "x0": 0.25,
        "y0": -0.3,
        "a": 0.15,
        "b": 0.3,
        "phi": -30.0,
    }
    object2 = {
        "Obj": Objects2D.RECTANGLE,
        "C0": 1.00,
        "x0": -0.2,
        "y0": 0.2,
        "a": 0.25,
        "b": 0.4,
        "phi": 60.0,
    }
    myObjects = [object1, object2]  # dictionary of objects
    phantom = TomoP2D.Object(phantom_size, myObjects) # which will generate 256^2 phantom

.. _ref_object_3d:

Parameters of 3D objects explained
----------------------------------

3D object is an extension of a 2D object to a third dimension. 
Most of the parameters remain as in :ref:`ref_object_2d`, but we will list them bellow for clarity as we have now 11 parameters in total:

.. dropdown:: 1. object's name, type: string

	For 3D phantoms the accepted objects names are :code:`gaussian`, :code:`paraboloid`, :code:`ellipsoid`, :code:`cone`, :code:`cuboid` and :code:`elliptical_cylinder`. See more in the TomoPhantom paper [SX2018]_ for analytical formulae. Please note that when working with :func:`Object` functions, the object names need to be defined in capital letters, see :ref:`ref_object_api_python3d`.

.. dropdown:: 2. :math:`C_{0}`: intensity level, type: float

	Parameter :math:`C_{0}` defines the grayscale intensity level, given as a floating point number. :math:`C_{0}` can be either negative or positive. Objects in the model are concatenated by summation, so one can do a subtraction of objects by defining negative intensities.  

.. dropdown:: 3. :math:`x_{0}` : object's center along the :math:`X` axis, type: float

	Parameter :math:`x_{0}`  defines the location of the object's center along the :math:`X` axis. The range of the parameter is normally within [-0.5, 0.5]. 

.. dropdown:: 4. :math:`y_{0}` : object's center along the :math:`Y` axis, type: float

	Parameter :math:`y_{0}` defines the location of the object's center along the :math:`Y` axis. The range of the parameter is normally within [-0.5, 0.5]. 

.. dropdown:: 5. :math:`z_{0}` : object's center along the :math:`Z` axis, type: float

	Parameter :math:`z_{0}` defines the location of the object's center along the :math:`Z` axis. The range of the parameter is normally within [-0.5, 0.5]. 

.. dropdown:: 6. :math:`a` : object's half-width, type: float

	Parameter :math:`a`  defines the half-width of the object. The range of the parameter is normally within (0, 2). 

.. dropdown:: 7. :math:`b` : object's half-length, type: float

	Parameter :math:`b`  defines the half-length of the object. The range of the parameter is normally within (0, 2).

.. dropdown:: 8. :math:`c` : object's half-depth, type: float

	Parameter :math:`c`  defines the half-depth of the object. The range of the parameter is normally within (0, 2).

.. dropdown:: 9. :math:`\phi_{1}` : object's rotation angle in degrees, type: float

	Parameter :math:`\phi_{1}` First Euler angle in degrees. Positive values lead to clockwise rotation. 

.. dropdown:: 10. :math:`\phi_{2}` : object's rotation angle in degrees, type: float

	Parameter :math:`\phi_{2}` Second Euler angle in degrees. For 3D objects one can use the Euler angles to describe object's orientation. Please note, that although you can use that angle to control the rotation of the object, for projection data it will be ignored.

.. dropdown:: 11. :math:`\phi_{3}` : object's rotation angle in degrees, type: float

	Parameter :math:`\phi_{3}` Third Euler angle in degrees. For 3D objects one can use the Euler angles to describe object's orientation. Please note, that although you can use that angle to control the rotation of the object, for projection data it will be ignored.


.. warning:: The phantom objects can be controlled using all three Euler angles :math:`\phi_{1-3}`. However, if one also would like to generate projection data analytically, please use :math:`\phi_{2-3} = 0`. The full support of all three angles is not working currently, you still can use :math:`\phi_{1}` for both phantom and data generation.


.. _ref_object_api3d:

Parameters of 3D objects in library file
----------------------------------------

More general information about what the library file can contain is given in :ref:`ref_library_files_api`. Here we look at the part where the object is described: 

.. code-block:: text

    Object : ellipsoid 1.0 0.0 0.0 0.0 0.69 0.92 0.81 0.0 0.0 0.0

TomoPhantom converts this string for a 3D model into the following behind the scenes:

.. code-block:: python

    objects_name="ellipsoid" # object's name
    C0=1.0 # object's intensity level
    x0=0.0 # object's center along the X-axis
    y0=0.0 # object's center along the Y-axis
    z0=0.0 # object's center along the Z-axis
    a=0.69 # object's half-width
    b=0.92 # object's half-length
    c=0.81 # object's half-depth
    phi1=0.0 # object's rotation angle
    phi2=0.0 # object's rotation angle
    phi3=0.0 # object's rotation angle

.. _ref_object_api_python3d:

Defining 3D objects in Python
-----------------------------

In order to create a 3D object in Python one should use :func:`Object` function of :func:`tomophantom.TomoP3D` to create a phantom-object and 
:func:`ObjectSino` function to generate projection data for that object. 

You will need to create a dictionary for one object or a list of dictionaries for multiple objects. For example, to create a phantom with one object one can do:

.. code-block:: python

    from tomophantom import TomoP3D
    from tomophantom.TomoP3D import Objects3D
    phantom_size = 256
    object1 = {
        "Obj": Objects3D.GAUSSIAN,
        "C0": 1.0,
        "x0": -0.25,
        "y0": -0.15,
        "z0": 0.0,
        "a": 0.3,
        "b": 0.2,
        "c": 0.3,
        "phi1": 35.0, # note that phi2,3 are not set here. They are equal to zero by default.
    }
    phantom = TomoP3D.Object(phantom_size, object1) # which will generate 256^3 phantom

.. note:: To define object's name in the enum class :code:`Objects3D` one need to use capital letters for objects: :code:`GAUSSIAN`, :code:`PARABOLOID`, :code:`ELLIPSOID`, :code:`CONE`, :code:`CUBOID`, :code:`ELLIPCYLINDER`.

And for multiple 3D objects stacked together, one can do 

.. code-block:: python

    from tomophantom import TomoP3D
    from tomophantom.TomoP3D import Objects3D
    phantom_size = 256
    object1 = {
        "Obj": Objects3D.GAUSSIAN,
        "C0": 1.0,
        "x0": -0.25,
        "y0": -0.15,
        "z0": 0.0,
        "a": 0.3,
        "b": 0.2,
        "c": 0.3,
        "phi1": 35.0,
    }

    object2 = {
        "Obj": Objects3D.CUBOID,
        "C0": 1.00,
        "x0": 0.1,
        "y0": 0.2,
        "z0": 0.0,
        "a": 0.15,
        "b": 0.35,
        "c": 0.6,
        "phi1": -60.0,
    }
    myObjects = [object1, object2]  # dictionary of objects
    phantom = TomoP3D.Object(phantom_size, myObjects) # which will generate 256^3 phantom