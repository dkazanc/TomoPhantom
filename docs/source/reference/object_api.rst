.. _ref_object_api:

Object API
==========

In TomoPhantom, we have a dedicated functionality that allows the user to build bespoke phantoms by stacking multiple objects together.
One should look into the :func:`Object` function of :func:`tomophantom.TomoP2D` or :func:`tomophantom.TomoP3D` modules if one wishes to build a phantom (2D or 3D, respectively).
If one needs also projection data to be generated, then please see the :func:`ObjectSino` function in the same modules. 

One can notice that one of the input parameters for :func:`Object` and :func:`ObjectSino` functions is :data:`obj_params` (*a dictionary or a list of dictionaries with the parameters for object(s)*).
This is where one needs the understanding of what those parameters represent. Bellow, we provide the details about these parameters and how these parameters 
should be passed to the chosen function, see :ref:`ref_object_api_python2d`.

.. _ref_object_2d:

Parameters of 2D objects explained
----------------------------------

There are 7 parameters in total to describe a 2D object. Please refer to :numref:`fig_objects2d_notations` when exploring parameters bellow.

.. dropdown:: 1. object's name, type: string

	For 2D phantoms the accepted objects names are :code:`gaussian`, :code:`parabola`, :code:`parabola1`, :code:`ellipse`, :code:`cone` and :code:`rectangle`. See more in the TomoPhantom paper [SX2018]_ for analytical formulae.

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

.. _fig_objects2d_notations:
.. figure::  ../_static/api_objects/objects2d_notations.png
    :scale: 50 %
    :alt: Notation for a 2D object

    The coordinate system in which the objects are defined and an ellipsoid as an example. Note the location of :math:`X` and :math:`Y` axes and the ranges for the visible field of view. Parameters describing a 2D object are: object's name, :math:`C_{0}`, :math:`x_{0}`, :math:`y_{0}`, :math:`a`, :math:`b` and the rotation angle :math:`\phi`.

.. _ref_object_api2d:

Parameters of 2D object in library file
---------------------------------------

Let us look at this line in detail: 

.. code-block:: text

    Object : ellipse 1.0 0.0 0.0 0.69 0.92 0.0;

TomoPhantom converts this string of 7 parameters for 2D model into the following behind the scenes:

.. code-block:: text

    objects_name="ellipse"
    C0=1.0
    x0=0.0
    y0=0.0
    a=0.69
    b=0.92
    angle=0.0


.. _ref_object_api_python2d:

Defining 2D object in Python
----------------------------



.. _ref_object_3d:

Parameters of 3D objects explained
----------------------------------


.. _ref_object_api3d:

Parameters of 3D object in library file
---------------------------------------

