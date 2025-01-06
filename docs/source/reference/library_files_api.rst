.. _ref_library_files_api:

Library files API
=================

There is a certain API protocol how the user should define objects in the :ref:`ref_glossary_library`. If one is interested to learn about parameters of the objects, it is better to skip to :ref:`ref_object_api` directly.


.. _ref_library_files_api2d:

Library files for 2D models
---------------------------

Let us consider a 2D model no. 1 (a classical Shepp-Logan phantom) as an example from the :data:`Phantom2DLibrary.dat` file.

.. code-block:: text

    #----------------------------------------------------
    # Classical Shepp-Logan phantom (10 ellipses)
    Model : 01;
    Components : 10;
    TimeSteps : 1;
    Object : ellipse 1.0 0.0 0.0 0.69 0.92 0.0;
    Object : ellipse -0.8 0.0184 0.0 0.6624 0.874 0.0;
    Object : ellipse -0.2 0.0 0.22 0.11 0.31 18.0;
    Object : ellipse -0.2 0.0 -0.22 0.16 0.41 -18.0;
    Object : ellipse 0.1 -0.35 0.0 0.21 0.25 0.0;
    Object : ellipse 0.1 -0.1 0.0 0.046 0.046 0.0;
    Object : ellipse 0.1 0.1 0.0 0.046 0.046 0.0;
    Object : ellipse 0.1 0.605 -0.08 0.046 0.023 0.0;
    Object : ellipse 0.1 0.605 0.0 0.023 0.023 0.0;
    Object : ellipse 0.1 0.605 0.06 0.023 0.046 0.0;
    #----------------------------------------------------

.. note:: All the lines starting with the :code:`#` symbol will be parsed into comments. Also each line must end with the :code:`;` symbol.

* 1. :code:`Model : 01;` this is how the unique model number should be specified. Please note that the numbers [0-99] are reserved for static phantoms while 100 and larger for dynamic extensions.
* 2. :code:`Components : 10;` Here one should specify the number of components present in that model.
* 3. :code:`TimeSteps : 1;` TimeSteps are only relevant for :ref:`ref_glossary_dynamic_model` and for a stationary object this should be set to 1. 
* 4. :code:`Object : ellipse 1.0 0.0 0.0 0.69 0.92 0.0;` This is where the name of the elementary function and its parameters should be set. Here :code:`Object :` is a necessary syntax to initiate an object build.

Let us look at this line in detail: 

.. code-block:: text

    Object : ellipse 1.0 0.0 0.0 0.69 0.92 0.0;

TomoPhantom converts this string [#f2]_ into the following behind the scenes:

.. code-block:: python

    objects_name="ellipse" # object's name
    C0=1.0 # object's intensity level
    x0=0.0 # object's center along the X-axis
    y0=0.0 # object's center along the Y-axis
    a=0.69 # object's half-width
    b=0.92 # object's half-length
    phi=0.0 # object's rotation angle

As we move here into the :ref:`ref_glossary_object` API territory, it is best to re-direct the reader to the  :ref:`ref_object_2d` section. The meaning of each parameter for a 2D object is explained in detail there. 


.. _ref_library_files_api3d:

Library files for 3D models
---------------------------

3D models and the dynamic phantoms from the :data:`Phantom3DLibrary.dat` file are similar to 2D, but with an addition of few parameters. For instance, the 3D model 13 is given as: 

.. code-block:: text

    #----------------------------------------------------
    #  3D Shepp-Logan
    Model : 13;
    Components : 10;
    TimeSteps : 1;
    Object : ellipsoid 1.0 0.0 0.0 0.0 0.69 0.92 0.81 0.0 0.0 0.0
    Object : ellipsoid -0.8 0.0184 0.0 0.0 0.6624 0.874 0.78 0.0 0.0 0.0
    Object : ellipsoid -0.2 0.0 0.22 0.0 0.11 0.31 0.22 18.0 0.0 0.0
    Object : ellipsoid -0.2 0.0 -0.22 0.0 0.16 0.41 0.28 -18.0 0.0 0.0
    Object : ellipsoid 0.1 -0.35 0.0 -0.15 0.21 0.250 0.41 0.0 0.0 0.0
    Object : ellipsoid 0.1 -0.1 0.0 0.0 0.046 0.046 0.05 0.0 0.0 0.0
    Object : ellipsoid 0.1 0.1 0.0 0.0 0.046 0.046 0.05 0.0 0.0 0.0
    Object : ellipsoid 0.1 0.605 -0.08 0.0 0.046 0.023 0.05 0.0 0.0 0.0
    Object : ellipsoid 0.1 0.606 0.0 0.0 0.023 0.023 0.02 0.0 0.0 0.0
    Object : ellipsoid 0.1 0.605 0.06 0.0 0.023 0.046 0.02 0.0 0.0 0.0
    #----------------------------------------------------

and the objects definition is getting longer:

.. code-block:: text

    Object : ellipsoid 1.0 0.0 0.0 0.0 0.69 0.92 0.81 0.0 0.0 0.0

Which is converted into the following:

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

Again, for more in-depth read about parameters of 3D objects see :ref:`ref_object_api3d`.

.. _ref_library_dynamic_extensions:

Dynamic phantoms API
---------------------

.. rubric:: Footnotes

.. [#f2] Arguably, the choice of the text format to represent configuration library files might not be the best choice. Potentially, choosing the YAML language instead would make parameters for objects more informative and readable. However, the choice of text file libraries was historically inherited from the initial implementation of software in C language. We understand that it would be great to refactor it at some point.   





