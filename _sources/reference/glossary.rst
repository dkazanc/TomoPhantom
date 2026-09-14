.. _ref_glossary:

Glossary
========

Some terms used in TomoPhantom documentation can be ambiguous, here we provide explanations on that terminology. 

.. _ref_glossary_model:

Model
-----

When using  **model** term, the authors normally refer to a composite phantom from the pre-defined phantom model
that is stored in :ref:`ref_glossary_library`. The term can be applied both to a phantom and projection data that can be
generated for it using the same model. Normally, the model consists of a number of elementary :ref:`ref_glossary_object` (s). 

For reproducibility sake, one can re-use the same model and refer to it while testing an algorithm. There is no need to save the phantom 
on the disk as it can be quickly generated. There is a variety of models exist, for instance, 2D/3D versions of Shepp-Logan, Defrise, and QRM phantoms. 
The reserved numbers for stationary models should be in the range of [0-99] as higher numbers are reserved for :ref:`ref_glossary_dynamic_model`.

In terms of the API provided, please look at :func:`Model` in :func:`tomophantom.TomoP2D` for 2D phantoms and :func:`tomophantom.TomoP3D` for 3D phantoms, respectively.

.. _ref_glossary_dynamic_model:

Dynamic model
-------------

Dynamic phantoms or `time-lapse` phantoms are the extensions over the existing 2D and 3D phantoms that are provided as models in :ref:`ref_glossary_library`.
The are implemented to test algorithms that involve motion through a some kind of structural or intensity change in time. So far only simple linear 
changes can be programmed, but, in theory, this feature can be extended to accept more complex motion patterns. In addition, if projection data is simulated,
the change happens in-between scans, i.e., the motion is NOT encoded into projections. This can be also a desirable feature if more realistic change needs to be
modelled and tested through reconstruction algorithms. More information can be found in :ref:`ref_library_dynamic_extensions`.


.. _ref_glossary_object:

Object
------

An **object** is an elementary function that defines a phantom or contributes to the projection data. TomoPhantom offers a variety of simple geometrical objects,
through which a phantom is composed. Using objects allows you to build customised phantoms quicker, which later can be stored as a :ref:`ref_glossary_model` in 
:ref:`ref_glossary_library`. It also offers a possibility to build :ref:`ref_glossary_random_phantoms`, i.e., where parameters are changing stochastically. 

TomoPhantom provides a dedicated :ref:`ref_object_api`, if one would like to build customised phantoms without the use of the library files.


.. _ref_glossary_random_phantoms:

Random Phantoms
---------------

**Random phantoms** are the phantoms that are built using :ref:`ref_object_api` where parameters of the objects changing by sampling from some 
probability distribution. In TomoPhantom, there are functions that can help you to build phantoms that are packed with 
circular objects of different types, e.g., a 2D or 3D foam phantom. Please see :mod:`tomophantom.generator`.

.. _ref_glossary_library:

Library files
-------------

The library files are the text files which are provided to users with the downloaded repository or an installation. 
Parameter for different objects are stored in the files resulting in the stacked phantom to be built reproducibly. 
Specifically, 2D phantoms and (2D + time) *dynamic* phantoms are stored in the :data:`Phantom2DLibrary.dat` file, while 3D phantoms and (3D + time) dynamic
phantoms are stored in the :data:`Phantom3DLibrary.dat` file. See more how to :ref:`howto_libraries`. 

To understand how object parameters stored so that one can edit or add a new model, please
see :ref:`ref_library_files_api`.

The library files are located in the :data:`tomophantom/phantomlib` in your installed location
or in the `Github repository <https://github.com/dkazanc/TomoPhantom/tree/master/tomophantom/phantomlib>`_.


