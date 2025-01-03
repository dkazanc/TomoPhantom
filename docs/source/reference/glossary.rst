.. _ref_glossary:

Glossary
========

As some terms used in TomoPhantom can be ambiguous, we provide explanations on the terminology used. 

.. _ref_glossary_model:

Model
-----

The **model** term here describes a phantom and also associated projection data that can be generated for it. Normally model 
consists of the number of elementary geometrical **objects**. Note that objects can be 
also build independently, please see :ref:`ref_api`. For instance, 2D/3D versions of Shepp-Logan, Defrise, and QRM phantoms. 
:mod:`tomophantom.TomoP2D.model`.

.. _ref_glossary_dynamic_model:

Dynamic model
-------------

.. _ref_glossary_object:

Object
------

Object blabla


.. _ref_glossary_library:

Library files
-------------

The library files are the text files which are provided to users with the downloaded repository or an installation. Different phantoms are hardcoded using parameters 
for objects in those files. 2D phantoms and (2D + time) *dynamic* phantoms are stored in the :data:`Phantom2DLibrary.dat` file, while 3D phantoms and (3D + time)  dynamic
phantoms are stored in the :data:`Phantom3DLibrary.dat` file. The files are located in the :data:`tomophantom/phantomlib` folder.


