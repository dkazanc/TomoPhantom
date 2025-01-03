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

Library files are the text files which are provided to users with the downloaded repository or installation. Different phantoms are hardcoded using parameters 
in the files. 2D phantoms and their dynamic extenstions are stored in the :data:`Phantom2DLibrary.dat` file, while 3D phantoms and their dynamic
extensions are stored in the :data:`Phantom3DLibrary.dat` file. The files are located in the `tomophantom/phantomlib` folder.


