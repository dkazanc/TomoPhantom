.. _howto_libraries:

Use Phantom Libraries
=====================

TomoPhantom provides a number of phantoms for users that can be used `out-of-the-box`.
There are static (2D, 3D) or dynamic (2D+time, 3D+time) phantoms with the parameters for them stored in :ref:`ref_glossary_library`. 

The users can generate phantoms by specifying a :ref:`ref_glossary_model` number from the associated library file. 
Alternatively, more experienced user can build a phantom using :ref:`ref_glossary_object` API.


.. _howto_2d_libs:

Library file for 2D phantoms
----------------------------

.. dropdown:: The library file with 2D models + dynamic extensions

    .. literalinclude:: ../../../tomophantom/phantomlib/Phantom2DLibrary.dat


.. _howto_3d_libs:

Library file for 3D phantoms
----------------------------

.. dropdown:: The library file with 3D models + dynamic extensions

    .. literalinclude:: ../../../tomophantom/phantomlib/Phantom3DLibrary.dat


See more on :ref:`ref_library_files_api` where library files syntax explained. It is also useful if one is interested in adding new or modifying the existing phantoms [#f1]_. 

.. rubric:: Footnotes

.. [#f1] The generation of new phantoms is highly encouraged and the resulting models can be added to the repo. Please submit a pull requests for that through the Github page of the repository.

