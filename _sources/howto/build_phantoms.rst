.. _howto_buildphantoms:

Build Phantoms
==============

This sections explains how to build phantoms by presenting two main ways of doing that. 

* 1. The easiest and the quickest way to build a phantom is to use provided library files, see more in :ref:`howto_libraries`.  This is the most suitable way to start for someone who is new to TomoPhantom. We recommend to try building already pre-existing models by following this guide on how to generate :ref:`tutorial_model`.

* 2. For more advanced approach in building bespoke phantoms consisting of multiple objects, we recommend to follow the guidelines of :ref:`ref_object_api`. This is how one can build a phantom directly, i.e., without the use of the library file. One can specify parameters as a dictionary in Python and pass into a function that works with the objects. This can be more versatile way of building phantoms when one need to update or manipulate with parameters of objects. We recommend to have a look at the tutorial how to generate :ref:`tutorial_object`.