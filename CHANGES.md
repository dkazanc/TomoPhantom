# ChangeLog

## v1.4.3
* Artifacts simulation module has been modified. It is based now on specifying the dictionaries where artifact types are
described. Stipes can be simulated to be partial and with variable intensity.
* 2D Model 15 added - DLS phantom


## v1.4.2
* Jupiter Notebook demo added which uses deep learning algorithms to employ TomoPhantom for data generation

## v1.4.1
* All demos renamed and placed into categories: 2D/3D/4D, Random
* Random generation of 2D/3D phantoms is initiated using methods from "generator.py" script
* Foam phantoms can be built
* Demo 'RandPhantGen.py' shows how phantoms can be generated

## v1.4

* Artifacts generation have been re-written from the class-based structure to function based structure
* Each method in Artifacts can be called independently or all artifacts+noise can be simulated using _Artifacts_ method
* Seeding added for noise module
* Demos changed to adapt new artifacts functions and ToMoBAR changes

## v1.3

* Objects3D (generation of the customised objects) function is fixed
