# ChangeLog

## v1.4.9
* Artifacts function syntax has changed. Now the dictionaries with keys can be be provided in a random order. Demos updated to make use of the updated function
* Partial volume effect generator added to the Artifacts script
* Fresnel propagator simulator added to the Artifacts script

## v1.4.7
* Flat field simulator has been modified to be based on a speckle generator of the background
* Flat fields and the background jitter has been incorporated
* Model 17 has been added to simulate i23 data

## v1.4.5
* Flat field simulation for 3D case has been re-written to conform the traditional imaging scenario. One can modify the
inensity of the X-ray source which leads to more imaging artefacts
* Demo ReconASTRA3D_realistic shows the incorporated changes
* normraw function has been deleted and replaced with the conventional normaliser from the ToMoBAR software

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
* Each method in Artifacts can be called independently or all artefacts+noise can be simulated using _Artifacts_ method
* Seeding added for noise module
* Demos changed to adapt new artefacts functions and ToMoBAR changes

## v1.3

* Objects3D (generation of the customised objects) function is fixed
