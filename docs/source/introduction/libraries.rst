Phantom Libraries
*******************

TomoPhantom provides a variety of phantoms that are static (2D,3D) and dynamic one (3D, 4D).
The users can start using them right away by specifying a **model** number from 
the associated library file. Model is usually reffered to a combined phantom which 
consists of simpler geometrical **objects**. Note that objects can be 
also build independently, please see :ref:`ref_api`.

`Phantom2DLibrary.dat` and `Phantom3DLibrary.dat` are editable 
text files located in `tomophantom/phantomlib` folder with parametrised models.
For instance, 2D/3D versions of Shepp-Logan, Defrise, and QRM phantoms. 

The generation of new phantoms is highly encouraged, 
please submit them through the pull requests to the repo.