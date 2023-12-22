===================================
How to build the html pages locally
===================================

Create a conda environment
==========================

Setup a new conda environment using the requirements file docs/source/doc-conda-requirements.yml

    >>> conda env create --prefix /path/to/env/doc-env --file /path/to/tomophantom/docs/source/doc-conda-requirements.yml


Update API documentation and build
==================================

While inside your virtual environment, run the sphinx-build.sh script.

    >>> conda activate /path/to/env/doc-env
    >>> source /path/to/tomophantom/docs/sphinx-build.sh

The script will:

1. Remove previous api and build directories.
2. Generate current tomophantom api files
3. Run a python script to add yaml file downloads for every function for every module rst file.
4. Run sphinx to create html documentation pages in docs/build.

You can view the completed pages by opening the tomophantom/docs/build/index.html page inside a browser.
