#!/bin/bash

PKG_NAME=libtomophantom
USER=httomo-team
OS=linux-64

mkdir ~/conda-bld
conda config --set anaconda_upload no
export CONDA_BLD_PATH=~/conda-bld

export CIL_VERSION=3.0.2
$CONDA/bin/conda build conda-recipe_library . -c httomo

# upload packages to conda
find $CONDA_BLD_PATH/$OS -name *.conda | while read file
do
    echo $file
    $CONDA/bin/anaconda -v --show-traceback -t $ANACONDA_API_TOKEN upload $file --force
done
