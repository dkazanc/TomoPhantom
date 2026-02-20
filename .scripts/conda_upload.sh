#!/bin/bash

PKG_NAME=tomophantom
USER=httomo-team
OS=noarch

mkdir ~/conda-bld
conda config --set anaconda_upload no
export CONDA_BLD_PATH=~/conda-bld

export CIL_VERSION=3.0.2
$CONDA/bin/conda build conda-recipe . -c httomo

# upload packages to conda
find $CONDA_BLD_PATH/$OS -name *.conda | while read file
do
    echo $file
    $CONDA/bin/anaconda -v --show-traceback -t $ANACONDA_API_TOKEN upload $file --force
done
