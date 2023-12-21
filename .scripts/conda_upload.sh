#!/bin/bash

PKG_NAME=tomophantom
USER=httomo-team
OS=noarch
CONDA_TOKEN=$(cat $HOME/.secrets/my_secret.json)

mkdir ~/conda-bld
conda config --set anaconda_upload no
export CONDA_BLD_PATH=~/conda-bld

export CIL_VERSION=3.0
$CONDA/bin/conda build conda-recipe . -c httomo

# upload packages to conda
find $CONDA_BLD_PATH/$OS -name *.tar.bz2 | while read file
do
    echo $file
    $CONDA/bin/anaconda -v --show-traceback --token $CONDA_TOKEN upload $file --force
done
