#!/bin/bash

is_pr=0
while getopts :p option; do
  case "${option}" in
  p) is_pr=1 ;; # pull request
  *) ;;
  esac
done

PKG_NAME=tomophantom
USER=httomo-team
OS=noarch

mkdir ~/conda-bld
conda config --set anaconda_upload no
export CONDA_BLD_PATH=~/conda-bld

export CIL_VERSION=3.1.0
$CONDA/bin/conda build conda-recipe . -c httomo
$CONDA/bin/conda install anaconda-client

# upload packages to conda
if [ $is_pr -eq 1 ]; then
    echo "Pull request detected, skipping conda upload."
    exit 0
else
    # upload packages to conda
    find $CONDA_BLD_PATH/$OS -name *.conda | while read file
    do
        echo $file
        $CONDA/bin/anaconda -v --show-traceback -t $ANACONDA_API_TOKEN upload $file --force
    done

fi
