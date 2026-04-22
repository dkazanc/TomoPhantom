set -xe
cd $SRC_DIR

cmake ${RECIPE_DIR}/../ -DLIBRARY_INC=$CONDA_PREFIX -DCMAKE_INSTALL_PREFIX=$PREFIX

cmake --build .
cmake --install .

#$PYTHON -m pip install . --no-build-isolation \
#    --no-deps --ignore-installed --no-index --no-cache-dir -vv
