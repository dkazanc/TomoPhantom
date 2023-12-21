set -xe 
cd $SRC_DIR

cmake -G "Unix Makefiles" ${RECIPE_DIR}/../ -DBUILD_PYTHON_WRAPPER=OFF -DCONDA_BUILD=ON -DCMAKE_BUILD_TYPE="Release" -DLIBRARY_LIB=$CONDA_PREFIX/lib -DLIBRARY_INC=$CONDA_PREFIX -DCMAKE_INSTALL_PREFIX=$PREFIX

cmake --build .
cmake --install .

#$PYTHON -m pip install . --no-build-isolation \
#    --no-deps --ignore-installed --no-index --no-cache-dir -vv
