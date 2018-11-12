

mkdir -p ${SRC_DIR}/tomophantom
#cp -r "${RECIPE_DIR}/.." ${SRC_DIR}/tomophantom
cd tomophantom

#issue cmake to create setup.py
cmake -G "Unix Makefiles" ${RECIPE_DIR}/../../../ -DBUILD_PYTHON_WRAPPER=ON -DCONDA_BUILD=ON -DCMAKE_BUILD_TYPE="Release" -DLIBRARY_LIB=$CONDA_PREFIX/lib -DLIBRARY_INC=$CONDA_PREFIX -DCMAKE_INSTALL_PREFIX=$PREFIX
#cp -rv "${RECIPE_DIR}/../../../PhantomLibrary" ${SRC_DIR}/tomophantom/Wrappers/Python/tomophantom
make install
