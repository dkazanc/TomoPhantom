
mkdir ${SRC_DIR}/tomophantom
cp -r "${RECIPE_DIR}/../../" ${SRC_DIR}/tomophantom

cd ${SRC_DIR}/tomophantom/python
$PYTHON setup.py install
