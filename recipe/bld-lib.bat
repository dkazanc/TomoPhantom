mkdir "%SRC_DIR%\tomophantom"
ROBOCOPY /E "%RECIPE_DIR%\.." "%SRC_DIR%\tomophantom"
cd tomophantom

:: issue cmake to create setup.py
cmake %RECIPE_DIR%\..\ -DLIBRARY_INC="%CONDA_PREFIX%" -DCMAKE_INSTALL_PREFIX="%PREFIX%\Library"

cmake --build .
cmake --install .

::%PYTHON% setup.py build_ext
if errorlevel 1 exit 1
::%PYTHON% setup.py install
::if errorlevel 1 exit 1
