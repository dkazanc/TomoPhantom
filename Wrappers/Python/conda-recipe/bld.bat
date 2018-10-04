mkdir "%SRC_DIR%\tomophantom"
ROBOCOPY /E "%RECIPE_DIR%\.." "%SRC_DIR%\tomophantom"
cd tomophantom

:: issue cmake to create setup.py
cmake -G "NMake Makefiles" %RECIPE_DIR%\..\..\..\ -DBUILD_PYTHON_WRAPPERS=ON -DCONDA_BUILD=ON -DLIBRARY_LIB="%CONDA_PREFIX%\lib" -DLIBRARY_INC="%CONDA_PREFIX%" -DCMAKE_INSTALL_PREFIX="%PREFIX%\Library"
nmake install
::%PYTHON% setup.py build_ext
if errorlevel 1 exit 1
::%PYTHON% setup.py install
::if errorlevel 1 exit 1
