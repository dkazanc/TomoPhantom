mkdir "%SRC_DIR%\tomophantom"
xcopy /e "%RECIPE_DIR%\..\.." "%SRC_DIR%\tomophantom"
cd tomophantom\python

%PYTHON% setup.py build_ext
if errorlevel 1 exit 1
%PYTHON% setup.py install
if errorlevel 1 exit 1
