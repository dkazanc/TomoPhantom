% compile mex files once from Matlab running this script

cd ../
mkdir compiled
cd ../functions/

fprintf('%s \n', 'Building functions...');
mex TomoP2DModel.c TomoP2DModel_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile TomoP2DModel.mex* ../matlab/compiled/
mex TomoP2DObject.c TomoP2DModel_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile TomoP2DObject.mex* ../matlab/compiled/
mex TomoP2DModelSino.c TomoP2DModelSino_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile TomoP2DModelSino.mex* ../matlab/compiled/
mex TomoP2DObjectSino.c TomoP2DModelSino_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile TomoP2DObjectSino.mex* ../matlab/compiled/
mex TomoP3DModel.c TomoP3DModel_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile TomoP3DModel.mex* ../matlab/compiled/
mex TomoP3DObject.c TomoP3DModel_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile TomoP3DObject.mex* ../matlab/compiled/
mex BackProjCPU.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile BackProjCPU.mex* ../matlab/compiled/
fprintf('%s \n', 'All compiled!');

cd ../
cd matlab
