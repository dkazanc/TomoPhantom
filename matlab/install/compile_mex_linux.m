% compile mex files once from Matlab running this script

fsep = '/';
UpPath = sprintf(['..' fsep], 1i);
cd(UpPath);

mkdir compiled

PathFunc = sprintf(['..' fsep 'functions' fsep], 1i);
cd(PathFunc);

Pathmove = sprintf(['..' fsep 'matlab' fsep 'compiled' fsep], 1i);

fprintf('%s \n', 'Building functions...');
mex TomoP2DModel.c TomoP2DModel_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile('TomoP2DModel.mex*',Pathmove);
mex TomoP2DObject.c TomoP2DModel_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile('TomoP2DObject.mex*',Pathmove);
mex TomoP2DModelSino.c TomoP2DModelSino_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile('TomoP2DModelSino.mex*',Pathmove);
mex TomoP2DObjectSino.c TomoP2DModelSino_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile('TomoP2DObjectSino.mex*',Pathmove);
mex TomoP3DModel.c TomoP3DModel_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile('TomoP3DModel.mex*',Pathmove);
mex TomoP3DObject.c TomoP3DModel_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile('TomoP3DObject.mex*',Pathmove);
fprintf('%s \n', 'All compiled!');

cd(UpPath);
cd matlab