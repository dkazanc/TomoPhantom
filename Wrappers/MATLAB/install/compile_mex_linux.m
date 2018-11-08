% compile mex files once from Matlab running this script

fsep = '/';
UpPath = sprintf(['..' fsep], 1i);
cd(UpPath);

mkdir compiled

PathFunc = sprintf(['..' fsep '..' fsep 'Core' fsep], 1i);

copyfile(PathFunc, 'src');

cd('src');

Pathmove = sprintf(['..' fsep 'compiled' fsep], 1i);

fprintf('%s \n', 'Building functions...');
mex TomoP2DModel.c TomoP2DModel_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile('TomoP2DModel.mex*',Pathmove);
mex TomoP2DObject.c TomoP2DModel_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile('TomoP2DObject.mex*',Pathmove);
mex TomoP2DModelSino.c TomoP2DModelSino_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile('TomoP2DModelSino.mex*',Pathmove);
mex TomoP2DObjectSino.c TomoP2DModelSino_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile('TomoP2DObjectSino.mex*',Pathmove);
mex TomoP2DSinoNum.c TomoP2DSinoNum_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile('TomoP2DSinoNum.mex*',Pathmove);
mex TomoP3DModel.c TomoP3DModel_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile('TomoP3DModel.mex*',Pathmove);
mex TomoP3DObject.c TomoP3DModel_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile('TomoP3DObject.mex*',Pathmove);
mex TomoP3DModelSino.c TomoP3DModelSino_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile('TomoP3DModelSino.mex*',Pathmove);
delete TomoP2DModel_core* TomoP2DModelSino_core* TomoP3DModel_core* TomoP2DSinoNum_core* TomoP3DModelSino_core* CCPiDefines.h utils* CMakeLists.txt;

fprintf('%s \n', 'All compiled!');

cd(UpPath);
