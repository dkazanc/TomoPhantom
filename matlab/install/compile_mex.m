% compile mex files once from Matlab running this script

fsep = '/';
UpPath = sprintf(['..' fsep], 1i);
cd(UpPath);

mkdir compiled

PathFunc = sprintf(['..' fsep 'functions' fsep], 1i);
cd(PathFunc);

Pathmove = sprintf(['..' fsep 'matlab' fsep 'compiled' fsep], 1i);

fprintf('%s \n', 'Building functions...');
% mex TomoP2DModel.c TomoP2DModel_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
% movefile('TomoP2DModel.mex*',Pathmove);
% mex TomoP2DObject.c TomoP2DModel_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
% movefile('TomoP2DObject.mex*',Pathmove);
% mex TomoP2DModelSino.c TomoP2DModelSino_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
% movefile('TomoP2DModelSino.mex*',Pathmove);
% mex TomoP2DObjectSino.c TomoP2DModelSino_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
% movefile('TomoP2DObjectSino.mex*',Pathmove);
% mex TomoP3DModel.c TomoP3DModel_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
% movefile('TomoP3DModel.mex*',Pathmove);
% mex TomoP3DObject.c TomoP3DModel_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
% movefile('TomoP3DObject.mex*',Pathmove);
mex TomoP2DModel.c TomoP2DModel_core.c utils.c COMPFLAGS="\$COMPFLAGS -fopenmp -Wall -std=c99"
movefile('TomoP2DModel.mex*',Pathmove);
mex TomoP2DObject.c TomoP2DModel_core.c utils.c COMPFLAGS="\$COMPFLAGS -fopenmp -Wall -std=c99"
movefile('TomoP2DObject.mex*',Pathmove);
mex TomoP2DModelSino.c TomoP2DModelSino_core.c utils.c COMPFLAGS="\$COMPFLAGS -fopenmp -Wall -std=c99"
movefile('TomoP2DModelSino.mex*',Pathmove);
mex TomoP2DObjectSino.c TomoP2DModelSino_core.c utils.c COMPFLAGS="\$COMPFLAGS -fopenmp -Wall -std=c99"
movefile('TomoP2DObjectSino.mex*',Pathmove);
mex TomoP3DModel.c TomoP3DModel_core.c utils.c COMPFLAGS="\$COMPFLAGS -fopenmp -Wall -std=c99"
movefile('TomoP3DModel.mex*',Pathmove);
mex TomoP3DObject.c TomoP3DModel_core.c utils.c COMPFLAGS="\$COMPFLAGS -fopenmp -Wall -std=c99"
movefile('TomoP3DObject.mex*',Pathmove);

%1. Install TDM-GCC independently from http://tdm-gcc.tdragon.net/ (I installed 5.1.0)
%2. Link til libgomp.a in that installation when compilling your mex file.
% https://uk.mathworks.com/matlabcentral/answers/279171-using-mingw-compiler-and-open-mp#comment_359122
% mex C:\'Program Files'\TDM\lib\gcc\x86_64-w64-mingw32\5.1.0\libgomp.a CXXFLAGS="$CXXFLAGS -std=c++11 -fopenmp" TomoP3DModel.c TomoP3DModel_core.c utils.c
% mex BackProjCPU.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
% movefile BackProjCPU.mex* ../matlab/compiled/
fprintf('%s \n', 'All compiled!');

cd(UpPath);
cd matlab