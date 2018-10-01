% compile mex files once from Matlab running this script

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% I've been able to compile on Windows 7 with MinGW and Matlab 2016b, however, 
% not sure if openmp is enabled after the compilation. 

% Here I present two ways how software can be compiled, if you have some
% other suggestions please contact me at dkazanc@hotmail.com 
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>

fsep = '/';
UpPath = sprintf(['..' fsep], 1i);
cd(UpPath);

mkdir compiled

PathFunc = sprintf(['..' fsep 'functions' fsep], 1i);
cd(PathFunc);
Pathmove = sprintf(['..' fsep 'matlab' fsep 'compiled' fsep], 1i);

% One way to compile using Matlab-installed MinGW: 
fprintf('%s \n', 'Building functions...');
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

%%% The second approach to compile using TDM-GCC which follows this
%%% discussion:
%%% https://uk.mathworks.com/matlabcentral/answers/279171-using-mingw-compiler-and-open-mp#comment_359122
%%% 1. Install TDM-GCC independently from http://tdm-gcc.tdragon.net/ (I installed 5.1.0)
%%% Install openmp version: http://sourceforge.net/projects/tdm-gcc/files/TDM-GCC%205%20series/5.1.0-tdm64-1/gcc-5.1.0-tdm64-1-openmp.zip/download
%%% 2. Link til libgomp.a in that installation when compilling your mex file.

%%% assuming you unzipped TDM GCC (OpenMp) in folder TDMGCC on C drive, uncomment
%%% bellow
% fprintf('%s \n', 'Building functions...');
% mex C:\TDMGCC\lib\gcc\x86_64-w64-mingw32\5.1.0\libgomp.a CXXFLAGS="$CXXFLAGS -std=c++11 -fopenmp" TomoP2DModel.c TomoP2DModel_core.c utils.c
% movefile('TomoP2DModel.mex*',Pathmove);
% mex C:\TDMGCC\lib\gcc\x86_64-w64-mingw32\5.1.0\libgomp.a CXXFLAGS="$CXXFLAGS -std=c++11 -fopenmp" TomoP2DObject.c TomoP2DModel_core.c utils.c
% movefile('TomoP2DObject.mex*',Pathmove);
% mex C:\TDMGCC\lib\gcc\x86_64-w64-mingw32\5.1.0\libgomp.a CXXFLAGS="$CXXFLAGS -std=c++11 -fopenmp" TomoP2DModelSino.c TomoP2DModelSino_core.c utils.c
% movefile('TomoP2DModelSino.mex*',Pathmove);
% mex C:\TDMGCC\lib\gcc\x86_64-w64-mingw32\5.1.0\libgomp.a CXXFLAGS="$CXXFLAGS -std=c++11 -fopenmp" TomoP2DObjectSino.c TomoP2DModelSino_core.c utils.c
% movefile('TomoP2DObjectSino.mex*',Pathmove);
% mex C:\TDMGCC\lib\gcc\x86_64-w64-mingw32\5.1.0\libgomp.a CXXFLAGS="$CXXFLAGS -std=c++11 -fopenmp" TomoP3DModel.c TomoP3DModel_core.c utils.c
% movefile('TomoP3DModel.mex*',Pathmove);
% mex C:\TDMGCC\lib\gcc\x86_64-w64-mingw32\5.1.0\libgomp.a CXXFLAGS="$CXXFLAGS -std=c++11 -fopenmp" TomoP3DObject.c TomoP3DModel_core.c utils.c
% movefile('TomoP3DObject.mex*',Pathmove);
% fprintf('%s \n', 'All compiled!');

cd(UpPath);
cd matlab
