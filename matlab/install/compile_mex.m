% compile mex files once from Matlab running this script

cd ../
cd ../functions/

fprintf('%s \n', 'Building functions...');
mex buildPhantom2D.c buildPhantom2D_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile buildPhantom2D.mexa64 ../matlab/compiled/
mex buildSino2D.c buildSino2D_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile buildSino2D.mexa64 ../matlab/compiled/
mex buildPhantom3D.c buildPhantom3D_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile buildPhantom3D.mexa64 ../matlab/compiled/
mex buildSino3D.c buildSino3D_core.c utils.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile buildSino3D.mexa64 ../matlab/compiled/
mex DeformObject_C.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile DeformObject_C.mexa64 ../matlab/compiled/
fprintf('%s \n', 'All compiled!');

cd ../
cd matlab
