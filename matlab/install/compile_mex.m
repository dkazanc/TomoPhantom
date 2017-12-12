% compile mex files once from Matlab running this script

cd ../../functions/

fprintf('%s \n', 'Building functions...');
mex buildPhantom2D.c buildPhantom2D_core.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile buildPhantom2D.mexa64 ../matlab/install/
mex buildSino2D.c buildSino2D_core.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile buildSino2D.mexa64 ../matlab/install/
mex buildPhantom3D.c buildPhantom3D_core.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile buildPhantom3D.mexa64 ../matlab/install/
mex buildSino3D.c buildSino3D_core.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile buildSino3D.mexa64 ../matlab/install/
mex DeformObject_C.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
movefile DeformObject_C.mexa64 ../matlab/install/
fprintf('%s \n', 'All done!');

cd ../
cd matlab
