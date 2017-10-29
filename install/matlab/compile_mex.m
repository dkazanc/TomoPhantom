% compile mex's in Matlab once
cd ../../functions/

fprintf('%s \n', 'Building functions...');
mex buildPhantom2D.c buildPhantom2D_core.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
mex buildSino2D.c buildSino2D_core.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
mex buildPhantom3D.c buildPhantom3D_core.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
mex buildSino3D.c buildSino3D_core.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
mex DeformObject_C.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
fprintf('%s \n', 'All done!');

cd ../
