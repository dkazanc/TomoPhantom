% compile mex's in Matlab once
cd ../../functions/

mex buildPhantom3D.c buildPhantom3D_core.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"

cd ../
