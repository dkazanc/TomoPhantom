#include <matrix.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
#include "omp.h"

/*
* License: Version 2.0
* Copyright {2017} {Daniil Kazantsev, The University of Manchester}
*/

float buildPhantom3D_core(float *A, int ModelSelected, int N, char *ModelParametersFilename);

