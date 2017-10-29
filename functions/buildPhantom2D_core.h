#include <matrix.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
#include "omp.h"

/*
* License: Apache Version 2.0
* Copyright {2017} {Daniil Kazantsev, The University of Manchester}
*/

float buildPhantom2D_core(float *A, int ModelSelected, int N, char *ModelParametersFilename);
float buildPhantom2D_core_single(float *A, int N,  int Object, float C0, float x0, float y0, float a, float b, float phi_rot);
