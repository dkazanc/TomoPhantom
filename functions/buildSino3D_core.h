#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
#include "omp.h"

/*
* License: Apache Version 2.0
* Copyright {2017} {Daniil Kazantsev, The University of Manchester}
*/

float buildSino3D_core(float *A, int ModelSelected, int N, int P, float *Th, int AngTot, int CenTypeIn,char* ModelParametersFilename);
float buildSino3D_core_single(float *A, int N, int P, float *Th, int AngTot, int CenTypeIn, int Obj, float C0, float x0, float y0, float z0, float a, float b, float c, float phi_rot);