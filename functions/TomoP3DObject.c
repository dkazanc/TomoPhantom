/*
 * Copyright 2017 Daniil Kazantsev
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
#include "omp.h"

#include "TomoP3DModel_core.h"
#include "utils.h"

#define M_PI 3.14159265358979323846

/* Function to build 3D objects without using the library file
 * (MATLAB mex-wrapper)
 *
 * Input Parameters:
 * 1. parameters for 3D object
 * 2. Size of the 3D object
 *
 * VolumeSize in voxels (N x N x N)
 *
 * Output:
 * 1. The analytical phantom size of [N x N x N]
 *
 */
#define MAXCHAR 1000

void mexFunction(
        int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
        
{
    int N;
    float *A;
    char *tmpstr2;
    float C0 = 0.0f, x0 = 0.0f, y0 = 0.0f, z0 = 0.0f, a = 0.0f, b = 0.0f, c = 0.0f, psi_gr1 = 0.0f, psi_gr2 = 0.0f, psi_gr3 = 0.0f;
    
    tmpstr2 = mxArrayToString(prhs[0]); /* name of the object */
    C0 = (float) mxGetScalar(prhs[1]); /* intensity */
    x0 = (float) mxGetScalar(prhs[2]); /* x0 position */
    y0 = (float) mxGetScalar(prhs[3]); /* y0 position */
    z0 = (float) mxGetScalar(prhs[4]); /* z0 position */
    a = (float) mxGetScalar(prhs[5]); /* a - size object */
    b = (float) mxGetScalar(prhs[6]); /* b - size object */
    c = (float) mxGetScalar(prhs[7]); /* c - size object */
    psi_gr1 = (float) mxGetScalar(prhs[8]); /* rotation angle 1*/
    psi_gr2 = (float) mxGetScalar(prhs[9]); /* rotation angle 2*/
    psi_gr3 = (float) mxGetScalar(prhs[10]); /* rotation angle 3*/
    N  = (int) mxGetScalar(prhs[11]); /* choosen dimension (N x N x N) */
    
    /*Handling Matlab input data*/
    if (nrhs != 12) mexErrMsgTxt("Input of 12 parameters is required");   
    
    const mwSize N_dims[3] = {N, N, N}; /* image dimensions */
    A = (float*)mxGetPr(plhs[0] = mxCreateNumericArray(3, N_dims, mxSINGLE_CLASS, mxREAL)); /*output array*/
    
    
    if ((strcmp("gaussian",tmpstr2) != 0) && (strcmp("paraboloid",tmpstr2) != 0) && (strcmp("ellipsoid",tmpstr2) != 0) && (strcmp("cone",tmpstr2) != 0) && (strcmp("cuboid",tmpstr2) != 0) && (strcmp("ellipticalcylinder",tmpstr2) != 0) ) {
        printf("%s %s\n", "Unknown name of the object, the given name is", tmpstr2);
        mexErrMsgTxt("Unknown name of the object");
        }
    if (C0 == 0) {
        printf("%s %f\n", "C0 should not be equal to zero, the given value is", C0);
        mexErrMsgTxt("C0 should not be equal to zero");
        }
    if ((x0 < -1) || (x0 > 1)) {
        printf("%s %f\n", "x0 (object position) must be in [-1,1] range, the given value is", x0);
        mexErrMsgTxt("x0 (object position) must be in [-1,1] range");
        }
    if ((y0 < -1) || (y0 > 1)) {
        printf("%s %f\n", "y0 (object position) must be in [-1,1] range, the given value is", y0);
        mexErrMsgTxt("y0 (object position) must be in [-1,1] range");
        }
    if ((z0 < -1) || (z0 > 1)) {
        printf("%s %f\n", "z0 (object position) must be in [-1,1] range, the given value is", z0);
        mexErrMsgTxt("z0 (object position) must be in [-1,1] range");
        }
    if ((a <= 0) || (a > 2)) {
        printf("%s %f\n", "a (object size) must be positive in [0,2] range, the given value is", a);
        mexErrMsgTxt("a (object size) must be positive in [0,2] range");
        }
    if ((b <= 0) || (b > 2)) {
        printf("%s %f\n", "b (object size) must be positive in [0,2] range, the given value is", b);
        mexErrMsgTxt("b (object size) must be positive in [0,2] range");
        }
    if ((c <= 0) || (c > 2)) {
        printf("%s %f\n", "c (object size) must be positive in [0,2] range, the given value is", c);
        mexErrMsgTxt("c (object size) must be positive in [0,2] range");
        }
    printf("\nObject : %s \nC0 : %f \nx0 : %f \ny0 : %f \nz0 : %f \na : %f \nb : %f \nc : %f \n", tmpstr2, C0, x0, y0, z0, a, b, c);
    
    TomoP3DObject_core(A, N, tmpstr2, C0, x0, y0, z0, b, a, c, -psi_gr1, psi_gr2, psi_gr3, 0); /* Matlab */    
    mxFree(tmpstr2);
}
