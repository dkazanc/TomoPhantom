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

#include "TomoP2DModelSino_core.h"
#include "utils.h"

/* Function to build 2D sinogram objects without using the library file
 * (MATLAB mex-wrapper)
 *
 * Input Parameters:
 * 1. parameters for 2D sinogram
 * 3. 
 * 2. Size of the 2D object
 *
 * VolumeSize in voxels (N x N )
 *
 * Output:
 * 1. The analytical phantom size of [N x N]
 *
 */
#define MAXCHAR 1000

void mexFunction(
        int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
        
{
    mwSize  NStructElems;
    int N, P;
    float *A, *Th;
    char *tmpstr2;
    float C0 = 0.0f, x0 = 0.0f, y0 = 0.0f, a = 0.0f, b = 0.0f, psi_gr1 = 0.0f;
    
    tmpstr2 = mxArrayToString(prhs[0]); /* name of the object */
    C0 = (float) mxGetScalar(prhs[1]); /* intensity */
    x0 = (float) mxGetScalar(prhs[2]); /* x0 position */
    y0 = (float) mxGetScalar(prhs[3]); /* y0 position */
    a = (float) mxGetScalar(prhs[4]); /* a - size object */
    b = (float) mxGetScalar(prhs[5]); /* b - size object */
    psi_gr1 = (float) mxGetScalar(prhs[6]); /* rotation angle 1*/        
    N  = (int) mxGetScalar(prhs[7]); /* choosen dimension (N x N) */    
    P  = (int) mxGetScalar(prhs[8]); /* detector size */
    Th  = (float*) mxGetData(prhs[9]); /* angles */
    NStructElems = mxGetNumberOfElements(prhs[9]);
    
    int CenTypeIn = 1; /* astra-type centering is the default one */       
    
    /*Handling Matlab input data*/
    if (nrhs != 10) mexErrMsgTxt("Input of 10 parameters is required");   
    
    const mwSize N_dims[2] = {NStructElems,P}; /*format: X-detectors, Y-angles dim*/
    A = (float*)mxGetPr(plhs[0] = mxCreateNumericArray(2, N_dims, mxSINGLE_CLASS, mxREAL)); 
    
    if ((strcmp("gaussian",tmpstr2) != 0) && (strcmp("parabola",tmpstr2) != 0) && (strcmp("ellipse",tmpstr2) != 0) && (strcmp("parabola1",tmpstr2) != 0) && (strcmp("cone",tmpstr2) != 0) && (strcmp("rectangle",tmpstr2) != 0) ) {
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
    if ((a <= 0) || (a > 2)) {
        printf("%s %f\n", "a (object size) must be positive in [0,2] range, the given value is", a);
        mexErrMsgTxt("a (object size) must be positive in [0,2] range");
        }
    if ((b <= 0) || (b > 2)) {
        printf("%s %f\n", "b (object size) must be positive in [0,2] range, the given value is", b);
        mexErrMsgTxt("b (object size) must be positive in [0,2] range");
        }    
    printf("\nObject : %s \nC0 : %f \nx0 : %f \ny0 : %f \na : %f \nb : %f \n", tmpstr2, C0, x0, y0, a, b);        

    TomoP2DObjectSino_core(A, N, P, Th, (int)NStructElems, CenTypeIn, tmpstr2, C0, x0, y0, b, a, -psi_gr1, 0); /* Matlab */
    
    mxFree(tmpstr2);
}
