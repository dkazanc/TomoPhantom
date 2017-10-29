#include <matrix.h>
#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
#include "omp.h"

#include "buildPhantom2D_core.h"

#define M_PI 3.14159265358979323846

/* Function to read from the file Phantom2DLibrary.dat the required parameters to build 2D analytical models
 * (MATLAB mex-wrapper)
 *
 * Input Parameters:
 * 1. ModelNo - a model number from Phantom2DLibrary file
 * 2. VolumeSize in voxels (N x N)
 * 3. An absolute path to the file Phantom2DLibrary.dat (see OS-specific syntax-differences)
 *
 * Output:
 * 1. The analytical phantom size of [N x N]
 *
 * License: Apache Version 2.0
 * Copyright {2017} {Daniil Kazantsev, The University of Manchester}
 */

void mexFunction(
        int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
        
{
    int ModelSelected, N;
    float *A;
    char *ModelParameters_PATH;
    
    /*Handling Matlab input data*/
    if (nrhs != 3) mexErrMsgTxt("Input of 3 parameters is required: model, dimension, PATH");
    
    ModelSelected  = (int) mxGetScalar(prhs[0]); /* selected model */
    N  = (int) mxGetScalar(prhs[1]); /* choosen dimension (N x N x N) */
    ModelParameters_PATH = mxArrayToString(prhs[2]); /* provide an absolute path to the file */      
    
    /*Handling Matlab output data*/
    int N_dims[] = {N, N};
    A = (float*)mxGetPr(plhs[0] = mxCreateNumericArray(2, N_dims, mxSINGLE_CLASS, mxREAL));   
    
    buildPhantom2D_core(A, ModelSelected, N, ModelParameters_PATH);
    
    mxFree(ModelParameters_PATH);
}

