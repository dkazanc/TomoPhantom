#include <matrix.h>
#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
#include "omp.h"

#include "buildPhantom3D_core.h"

#define M_PI 3.14159265358979323846

/* Function to read from a file the required parameters to build 3D analytical model, see Phantom3DLibrary.dat to modify parameters
 * (MATLAB wrapper)
 *
 * Input Parameters:
 * 1. ModelNo - the model number from Phantom3DLibrary file
 * 2. VolumeSize in voxels (N x N x N)
 *
 * Output:
 * 1. The analytical phantom size of [N x N x N]
 *
 * License: Version 2.0
 * Copyright {2017} {Daniil Kazantsev, The University of Manchester}
 */

void mexFunction(
        int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
        
{
    int ModelSelected, N;
    float *A;
    
    /*Handling Matlab input data*/
    if (nrhs != 2) mexErrMsgTxt("Input of 2 parameters is required: model and dimension");
    
    ModelSelected  = (int) mxGetScalar(prhs[0]); /* selected model */
    N  = (int) mxGetScalar(prhs[1]); /* choosen dimension (N x N x N) */
    
    /*Handling Matlab output data*/
    int N_dims[] = {N, N, N};
    A = (float*)mxGetPr(plhs[0] = mxCreateNumericArray(3, N_dims, mxSINGLE_CLASS, mxREAL));   
    
    buildPhantom3D_core(A, ModelSelected, N, "models/Phantom3DLibrary.dat");
}

