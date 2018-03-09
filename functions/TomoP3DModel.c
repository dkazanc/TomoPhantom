/*
Copyright 2017 Daniil Kazantsev

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
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

/* Function to read from the file Phantom3DLibrary.dat the required parameters to build 3D or 3D + time analytical models
 * (MATLAB mex-wrapper)
 *
 * Input Parameters:
 * 1. ModelNo - a model number from Phantom3DLibrary file
 * 2. VolumeSize in voxels (N x N x N) or (N x N x N x time-frames)
 * 3. An absolute path to the file Phantom3DLibrary.dat (ssee OS-specific syntax-differences)
 *
 * Output:
 * 1. The analytical phantom size of [N x N x N] or a temporal 4D phantom (N x N x N x time-frames)
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
    
    /* first to check if the model is stationary or temporal */
    int *timeFrames = (int *) mxGetData(mxCreateNumericMatrix(1, 1, mxINT16_CLASS, mxREAL));
    /* extract the number of time-step */
    extractTimeFrames(timeFrames, ModelSelected, ModelParameters_PATH);    
    
    int *params_switch = (int *) mxGetData(mxCreateNumericMatrix(1, 12, mxINT16_CLASS, mxREAL));
    /* run checks for parameters (only integer switches and the number of time-frames) */
    checkParams3D(params_switch, ModelSelected, ModelParameters_PATH, timeFrames[0]);    
    
    if (params_switch[0] == 0) mexErrMsgTxt("Parameters file (Phantom3DLibrary.dat) cannot be read, check that your PATH to the file is a valid one and the file exists");
    if (params_switch[1] == 0) {
        printf("%s %i\n", ">>>>> No such model available in the library, the given model is", ModelSelected);
        mexErrMsgTxt("Check the syntax in Phantom3DLibrary.dat");
    }
    if (params_switch[2] == 0) {
        printf("%s %i\n", ">>>>> Components number cannot be negative for the given model", ModelSelected);
        mexErrMsgTxt("Check the syntax in Phantom3DLibrary.dat");
    }
    if (params_switch[3] == 0) {
        printf("%s %i\n", ">>>>> TimeSteps value cannot be negative for the given model", ModelSelected);
        mexErrMsgTxt("Check the syntax in Phantom3DLibrary.dat");
    }
    if (params_switch[4] == 0) {
        printf("%s %i\n", ">>>>> Unknown name of the object for the given model", ModelSelected);
        mexErrMsgTxt("Check the syntax in Phantom3DLibrary.dat");
    }
    if (params_switch[5] == 0) {
        printf("%s %i\n", ">>>>> C0 should not be equal to zero for the given model", ModelSelected);
        mexErrMsgTxt("Check the syntax in Phantom3DLibrary.dat");
    }
    if (params_switch[6] == 0) {
        printf("%s %i\n", ">>>>> x0 (object position) must be in [-1,1] range for the given model", ModelSelected);
        mexErrMsgTxt("Check the syntax in Phantom3DLibrary.dat");
    }
    if (params_switch[7] == 0) {
        printf("%s %i\n", ">>>>> y0 (object position) must be in [-1,1] range for the given model", ModelSelected);
        mexErrMsgTxt("Check the syntax in Phantom3DLibrary.dat");
    }
    if (params_switch[8] == 0) {
        printf("%s %i\n", ">>>>> z0 (object position) must be in [-1,1] range for the given model", ModelSelected);
        mexErrMsgTxt("Check the syntax in Phantom3DLibrary.dat");
    }    
    if (params_switch[9] == 0) {
        printf("%s %i\n", ">>>>> a (object size) must be positive in [0,2] range", ModelSelected);
        mexErrMsgTxt("Check the syntax in Phantom3DLibrary.dat");
    }
    if (params_switch[10] == 0) {
        printf("%s %i\n", ">>>>> b (object size) must be positive in [0,2] range", ModelSelected);
        mexErrMsgTxt("Check the syntax in Phantom3DLibrary.dat");
    }
    if (params_switch[11] == 0) {
        printf("%s %i\n", ">>>>> c (object size) must be positive in [0,2] range", ModelSelected);
        mexErrMsgTxt("Check the syntax in Phantom3DLibrary.dat");
    }     
    
    /*Handling Matlab output data*/
    if (params_switch[3] == 1) {    
    /* static phantom case */
    int N_dims[] = {N, N, N};
    A = (float*)mxGetPr(plhs[0] = mxCreateNumericArray(3, N_dims, mxSINGLE_CLASS, mxREAL));}
    else {
    /* temporal phantom 3D + time (4D) */   
    int N_dims[] = {N, N, N, params_switch[3]};
    A = (float*)mxGetPr(plhs[0] = mxCreateNumericArray(4, N_dims, mxSINGLE_CLASS, mxREAL));}   
    //TomoP3DModel_core(A, ModelSelected, N, ModelParameters_PATH, 0);    
    mxFree(ModelParameters_PATH);
}
