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

#define M_PI 3.14159265358979323846

/* MATLAB mex-wrapper to generate parallel beam sinograms using TomoP2DModelSino_core.c
 *
 * Input Parameters:
 * 1. Model number (see Phantom2DLibrary.dat) [required]
 * 2. ImageSize in pixels (N x N) [required]
 * 3. Detector array size P (in pixels) [required]
 * 4. Projection angles Theta (in degrees) [required]
 * 5. An absolute path to the file Phantom2DLibrary.dat (see OS-specific syntax-differences) [required]
 * 6. ImageCentring, choose 'radon' or 'astra' (default) [optional]
 *
 * Output:
 * 1. 2D sinogram size of [length(angles), P] or a temporal sinogram size of [length(angles), P, Time-Frames]
 */

void mexFunction(
        int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
        
{
    int ModelSelected, N, CenTypeIn, P;
    float *A, *Th;
    mwSize  NStructElems;
    char *ModelParameters_PATH;
    
    /*Handling Matlab input data*/
    if ((nrhs < 5) && (nrhs > 6)) mexErrMsgTxt("Input of 5 or 6 parameters is required: Model, ImageSize, Detector Size, Projection angles, PATH, Centering");
    if (mxGetClassID(prhs[3]) != mxSINGLE_CLASS) {mexErrMsgTxt("The vector of angles must be in a single precision"); }
    
    ModelSelected  = (int) mxGetScalar(prhs[0]); /* selected model */
    N  = (int) mxGetScalar(prhs[1]); /* choosen dimension (N x N) */
    P  = (int) mxGetScalar(prhs[2]); /* detector size */
    Th  = (float*) mxGetData(prhs[3]); /* angles */
    ModelParameters_PATH = mxArrayToString(prhs[4]); /* provide an absolute path to the file */
    CenTypeIn = 1; /* astra-type centering is the default one */
    
    if (nrhs == 6)  {
        char *CenType;
        CenType = mxArrayToString(prhs[5]); /* 'radon' or 'astra' (default) */
        if ((strcmp(CenType, "radon") != 0) && (strcmp(CenType, "astra") != 0)) mexErrMsgTxt("Choose 'radon' or 'astra''");
        if (strcmp(CenType, "radon") == 0)  CenTypeIn = 0;  /* enable 'radon'-type centaering */
        mxFree(CenType);
    }
    NStructElems = mxGetNumberOfElements(prhs[3]);
    
    /* first to check if the model is stationary or temporal */
    int *timeFrames = (int *) mxGetData(mxCreateNumericMatrix(1, 1, mxINT16_CLASS, mxREAL));
    /* extract the number of time-step */
    extractTimeFrames(timeFrames, ModelSelected, ModelParameters_PATH);    
    
    int *params_switch = (int *) mxGetData(mxCreateNumericMatrix(1, 10, mxINT16_CLASS, mxREAL));
    /* run checks for parameters (only integer switches and the number of time-frames) */
    checkParams2D(params_switch, ModelSelected, ModelParameters_PATH, timeFrames[0]);
    
    if (params_switch[0] == 0) mexErrMsgTxt("Parameters file (Phantom2DLibrary.dat) cannot be read, check that your PATH to the file is a valid one and the file exists");
    if (params_switch[1] == 0) {
        printf("%s %i\n", ">>>>> No such model available in the library, the given model is", ModelSelected);
        mexErrMsgTxt("Check the syntax in Phantom2DLibrary.dat");
    }
    if (params_switch[2] == 0) {
        printf("%s %i\n", ">>>>> Components number cannot be negative for the given model", ModelSelected);
        mexErrMsgTxt("Check the syntax in Phantom2DLibrary.dat");
    }
    if (params_switch[3] == 0) {
        printf("%s %i\n", ">>>>> TimeSteps value cannot be negative for the given model", ModelSelected);
        mexErrMsgTxt("Check the syntax in Phantom2DLibrary.dat");
    }
    if (params_switch[4] == 0) {
        printf("%s %i\n", ">>>>> Unknown name of the object for the given model", ModelSelected);
        mexErrMsgTxt("Check the syntax in Phantom2DLibrary.dat");
    }
    if (params_switch[5] == 0) {
        printf("%s %i\n", ">>>>> C0 should not be equal to zero for the given model", ModelSelected);
        mexErrMsgTxt("Check the syntax in Phantom2DLibrary.dat");
    }
    if (params_switch[6] == 0) {
        printf("%s %i\n", ">>>>> x0 (object position) must be in [-1,1] range for the given model", ModelSelected);
        mexErrMsgTxt("Check the syntax in Phantom2DLibrary.dat");
    }
    if (params_switch[7] == 0) {
        printf("%s %i\n", ">>>>> y0 (object position) must be in [-1,1] range for the given model", ModelSelected);
        mexErrMsgTxt("Check the syntax in Phantom2DLibrary.dat");
    }
    if (params_switch[8] == 0) {
        printf("%s %i\n", ">>>>> a (object size) must be positive in [0,2] range", ModelSelected);
        mexErrMsgTxt("Check the syntax in Phantom2DLibrary.dat");
    }
    if (params_switch[9] == 0) {
        printf("%s %i\n", ">>>>> b (object size) must be positive in [0,2] range", ModelSelected);
        mexErrMsgTxt("Check the syntax in Phantom2DLibrary.dat");
    }
    
    /*Handling Matlab output data*/
    if (params_switch[3] == 1) {
        /* static phantom case */
        int N_dims[] = {NStructElems,P}; /*format: X-detectors, Y-angles dim*/
        A = (float*)mxGetPr(plhs[0] = mxCreateNumericArray(2, N_dims, mxSINGLE_CLASS, mxREAL)); }
    else {
        int N_dims[] = {NStructElems, P, params_switch[3]}; /*format: X-detectors, Y-angles dim, Time-Frames*/
        A = (float*)mxGetPr(plhs[0] = mxCreateNumericArray(3, N_dims, mxSINGLE_CLASS, mxREAL)); }
    TomoP2DModelSino_core(A, ModelSelected, N, P, Th, (int)NStructElems, CenTypeIn, ModelParameters_PATH, 0);
    mxFree(ModelParameters_PATH);
}