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
#include <matrix.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
#include <malloc.h>
#include "omp.h"

#include "TomoP2DSinoNum_core.h"
#include "utils.h"


/* C OMP implementation of the forward projection (The Radon Transform)
 * by rotating a padded image and summing over columns (rotation-based projector 
 * for parallel beam)
 *
 * Input Parameters:
 * 1. Phantom to calculate projections from [required]
 * 2. Detector array size P (in pixels) [required]
 * 3. Projection angles Theta (in degrees) [required]
 *
 * Output:
 * Sinogram [No. angles] x [No. detectors]
 *
 * mex ForwardProjNum.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
 * sinogram = TomoP2DSinoNum(single(G), P, single(angles));
 */

void mexFunction(
        int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    int dimX, dimY,ThetaLength, DetSize, sys=1;
    const int  *dim_array, *dimsTh;
    float *Phantom, *Sinogram, *Theta;    
        
    /*Handling Matlab input data*/
    Phantom  = (float *) mxGetData(prhs[0]);
    DetSize  = (int) mxGetScalar(prhs[1]);    
    Theta  = (float *) mxGetData(prhs[2]);
    
    dim_array = mxGetDimensions(prhs[0]);
    dimsTh = mxGetDimensions(prhs[2]);
    
    dimX = dim_array[0]; dimY = dim_array[1]; /*Image Dimensions (must be squared) */
    if( dimX != dimY) mexErrMsgTxt("The input image must be squared");
    if( DetSize < dimX) mexErrMsgTxt("The detector dimension is smaller than the phantom dimension! ");
    
    ThetaLength = dimsTh[1]; /* angles dimension */           
        
    /*Handling Matlab output data*/
    const mwSize N_dims[2] = {ThetaLength, DetSize}; /*format: Y-angles x X-detectors dim*/
    Sinogram = (float*)mxGetPr(plhs[0] = mxCreateNumericArray(2, N_dims, mxSINGLE_CLASS, mxREAL));
    
    TomoP2DSinoNum_core(Sinogram, Phantom, dimX, DetSize, Theta, ThetaLength, sys);
}
