#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
#include "omp.h"

#include "buildSino3D_core.h"

#define M_PI 3.14159265358979323846

/* Function to create 3D analytical sinograms (parallel beam geometry) to 3D phantoms using Phantom3DLibrary.dat
 * (MATLAB wrapper)
 *
 * Input Parameters:
 * 1. Model number (see Phantom3DLibrary.dat) [required]
 * 2. VolumeSize in voxels (N x N x N) [required]
 * 3. Detector array size P (in pixels) [required]
 * 4. Projection angles Th (in degrees) [required]
 * 5. VolumeCentring, choose 'radon' or 'astra' (default) [optional]
 *
 * Output:
 * 1. 3D sinogram size of [P, length(Th), N]
 *
 * License: Apache Version 2.0
 * Copyright {2017} {Daniil Kazantsev, The University of Manchester}
 */

void mexFunction(
        int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
        
{
    int ModelSelected, N, CenTypeIn, P;
    float *A, *Th;
    mwSize  NStructElems;
    
    /*Handling Matlab input data*/
    if ((nrhs < 4) && (nrhs > 5)) mexErrMsgTxt("Input of 4 or 5 parameters is required: Model, VolumeSize, Detector Size, Projection angles");
    if (mxGetClassID(prhs[3]) != mxSINGLE_CLASS) {mexErrMsgTxt("The vector of angles must be in a single precision"); }
    
    ModelSelected  = (int) mxGetScalar(prhs[0]); /* selected model */
    N  = (int) mxGetScalar(prhs[1]); /* choosen dimension (N x N x N) */
    P  = (int) mxGetScalar(prhs[2]); /* detector size */
    Th  = (float*) mxGetData(prhs[3]); /* angles */
    CenTypeIn = 1; /* astra-type centering is the default one */
    
   if (nrhs == 5)  {
        char *CenType;
        CenType = mxArrayToString(prhs[4]); /* 'radon' or 'astra' (default) */
        if ((strcmp(CenType, "radon") != 0) && (strcmp(CenType, "astra") != 0)) mexErrMsgTxt("Choose 'radon' or 'astra''");
        if (strcmp(CenType, "radon") == 0)  CenTypeIn = 0;  /* enable 'radon'-type centaering */
        mxFree(CenType);
    }
    NStructElems = mxGetNumberOfElements(prhs[3]);
    
    /*Handling Matlab output data*/
    int N_dims[] = {NStructElems, P, N};
    A = (float*)mxGetPr(plhs[0] = mxCreateNumericArray(3, N_dims, mxSINGLE_CLASS, mxREAL));    
        
    buildSino3D_core(A, ModelSelected, N, P, Th, (int)NStructElems, CenTypeIn);
}

