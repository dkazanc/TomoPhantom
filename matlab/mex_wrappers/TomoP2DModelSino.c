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
#define MAXCHAR 1000

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
    char *filename;
    FILE * fp;
    
  /*Handling Matlab input data*/
    if ((nrhs < 5) && (nrhs > 6)) mexErrMsgTxt("Input of 5 or 6 parameters is required: Model, ImageSize, Detector Size, Projection angles, PATH, Centering");
    if (mxGetClassID(prhs[3]) != mxSINGLE_CLASS) {mexErrMsgTxt("The vector of angles must be in a single precision"); }
    
    ModelSelected  = (int) mxGetScalar(prhs[0]); /* selected model */
    N  = (int) mxGetScalar(prhs[1]); /* choosen dimension (N x N) */
    P  = (int) mxGetScalar(prhs[2]); /* detector size */
    Th  = (float*) mxGetData(prhs[3]); /* angles */
    CenTypeIn = 1; /* astra-type centering is the default one */
    
    if (nrhs == 6)  {
        char *CenType;
        CenType = mxArrayToString(prhs[5]); /* 'radon' or 'astra' (default) */
        if ((strcmp(CenType, "radon") != 0) && (strcmp(CenType, "astra") != 0)) mexErrMsgTxt("Choose 'radon' or 'astra''");
        if (strcmp(CenType, "radon") == 0)  CenTypeIn = 0;  /* enable 'radon'-type centaering */
        mxFree(CenType);
    }
    NStructElems = mxGetNumberOfElements(prhs[3]);
    
    int Model=0, Components=0, steps = 0, counter=0, ii;
    float C0 = 0.0f, x0 = 0.0f, y0 = 0.0f, a = 0.0f, b = 0.0f, psi_gr1;
    
    if(!mxIsChar(prhs[4]) ) {
        mexErrMsgTxt("Need filename absolute path string input.");
    }
    filename = mxArrayToString(prhs[4]);
    fp= fopen(filename, "r");
    mxFree(filename);
    if( fp == NULL ) {
        mexErrMsgTxt("Cannot open the model library file (Phantom3DLibrary.dat)");
    }
    else {
        
        char str[MAXCHAR];
        char tmpstr1[16];
        char tmpstr2[22];
        char tmpstr3[16];
        char tmpstr4[16];
        char tmpstr5[16];
        char tmpstr6[16];
        char tmpstr7[16];
        char tmpstr8[16];
        
        while (fgets(str, MAXCHAR, fp) != NULL)
        {
            /* work with non-# commented lines */
            if(str[0] != '#') {
                sscanf(str, "%15s : %21[^;];", tmpstr1, tmpstr2);
                if (strcmp(tmpstr1,"Model")==0)
                {
                    Model = atoi(tmpstr2);
                    if ((ModelSelected == Model) && (counter == 0)) {
                        /* check if we have a right model */
                        if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %21[^;];", tmpstr1, tmpstr2);
                        else {
                            mexErrMsgTxt("Unexpected the end of the line (Components) in parameters file");
                            break; }
                        if  (strcmp(tmpstr1,"Components") == 0) Components = atoi(tmpstr2);
                        printf("%s %i\n", "Components:", Components);
                        if (Components <= 0) {
                            printf("%s %i\n", "Components cannot be negative, the given value is", Components);
                            mexErrMsgTxt("Components cannot be negative");
                            break; }
                        if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %21[^;];", tmpstr1, tmpstr2);
                        else {
                            mexErrMsgTxt("Unexpected the end of the line (TimeSteps) in parameters file");
                            break; }
                        if  (strcmp(tmpstr1,"TimeSteps") == 0) steps = atoi(tmpstr2);
                        if (steps <= 0) {
                            printf("%s %i\n", "TimeSteps cannot be negative, the given value is", steps);
                            mexErrMsgTxt("TimeSteps cannot be negative");
                            break; }
                        printf("%s %i\n", "TimeSteps:", steps);
                        if (steps == 1) {
                            /**************************************************/
                            printf("\n %s %i %s \n", "Stationary 2D model", ModelSelected, " is selected");

                            const mwSize N_dims[2] = {NStructElems,P}; /*format: X-detectors, Y-angles dim*/
                             A = (float*)mxGetPr(plhs[0] = mxCreateNumericArray(2, N_dims, mxSINGLE_CLASS, mxREAL)); 
                            
                            /* loop over all components */
                            for(ii=0; ii<Components; ii++) {
                                if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %21s %15s %15s %15s %15s %15s %15[^;];", tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5, tmpstr6, tmpstr7, tmpstr8);
                                else {
                                    mexErrMsgTxt("Unexpected the end of the line (objects loop) in parameters file");
                                    break; }
                                
                                if  (strcmp(tmpstr1,"Object") == 0) {
                                    C0 = (float)atof(tmpstr3); /* intensity */
                                    x0 = (float)atof(tmpstr4); /* x0 position */
                                    y0 = (float)atof(tmpstr5); /* y0 position */
                                    a = (float)atof(tmpstr6); /* a - size object */
                                    b = (float)atof(tmpstr7); /* b - size object */
                                    psi_gr1 = (float)atof(tmpstr8); /* rotation angle 1*/
                                }
                                else {
                                    mexErrMsgTxt("Cannot find 'Object' string in parameters file");
                                    break; }
                                
                                if ((strcmp("gaussian",tmpstr2) != 0) && (strcmp("parabola",tmpstr2) != 0) && (strcmp("ellipse",tmpstr2) != 0) && (strcmp("parabola1",tmpstr2) != 0) && (strcmp("cone",tmpstr2) != 0) && (strcmp("rectangle",tmpstr2) != 0) ) {
                                    printf("%s %s\n", "Unknown name of the object, the given name is", tmpstr2);
                                    mexErrMsgTxt("Unknown name of the object");
                                    break; }
                                if (C0 == 0) {
                                    printf("%s %f\n", "C0 should not be equal to zero, the given value is", C0);
                                    mexErrMsgTxt("C0 should not be equal to zero");
                                    break; }
                                if ((x0 < -1) || (x0 > 1)) {
                                    printf("%s %f\n", "x0 (object position) must be in [-1,1] range, the given value is", x0);
                                    mexErrMsgTxt("x0 (object position) must be in [-1,1] range");
                                    break; }
                                if ((y0 < -1) || (y0 > 1)) {
                                    printf("%s %f\n", "y0 (object position) must be in [-1,1] range, the given value is", y0);
                                    mexErrMsgTxt("y0 (object position) must be in [-1,1] range");
                                    break; }
                                if ((a <= 0) || (a > 2)) {
                                    printf("%s %f\n", "a (object size) must be positive in [0,2] range, the given value is", a);
                                    mexErrMsgTxt("a (object size) must be positive in [0,2] range");
                                    break; }
                                if ((b <= 0) || (b > 2)) {
                                    printf("%s %f\n", "b (object size) must be positive in [0,2] range, the given value is", b);
                                    mexErrMsgTxt("b (object size) must be positive in [0,2] range");
                                    break; }
                                printf("\nObject : %s \nC0 : %f \nx0 : %f \ny0 : %f \na : %f \nb : %f \n", tmpstr2, C0, x0, y0, a, b);

                                TomoP2DObjectSino_core(A, N, P, Th, (int)NStructElems, CenTypeIn, tmpstr2, C0, x0, y0, b, a, -psi_gr1, 0); /* Matlab */
                            }
                        }
                        else {
                            /**************************************************/
                            printf("\n %s %i %s \n", "Temporal 2D+time model", ModelSelected, " is selected");
                            /* temporal phantom 2D + time (3D) */ 
                            const mwSize N_dims[3] = {NStructElems, P, steps}; /*format: X-detectors, Y-angles dim, Time-Frames*/
                            A = (float*)mxGetPr(plhs[0] = mxCreateNumericArray(3, N_dims, mxSINGLE_CLASS, mxREAL)); 
                            
                            float C1 = 0.0f, x1 = 0.0f, y1 = 0.0f, a1 = 0.0f, b1 = 0.0f, psi_gr1_1 = 0.0f;
                            /* loop over all components */
                            for(ii=0; ii<Components; ii++) {
                                
                                if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %15s %15s %15s %15s %15s %15s %15[^;];", tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5, tmpstr6, tmpstr7, tmpstr8);
                                else {
                                    mexErrMsgTxt("Unexpected the end of the line (objects loop) in parameters file");
                                    break; }
                                
                                if  (strcmp(tmpstr1,"Object") == 0) {
                                    C0 = (float)atof(tmpstr3); /* intensity */
                                    x0 = (float)atof(tmpstr4); /* x0 position */
                                    y0 = (float)atof(tmpstr5); /* y0 position */
                                    a = (float)atof(tmpstr6); /* a - size object */
                                    b = (float)atof(tmpstr7); /* b - size object */
                                    psi_gr1 = (float)atof(tmpstr8); /* rotation angle 1*/
                                }
                                else {
                                    mexErrMsgTxt("Cannot find 'Object' string in parameters file");
                                    break; }
                                
                                if ((strcmp("gaussian",tmpstr2) != 0) && (strcmp("parabola",tmpstr2) != 0) && (strcmp("ellipse",tmpstr2) != 0) && (strcmp("parabola1",tmpstr2) != 0) && (strcmp("cone",tmpstr2) != 0) && (strcmp("rectangle",tmpstr2) != 0) ) {
                                    printf("%s %s\n", "Unknown name of the object, the given name is", tmpstr2);
                                    mexErrMsgTxt("Unknown name of the object");
                                    break; }
                                if (C0 == 0) {
                                    printf("%s %f\n", "C0 should not be equal to zero, the given value is", C0);
                                    mexErrMsgTxt("C0 should not be equal to zero");
                                    break; }
                                if ((x0 < -1) || (x0 > 1)) {
                                    printf("%s %f\n", "x0 (object position) must be in [-1,1] range, the given value is", x0);
                                    mexErrMsgTxt("x0 (object position) must be in [-1,1] range");
                                    break; }
                                if ((y0 < -1) || (y0 > 1)) {
                                    printf("%s %f\n", "y0 (object position) must be in [-1,1] range, the given value is", y0);
                                    mexErrMsgTxt("y0 (object position) must be in [-1,1] range");
                                    break; }
                                if ((a <= 0) || (a > 2)) {
                                    printf("%s %f\n", "a (object size) must be positive in [0,2] range, the given value is", a);
                                    mexErrMsgTxt("a (object size) must be positive in [0,2] range");
                                    break; }
                                if ((b <= 0) || (b > 2)) {
                                    printf("%s %f\n", "b (object size) must be positive in [0,2] range, the given value is", b);
                                    mexErrMsgTxt("b (object size) must be positive in [0,2] range");
                                    break; }
                                // printf("\nObject : %s \nC0 : %f \nx0 : %f \ny0 : %f \nz0 : %f \na : %f \nb : %f \n", tmpstr2, C0, x0, y0, z0, a, b, c);
                                
                                /* check Endvar relatedparameters */
                                if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %15s %15s %15s %15s %15s %15[^;];", tmpstr1, tmpstr3, tmpstr4, tmpstr5, tmpstr6, tmpstr7, tmpstr8);
                                else {
                                    mexErrMsgTxt("Unexpected the end of the line (Endvar loop) in parameters file");
                                    break; }
                                
                                if  (strcmp(tmpstr1,"Endvar") == 0) {
                                    C1 = (float)atof(tmpstr3); /* intensity */
                                    x1 = (float)atof(tmpstr4); /* x0 position */
                                    y1 = (float)atof(tmpstr5); /* y0 position */
                                    a1 = (float)atof(tmpstr6); /* a - size object */
                                    b1 = (float)atof(tmpstr7); /* b - size object */
                                    psi_gr1_1 = (float)atof(tmpstr8); /* rotation angle 1*/
                                }
                                else {
                                    printf("%s\n", "Cannot find 'Endvar' string in parameters file");
                                    break; }
                                
                                if (C1 == 0) {
                                    printf("%s %f\n", "Endvar C1 should not be equal to zero, the given value is", C1);
                                    mexErrMsgTxt("Endvar C0 should not be equal to zero");
                                    break; }
                                if ((x1 < -1) || (x1 > 1)) {
                                    printf("%s %f\n", "Endvar x1 (object position) must be in [-1,1] range, the given value is", x1);
                                    mexErrMsgTxt("Endvar x0 (object position) must be in [-1,1] range");
                                    break; }
                                if ((y1 < -1) || (y1 > 1)) {
                                    printf("%s %f\n", "Endvar y1 (object position) must be in [-1,1] range, the given value is", y1);
                                    mexErrMsgTxt("Endvar y0 (object position) must be in [-1,1] range");
                                    break; }
                                if ((a1 <= 0) || (a1 > 2)) {
                                    printf("%s %f\n", "Endvar a1 (object size) must be positive in [0,2] range, the given value is", a1);
                                    mexErrMsgTxt("Endvar a (object size) must be positive in [0,2] range");
                                    break; }
                                if ((b1 <= 0) || (b1 > 2)) {
                                    printf("%s %f\n", "Endvar b1 (object size) must be positive in [0,2] range, the given value is", b1);
                                    mexErrMsgTxt("Endvar b (object size) must be positive in [0,2] range");
                                    break; }
                                //printf("\nObject : %s \nC0 : %f \nx0 : %f \ny0 : %f \nz0 : %f \na : %f \nb : %f \nc : %f \n", tmpstr2, C0, x0, y0, z0, a1, b1, c1);
                                
                                /*now we know the initial parameters of the object and the final ones. We linearly extrapolate to establish steps and coordinates. */
                                /* calculating the full distance berween the start and the end points */                                
                                float distance = sqrtf(pow((x1 - x0),2) + pow((y1 - y0),2));
                                float d_dist = distance/(steps-1); /*a step over line */
                                float C_step = (C1 - C0)/(steps-1);
                                float a_step = (a1 - a)/(steps-1);
                                float b_step = (b1 - b)/(steps-1);
                                float phi_rot_step = (psi_gr1_1 - psi_gr1)/(steps-1);
                                
                                int tt;
                                float x_t, y_t, a_t, b_t, C_t, phi_t, d_step;
                                /* initialize */
                                x_t = x0; y_t = y0; a_t = a; b_t = b; C_t = C0; phi_t = psi_gr1; d_step = d_dist;
                                /*loop over time frames*/
                                for(tt=0; tt < steps; tt++) {
                                    
                                    TomoP2DObjectSino_core(A, N, P, Th, (int)NStructElems, CenTypeIn, tmpstr2, C_t, x_t, y_t, b_t, a_t, -phi_t, tt);
                                    
                                    /* calculating new coordinates of an object */
                                    if (distance != 0.0f) {
                                        float t = d_step/distance;
                                        x_t = (1-t)*x0 + t*x1;
                                        y_t = (1-t)*y0 + t*y1; }
                                    else {
                                        x_t = x0;
                                        y_t = y0;   }
                                   
                                    d_step += d_dist;
                                    a_t += a_step;
                                    b_t += b_step;
                                    C_t += C_step;
                                    phi_t += phi_rot_step;
                                } /*time steps*/
                                
                            } /*components loop*/
                        }
                        counter++;
                    }
                }
            }
        }
    }
    fclose(fp);
    if (counter == 0) {
        printf("%s %i %s \n", "Model no. ", ModelSelected, "is not found!");
        mexErrMsgTxt("No object found, check models file"); 
    }    
}    
    
    
    
//     /* first to check if the model is stationary or temporal */
//     int *timeFrames = (int *) mxGetData(mxCreateNumericMatrix(1, 1, mxINT16_CLASS, mxREAL));
//     /* extract the number of time-step */
//     extractTimeFrames(timeFrames, ModelSelected, ModelParameters_PATH);    
//     
//     int *params_switch = (int *) mxGetData(mxCreateNumericMatrix(1, 10, mxINT16_CLASS, mxREAL));
//     /* run checks for parameters (only integer switches and the number of time-frames) */
//     checkParams2D(params_switch, ModelSelected, ModelParameters_PATH, timeFrames[0]);
//     
//     if (params_switch[0] == 0) mexErrMsgTxt("Parameters file (Phantom2DLibrary.dat) cannot be read, check that your PATH to the file is a valid one and the file exists");
//     if (params_switch[1] == 0) {
//         printf("%s %i\n", ">>>>> No such model available in the library, the given model is", ModelSelected);
//         mexErrMsgTxt("Check the syntax in Phantom2DLibrary.dat");
//     }
//     if (params_switch[2] == 0) {
//         printf("%s %i\n", ">>>>> Components number cannot be negative for the given model", ModelSelected);
//         mexErrMsgTxt("Check the syntax in Phantom2DLibrary.dat");
//     }
//     if (params_switch[3] == 0) {
//         printf("%s %i\n", ">>>>> TimeSteps value cannot be negative for the given model", ModelSelected);
//         mexErrMsgTxt("Check the syntax in Phantom2DLibrary.dat");
//     }
//     if (params_switch[4] == 0) {
//         printf("%s %i\n", ">>>>> Unknown name of the object for the given model", ModelSelected);
//         mexErrMsgTxt("Check the syntax in Phantom2DLibrary.dat");
//     }
//     if (params_switch[5] == 0) {
//         printf("%s %i\n", ">>>>> C0 should not be equal to zero for the given model", ModelSelected);
//         mexErrMsgTxt("Check the syntax in Phantom2DLibrary.dat");
//     }
//     if (params_switch[6] == 0) {
//         printf("%s %i\n", ">>>>> x0 (object position) must be in [-1,1] range for the given model", ModelSelected);
//         mexErrMsgTxt("Check the syntax in Phantom2DLibrary.dat");
//     }
//     if (params_switch[7] == 0) {
//         printf("%s %i\n", ">>>>> y0 (object position) must be in [-1,1] range for the given model", ModelSelected);
//         mexErrMsgTxt("Check the syntax in Phantom2DLibrary.dat");
//     }
//     if (params_switch[8] == 0) {
//         printf("%s %i\n", ">>>>> a (object size) must be positive in [0,2] range", ModelSelected);
//         mexErrMsgTxt("Check the syntax in Phantom2DLibrary.dat");
//     }
//     if (params_switch[9] == 0) {
//         printf("%s %i\n", ">>>>> b (object size) must be positive in [0,2] range", ModelSelected);
//         mexErrMsgTxt("Check the syntax in Phantom2DLibrary.dat");
//     }
//     
//     /*Handling Matlab output data*/
//     if (params_switch[3] == 1) {
//         /* static phantom case */
//         int N_dims[] = {NStructElems,P}; /*format: X-detectors, Y-angles dim*/
//         A = (float*)mxGetPr(plhs[0] = mxCreateNumericArray(2, N_dims, mxSINGLE_CLASS, mxREAL)); }
//     else {
//         int N_dims[] = {NStructElems, P, params_switch[3]}; /*format: X-detectors, Y-angles dim, Time-Frames*/
//         A = (float*)mxGetPr(plhs[0] = mxCreateNumericArray(3, N_dims, mxSINGLE_CLASS, mxREAL)); }
//     TomoP2DModelSino_core(A, ModelSelected, N, P, Th, (int)NStructElems, CenTypeIn, ModelParameters_PATH, 0);
//     mxFree(ModelParameters_PATH);
// }
