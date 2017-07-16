#include "mex.h"
#include <matrix.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
#include "omp.h"

#define M_PI 3.14159265358979323846

/* C-OMP function to perform forward/inverse deformation
 *
 * Input Parameters:
 * 1. A - image to deform
 * 2. RFP - propotional to focal point distance
 * 3. angle - deformation angle in degrees
 * 4. DeformType - deformation type, 0 - forward, 1 - inverse
 *
 * Output:
 * 1. Deformed image
 *
 * to compile with OMP support: mex DeformObject_C.c CFLAGS="\$CFLAGS -fopenmp -Wall -std=c99" LDFLAGS="\$LDFLAGS -fopenmp"
 * References:
 * [1] D. Kazantsev & V. Pickalov, "New iterative reconstruction methods for fan-beam tomography" IPSE, 2017
 * D. Kazantsev, 2016-17
 */

double Deform_func(double *A, double *B, double *Tomorange_X_Ar, double H_x, double RFP, double angleRad, int DeformType, int dimX, int dimY);
double pad_crop(double *A, double *Ap, int OldSizeX, int OldSizeY, int NewSizeX, int NewSizeY, int padXY, int switchpad_crop);

void mexFunction(
        int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
        
{
    int number_of_dims, dimX, dimY, newsizeX, newsizeY, DeformType, i, padXY;
    const int  *dim_array;
    double *A, *Ap, *B, *Tomorange_X_Ar=NULL, RFP, angle, angleRad, Tomorange_Xmin, Tomorange_Xmax, H_x;
    int N_dims2D[2];
    
    number_of_dims = mxGetNumberOfDimensions(prhs[0]);
    dim_array = mxGetDimensions(prhs[0]);
    
    /*Handling Matlab input data*/
    if (nrhs != 4) mexErrMsgTxt("Input of 4 parameters is required");
    
    A  = (double *) mxGetData(prhs[0]); /* image to deform */
    RFP =  (double) mxGetScalar(prhs[1]); /*  propotional to focal point distance*/
    angle =  (double) mxGetScalar(prhs[2]); /*  deformation angle in degrees */
    DeformType =  (int) mxGetScalar(prhs[3]); /* deformation type, 0 - forward, 1 - inverse  */
    
    /*Handling Matlab output data*/
    dimX = dim_array[0]; dimY = dim_array[1];
    
    padXY = 150;
    if (DeformType == 0) {    
    newsizeX = dimX + 2*(padXY); /* the X size of the padded array */
    newsizeY = dimY + 2*(padXY); /* the Y size of the padded array */      
    N_dims2D[0] = newsizeX; N_dims2D[1] = newsizeY ;  }
    else {
    newsizeX = dimX;
    newsizeY = dimY;
    dimX = newsizeX - 2*(padXY);
    dimY = newsizeY - 2*(padXY);            
    N_dims2D[0] = dimX; N_dims2D[1] = dimY;  }     
    
    if (number_of_dims == 2) {
        
        B = (double*)mxGetPr(plhs[0] = mxCreateNumericArray(2, N_dims2D, mxDOUBLE_CLASS, mxREAL));
        Ap = (double*)mxGetPr(mxCreateNumericArray(2, N_dims2D, mxDOUBLE_CLASS, mxREAL));
        
        Tomorange_X_Ar = malloc(newsizeX*sizeof(double));
        
        angleRad = angle*(M_PI/180.0f);
        Tomorange_Xmin = -1.0f;
        Tomorange_Xmax = 1.0f;
        H_x = (Tomorange_Xmax - Tomorange_Xmin)/(newsizeX);
        for(i=0; i<newsizeX; i++)  {Tomorange_X_Ar[i] = Tomorange_Xmin + (double)i*H_x;}
        
        /* do image padding */        
        pad_crop(A, Ap, dimX, dimY, newsizeX, newsizeY, padXY, DeformType);
        
        /* perform deformation */
        Deform_func(Ap, B, Tomorange_X_Ar, H_x, RFP, angleRad, DeformType, newsizeX, newsizeY);
        free(Tomorange_X_Ar);
    }
}

/* Related Functions */
/*****************************************************************/
double Deform_func(double *A, double *B, double *Tomorange_X_Ar, double H_x, double RFP, double angleRad, int DeformType, int dimX, int dimY)
{
    double i0,j0,xx, yy, xPersp, yPersp, xPersp1, yPersp1, u, v, ll, mm, Ar1step_inv, a, b, c, d;
    int i,j,i1,j1,i2,j2;
    Ar1step_inv = 1.0f/H_x;
#pragma omp parallel for shared(A,B,Tomorange_X_Ar,Ar1step_inv) private(i,j,i1,j1,i2,j2,i0,j0,xx, yy, xPersp, yPersp, xPersp1, yPersp1, u, v, ll, mm, a, b, c, d )
    for(i=0; i<dimX; i++) {
        for(j=0; j<dimY; j++) {
            xx = Tomorange_X_Ar[i]*cos(angleRad) + Tomorange_X_Ar[j]*sin(angleRad);
            yy = -Tomorange_X_Ar[i]*sin(angleRad) + Tomorange_X_Ar[j]*cos(angleRad);
            
            if (DeformType == 0) {
                /* do forward transform*/
                xPersp1 = xx*(1.0f - yy*RFP);
                yPersp1 = yy;
                xPersp = xPersp1*(1.0f - yPersp1*RFP)*cos(angleRad) - yPersp1*sin(angleRad);
                yPersp = xPersp1*(1.0f - yPersp1*RFP)*sin(angleRad) + yPersp1*cos(angleRad);
            }
            else {
                /* do inverse transform */
                xPersp1 = xx/(1.0f - yy*RFP);
                yPersp1 = yy;
                xPersp = xPersp1/(1.0f - yPersp1*RFP)*cos(angleRad) - yPersp1*sin(angleRad);
                yPersp = xPersp1/(1.0f - yPersp1*RFP)*sin(angleRad) + yPersp1*cos(angleRad);
            }
//             printf("%f %f \n", xPersp, yPersp);
            /*Bilinear 2D Interpolation */
            ll = (xPersp - (-1.0f))*Ar1step_inv;
            mm = (yPersp - (-1.0f))*Ar1step_inv;
            
            i0 = (double)floor((double)ll);
            j0 = (double)floor((double)mm);
            u = ll - i0;
            v = mm - j0;
            
            i2 = (int)i0;
            j2 = (int)j0;
            
            i1 = i2+1;
            j1 = j2+1;
            
            a = 0.0f; b = 0.0f; c = 0.0f; d = 0.0f;
            
            if ((i2 >= 0) && (i2 < dimX) && (j2 >= 0) && (j2 < dimX))   a = A[(i2)*dimY + (j2)];
            if ((i2 >= 0) && (i2 < dimX-1) && (j2 >= 0) && (j2 < dimX)) b = A[(i1)*dimY + (j2)];
            if ((i2 >= 0) && (i2 < dimX) && (j2 >= 0) && (j2 < dimX-1)) c = A[(i2)*dimY + (j1)];
            if ((i2 >= 0) && (i2 < dimX-1) && (j2 >= 0) && (j2 < dimX-1))  d = A[(i1)*dimY + (j1)];
            
            B[(i)*dimY + (j)] = (1.0f - u)*(1.0f - v)*a + u*(1.0f - v)*b + (1.0f - u)*v*c+ u*v*d;
            
        }}
    return *B;
}

double pad_crop(double *A, double *Ap, int OldSizeX, int OldSizeY, int NewSizeX, int NewSizeY, int padXY, int switchpad_crop)
{
    /* padding-cropping function */
    int i, j;
    for (i=0; i < NewSizeX; i++) {
        for (j=0; j < NewSizeY; j++) {
                if (((i >= padXY) && (i < NewSizeX-padXY)) &&  ((j >= padXY) && (j < NewSizeY-padXY))) {
                    if (switchpad_crop == 0)  Ap[i*NewSizeY+j] = A[(i-padXY)*OldSizeY+(j-padXY)];
                    else  Ap[(i-padXY)*OldSizeY+(j-padXY)] = A[i*NewSizeY+j];
                }
            }}
return *Ap;
}