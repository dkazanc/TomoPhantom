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

#include "TomoP2DModelSino_core.h"

#define M_PI 3.14159265358979323846
#define EPS 0.000000001
#define MAXCHAR 1000

/* Function to create 2D analytical sinograms (parallel beam geometry) to from models using Phantom2DLibrary.dat
 *
 * Input Parameters:
 * 1. Model number (see Phantom2DLibrary.dat) [required]
 * 2. VolumeSize in voxels (N x N ) [required]
 * 3. Detector array size P (in pixels) [required]
 * 4. Projection angles Th (in degrees) [required]
 * 5. An absolute path to the file Phantom2DLibrary.dat (see OS-specific syntax-differences) [required]
 * 6. VolumeCentring, choose 'radon' or 'astra' (default) [optional]
 *
 * Output:
 * 1.2D sinogram size of [length(angles), P] or a temporal sinogram size of [length(angles), P, Time-Frames]
 */

/* function to build a single sinogram - object */
float TomoP2DObjectSino_core(float *A, int N, int P, float *Th, int AngTot, int CenTypeIn, char *Object, float C0, float x0, float y0, float a, float b, float phi_rot, int tt)
{
    int i, j;
    float *Tomorange_X_Ar=NULL, Tomorange_Xmin, Tomorange_Xmax, Sinorange_Pmax, Sinorange_Pmin, H_p, H_x, C1, a22, b22, phi_rot_radian;
    float *Sinorange_P_Ar=NULL, *AnglesRad=NULL;
    float AA5, sin_2, cos_2, delta1, delta_sq, first_dr, AA2, AA3, AA6, under_exp, x00, y00;
    
    Sinorange_Pmax = (float)(P)/(float)(N+1);
    Sinorange_Pmin = -Sinorange_Pmax;
    
    Sinorange_P_Ar = malloc(P*sizeof(float));
    H_p = (Sinorange_Pmax - Sinorange_Pmin)/(P-1);
    for(i=0; i<P; i++) {Sinorange_P_Ar[i] = (Sinorange_Pmax) - (float)i*H_p;}
    
    Tomorange_X_Ar = malloc(N*sizeof(float));
    Tomorange_Xmin = -1.0f;
    Tomorange_Xmax = 1.0f;
    H_x = (Tomorange_Xmax - Tomorange_Xmin)/(N);
    for(i=0; i<N; i++)  {Tomorange_X_Ar[i] = Tomorange_Xmin + (float)i*H_x;}
    AnglesRad = malloc(AngTot*sizeof(float));
    for(i=0; i<AngTot; i++)  AnglesRad[i] = (Th[i])*((float)M_PI/180.0f) + M_PI;
    
    C1 = -4.0f*logf(2.0f);
    
    if (CenTypeIn == 0) {
        /* matlab radon-iradon settings */
        x00 = x0 + H_x;
        y00 = y0 + H_x;
    }
    else {
        /* astra-toolbox settings */
        /*2D parallel beam*/
        x00 = x0 + 0.5f*H_x;
        y00 = y0 + 0.5f*H_x;
    }
    
    /************************************************/
    phi_rot_radian = (phi_rot)*((float)M_PI/180.0f);
    a22 = a*a;
    b22 = b*b;
    
    /* parameters of an object have been extracted, now run the building module */
    if (strcmp("gaussian",Object) == 0) {
        /* The object is a gaussian */
        AA5 = (N/2.0f)*(C0*(a)*(b)/2.0f)*sqrtf((float)M_PI/logf(2.0f));
#pragma omp parallel for shared(A) private(i,j,sin_2,cos_2,delta1,delta_sq,first_dr,under_exp,AA2,AA3)
        for(i=0; i<AngTot; i++) {
            sin_2 = powf((sinf((AnglesRad[i]) - phi_rot_radian)),2);
            cos_2 = powf((cosf((AnglesRad[i]) - phi_rot_radian)),2);
            delta1 = 1.0f/(a22*sin_2+b22*cos_2);
            delta_sq = sqrtf(delta1);
            first_dr = AA5*delta_sq;
            AA2 = -x00*sinf(AnglesRad[i])+y00*cosf(AnglesRad[i]); /*p0*/
            for(j=0; j<P; j++) {
                AA3 = powf((Sinorange_P_Ar[j] - AA2),2); /*(p-p0)^2*/
                under_exp = (C1*AA3)*delta1;
                A[tt*AngTot*P + j*AngTot+i] += first_dr*expf(under_exp);
            }}
    }
    else if (strcmp("parabola",Object) == 0) {
        /* the object is a parabola Lambda = 1/2 */
        AA5 = (N/2.0f)*(((float)M_PI/2.0f)*C0*((a))*((b)));
#pragma omp parallel for shared(A) private(i,j,sin_2,cos_2,delta1,delta_sq,first_dr,AA2,AA3,AA6)
        for(i=0; i<AngTot; i++) {
            sin_2 = powf((sinf((AnglesRad[i]) - phi_rot_radian)),2);
            cos_2 = powf((cosf((AnglesRad[i]) - phi_rot_radian)),2);
            delta1 = 1.0f/(a22*sin_2+b22*cos_2);
            delta_sq = sqrtf(delta1);
            first_dr = AA5*delta_sq;
            AA2 = -x00*sinf(AnglesRad[i])+y00*cosf(AnglesRad[i]); /*p0*/
            for(j=0; j<P; j++) {
                AA3 = powf((Sinorange_P_Ar[j] - AA2),2); /*(p-p0)^2*/
                AA6 = AA3*delta1;
                if (AA6 < 1.0f) {
                    A[tt*AngTot*P + j*AngTot+i] += first_dr*(1.0f - AA6);
                }
            }}
    }
    else if (strcmp("ellipse",Object) == 0) {
        /* the object is an elliptical disk */
        AA5 = (N*C0*a*b);
#pragma omp parallel for shared(A) private(i,j,sin_2,cos_2,delta1,delta_sq,first_dr,AA2,AA3,AA6)
        for(i=0; i<AngTot; i++) {
            sin_2 = powf((sinf((AnglesRad[i]) - phi_rot_radian)),2);
            cos_2 = powf((cosf((AnglesRad[i]) - phi_rot_radian)),2);
            delta1 = 1.0f/(a22*sin_2 + b22*cos_2);
            delta_sq = sqrtf(delta1);
            first_dr = AA5*delta_sq;
            AA2 = -x00*sinf(AnglesRad[i]) + y00*cosf(AnglesRad[i]); /*p0*/
            for(j=0; j<P; j++) {
                AA3 = powf((Sinorange_P_Ar[j] - AA2),2); /*(p-p0)^2*/
                AA6 = (AA3)*delta1;
                if (AA6 < 1.0f) {
                    A[tt*AngTot*P + j*AngTot+i] += first_dr*sqrtf(1.0f - AA6);
                }
            }}
    }
    else if (strcmp("parabola1",Object) == 0) {
        /* the object is a parabola Lambda = 1 (12)*/
        AA5 = (N/2.0f)*(4.0f*((0.25f*(a)*(b)*C0)/2.5f));
#pragma omp parallel for shared(A) private(i,j,sin_2,cos_2,delta1,delta_sq,first_dr,AA2,AA3,AA6)
        for(i=0; i<AngTot; i++) {
            sin_2 = powf((sinf((AnglesRad[i]) - phi_rot_radian)),2);
            cos_2 = powf((cosf((AnglesRad[i]) - phi_rot_radian)),2);
            delta1 = 1.0f/(0.25f*(a22)*sin_2 + 0.25f*b22*cos_2);
            delta_sq = sqrtf(delta1);
            first_dr = AA5*delta_sq;
            AA2 = -x00*sinf(AnglesRad[i]) + y00*cosf(AnglesRad[i]); /*p0*/
            for(j=0; j<P; j++) {
                AA3 = powf((Sinorange_P_Ar[j] - AA2),2); /*(p-p0)^2*/
                AA6 = AA3*delta1;
                if (AA6 < 1.0f) {
                    A[tt*AngTot*P + j*AngTot+i] += first_dr*(1.0f - AA6);
                }
            }}
    }
    else if (strcmp("cone",Object) == 0) {
        /* the object is a cone */
        float pps2,rlogi,ty1;
        AA5 = (N/2.0f)*(a*b*C0);
#pragma omp parallel for shared(A) private(i,j,sin_2,cos_2,delta1,delta_sq,first_dr,AA2,AA3,AA6,pps2,rlogi,ty1)
        for(i=0; i<AngTot; i++) {
            sin_2 = powf((sinf((AnglesRad[i]) - phi_rot_radian)),2);
            cos_2 = powf((cosf((AnglesRad[i]) - phi_rot_radian)),2);
            delta1 = 1.0f/(a22*sin_2 + b22*cos_2);
            delta_sq = sqrtf(delta1);
            first_dr = AA5*delta_sq;
            AA2 = -x00*sinf(AnglesRad[i]) + y00*cosf(AnglesRad[i]); /*p0*/
            for(j=0; j<P; j++) {
                AA3 = powf((Sinorange_P_Ar[j] - AA2),2); /*(p-p0)^2*/
                AA6 = AA3*delta1;
                pps2 = 0.0f; rlogi=0.0f;
                if (AA6 < 1.0f)
                    if (AA6 < (1.0f - EPS)) {
                        pps2 = sqrtf(fabs(1.0f - AA6));
                        rlogi=0.0f;
                    }
                if ((AA6 > EPS) && (pps2 != 1.0f)) {
                    ty1 = (1.0f + pps2)/(1.0f - pps2);
                    if (ty1 > 0.0f) rlogi = 0.5f*AA6*logf(ty1);
                }
                A[tt*AngTot*P+ j*AngTot+i] += first_dr*(pps2 - rlogi);
            }}
    }
    else if (strcmp("rectangle",Object) == 0) {
        /* the object is a rectangle */
        float xwid,ywid,p00,ksi00,ksi1;
        float PI2,p,ksi,C,S,A2,B2,FI,CF,SF,P0,TF,PC,QM,DEL,XSYC,QP,SS,x11,y11;
        
        if (CenTypeIn == 0) {
            /* matlab radon-iradon settings */
            x11 = -2.0f*y0 - H_x;
            y11 = 2.0f*x0 + H_x;
        }
        else {
            /* astra-toolbox settings */
            x11 = -2.0f*y0 - 0.5f*H_x;
            y11 = 2.0f*x0 + 0.5f*H_x;
        }
        
        xwid = b;
        ywid = a;
        
        if (phi_rot_radian < 0)  {ksi1 = (float)M_PI + phi_rot_radian;}
        else ksi1 = phi_rot_radian;
        
#pragma omp parallel for shared(A) private(i,j,PI2,p,ksi,C,S,A2,B2,FI,CF,SF,P0,TF,PC,QM,DEL,XSYC,QP,SS,p00,ksi00)
        for(i=0; i<AngTot; i++) {
            ksi00 = AnglesRad[(AngTot-1)-i] - M_PI;
            for(j=0; j<P; j++) {
                p00 = Sinorange_P_Ar[j];
                
                PI2 = (float)M_PI*0.5f;
                p = p00;
                ksi=ksi00;
                
                if (ksi > (float)M_PI) {
                    ksi = ksi - (float)M_PI;
                    p = -p00; }
                
                S = sinf(ksi); C = cosf(ksi);
                XSYC = -x11*C + y11*S;
                A2 = xwid*0.5f;
                B2 = ywid*0.5f;
                
                if ((ksi - ksi1) < 0.0f)  FI = (float)M_PI + ksi - ksi1;
                else FI = ksi - ksi1;
                
                if (FI > PI2) FI = (float)M_PI - FI;
                
                CF = cosf(FI);
                SF = sinf(FI);
                P0 = fabs(p-XSYC);
                
                SS = xwid/CF*C0;
                
                if (fabs(CF) <= (float)EPS) {
                    SS = ywid*C0;
                    if ((P0 - A2) > (float)EPS) {
                        SS=0.0f;
                    }
                }
                if (fabs(SF) <= (float)EPS) {
                    SS = xwid*C0;
                    if ((P0 - B2) > (float)EPS) {
                        SS=0.0f;
                    }
                }
                TF = SF/CF;
                PC = P0/CF;
                QP = B2+A2*TF;
                QM = QP+PC;
                if (QM > ywid) {
                    DEL = P0+B2*CF;
                    SS = ywid/SF*C0;
                    if (DEL > (A2*SF)) {
                        SS = (QP-PC)/SF*C0;
                    }}
                if (QM > ywid) {
                    DEL = P0+B2*CF;
                    if (DEL > A2*SF) SS = (QP-PC)/SF*C0;
                    else SS = ywid/SF*C0;
                }
                else SS = xwid/CF*C0;
                if (PC >= QP) SS=0.0f;
                A[tt*AngTot*P+ j*AngTot+i] += (N/2.0f)*SS;
            }}
    }
    else {        
        return 0;
    }
    /************************************************/
    free(Tomorange_X_Ar); free(Sinorange_P_Ar); free(AnglesRad);
    return *A;
}

float TomoP2DModelSino_core(float *A, int ModelSelected, int N, int P, float *Th, int AngTot, int CenTypeIn, char *ModelParametersFilename)
{    
    int Model=0, Components=0, steps = 0, counter=0, ii;
    float C0 = 0.0f, x0 = 0.0f, y0 = 0.0f, a = 0.0f, b = 0.0f, psi_gr1 = 0.0f;
    
    FILE *fp = fopen(ModelParametersFilename, "r"); // read parameters file    
    
    if( fp == NULL ) {        
          printf("%s \n","Cannot open the model library file (Phantom2DLibrary.dat)");
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
                            //mexErrMsgTxt("Unexpected the end of the line (Components) in parameters file");
                            break; }
                        if  (strcmp(tmpstr1,"Components") == 0) Components = atoi(tmpstr2);
                        printf("%s %i\n", "Components:", Components);
                        if (Components <= 0) {
                            printf("%s %i\n", "Components cannot be negative, the given value is", Components);
                            //mexErrMsgTxt("Components cannot be negative");
                            break; }
                        if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %21[^;];", tmpstr1, tmpstr2);
                        else {
                            //mexErrMsgTxt("Unexpected the end of the line (TimeSteps) in parameters file");
                            break; }
                        if  (strcmp(tmpstr1,"TimeSteps") == 0) steps = atoi(tmpstr2);
                        if (steps <= 0) {
                            printf("%s %i\n", "TimeSteps cannot be negative, the given value is", steps);
                            //mexErrMsgTxt("TimeSteps cannot be negative");
                            break; }
                        printf("%s %i\n", "TimeSteps:", steps);
                        if (steps == 1) {
                            /**************************************************/
                            printf("\n %s %i %s \n", "Stationary 2D sinogram for model", ModelSelected, " is selected");
                            
                            /* loop over all components */
                            for(ii=0; ii<Components; ii++) {
                                if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %21s %15s %15s %15s %15s %15s %15[^;];", tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5, tmpstr6, tmpstr7, tmpstr8);
                                else {
                                    //mexErrMsgTxt("Unexpected the end of the line (objects loop) in parameters file");
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
                                    //mexErrMsgTxt("Cannot find 'Object' string in parameters file");
                                    break; }
                                
                                if ((strcmp("gaussian",tmpstr2) != 0) && (strcmp("parabola",tmpstr2) != 0) && (strcmp("ellipse",tmpstr2) != 0) && (strcmp("parabola1",tmpstr2) != 0) && (strcmp("cone",tmpstr2) != 0) && (strcmp("rectangle",tmpstr2) != 0) ) {
                                    printf("%s %s\n", "Unknown name of the object, the given name is", tmpstr2);
                                    //mexErrMsgTxt("Unknown name of the object");
                                    break; }
                                if (C0 == 0) {
                                    printf("%s %f\n", "C0 should not be equal to zero, the given value is", C0);
                                    //mexErrMsgTxt("C0 should not be equal to zero");
                                    break; }
                                if ((x0 < -1) || (x0 > 1)) {
                                    printf("%s %f\n", "x0 (object position) must be in [-1,1] range, the given value is", x0);
                                    //mexErrMsgTxt("x0 (object position) must be in [-1,1] range");
                                    break; }
                                if ((y0 < -1) || (y0 > 1)) {
                                    printf("%s %f\n", "y0 (object position) must be in [-1,1] range, the given value is", y0);
                                    //mexErrMsgTxt("y0 (object position) must be in [-1,1] range");
                                    break; }
                                if ((a <= 0) || (a > 2)) {
                                    printf("%s %f\n", "a (object size) must be positive in [0,2] range, the given value is", a);
                                    //mexErrMsgTxt("a (object size) must be positive in [0,2] range");
                                    break; }
                                if ((b <= 0) || (b > 2)) {
                                    printf("%s %f\n", "b (object size) must be positive in [0,2] range, the given value is", b);
                                    //mexErrMsgTxt("b (object size) must be positive in [0,2] range");
                                    break; }
                                //printf("\nObject : %s \nC0 : %f \nx0 : %f \ny0 : %f \na : %f \nb : %f \n", tmpstr2, C0, x0, y0,  a, b);
                               TomoP2DObjectSino_core(A, N, P, Th, AngTot, CenTypeIn, tmpstr2, C0, x0, y0, b, a, -psi_gr1, 0);
                            }
                        }
                        else {
                            /**************************************************/                            
                            printf("\n %s %i %s \n", "Temporal 2D+time sinogram for model", ModelSelected, " is selected");
                            /* temporal phantom 2D + time (3D) */                            
                            
                            float C1 = 0.0f, x1 = 0.0f, y1 = 0.0f, a1 = 0.0f, b1 = 0.0f, psi_gr1_1 = 0.0f;
                            /* loop over all components */
                            for(ii=0; ii<Components; ii++) {
                                
                                if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %15s %15s %15s %15s %15s %15s %15[^;];", tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5, tmpstr6, tmpstr7, tmpstr8);
                                else {
                                    //mexErrMsgTxt("Unexpected the end of the line (objects loop) in parameters file");
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
                                    //mexErrMsgTxt("Cannot find 'Object' string in parameters file");
                                    break; }                                
                               
                                // printf("\nObject : %s \nC0 : %f \nx0 : %f \ny0 : %f \nz0 : %f \na : %f \nb : %f \n", tmpstr2, C0, x0, y0, z0, a, b, c);
                                
                                /* check Endvar relatedparameters */
                                if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %15s %15s %15s %15s %15s %15[^;];", tmpstr1, tmpstr3, tmpstr4, tmpstr5, tmpstr6, tmpstr7, tmpstr8);
                                else {
                                    //mexErrMsgTxt("Unexpected the end of the line (Endvar loop) in parameters file");
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
                                    
                                    TomoP2DObjectSino_core(A, N, P, Th, AngTot, CenTypeIn, tmpstr2, C_t, -y_t, x_t, b_t, a_t, -phi_t, tt);                                                         
                                    
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
        //mexErrMsgTxt("No object found, check models file"); 
    }
    return *A;
}    
