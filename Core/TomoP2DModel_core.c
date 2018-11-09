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

#include "TomoP2DModel_core.h"

#define M_PI 3.14159265358979323846
#define MAXCHAR 1000

/* Functions to build spatial (2D) and temporal (2D +time) phantoms from the library of models: Phantom2DLibrary.dat
 *
 * Input Parameters:
 * 1. ModelNo - the model number from Phantom3DLibrary file
 * 2. VolumeSize in voxels (N x N)
 * 3. Object - Analytical Model
 * 4. C0 - intensity
 * 5. x0 - x0 position
 * 6. y0 - y0 position
 * 7. a  - size object
 * 8. b  - size object
 * 9. phi_rot - rotation angle
 *
 * Output:
 * 1. The analytical phantom size of [N x N] or  a temporal phantom size of [N xN x Time-frames]
 */

/* function to build a single object */
float TomoP2DObject_core(float *A, int N, char *Object,
        float C0, /* intensity */
        float x0, /* x0 position */
        float y0, /* y0 position */
        float a , /* a - size object */
        float b , /* b - size object */
        float phi_rot, /* phi - rotation angle */
        int tt /* time frame loop */)
{
    // printf ("Base C0 %.2e x0 %.2e y0 %.2e a %.2e b %.2e phi %.2e\n" , C0, x0, y0, a, b, phi_rot);                            
                                
    int i, j;
    float *Tomorange_X_Ar=NULL, Tomorange_Xmin, Tomorange_Xmax, H_x, C1, a2, b2, phi_rot_radian, sin_phi, cos_phi;
    float *Xdel = NULL, *Ydel = NULL, T;
    Tomorange_X_Ar = malloc(N*sizeof(float));
    Tomorange_Xmin = -1.0f;
    Tomorange_Xmax = 1.0f;
    H_x = (Tomorange_Xmax - Tomorange_Xmin)/(N);
    for(i=0; i<N; i++)  {Tomorange_X_Ar[i] = Tomorange_Xmin + (float)i*H_x;}
    C1 = -4.0f*logf(2.0f);
    
    /************************************************/
    phi_rot_radian = phi_rot*((float)M_PI/180.0f);
    sin_phi=sinf(phi_rot_radian); cos_phi=cosf(phi_rot_radian);
    
    Xdel = malloc(N*sizeof(float));
    Ydel = malloc(N*sizeof(float));
    for(i=0; i<N; i++)  {
        Xdel[i] = Tomorange_X_Ar[i] - x0;
        Ydel[i] = Tomorange_X_Ar[i] - y0;
    }
    
    a2 = 1.0f/(a*a);
    b2 = 1.0f/(b*b);
    
    /* all parameters of an object have been extracted, now run the building modules */
    if (strcmp("gaussian",Object) == 0) {
        /* The object is a gaussian */
#pragma omp parallel for shared(A) private(i,j,T)
        for(i=0; i<N; i++) {
            for(j=0; j<N; j++) {
                T = C1*(a2*powf((Xdel[i]*cos_phi + Ydel[j]*sin_phi),2) + b2*powf((-Xdel[i]*sin_phi + Ydel[j]*cos_phi),2));
                A[tt*N*N + j*N+i] += C0*expf(T);
            }}
    }
    else if (strcmp("parabola",Object) == 0) {
        /* the object is a parabola Lambda = 1/2 */
#pragma omp parallel for shared(A) private(i,j,T)
        for(i=0; i<N; i++) {
            for(j=0; j<N; j++) {
                T = a2*powf((Xdel[i]*cos_phi + Ydel[j]*sin_phi),2) + b2*powf((-Xdel[i]*sin_phi + Ydel[j]*cos_phi),2);
                if (T <= 1) T = C0*sqrtf(1.0f - T);
                else T = 0.0f;
                A[tt*N*N + j*N+i] += T;
            }}
    }
    else if (strcmp("ellipse",Object) == 0) {
        /* the object is an elliptical disk */
#pragma omp parallel for shared(A) private(i,j,T)
        for(i=0; i<N; i++) {
            for(j=0; j<N; j++) {
                T = a2*powf((Xdel[i]*cos_phi + Ydel[j]*sin_phi),2) + b2*powf((-Xdel[i]*sin_phi + Ydel[j]*cos_phi),2);
                if (T <= 1) T = C0;
                else T = 0.0f;
                A[tt*N*N + j*N+i] += T;
            }}
    }
    else if (strcmp("parabola1",Object) == 0) {
        /* the object is a parabola Lambda = 1*/
#pragma omp parallel for shared(A) private(i,j,T)
        for(i=0; i<N; i++) {
            for(j=0; j<N; j++) {
                T = (4.0f*a2)*powf((Xdel[i]*cos_phi + Ydel[j]*sin_phi),2) + (4.0f*b2)*powf((-Xdel[i]*sin_phi + Ydel[j]*cos_phi),2);
                if (T <= 1) T = C0*sqrtf(1.0f - T);
                else T = 0.0f;
                A[tt*N*N + j*N+i] += T;
            }}
    }
    else if (strcmp("cone",Object) == 0) {
        /*the object is a cone*/
#pragma omp parallel for shared(A) private(i,j,T)
        for(i=0; i<N; i++) {
            for(j=0; j<N; j++) {
                T = a2*powf((Xdel[i]*cos_phi + Ydel[j]*sin_phi),2) + b2*powf((-Xdel[i]*sin_phi + Ydel[j]*cos_phi),2);
                if (T <= 1) T = C0*(1.0f - sqrtf(T));
                else T = 0.0f;
                A[tt*N*N + j*N+i] += T;
            }}
    }
    else if (strcmp("rectangle",Object) == 0) {
        /* the object is a rectangle */
        float x0r, y0r, HX, HY;
        a2 = 0.5f*a;
        b2 = 0.5f*b;
        x0r=x0*cosf(0.0f) + y0*sinf(0.0f);
        y0r=-x0*sinf(0.0f) + y0*cosf(0.0f);
        if (phi_rot_radian < 0.0f) {
            phi_rot_radian = (float)M_PI + phi_rot_radian;
            sin_phi=sinf(phi_rot_radian);
            cos_phi=cosf(phi_rot_radian);
        }
#pragma omp parallel for shared(A) private(i,j,HX,HY,T)
        for(i=0; i<N; i++) {
            for(j=0; j<N; j++) {
                HX = fabsf((Xdel[i] - x0r)*sin_phi + (Ydel[j] - y0r)*cos_phi);
                T = 0.0f;
                if (HX <= a2) {
                    HY = fabsf((Ydel[j] - y0r)*sin_phi - (Xdel[i] - x0r)*cos_phi);
                    if (HY <= b2) {T = C0;}
                }
                A[tt*N*N + j*N+i] += T;
            }}
    }
    else {        
          return 0;
    }
    free(Xdel); free(Ydel);
    /************************************************/
    free(Tomorange_X_Ar);
    return *A;
}

float TomoP2DModel_core(float *A, int ModelSelected, int N, char *ModelParametersFilename)
{
    FILE *fp = fopen(ModelParametersFilename, "r"); // read parameters file    
    int Model=0, Components=0, steps = 0, counter=0, ii;
    float C0 = 0.0f, x0 = 0.0f, y0 = 0.0f, a = 0.0f, b = 0.0f, psi_gr1 = 0.0f;
    
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
                        //printf("%s %i\n", "Components:", Components);
                        if (Components <= 0) {
                            // printf("%s %i\n", "Components cannot be negative, the given value is", Components);
                           // mexErrMsgTxt("Components cannot be negative");
                            break; }
                        if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %21[^;];", tmpstr1, tmpstr2);
                        else {
                            //mexErrMsgTxt("Unexpected the end of the line (TimeSteps) in parameters file");
                            break; }
                        if  (strcmp(tmpstr1,"TimeSteps") == 0) steps = atoi(tmpstr2);
                        if (steps <= 0) {
                            // printf("%s %i\n", "TimeSteps cannot be negative, the given value is", steps);
                            //mexErrMsgTxt("TimeSteps cannot be negative");
                            break; }
                        //printf("%s %i\n", "TimeSteps:", steps);
                        if (steps == 1) {
                            /**************************************************/                            
                            printf("\n %s %i %s \n", "Stationary 2D model", ModelSelected, " is selected");
                            
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
                                // printf ("C0 %.2e x0 %.2e y0 %.2e a %.2e b %.2e phi %.2e\n" , C0, x0, y0, a, b, psi_gr1);                            
                                TomoP2DObject_core(A, N, tmpstr2, C0, y0, x0, a, b, psi_gr1, 0); /* python */
                            }
                        }
                        else {
                            /**************************************************/                            
                            printf("\n %s %i %s \n", "Temporal 2D+time model", ModelSelected, " is selected");
                            /* temporal phantom 2D + time (3D) */
                            
                            float C1 = 0.0f, x1 = 0.0f, y1 = 0.0f, a1 = 0.0f, b1 = 0.0f, psi_gr1_1 = 0.0f;
                            /* loop over all components */
                            for(ii=0; ii<Components; ii++) {
                                
                                if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %15s %15s %15s %15s %15s %15s %15[^;];", tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5, tmpstr6, tmpstr7, tmpstr8);
                                else {
                                   // mexErrMsgTxt("Unexpected the end of the line (objects loop) in parameters file");
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
                                   // mexErrMsgTxt("Cannot find 'Object' string in parameters file");
                                    break; }                               

                                /* check Endvar relatedparameters */
                                if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %15s %15s %15s %15s %15s %15[^;];", tmpstr1, tmpstr3, tmpstr4, tmpstr5, tmpstr6, tmpstr7, tmpstr8);
                                else {
                                   // mexErrMsgTxt("Unexpected the end of the line (Endvar loop) in parameters file");
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
                                    
                                     TomoP2DObject_core(A, N, tmpstr2, C_t, x_t, -y_t, a_t, b_t, phi_t, tt); /* python */
                                    
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
    return *A;
}
