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

#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
#include "omp.h"
#include "utils.h"

#define M_PI 3.14159265358979323846
#define MAXCHAR 1000

/* Function to read parameters from the file Phantom3DLibrary.dat to build 3D analytical models
 *
 * Input Parameters:
 * 1. ModelNo - the model number from Phantom3DLibrary file
 * 2. VolumeSize in voxels (N x N x N)
 * 3. Object - Analytical Model
 * 4. C0 - intensity
 * 5. x0 - x0 position
 * 6. y0 - y0 position
 * 7. z0 - z0 position
 * 8. a  - size object
 * 9. b  - size object
 * 10. c - size object
 * 11. psi_gr1 - rotation angle1
 * 12. psi_gr2 - rotation angle2
 * 12. psi_gr3 - rotation angle3
 *
 * Output:
 * 1. The analytical phantom size of [N x N x N] or a temporal 4D phantom (N x N x N x time-frames)
 */

/* function to build a single (stationary) object */
float TomoP3DObject_core(float *A, int N, char *Object,
        float C0, /* intensity */
        float x0, /* x0 position */
        float y0, /* y0 position */
        float z0, /* z0 position */
        float a , /* a - size object */
        float b , /* b - size object */
        float c , /* c - size object */
        float psi_gr1, /* rotation angle1 */
        float psi_gr2, /* rotation angle2 */
        float psi_gr3, /* rotation angle3 */
        int tt /*temporal index, 0 - for stationary */)
{
    int i, j, k;
    float Tomorange_Xmin, Tomorange_Xmax, H_x, C1, a2, b2, c2, phi_rot_radian, sin_phi, cos_phi, aa,bb,cc, psi1, psi2, psi3, T;
    float *Tomorange_X_Ar=NULL, *Xdel = NULL, *Ydel = NULL, *Zdel = NULL;
    Tomorange_X_Ar = malloc(N*sizeof(float));
    if(Tomorange_X_Ar == NULL ) printf("Allocation of 'Tomorange_X_Ar' failed");
    Tomorange_Xmin = -1.0f;
    Tomorange_Xmax = 1.0f;
    H_x = (Tomorange_Xmax - Tomorange_Xmin)/(N);
    for(i=0; i<N; i++)  {Tomorange_X_Ar[i] = Tomorange_Xmin + (float)i*H_x;}
    C1 = -4.0f*logf(2.0f);
    
    /* parameters of a model have been extracted, now run the building module */
    /************************************************/
    phi_rot_radian = psi_gr1*((float)M_PI/180.0f);
    sin_phi=sinf(phi_rot_radian); cos_phi=cosf(phi_rot_radian);
    
    Xdel = malloc(N*sizeof(float));
    if(Xdel == NULL ) printf("Allocation of 'Xdel' failed");
    Ydel = malloc(N*sizeof(float));
    if(Ydel == NULL ) printf("Allocation of 'Ydel' failed");
    Zdel = malloc(N*sizeof(float));
    if(Zdel == NULL ) printf("Allocation of 'Zdel' failed");
    for(i=0; i<N; i++)  {
        Xdel[i] = Tomorange_X_Ar[i] - x0;
        Ydel[i] = Tomorange_X_Ar[i] - y0;
        Zdel[i] = Tomorange_X_Ar[i] - z0;
    }
    
    psi1 = psi_gr1*((float)M_PI/180.0f);
    psi2 = psi_gr2*((float)M_PI/180.0f);
    psi3 = psi_gr3*((float)M_PI/180.0f);
    
    float bs[9];
    float xh[3];
    float xh3[3];
    
    a2 = 1.0f/(a*a);
    b2 = 1.0f/(b*b);
    c2 = 1.0f/(c*c);
    matrot3(bs,psi1,psi2,psi3); /* rotation of 3x3 matrix */
    
    xh3[0] = x0; xh3[1] = y0; xh3[2] = z0;
    matvet3(bs,xh3,xh);  /* matrix-vector multiplication */
    
    float xh1[3];
    float xh2[3];
    
    if ((strcmp("gaussian",Object) == 0) ||  (strcmp("paraboloid",Object) == 0) || (strcmp("ellipsoid",Object) == 0) || (strcmp("cone",Object) == 0)) {
 #pragma omp parallel for shared(A,bs) private(k,i,j,aa,bb,cc,T,xh2,xh1)
        for(k=0; k<N; k++) {
            for(i=0; i<N; i++) {
                for(j=0; j<N; j++) {
                    
                    if ((psi1 != 0.0f) || (psi2 != 0.0f) || (psi3 != 0.0f)) {
                        xh1[0]=Tomorange_X_Ar[i];
                        xh1[1]=Tomorange_X_Ar[j];
                        xh1[2]=Tomorange_X_Ar[k];
                        matvet3(bs,xh1,xh2);
                        aa = a2*powf((xh2[0]-xh[0]),2);
                        bb = b2*powf((xh2[1]-xh[1]),2);
                        cc = c2*powf((xh2[2]-xh[2]),2);
                    }
                    else {
                        aa = a2*powf(Xdel[i],2);
                        bb = b2*powf(Ydel[j],2);
                        cc = c2*powf(Zdel[k],2);
                    }
                    T = (aa + bb + cc);
                    if (strcmp("gaussian",Object) == 0) {
                        /* The object is a volumetric gaussian */
                        T = C0*expf(C1*T);
                    }
                    if (strcmp("paraboloid",Object) == 0) {
                        /* the object is a parabola Lambda = 1/2 */
                        if (T <= 1.0f) T = C0*sqrtf(1.0f - T);
                        else T = 0.0f;
                    }
                    if (strcmp("ellipsoid",Object) == 0) {
                        /* the object is en ellipsoid */
                        if (T <= 1.0f) T = C0;
                        else T = 0.0f;
                    }
                    if (strcmp("cone",Object) == 0) {
                        /* the object is a cone */
                        if (T <= 1.0f) T = C0*(1.0f - sqrtf(T));
                        else T = 0.0f;
                    }
                    A[tt*N*N*N + (k)*N*N + j*N+i] += T;
                }}}
    }
    if (strcmp("cuboid",Object) == 0) {
        /* the object is a cuboid */
        float x0r, y0r, HX, HY;
        a2 = 0.5f*a;
        b2 = 0.5f*b;
        c2 = 0.5f*c;
        x0r=x0*cosf(0.0f) + y0*sinf(0.0f);
        y0r=-x0*sinf(0.0f) + y0*cosf(0.0f);
        if (phi_rot_radian < 0.0f) {
            phi_rot_radian = (float)M_PI + phi_rot_radian;
            sin_phi=sinf(phi_rot_radian);
            cos_phi=cosf(phi_rot_radian);
        }
#pragma omp parallel for shared(A,Zdel) private(k,i,j,HX,HY,T)
        for(k=0; k<N; k++) {
            if  (fabs(Zdel[k]) < c2) {
                
                for(i=0; i<N; i++) {
                    for(j=0; j<N; j++) {
                        HX = fabsf((Xdel[i] - x0r)*cos_phi + (Ydel[j] - y0r)*sin_phi);
                        T = 0.0f;
                        if (HX <= a2) {
                            HY = fabsf((Ydel[j] - y0r)*cos_phi - (Xdel[i] - x0r)*sin_phi);
                            if (HY <= b2) {T = C0;}
                        }
                        A[tt*N*N*N + (k)*N*N + j*N+i] += T;
                    }
                }
            }
        }
    }
    if (strcmp("elliptical_cylinder",Object) == 0) {
        /* the object is an elliptical cylinder  */
#pragma omp parallel for shared(A) private(k,i,j,T)
        for(k=0; k<N; k++) {
            if  (fabs(Zdel[k]) < c) {
                for(i=0; i<N; i++) {
                    for(j=0; j<N; j++) {
                        T = a2*powf((Xdel[i]*cos_phi + Ydel[j]*sin_phi),2) + b2*powf((-Xdel[i]*sin_phi + Ydel[j]*cos_phi),2);
                        if (T <= 1) T = C0;
                        else T = 0.0f;
                        A[tt*N*N*N + (k)*N*N + j*N+i] += T;
                    }}
            }
        } /*k-loop*/
    }
    /****************************************************/
    free(Xdel); free(Ydel); free(Zdel);
    free(Tomorange_X_Ar);
    return *A;
}


/********************Core Function*****************************/
float TomoP3DModel_core(float *A, int ModelSelected, int N, char *ModelParametersFilename)
{
   
    int Model=0, Components=0, steps = 0, counter=0, ii;
    float C0 = 0.0f, x0 = 0.0f, y0 = 0.0f, z0 = 0.0f, a = 0.0f, b = 0.0f, c = 0.0f, psi_gr1 = 0.0f, psi_gr2 = 0.0f, psi_gr3 = 0.0f;    

    FILE *fp = fopen(ModelParametersFilename, "r"); // read parameters file
    if( fp == NULL ) {
        printf("%s \n", "Cannot open the file");
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
        char tmpstr9[16];
        char tmpstr10[16];
        char tmpstr11[16];
        char tmpstr12[16];
        
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
                            break; }
                        if  (strcmp(tmpstr1,"Components") == 0) Components = atoi(tmpstr2);
                        printf("%s %i\n", "Components:", Components);
                        if (Components <= 0) {
                            printf("%s %i\n", "Components cannot be negative, the given value is", Components);
                            break; }
                        if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %21[^;];", tmpstr1, tmpstr2);
                        else {
                            break; }
                        if  (strcmp(tmpstr1,"TimeSteps") == 0) steps = atoi(tmpstr2);
                        if (steps <= 0) {
                            printf("%s %i\n", "TimeSteps cannot be negative, the given value is", steps);
                            break; }
                        printf("%s %i\n", "TimeSteps:", steps);
                        if (steps == 1) {
                            /**************************************************/
                            //printf("\n %s %i %s \n", "Stationary 3D model", ModelSelected, " is selected");
                            /* loop over all components */
                            for(ii=0; ii<Components; ii++) {
                                if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %21s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15[^;];", tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5, tmpstr6, tmpstr7, tmpstr8, tmpstr9, tmpstr10, tmpstr11, tmpstr12);
                                else {
                                    break; }
                                
                                if  (strcmp(tmpstr1,"Object") == 0) {
                                    C0 = (float)atof(tmpstr3); /* intensity */
                                    x0 = (float)atof(tmpstr4); /* x0 position */
                                    y0 = (float)atof(tmpstr5); /* y0 position */
                                    z0 = (float)atof(tmpstr6); /* z0 position */
                                    a = (float)atof(tmpstr7); /* a - size object */
                                    b = (float)atof(tmpstr8); /* b - size object */
                                    c = (float)atof(tmpstr9); /* c - size object */
                                    psi_gr1 = (float)atof(tmpstr10); /* rotation angle 1*/
                                    psi_gr2 = (float)atof(tmpstr11); /* rotation angle 2*/
                                    psi_gr3 = (float)atof(tmpstr12); /* rotation angle 3*/
                                }
                                else {
                                    break; }                               
                                // printf("\nObject : %s \nC0 : %f \nx0 : %f \ny0 : %f \nz0 : %f \na : %f \nb : %f \nc : %f \n", tmpstr2, C0, x0, y0, z0, a, b, c);                                                          

                                 TomoP3DObject_core(A, N, tmpstr2, C0, y0, x0, z0, a, b, c, psi_gr1, psi_gr2, psi_gr3, 0); /* python */
                            }
                        }
                        else {
                            /**************************************************/
                            //printf("\n %s \n", "Temporal model is selected");
                            /* temporal phantom 3D + time (4D) */
                            
                            float C1 = 0.0f, x1 = 0.0f, y1 = 0.0f, z1 = 0.0f, a1 = 0.0f, b1 = 0.0f, c1 = 0.0f, psi_gr1_1 = 0.0f, psi_gr2_1 = 0.0f, psi_gr3_1 = 0.0f;
                            /* loop over all components */
                            for(ii=0; ii<Components; ii++) {
                                
                                if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %21s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15[^;];", tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5, tmpstr6, tmpstr7, tmpstr8, tmpstr9, tmpstr10, tmpstr11, tmpstr12);
                                else {
                                    break; }
                                
                                if  (strcmp(tmpstr1,"Object") == 0) {
                                    C0 = (float)atof(tmpstr3); /* intensity */
                                    x0 = (float)atof(tmpstr4); /* x0 position */
                                    y0 = (float)atof(tmpstr5); /* y0 position */
                                    z0 = (float)atof(tmpstr6); /* y0 position */
                                    a = (float)atof(tmpstr7); /* a - size object */
                                    b = (float)atof(tmpstr8); /* b - size object */
                                    c = (float)atof(tmpstr9); /* b - size object */
                                    psi_gr1 = (float)atof(tmpstr10); /* rotation angle 1*/
                                    psi_gr2 = (float)atof(tmpstr11); /* rotation angle 2*/
                                    psi_gr3 = (float)atof(tmpstr12); /* rotation angle 3*/
                                }
                                else {
                                    break; }                                
                                
                                // printf("\nObject : %s \nC0 : %f \nx0 : %f \ny0 : %f \nz0 : %f \na : %f \nb : %f \n", tmpstr2, C0, x0, y0, z0, a, b, c);
                                
                                /* check Endvar relatedparameters */
                                if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %15s %15s %15s %15s %15s %15s %15s %15s %15s %15[^;];", tmpstr1, tmpstr3, tmpstr4, tmpstr5, tmpstr6, tmpstr7, tmpstr8, tmpstr9, tmpstr10, tmpstr11, tmpstr12);
                                else break; 
                                
                                if  (strcmp(tmpstr1,"Endvar") == 0) {
                                    C1 = (float)atof(tmpstr3); /* intensity */
                                    x1 = (float)atof(tmpstr4); /* x0 position */
                                    y1 = (float)atof(tmpstr5); /* y0 position */
                                    z1 = (float)atof(tmpstr6); /* z0 position */
                                    a1 = (float)atof(tmpstr7); /* a - size object */
                                    b1 = (float)atof(tmpstr8); /* b - size object */
                                    c1 = (float)atof(tmpstr9); /* c - size object */
                                    psi_gr1_1 = (float)atof(tmpstr10); /* rotation angle 1*/
                                    psi_gr2_1 = (float)atof(tmpstr11); /* rotation angle 2*/
                                    psi_gr3_1 = (float)atof(tmpstr12); /* rotation angle 3*/
                                }
                                else {
                                    printf("%s\n", "Cannot find 'Endvar' string in parameters file");
                                    break; }                        
                                
                                //printf("\nObject : %s \nC0 : %f \nx0 : %f \ny0 : %f \nz0 : %f \na : %f \nb : %f \nc : %f \n", tmpstr2, C0, x0, y0, z0, a1, b1, c1);
                                
                                /*now we know the initial parameters of the object and the final ones. We linearly extrapolate to establish steps and coordinates. */
                                /* calculating the full distance berween the start and the end points */
                                float distance = sqrtf(pow((x1 - x0),2) + pow((y1 - y0),2) + pow((z1 - z0),2));
                                float d_dist = distance/(steps-1); /*a step over line */
                                float C_step = (C1 - C0)/(steps-1);
                                float a_step = (a1 - a)/(steps-1);
                                float b_step = (b1 - b)/(steps-1);
                                float c_step = (c1 - c)/(steps-1);
                                float phi_rot_step1 = (psi_gr1_1 - psi_gr1)/(steps-1);
                                float phi_rot_step2 = (psi_gr2_1 - psi_gr2)/(steps-1);
                                float phi_rot_step3 = (psi_gr3_1 - psi_gr3)/(steps-1);
                                
                                int tt;
                                float x_t, y_t, z_t, a_t, b_t, c_t, C_t, phi1_t, phi2_t, phi3_t, d_step;
                                /* initialize */
                                x_t = x0; y_t = y0; z_t = z0; a_t = a; b_t = b; c_t = c; C_t = C0; phi1_t = psi_gr1; phi2_t = psi_gr2; phi3_t = psi_gr3; d_step = d_dist;
                                
                                /*loop over time frames*/
                                for(tt=0; tt < steps; tt++) {
                                    
                                    TomoP3DObject_core(A, N, tmpstr2, C_t, y_t, x_t, z_t, a_t, b_t, c_t, phi1_t, phi2_t, phi3_t, tt); /* python */
                                    
                                    /* calculating new coordinates of an object */
                                    if (distance != 0.0f) {
                                        float t = d_step/distance;
                                        x_t = (1-t)*x0 + t*x1;
                                        y_t = (1-t)*y0 + t*y1;
                                        z_t = (1-t)*z0 + t*z1;}
                                    else {
                                        x_t = x0;
                                        y_t = y0;
                                        z_t = z0;  }
                                    
                                    d_step += d_dist;
                                    a_t += a_step;
                                    b_t += b_step;
                                    c_t += c_step;
                                    C_t += C_step;
                                    phi1_t += phi_rot_step1;
                                    phi2_t += phi_rot_step2;
                                    phi3_t += phi_rot_step3;
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
