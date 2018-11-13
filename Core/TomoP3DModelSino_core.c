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


#include "TomoP3DModelSino_core.h"

#define M_PI 3.14159265358979323846
#define M_PI2 1.57079632679
#define EPS 0.000000001
#define MAXCHAR 1000

/* Function to create 3D analytical projection data (parallel beam geometry) for 3D models
 *
 * Input Parameters:
 * - ModelNo - the model number from Phantom3DLibrary file
 * - Horiz_det - the number of horizontal detectors
 * - Vert_det - - the number of vertical detectors
 * - N - the size of the Phantom (N x N x N), currently  Vert_det must be set to N
 * - Theta_proj - a vector of projection anges in degrees
 * - C0 - intensity
 * - x0 - x0 position
 * - y0 - y0 position
 * - z0 - z0 position
 * - a  - size object
 * - b  - size object
 * - c - size object
 * - psi_gr1 - rotation angle1
 *
 * Output:
 * 1. The analytical phantom size of [N1 x N2 x N3] or temporal 4D phantom (N1 x N2 x N3 x time-frames)
 *    Note if Z1, Z2 indeces selected then the size can be [N1 x N2 x Z2-Z1]
 *
 */

float TomoP3DObjectSino_core(float *A, long Horiz_det, long Vert_det, long Z1, long Z2, long N, float *Theta_proj, int AngTot, char *Object,
        float C0, /* intensity */
        float x0, /* x0 position */
        float y0, /* y0 position */
        float z0, /* z0 position */
        float a, /* a - size object */
        float b, /* b - size object */
        float c, /* c - size object */
        float psi_gr1, /* rotation angle1 */
        float psi_gr2, /* rotation angle2 */
        float psi_gr3, /* rotation angle3 */
        long tt /*temporal index, 0 - for stationary */)
        
{
    int ll;
    long i, j, k, index, sub_vol_size;
    float *DetectorRange_Horiz_ar=NULL, DetectorRange_Horiz_min, DetectorRange_Horiz_max, *DetectorRange_Vert_ar=NULL, DetectorRange_Vert_min, DetectorRange_Vert_max, U_step, V_step, a22, b22, c2, *Tomorange_Z_Ar = NULL, *Zdel = NULL;
    
    float *AnglesRad=NULL;
    float Tomorange_Xmin, Tomorange_Xmax, H_x, multiplier, x00, y00, z00, a2, b2;
    
    DetectorRange_Horiz_max = (float)(Horiz_det)/(float)(N+1); /* horizontal detector range */
    DetectorRange_Horiz_min = -DetectorRange_Horiz_max;
    
    /* Here assuming that the size of the vertical detector array is always equal to Z-dim of the phantom */
    DetectorRange_Vert_max = (float)(Vert_det)/(float)(N+1); /* vertical detector range */
    DetectorRange_Vert_min = -DetectorRange_Vert_max;
    
    DetectorRange_Horiz_ar = malloc(Horiz_det*sizeof(float)); /* horizontal array */
    DetectorRange_Vert_ar = malloc(Vert_det*sizeof(float));   /* vertical array */
    
    U_step = (DetectorRange_Horiz_max - DetectorRange_Horiz_min)/(float)(Horiz_det-1);
    V_step = (DetectorRange_Vert_max - DetectorRange_Vert_min)/(float)(Vert_det-1);
    
    for(i=0; i<Horiz_det; i++) {DetectorRange_Horiz_ar[i] = (DetectorRange_Horiz_max) - (float)i*U_step;}
    for(i=0; i<Vert_det; i++) {DetectorRange_Vert_ar[i] = (DetectorRange_Vert_max) - (float)i*V_step;}
    
    sub_vol_size = Z2 - Z1;
    Tomorange_Xmin = -1.0f;
    Tomorange_Xmax = 1.0f;
    H_x = (Tomorange_Xmax - Tomorange_Xmin)/(float)(N);
    
    Tomorange_Z_Ar = malloc(N * sizeof(float));
    for (i = 0; i<N; i++) { Tomorange_Z_Ar[i] = Tomorange_Xmin + (float)i*H_x; }
    Zdel = malloc(N * sizeof(float));
    for (i = 0; i<N; i++) Zdel[i] = Tomorange_Z_Ar[i] - z0;
    
    /* convert to radians */
    AnglesRad = malloc(AngTot*sizeof(float));
    for(ll=0; ll<AngTot; ll++)  AnglesRad[ll] = (Theta_proj[ll])*((float)M_PI/180.0f);
    
    float alog2 = logf(2.0f);
    multiplier = (C0*(N/2.0f));
    
    float RS = 1.0f;
    /*Angles: TETA1, PSIs, FI1 ? */
    float TETAs, TETA1, FIs, FI1, PSI1, PSIs;
    float psi1, psi2, psi3;
    
    psi_gr1 = psi_gr1 + 90.0f;
    psi_gr2 = psi_gr2 + 90.0f;
    
    psi1 = psi_gr1*((float)M_PI/180.0f);
    psi2 = psi_gr2*((float)M_PI/180.0f);
    psi3 = psi_gr3*((float)M_PI/180.0f);
    
    float xwid=0.0f,ywid=0.0f,p00,ksi00,ksi1=0.0f; /* cuboid-related */
    float p,ksi,C,S,A2,B2,FI,CF,SF,P0,TF,PC,QM,DEL,XSYC,QP,SS,x11=0.0f,y11=0.0f;
    
    /* fix for centering */
    x00 = x0 + 0.5f*H_x;
    y00 = y0 - 0.5f*H_x;
    z00 = z0 - 0.5f*H_x;    
    
    /* parameters of an object have been extracted, now run the building module */
    /************************************************/
    if (strcmp("gaussian",Object) == 0) {
        a = 0.5f*a;
        b = 0.5f*b;
        c = 0.5f*c;
    }
    a22 = a*a;
    b22 = b*b;
    
    a2 = 1.0f/(a22);
    b2 = 1.0f/(b22);
    c2 = 1.0f/(c*c);
    
    float xh[3] = {0.0f, 0.0f, 0.0f};
    float xh1[3] = {0.0f, 0.0f, 0.0f};
    float vh1[3] = {0.0f, 0.0f, 0.0f};
    float aa[3] = {0.0f, 0.0f, 0.0f};
    float aa1[3] = {0.0f, 0.0f, 0.0f};
    float al[3] = {0.0f, 0.0f, 0.0f};
    
    float ai[3][3] = {
        {0.0f,0.0f,0.0f},
        {0.0f,0.0f,0.0f},
        {0.0f,0.0f,0.0f} };
        float bs[3][3] = {
            {0.0f,0.0f,0.0f},
            {0.0f,0.0f,0.0f},
            {0.0f,0.0f,0.0f} };
            float bsai[3][3] = {
                {0.0f,0.0f,0.0f},
                {0.0f,0.0f,0.0f},
                {0.0f,0.0f,0.0f} };
                
                matrot3(bs,psi1,psi2,psi3); /* rotation of 3x3 matrix */
                
                xh1[0] = x00; xh1[1] = y00; xh1[2] = z00;
                matvet3(bs,xh1,xh);  /* matrix-vector multiplication */
                
                float a_v, b_v, c_v, d_v, p1, p2, alh, bth, gmh;
                
                if (strcmp("cuboid",Object) == 0) {
                    /* the object is a cuboid */
                    x11 = 2.0f*x0;
                    y11 = 2.0f*y0 + 1.0f*H_x;
                    
                    xwid = b;
                    ywid = a;
                    c2 = 0.5f*c;
                    
                    if (psi3 < 0)  {ksi1 = (float)M_PI + psi3;}
                    else ksi1 = psi3;
                }
                
                float AA5, sin_2, cos_2, delta1, delta_sq, first_dr, AA2, AA3, AA6;
                AA5 = (N*C0*a*b);
                
#pragma omp parallel for shared(A) private(index,k,j,ll,TETAs,FIs,PSIs,aa1,aa,FI1,TETA1,PSI1,ai,bsai,vh1,al,a_v,b_v,c_v,d_v,p1,p2,alh,bth,gmh,sin_2, cos_2,delta1,delta_sq,first_dr, AA2, AA3, AA6, p00,ksi00,p,ksi,C,S,A2,B2,FI,CF,SF,P0,TF,PC,QM,DEL,XSYC,QP,SS)
                for(ll=0; ll<AngTot; ll++) {
                    
                    TETAs = AnglesRad[ll]; /* the variable projection angle (AnglesRad) */
                    
                    TETA1 = TETAs - M_PI2;
                    FIs =  0.0f;  /* always zero for the fixed source? */
                    PSIs = 0.0f;  /* always zero for the fixed source? */
                    
                    aa1[0]=-RS*sinf(TETAs)*cosf(FIs);
                    aa1[1]=-RS*sinf(TETAs)*sinf(FIs);
                    aa1[2]=-RS*cosf(TETAs);
                    
                    matvet3(bs,aa1,aa);  /* matrix-vector multiplication */
                    
                    /* calculation of inverse matrix */
                    FI1=-FIs;
                    TETA1=-TETA1;
                    PSI1=-PSIs;
                    matrot3(ai,PSI1,TETA1,FI1); /* rotation of 3x3 matrix */
                    
                    /* A transform matrix from projection space to object space */
                    matmat3(bs,ai,bsai);
                    vh1[0]=0.0f;
                    
                    /* the object is an ellipsoid */
                    for(j=0; j<Horiz_det; j++) {
                        // for(k=0; k<Vert_det; k++) {
                        for(k=Z1; k<Z2; k++) {
                            //index = tt*Vert_det*Horiz_det*AngTot + ll*Vert_det*Horiz_det + k*Horiz_det + j;
                            index = tt*Horiz_det*AngTot*sub_vol_size + ll*sub_vol_size*Horiz_det + (k - Z1)*Horiz_det + j;
                            
                            vh1[2]=DetectorRange_Horiz_ar[j];
                            vh1[1]=DetectorRange_Vert_ar[k];
                            
                            matvet3(bsai,vh1,al); /*matrix-vector multiplication */
                            
                            if (strcmp("ellipsoid",Object) == 0) {
                                a_v = powf((aa[0]/a),2) + powf((aa[1]/b),2) + powf((aa[2]/c),2);
                                b_v = aa[0]*(al[0]-xh[0])*a2 + aa[1]*(al[1]-xh[1])*b2 + aa[2]*(al[2]-xh[2])*c2;
                                c_v = powf(((al[0]-xh[0])/a),2) + powf(((al[1]-xh[1])/b), 2) + powf(((al[2]-xh[2])/c),2) - 1.0f;
                                d_v = b_v*b_v - a_v*c_v;
                                
                                if(d_v > 0) {
                                    p1 = -(sqrtf(d_v)+b_v)/a_v;
                                    p2 = (sqrtf(d_v)-b_v)/a_v;
                                    A[index] += (p2-p1)*multiplier;
                                }
                            }
                            if (strcmp("paraboloid",Object) == 0) {
                                /* the object is a parabola Lambda = 1 */
                                a_v = powf((aa[0]/a),2) + powf((aa[1]/b),2) + powf((aa[2]/c),2);
                                b_v = aa[0]*(al[0]-xh[0])*a2 + aa[1]*(al[1]-xh[1])*b2 + aa[2]*(al[2]-xh[2])*c2;
                                c_v = powf(((al[0]-xh[0])/a),2) + powf(((al[1]-xh[1])/b), 2) + powf(((al[2]-xh[2])/c),2) - 1.0f;
                                d_v = b_v*b_v - a_v*c_v;
                                
                                if(d_v > 0) {
                                    p1 = -(sqrtf(d_v)+b_v)/a_v;
                                    p2 = (sqrtf(d_v)-b_v)/a_v;
                                    A[index] += multiplier*(a_v/3.0f*(pow(p1,3.0f) - pow(p2,3.0f)) +  b_v*(pow(p1,2.0f) - pow(p2,2.0f)) + c_v*(p1-p2));
                                }
                            }
                            if (strcmp("gaussian",Object) == 0) {
                                /* The object is a volumetric gaussian */
                                alh=alog2*a2;
                                bth=alog2*b2;
                                gmh=alog2*c2;
                                
                                a_v = 1.0f/(alh*powf((aa[0]),2) + bth*powf((aa[1]),2) + gmh*powf((aa[2]),2));
                                b_v = aa[0]*alh*(al[0]-xh[0]) + aa[1]*bth*(al[1]-xh[1]) + aa[2]*gmh*(al[2]-xh[2]);
                                c_v = alh*powf(((al[0]-xh[0])),2) + bth*powf(((al[1]-xh[1])),2) + gmh*powf(((al[2]-xh[2])),2);
                                
                                A[index] += multiplier*sqrtf(M_PI*a_v)*expf((pow(b_v,2))*a_v-c_v);
                            }
                        }} /*main for j-k loop*/
                    
                    if (strcmp("elliptical_cylinder",Object) == 0)  {
                        sin_2 = powf((sinf(TETAs - psi3)),2);
                        cos_2 = powf((cosf(TETAs - psi3)),2);
                        
                        delta1 = 1.0f/(a22*sin_2 + b22*cos_2);
                        delta_sq = sqrtf(delta1);
                        first_dr = AA5*delta_sq;
                        AA2 = -x00*sinf(TETAs) + y00*cosf(TETAs);
                        for(k=Z1; k<Z2; k++) {
                            if (fabs(Zdel[k]) < c)  {
                                for(j=0; j<Horiz_det; j++) {
                                    AA3 = powf((DetectorRange_Horiz_ar[j] - AA2),2);
                                    AA6 = (AA3)*delta1;
                                    //index = tt*Vert_det*Horiz_det*AngTot + ll*Vert_det*Horiz_det + k*Horiz_det + j;
                                    index = tt*sub_vol_size*Horiz_det*AngTot + ll*sub_vol_size*Horiz_det + (k - Z1)*Horiz_det + j;
                                    if (AA6 < 1.0f) A[index] += first_dr*sqrtf(1.0f - AA6);
                                }
                            }
                        }
                    }
                    if (strcmp("cuboid",Object) == 0) {
                        /* the object is a cuboid */
                        ksi00 = AnglesRad[(AngTot-1)-ll];
                        for(k=Z1; k<Z2; k++) {
                            if (fabs(Zdel[k]) < c2) {
                                for(j=0; j< Horiz_det; j++) {
                                    p00 = DetectorRange_Horiz_ar[j];
                                    
                                    index = tt*sub_vol_size*Horiz_det*AngTot + ll*sub_vol_size*Horiz_det + (k - Z1)*Horiz_det + j;
                                    
                                    p = p00;
                                    ksi=ksi00;
                                    
                                    if (ksi > (float)M_PI) {
                                        ksi = ksi - (float)M_PI;
                                        p = -p00; }
                                    
                                    C = cosf(ksi); S = sinf(ksi);
                                    XSYC = -x11*S + y11*C;
                                    A2 = xwid*0.5f;
                                    B2 = ywid*0.5f;
                                    
                                    if ((ksi - ksi1) < 0.0f)  FI = (float)M_PI + ksi - ksi1;
                                    else FI = ksi - ksi1;
                                    
                                    if (FI > M_PI2) FI = (float)M_PI - FI;
                                    
                                    CF = cosf(FI);
                                    SF = sinf(FI);
                                    P0 = fabs(p-XSYC);
                                    
                                    SS = xwid/CF*C0;
                                    
                                    if (fabs(CF) <= (float)EPS) {
                                        SS = ywid*C0;
                                        if ((P0 - A2) > (float)EPS) SS=0.0f;
                                    }
                                    if (fabs(SF) <= (float)EPS) {
                                        SS = xwid*C0;
                                        if ((P0 - B2) > (float)EPS) SS=0.0f;
                                    }
                                    
                                    TF = SF/CF;
                                    PC = P0/CF;
                                    QP = B2+A2*TF;
                                    QM = QP+PC;
                                    if (QM > ywid) {
                                        DEL = P0+B2*CF;
                                        SS = ywid/SF*C0;
                                        if (DEL > (A2*SF)) SS = (QP-PC)/SF*C0;
                                    }
                                    if (QM > ywid) {
                                        DEL = P0+B2*CF;
                                        if (DEL > A2*SF) SS = (QP-PC)/SF*C0;
                                        else SS = ywid/SF*C0;
                                    }
                                    else SS = xwid/CF*C0;
                                    if (PC >= QP) SS=0.0f;
                                    
                                    A[index] += (N/2.0f)*SS;
                                } /*j-loop*/
                            }
                        } /*k-loop*/
                    }
                }
                /************************************************/
                free(AnglesRad);
                free(DetectorRange_Horiz_ar);
                free(DetectorRange_Vert_ar);
                free(Tomorange_Z_Ar);
                free(Zdel);
                return *A;
}

/********************Core Function*****************************/
float TomoP3DModelSino_core(float *A, int ModelSelected, long Horiz_det, long Vert_det, long Z1, long Z2, long N, float *Angl_vector, int AngTot, char* ModelParametersFilename)
{
    
    int Model = 0, Components = 0, steps = 0, counter = 0, ii;
    float C0 = 0.0f, x0 = 0.0f, y0 = 0.0f, z0 = 0.0f, a = 0.0f, b = 0.0f, c = 0.0f, psi_gr1 = 0.0f, psi_gr2 = 0.0f, psi_gr3 = 0.0f;
    
    FILE *fp = fopen(ModelParametersFilename, "r"); // read parameters file
    if (fp == NULL) {
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
            if (str[0] != '#') {
                sscanf(str, "%15s : %21[^;];", tmpstr1, tmpstr2);
                if (strcmp(tmpstr1, "Model") == 0)
                {
                    Model = atoi(tmpstr2);
                    if ((ModelSelected == Model) && (counter == 0)) {
                        /* check if we have a right model */
                        if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %21[^;];", tmpstr1, tmpstr2);
                        else {
                            break;
                        }
                        if (strcmp(tmpstr1, "Components") == 0) Components = atoi(tmpstr2);
                        printf("%s %i\n", "Components:", Components);
                        if (Components <= 0) {
                            printf("%s %i\n", "Components cannot be negative, the given value is", Components);
                            break;
                        }
                        if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %21[^;];", tmpstr1, tmpstr2);
                        else {
                            break;
                        }
                        if (strcmp(tmpstr1, "TimeSteps") == 0) steps = atoi(tmpstr2);
                        if (steps <= 0) {
                            printf("%s %i\n", "TimeSteps cannot be negative, the given value is", steps);
                            break;
                        }
                        printf("%s %i\n", "TimeSteps:", steps);
                        if (steps == 1) {
                            /**************************************************/
                            //printf("\n %s %i %s \n", "Stationary 3D model", ModelSelected, " is selected");
                            /* loop over all components */
                            for (ii = 0; ii<Components; ii++) {
                                if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %21s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15[^;];", tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5, tmpstr6, tmpstr7, tmpstr8, tmpstr9, tmpstr10, tmpstr11, tmpstr12);
                                else {
                                    break;
                                }
                                
                                if (strcmp(tmpstr1, "Object") == 0) {
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
                                    break;
                                }
                                // printf("\nObject : %s \nC0 : %f \nx0 : %f \ny0 : %f \nz0 : %f \na : %f \nb : %f \nc : %f \n", tmpstr2, C0, x0, y0, z0, a, b, c);
                                if ((strcmp("gaussian",tmpstr2) == 0) || (strcmp("paraboloid",tmpstr2) == 0) || (strcmp("ellipsoid",tmpstr2) == 0)) {
                                    TomoP3DObjectSino_core(A, Horiz_det, Vert_det, Z1, Z2, N, Angl_vector, AngTot, tmpstr2, C0, y0, -z0, -x0, b, a, c, psi_gr3, psi_gr2, psi_gr1, 0l); //python
                                }
                                else if (strcmp("elliptical_cylinder",tmpstr2) == 0) {
                                    TomoP3DObjectSino_core(A, Horiz_det, Vert_det, Z1, Z2, N, Angl_vector, AngTot, tmpstr2, C0, x0, -y0, z0, b, a, c, psi_gr3, psi_gr2, psi_gr1, 0l); //python
                                }
                                else {
                                    TomoP3DObjectSino_core(A, Horiz_det, Vert_det,  Z1, Z2, N, Angl_vector, AngTot, tmpstr2, C0, x0, y0, z0, a, b, c, psi_gr3, psi_gr2, -psi_gr1, 0l); //python
                                }
                            }
                        }
                        else {
                            /**************************************************/
                            //printf("\n %s \n", "Temporal model is selected");
                            /* temporal phantom 3D + time (4D) */
                            float C1 = 0.0f, x1 = 0.0f, y1 = 0.0f, z1 = 0.0f, a1 = 0.0f, b1 = 0.0f, c1 = 0.0f, psi_gr1_1 = 0.0f, psi_gr2_1 = 0.0f, psi_gr3_1 = 0.0f;
                            /* loop over all components */
                            for (ii = 0; ii<Components; ii++) {
                                
                                if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %21s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15[^;];", tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5, tmpstr6, tmpstr7, tmpstr8, tmpstr9, tmpstr10, tmpstr11, tmpstr12);
                                else {
                                    break;
                                }
                                
                                if (strcmp(tmpstr1, "Object") == 0) {
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
                                    break;
                                }
                                
                                // printf("\nObject : %s \nC0 : %f \nx0 : %f \ny0 : %f \nz0 : %f \na : %f \nb : %f \n", tmpstr2, C0, x0, y0, z0, a, b, c);
                                
                                /* check Endvar relatedparameters */
                                if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %15s %15s %15s %15s %15s %15s %15s %15s %15s %15[^;];", tmpstr1, tmpstr3, tmpstr4, tmpstr5, tmpstr6, tmpstr7, tmpstr8, tmpstr9, tmpstr10, tmpstr11, tmpstr12);
                                else break;
                                
                                if (strcmp(tmpstr1, "Endvar") == 0) {
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
                                    break;
                                }
                                
                                //printf("\nObject : %s \nC0 : %f \nx0 : %f \ny0 : %f \nz0 : %f \na : %f \nb : %f \nc : %f \n", tmpstr2, C0, x0, y0, z0, a1, b1, c1);
                                
                                /*now we know the initial parameters of the object and the final ones. We linearly extrapolate to establish steps and coordinates. */
                                /* calculating the full distance berween the start and the end points */
                                float distance = sqrtf(pow((x1 - x0), 2) + pow((y1 - y0), 2) + pow((z1 - z0), 2));
                                float d_dist = distance / (steps - 1); /*a step over line */
                                float C_step = (C1 - C0) / (steps - 1);
                                float a_step = (a1 - a) / (steps - 1);
                                float b_step = (b1 - b) / (steps - 1);
                                float c_step = (c1 - c) / (steps - 1);
                                float phi_rot_step1 = (psi_gr1_1 - psi_gr1) / (steps - 1);
                                float phi_rot_step2 = (psi_gr2_1 - psi_gr2) / (steps - 1);
                                float phi_rot_step3 = (psi_gr3_1 - psi_gr3) / (steps - 1);
                                
                                long tt;
                                float x_t, y_t, z_t, a_t, b_t, c_t, C_t, phi1_t, phi2_t, phi3_t, d_step;
                                /* initialize */
                                x_t = x0; y_t = y0; z_t = z0; a_t = a; b_t = b; c_t = c; C_t = C0; phi1_t = psi_gr1; phi2_t = psi_gr2; phi3_t = psi_gr3; d_step = d_dist;
                                
                                /*loop over time frames*/
                                for (tt = 0; tt < (long)steps; tt++) {
                                    
                                    if ((strcmp("gaussian",tmpstr2) == 0) || (strcmp("paraboloid",tmpstr2) == 0) || (strcmp("ellipsoid",tmpstr2) == 0)) {
                                        TomoP3DObjectSino_core(A, Horiz_det, Vert_det, Z1, Z2,  N, Angl_vector, AngTot, tmpstr2, C_t, y_t, -z_t, -x_t, b_t, a_t, c_t, phi3_t, phi2_t, phi1_t, tt); //python
                                    }
                                    else if (strcmp("elliptical_cylinder",tmpstr2) == 0) {
                                        TomoP3DObjectSino_core(A, Horiz_det, Vert_det,  Z1, Z2, N, Angl_vector, AngTot, tmpstr2, C_t, x_t, -y_t, z_t, b_t, a_t, c_t, phi3_t, phi2_t, phi1_t, tt); //python
                                    }
                                    else {
                                        TomoP3DObjectSino_core(A, Horiz_det, Vert_det,  Z1, Z2, N, Angl_vector, AngTot, tmpstr2, C_t, x_t, y_t, z_t, a_t, b_t, c_t, phi3_t, phi2_t, -phi1_t, tt); //python
                                    }
                                    /* calculating new coordinates of an object */
                                    if (distance != 0.0f) {
                                        float t = d_step / distance;
                                        x_t = (1 - t)*x0 + t*x1;
                                        y_t = (1 - t)*y0 + t*y1;
                                        z_t = (1 - t)*z0 + t*z1;
                                    }
                                    else {
                                        x_t = x0;
                                        y_t = y0;
                                        z_t = z0;
                                    }
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


