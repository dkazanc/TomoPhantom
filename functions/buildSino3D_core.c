#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
#include "omp.h"

#define M_PI 3.14159265358979323846
#define EPS 0.0000001




float buildSino3D_core_single(float *A, int N, int P, float *Th, int AngTot, int CenTypeIn, int Object, float C0, float x0, float y0, float z0, float a, float b, float c, float phi_rot)
{
    int i, j, k;
    float *Tomorange_X_Ar=NULL, Tomorange_Xmin, Tomorange_Xmax, Sinorange_Pmax, Sinorange_Pmin, H_p, H_x, C1, C00, a1, b1, a22, b22, c22, c2, phi_rot_radian;
    float *Zdel = NULL, *Zdel2 = NULL, *Sinorange_P_Ar=NULL, *AnglesRad=NULL;
    float AA5, sin_2, cos_2, delta1, delta_sq, first_dr, AA2, AA3, AA6, under_exp;
    
    Sinorange_Pmax = (float)(P)/(float)(N);
    Sinorange_Pmin = -Sinorange_Pmax;
    
    Sinorange_P_Ar = malloc(P*sizeof(float));
    H_p = (Sinorange_Pmax - Sinorange_Pmin)/(P-1);
    for(i=0; i<P; i++) {Sinorange_P_Ar[i] = Sinorange_Pmax - (float)i*H_p;}
    
    Tomorange_X_Ar = malloc(N*sizeof(float));
    Tomorange_Xmin = -1.0f;
    Tomorange_Xmax = 1.0f;
    H_x = (Tomorange_Xmax - Tomorange_Xmin)/(N);
    for(i=0; i<N; i++)  {Tomorange_X_Ar[i] = Tomorange_Xmin + (float)i*H_x;}
    AnglesRad = malloc(AngTot*sizeof(float));
    for(i=0; i<AngTot; i++)  AnglesRad[i] = (Th[i])*((float)M_PI/180.0f);
    
    C1 = -4.0f*logf(2.0f);
    
	if (CenTypeIn == 0) {
		/* matlab radon-iradon settings */
		x0 = x0 + H_x;
		y0 = y0 + H_x;
	}
	else {
		/* astra-toolbox settings */
		/*2D parallel beam*/
		x0 = x0 + 0.5f*H_x;
		y0 = y0 + 0.5f*H_x;
	}
	
	/* parameters of an object have been extracted, now run the building module */
	/************************************************/
	c22 = c*c;
	c2 = 1.0f/c22;
	if (Object == 4) c2 =4.0f*c2;
	phi_rot_radian = (phi_rot)*((float)M_PI/180.0f);
	
	Zdel = malloc(N*sizeof(float));
	Zdel2 = malloc(N*sizeof(float));
	for(i=0; i<N; i++)  {
		Zdel[i] = Tomorange_X_Ar[i] - z0;
	}
	for(i=0; i<N; i++)  {Zdel2[i] = c2*powf(Zdel[i],2);}
	
	if (Object == 1) {
		/* The object is a volumetric gaussian */
#pragma omp parallel for shared(A,Zdel2) private(k,i,j,a1,b1,C00,sin_2,cos_2,delta1,delta_sq,first_dr,under_exp,AA2,AA3,AA5)
		for(k=0; k<N; k++) {
			if (Zdel2[k] <= 1) {
				a1 = a*powf((1.0f - Zdel2[k]),2);
				if (a1 == 0.0f) a1 = (float)EPS;
				b1 = b*powf((1.0f - Zdel2[k]),2);
				if (b1 == 0.0f) b1 = (float)EPS;
				//C00 = C0*((exp(pow((1.0f - Zdel2[k]),1/2)))/2.7183f);
				C00 = C0*(powf((1.0f - Zdel2[k]),2));
				if (C00 == 0.0f) C00 = (float)EPS;
				
				AA5 = (N/2.0f)*(C00*sqrtf(a1)*sqrtf(b1)/2.0f)*sqrtf((float)M_PI/logf(2.0f));
				for(i=0; i<AngTot; i++) {
					sin_2 = powf((sinf((AnglesRad[i]) + phi_rot_radian)),2);
					cos_2 = powf((cosf((AnglesRad[i]) + phi_rot_radian)),2);
					delta1 = 1.0f/(a1*cos_2+b1*sin_2);
					delta_sq = sqrtf(delta1);
					first_dr = AA5*delta_sq;
					AA2 = -x0*cosf(AnglesRad[i])+y0*sinf(AnglesRad[i]); /*p0*/
					for(j=0; j<P; j++) {
						AA3 = powf((Sinorange_P_Ar[j] - AA2),2); /*(p-p0)^2*/
						under_exp = (C1*AA3)*delta1;
						A[(k)*P*AngTot + (j)*AngTot + (i)] = A[(k)*P*AngTot + (j)*AngTot + (i)] + first_dr*expf(under_exp);
					}}
			}
		} /*k-loop*/
	}
	else if (Object == 2) {
		/* the object is a parabola Lambda = 1/2 */
#pragma omp parallel for shared(A,Zdel2) private(k,i,j,a1,b1,C00,sin_2,cos_2,delta1,delta_sq,first_dr,AA2,AA3,AA5,AA6)
		for(k=0; k<N; k++) {
			if (Zdel2[k] <= 1) {
				a1 = a*powf((1.0f - Zdel2[k]),2);
				b1 = b*powf((1.0f - Zdel2[k]),2);
				C00 = C0*(powf((1.0f - Zdel2[k]),2));
				//printf("%f\n",C00);
				
				AA5 = (N/2.0f)*(((float)M_PI/2.0f)*C00*(sqrtf(a1))*(sqrtf(b1)));
				
				for(i=0; i<AngTot; i++) {
					sin_2 = powf((sinf((AnglesRad[i]) + phi_rot_radian)),2);
					cos_2 = powf((cosf((AnglesRad[i]) + phi_rot_radian)),2);
					delta1 = 1.0f/(a1*cos_2+b1*sin_2);
					delta_sq = sqrtf(delta1);
					first_dr = AA5*delta_sq;
					AA2 = -x0*cosf(AnglesRad[i])+y0*sinf(AnglesRad[i]); /*p0*/
					for(j=0; j<P; j++) {
						AA3 = powf((Sinorange_P_Ar[j] - AA2),2); /*(p-p0)^2*/
						AA6 = AA3*delta1;
						if (AA6 < 1.0f) {
							A[(k)*P*AngTot + (j)*AngTot + (i)] = A[(k)*P*AngTot + (j)*AngTot + (i)] + first_dr*(1.0f - AA6);
						}
					}}
			}
		} /*k-loop*/
	}
	else if (Object == 3) {
		/* the object is an elliptical disk */
		a22 = a*a;
		b22 = b*b;
		AA5 = (N*C0*a*b);
#pragma omp parallel for shared(A,Zdel2) private(k,i,j,sin_2,cos_2,delta1,delta_sq,first_dr,AA2,AA3,AA6)
		for(k=0; k<N; k++) {
			if (Zdel2[k] <= 1) {
				/* round objects case
				 * a22 = a*pow((1.0f - Zdel2[k]),2);
				 * a2 = 1.0f/a22;
				 * b22 = b*pow((1.0f - Zdel2[k]),2);
				 * b2 = 1.0f/b22;
				 */
				for(i=0; i<AngTot; i++) {
					sin_2 = powf((sinf((AnglesRad[i]) + phi_rot_radian)),2);
					cos_2 = powf((cosf((AnglesRad[i]) + phi_rot_radian)),2);
					delta1 = 1.0f/(a22*cos_2+b22*sin_2);
					delta_sq = sqrtf(delta1);
					first_dr = AA5*delta_sq;
					AA2 = -x0*cosf(AnglesRad[i])+y0*sinf(AnglesRad[i]); /*p0*/
					for(j=0; j<P; j++) {
						AA3 = powf((Sinorange_P_Ar[j] - AA2),2); /*(p-p0)^2*/
						AA6 = (AA3)*delta1;
						if (AA6 < 1.0f) {
							A[(k)*P*AngTot + (j)*AngTot + (i)] = A[(k)*P*AngTot + (j)*AngTot + (i)] + first_dr*sqrtf(1.0f - AA6);
						}
					}}
			}
		} /*k-loop*/
	}
	else if (Object == 4) {
		/* the object is a parabola Lambda = 1 (12)*/
#pragma omp parallel for shared(A,Zdel2) private(k,i,j,a1,b1,C00,sin_2,cos_2,delta1,delta_sq,first_dr,AA2,AA3,AA5,AA6)
		for(k=0; k<N; k++) {
			if (Zdel2[k] <= 1) {
				a1 = a*powf((1.0f - Zdel2[k]),2);
				b1 = b*powf((1.0f - Zdel2[k]),2);
				C00 = C0*(powf((1.0f - Zdel2[k]),2));
				
				AA5 = (N/2.0f)*(4.0f*((0.25f*sqrtf(a1)*sqrtf(b1)*C00)/2.5f));
				
				for(i=0; i<AngTot; i++) {
					sin_2 = powf((sinf((AnglesRad[i]) + phi_rot_radian)),2);
					cos_2 = powf((cosf((AnglesRad[i]) + phi_rot_radian)),2);
					delta1 = 1.0f/(0.25f*(a1)*cos_2+0.25f*b1*sin_2);
					delta_sq = sqrtf(delta1);
					first_dr = AA5*delta_sq;
					AA2 = -x0*cosf(AnglesRad[i])+y0*sinf(AnglesRad[i]); /*p0*/
					for(j=0; j<P; j++) {
						AA3 = powf((Sinorange_P_Ar[j] - AA2),2); /*(p-p0)^2*/
						AA6 = AA3*delta1;
						if (AA6 < 1.0f) {
							A[(k)*P*AngTot + (j)*AngTot + (i)] = A[(k)*P*AngTot + (j)*AngTot + (i)] + first_dr*(1.0f - AA6);
						}
					}}
			}
		} /*k-loop*/
	}
	else if (Object == 5) {
		/* the object is a rectangle (18) */
		// Rectang2DSino(Sinorange_P_Ar(j),AnglesRad(AnglesTot-ll+1),C00,-2*y0(i)-H_x,2*x0(i)+H_x,phi_rot_radian,b(i),a(i));
		float x00, y00,xwid,ywid,p00,ksi00,ksi1;
		float PI2,p,ksi,C,S,A2,B2,FI,CF,SF,P0,TF,PC,QM,DEL,XSYC,QP,SS;
		x00 = -2.0f*y0-H_x;
		y00 = 2.0f*x0+H_x;
		xwid = b;
		ywid = a;
		if (phi_rot_radian < 0)  {ksi1 = (float)M_PI + phi_rot_radian;}
		else ksi1 = phi_rot_radian;
		
// #pragma omp parallel for shared(A,Zdel) private(k,i,j,HX,HY,T)
		for(k=0; k<N; k++) {
			if (Zdel2[k] <= 1) {
				for(i=0; i<AngTot; i++) {
					ksi00 = AnglesRad[AngTot-i];
					for(j=0; j<P; j++) {
						p00 = Sinorange_P_Ar[j];
						
						PI2 = (float)M_PI*0.5f;
						p = p00;
						ksi=ksi00;
						
						if (ksi > M_PI) {
							ksi = ksi - (float)M_PI;
							p = -p00; }
						
						C = cosf(ksi); S = sinf(ksi);
						XSYC = -x00*S + y00*C;
						A2 = xwid*0.5f;
						B2 = ywid*0.5f;
						
						if ((ksi - ksi1) < 0.0f)  FI = (float)M_PI + ksi - ksi1;
						else FI = ksi - ksi1;
						
						if (FI > PI2) FI = (float)M_PI - FI;
						
						CF = cosf(FI);
						SF = sinf(FI);
						P0 = fabs(p-XSYC);
						
						SS = xwid/CF*C0;
						
						if (fabs(CF) <= EPS) {
							if ((P0 - A2) > EPS) {
								SS=0.0f;
								break;
							}
							SS = ywid*C0;
							break;
						}                                   
						if (fabs(SF) <= EPS) {
							if ((P0 - B2) > EPS) {
								SS=0.0f;
								break;
							}
							SS = xwid*C0;
							break;
						}
						
						TF = SF/CF;
						PC = P0/CF;
						QP = B2+A2*TF;
						
						if (PC >= QP) {
							SS=0.0f;
							break;
						}
						
						QM = QP+PC;
						if (QM > ywid) {
							DEL = P0+B2*CF;
							if (DEL > (A2*SF)) {
								SS = (QP-PC)/SF*C0;
								break;
							}
							
							SS = ywid/SF*C0;
							break;
						}                                        
						A[(k)*P*AngTot + (j)*AngTot + (i)] = A[(k)*P*AngTot + (j)*AngTot + (i)] + SS;
					}}
			}
		} /*k-loop*/
	}
	else {
		printf("%s\n", "No such object exist!");
		return 0;
	}
	free(Zdel); free(Zdel2);
	/************************************************/

    
    free(Tomorange_X_Ar); free(Sinorange_P_Ar); free(AnglesRad);
    return *A;	
}


/* Function to create 3D analytical sinograms (parallel beam geometry) to 3D phantoms using Phantom3DLibrary.dat
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

float buildSino3D_core(float *A, int ModelSelected, int N, int P, float *Th, int AngTot, int CenTypeIn, char *ModelParametersFilename)
{
    int i, ii;
    float *Tomorange_X_Ar=NULL, Tomorange_Xmin, Tomorange_Xmax, Sinorange_Pmax, Sinorange_Pmin, H_p, H_x, C1;
    float *Zdel = NULL, *Zdel2 = NULL, *Sinorange_P_Ar=NULL, *AnglesRad=NULL;
    
    Sinorange_Pmax = (float)(P)/(float)(N);
    Sinorange_Pmin = -Sinorange_Pmax;
    
    Sinorange_P_Ar = malloc(P*sizeof(float));
    H_p = (Sinorange_Pmax - Sinorange_Pmin)/(P-1);
    for(i=0; i<P; i++) {Sinorange_P_Ar[i] = Sinorange_Pmax - (float)i*H_p;}
    
    Tomorange_X_Ar = malloc(N*sizeof(float));
    Tomorange_Xmin = -1.0f;
    Tomorange_Xmax = 1.0f;
    H_x = (Tomorange_Xmax - Tomorange_Xmin)/(N);
    for(i=0; i<N; i++)  {Tomorange_X_Ar[i] = Tomorange_Xmin + (float)i*H_x;}
    AnglesRad = malloc(AngTot*sizeof(float));
    for(i=0; i<AngTot; i++)  AnglesRad[i] = (Th[i])*((float)M_PI/180.0f);
    
    C1 = -4.0f*logf(2.0f);
    FILE *in_file = fopen(ModelParametersFilename, "r"); // read parameters file
    
    if (! in_file )
    {
        printf("%s\n", "Parameters file does not exist or cannot be read!");
    }
    char tmpstr1[16];
    char tmpstr2[16];
    char tmpstr3[16];
    char tmpstr4[16];
    char tmpstr5[16];
    char tmpstr6[16];
    char tmpstr7[16];
    char tmpstr8[16];
    char tmpstr9[16];
    char tmpstr10[16];
    
    int Model = 0, Components = 0, Object = 0;
    float C0 = 0.0f, x0 = 0.0f, y0 = 0.0f, z0 = 0.0f, a = 0.0f, b = 0.0f, c = 0.0f, phi_rot = 0.0f;
    
    char tempbuff[100];
    while(!feof(in_file))
    {
        if (fgets(tempbuff,100,in_file)) {
            
            if(tempbuff[0] == '#') continue;
            
            sscanf(tempbuff, "%15s : %15[^;];", tmpstr1, tmpstr2);
            /*printf("<<%s>>\n",  tmpstr1);*/
            
            if (strcmp(tmpstr1,"Model")==0) {
                Model = atoi(tmpstr2);
            }
            
            /*check if we got the right model */
            if (ModelSelected == Model) {
                /* read the model parameters */
                printf("\nThe selected Model : %i \n", Model);
                if (fgets(tempbuff,100,in_file)) {
                    sscanf(tempbuff, "%15s : %15[^;];", tmpstr1, tmpstr2); }
                if  (strcmp(tmpstr1,"Components") == 0) {
                    Components = atoi(tmpstr2);
                }
                else {
                    printf("%s\n", "The number of components is unknown!");
                    return 0;
                }
                
                /* loop over all components */
                for(ii=0; ii<Components; ii++) {
                    
                    if (fgets(tempbuff,100,in_file)) {
                        sscanf(tempbuff, "%15s : %15s %15s %15s %15s %15s %15s %15s %15s %15[^;];", tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5, tmpstr6, tmpstr7, tmpstr8, tmpstr9, tmpstr10);
                    }
                    if  (strcmp(tmpstr1,"Object") == 0) {
                        Object = atoi(tmpstr2); /* analytical model */
                        C0 = (float)atof(tmpstr3); /* intensity */
                        x0 = (float)atof(tmpstr4); /* x0 position */
                        y0 = (float)atof(tmpstr5); /* y0 position */
                        z0 = (float)atof(tmpstr6); /* z0 position */
                        a = (float)atof(tmpstr7); /* a - size object */
                        b = (float)atof(tmpstr8); /* b - size object */
                        c = (float)atof(tmpstr9); /* c - size object */
                        phi_rot = (float)atof(tmpstr10); /* phi - rotation angle */
                        /*printf("\nObject : %i \nC0 : %f \nx0 : %f \nc : %f \n", Object, C0, x0, c);*/
						buildSino3D_core_single(A, N, P, Th, AngTot, CenTypeIn, Object, C0, x0,y0,z0,a,b,c,phi_rot);
                    }
                }
                break;
            }
        }
    }
    fclose(in_file);
    free(Tomorange_X_Ar); free(Sinorange_P_Ar); free(AnglesRad);
    return *A;
}