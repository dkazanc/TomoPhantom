#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
#include "omp.h"

#define M_PI 3.14159265358979323846

/* Function to read from a file the required parameters to build 3D analytical model, see Phantom3DLibrary.dat to modify parameters
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
 * 11. phi_rot - rotation angle
 *
 * Output:
 * 1. The analytical phantom size of [N x N x N]
 *
 * License: Version 2.0
 * Copyright {2017} {Daniil Kazantsev, The University of Manchester}
 */

float buildPhantom3D_core_single(float *A, int N,  int Object,
                        float C0, /* intensity */
                        float x0, /* x0 position */
                        float y0, /* y0 position */
                        float z0, /* z0 position */
                        float a , /* a - size object */
                        float b , /* b - size object */
                        float c , /* c - size object */
                        float phi_rot /* phi - rotation angle */)
{
    int i, j, k;
    float *Tomorange_X_Ar=NULL, Tomorange_Xmin, Tomorange_Xmax, H_x, C1, C00, a22, a2, b22, b2, c22, c2, phi_rot_radian, sin_phi, cos_phi;
    float *Xdel = NULL, *Ydel = NULL, *Zdel = NULL, *Zdel2 = NULL, T;
    Tomorange_X_Ar = malloc(N*sizeof(float));
    Tomorange_Xmin = -1.0f;
    Tomorange_Xmax = 1.0f;
    H_x = (Tomorange_Xmax - Tomorange_Xmin)/(N);
    for(i=0; i<N; i++)  {Tomorange_X_Ar[i] = Tomorange_Xmin + (float)i*H_x;}
    C1 = -4.0f*logf(2.0f);
    


                          
                
				/* parameters of an object have been extracted, now run the building module */
				/************************************************/
				c22 = c*c;
				c2 = 1.0f/c22;
				if (Object == 4) c2 =4.0f*c2;
				
				phi_rot_radian = phi_rot*((float)M_PI/180.0f);
				sin_phi=sinf(phi_rot_radian); cos_phi=cosf(phi_rot_radian);
				
				Xdel = malloc(N*sizeof(float));
				Ydel = malloc(N*sizeof(float));
				Zdel = malloc(N*sizeof(float));
				Zdel2 = malloc(N*sizeof(float));
				for(i=0; i<N; i++)  {
					Xdel[i] = Tomorange_X_Ar[i] - x0;
					Ydel[i] = Tomorange_X_Ar[i] - y0;
					Zdel[i] = Tomorange_X_Ar[i] - z0;
				}
				for(i=0; i<N; i++)  {Zdel2[i] = c2*powf(Zdel[i],2);}
				
				if (Object == 1) {
					/* The object is a volumetric gaussian */
#pragma omp parallel for shared(A,Zdel2) private(k,i,j,a22,a2,b22,b2,T,C00)
					for(k=0; k<N; k++) {
						if (Zdel2[k] <= 1) {
							a22 = (a)*powf((1.0f - Zdel2[k]),2);
							if (a22 == 0.0f) a22 = 0.000001f;
							a2 = 1.0f/a22;
							b22 = (b)*powf((1.0f - Zdel2[k]),2);
							if (b22 == 0.0f) b22 = 0.000001f;
							b2 = 1.0f/b22;
							//C00 = C0*((exp(pow((1.0f - Zdel2[k]),2)))/2.7183f);
							C00 = C0*(powf((1.0f - Zdel2[k]),2));
							
							for(i=0; i<N; i++) {
								for(j=0; j<N; j++) {
									T = C1*(a2*powf((Xdel[i]*cos_phi + Ydel[j]*sin_phi),2) + b2*powf((-Xdel[i]*sin_phi + Ydel[j]*cos_phi),2));
									A[(k)*N*N + (i)*N + (j)] = A[(k)*N*N + (i)*N + (j)] + C00*expf(T);                     
								}}
						}
					} /*k-loop*/
				}
				else if (Object == 2) {
					/* the object is a parabola Lambda = 1/2 */
#pragma omp parallel for shared(A,Zdel2) private(k,i,j,a22,a2,b22,b2,T,C00)
					for(k=0; k<N; k++) {
						if (Zdel2[k] <= 1) {
							a22 = a*powf((1.0f - Zdel2[k]),2);
							a2 = 1.0f/a22;
							b22 = b*powf((1.0f - Zdel2[k]),2);
							b2 = 1.0f/b22;
							C00 = C0*(powf((1.0f - Zdel2[k]),2));
							
							for(i=0; i<N; i++) {
								for(j=0; j<N; j++) {
									T = a2*powf((Xdel[i]*cos_phi + Ydel[j]*sin_phi),2) + b2*powf((-Xdel[i]*sin_phi + Ydel[j]*cos_phi),2);
									if (T <= 1) T = C00*sqrtf(1.0f - T);
									else T = 0.0f;
									A[(k)*N*N + (i)*N + (j)] = A[(k)*N*N + (i)*N + (j)] + T;
								}}
						}
					} /*k-loop*/
				}
				else if (Object == 3) {
					/* the object is an elliptical disk */
					a22 = a*a;
					a2 = 1.0f/a22;
					b22 = b*b;
					b2 = 1.0f/b22;
#pragma omp parallel for shared(A,Zdel2) private(k,i,j,T)
					for(k=0; k<N; k++) {
						if (Zdel2[k] <= 1) {
							/* round objects case
							 * a22 = a*pow((1.0f - Zdel2[k]),2);
							 * a2 = 1.0f/a22;
							 * b22 = b*pow((1.0f - Zdel2[k]),2);
							 * b2 = 1.0f/b22;
							 */
							for(i=0; i<N; i++) {
								for(j=0; j<N; j++) {
									T = a2*powf((Xdel[i]*cos_phi + Ydel[j]*sin_phi),2) + b2*powf((-Xdel[i]*sin_phi + Ydel[j]*cos_phi),2);
									if (T <= 1) T = C0;
									else T = 0.0f;
									A[(k)*N*N + (i)*N + (j)] = A[(k)*N*N + (i)*N + (j)] + T;
								}}
						}
					} /*k-loop*/
				}
				else if (Object == 4) {
					/* the object is a parabola Lambda = 1 (12)*/
#pragma omp parallel for shared(A,Zdel2) private(k,i,j,a22,a2,b22,b2,T,C00)
					for(k=0; k<N; k++) {
						if (Zdel2[k] <= 1) {
							a22 = a*powf((1.0f - Zdel2[k]),2);
							a2 = 1.0f/a22;
							b22 = b*powf((1.0f - Zdel2[k]),2);
							b2 = 1.0f/b22;
							C00 = C0*(powf((1.0f - Zdel2[k]),2));
							
							for(i=0; i<N; i++) {
								for(j=0; j<N; j++) {
									T = (4.0f*a2)*powf((Xdel[i]*cos_phi + Ydel[j]*sin_phi),2) + (4.0f*b2)*powf((-Xdel[i]*sin_phi + Ydel[j]*cos_phi),2);
									if (T <= 1) T = C00*sqrtf(1.0f - T);
									else T = 0.0f;
									A[(k)*N*N + (i)*N + (j)] = A[(k)*N*N + (i)*N + (j)] + T;
								}}
						}
					} /*k-loop*/

				}
				else if (Object == 5) {
					/*the object is a cone (13)*/
#pragma omp parallel for shared(A,Zdel2) private(k,i,j,a22,a2,b22,b2,T,C00)
					for(k=0; k<N; k++) {
						if (Zdel2[k] <= 1) {
							a22 = a*powf((1.0f - Zdel2[k]),2);
							a2 = 1.0f/a22;
							b22 = b*powf((1.0f - Zdel2[k]),2);
							b2 = 1.0f/b22;
							C00 = C0*(powf((1.0f - Zdel2[k]),2));
							
							for(i=0; i<N; i++) {
								for(j=0; j<N; j++) {
									T = a2*powf((Xdel[i]*cos_phi + Ydel[j]*sin_phi),2) + b2*powf((-Xdel[i]*sin_phi + Ydel[j]*cos_phi),2);
									if (T <= 1) T = C00*(1.0f - sqrtf(T));
									else T = 0.0f;
									A[(k)*N*N + (i)*N + (j)] = A[(k)*N*N + (i)*N + (j)] + T;
								}}
						}
					} /*k-loop*/
				}
				else if (Object == 6) {
					/* the object is a parabola Lambda = 3/2 (14)*/
#pragma omp parallel for shared(A,Zdel2) private(k,i,j,a22,a2,b22,b2,T,C00)
					for(k=0; k<N; k++) {
						if (Zdel2[k] <= 1) {
							a22 = a*powf((1.0f - Zdel2[k]),2);
							a2 = 1.0f/a22;
							b22 = b*powf((1.0f - Zdel2[k]),2);
							b2 = 1.0f/b22;
							C00 = C0*(powf((1.0f - Zdel2[k]),2));
							
							for(i=0; i<N; i++) {
								for(j=0; j<N; j++) {
									T = a2*powf((Xdel[i]*cos_phi + Ydel[j]*sin_phi),2) + b2*powf((-Xdel[i]*sin_phi + Ydel[j]*cos_phi),2);
									if (T <= 1) T = C00*(powf((1.0f - T), 1.5f));
									else T = 0.0f;
									A[(k)*N*N + (i)*N + (j)] = A[(k)*N*N + (i)*N + (j)] + T;
								}}
						}
					} /*k-loop*/
				}
				else if (Object == 7) {
					/* the object is a rectangle (18) */
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
									A[(k)*N*N + (i)*N + (j)] = A[(k)*N*N + (i)*N + (j)] + T;
								}
							}
						}
                    }
				}
                else {
                    printf("%s\n", "No such object exist!");
                    return 0;
                }
                free(Xdel); free(Ydel); free(Zdel); free(Zdel2);
                    /************************************************/
    free(Tomorange_X_Ar);
    return *A;
}

float buildPhantom3D_core(float *A, int ModelSelected, int N, char *ModelParametersFilename)
{
    FILE *in_file = fopen(ModelParametersFilename, "r"); // read parameters file
    int ii;
    if (! in_file )
    {
        printf("%s %s\n", "Parameters file does not exist or cannot be read!", ModelParametersFilename);
		printf("Trying models/Phantom3DLibrary.dat");
		in_file = fopen("models/Phantom3DLibrary.dat","r");
		if(! in_file)
		{
			printf("models/Phantom3DLibrary.dat is not found");
			return 0;
		}
    }    	
    char tempbuff[100];
    while(!feof(in_file))
    {

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
    		
        if (fgets(tempbuff,100,in_file)) {
            
            if(tempbuff[0] == '#') continue;
            
            sscanf(tempbuff, "%15s : %15[^;];", tmpstr1, tmpstr2);
            /*printf("<<%s>>\n",  tmpstr1);*/
			int Model = 0, Components = 0;			
            
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
					int Object = 0;
                    float C0 = 0.0f, x0 = 0.0f, y0 = 0.0f, z0 = 0.0f, a = 0.0f, b = 0.0f, c = 0.0f, phi_rot = 0.0f;
					
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
                    }
					buildPhantom3D_core_single(A, N, Object, C0, x0, y0, z0, a, b, c, phi_rot);
				}				
			}
			
		}
	}		
}