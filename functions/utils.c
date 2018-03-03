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

#include "utils.h"
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>

float extractSteps(int *steps, int ModelSelected, char *ModelParametersFilename)
{
    FILE *in_file = fopen(ModelParametersFilename, "r"); // read parameters file
    if (! in_file )
    {
        printf("%s %s\n", "Parameters file does not exist or cannot be read!", ModelParametersFilename);
        return 0;
    }		
        char tmpstr1[16];
        char tmpstr2[22];
    char tempbuff[100];
    while(!feof(in_file))
    {
        if (fgets(tempbuff,100,in_file)) {
            
            if(tempbuff[0] == '#') continue;
            
            sscanf(tempbuff, "%15s : %21[^;];", tmpstr1, tmpstr2);
            int Model = 0;
            
            if (strcmp(tmpstr1,"Model")==0) {
                Model = atoi(tmpstr2);
            }
            
            /*check if we have got the right model */
            if (ModelSelected == Model) {
                /* read the model parameters */
                printf("\nThe selected Model : %i \n", Model);
                if (fgets(tempbuff,100,in_file)) {
                    sscanf(tempbuff, "%15s : %21[^;];", tmpstr1, tmpstr2); }
                if (fgets(tempbuff,100,in_file)) {
                    sscanf(tempbuff, "%15s : %21[^;];", tmpstr1, tmpstr2); }                
                    // printf("<<%s>>\n",  tmpstr1);
                if  (strcmp(tmpstr1,"TimeSteps") == 0) {
                    steps[0] = atoi(tmpstr2);
                }
                else {printf("%s\n", "Number of steps should be >= 1");}              
            }            
        }
    }
    return *steps;
}
float parameters_check2D(float C0, float x0, float y0, float a, float b, float phi_rot)
{
    if (C0 <= 0) {
        printf("%s %f\n", "C0 (intensity) cannot be negative or equal to zero, the given value is", C0);
    }
    if ((x0 < -1) || (x0 > 1)) {
        printf("%s %f\n", "x0 (object position) must be in [-1,1] range, the given value is", x0);
        return -1;
    }
    if ((y0 < -1) || (y0 > 1)) {
        printf("%s %f\n", "y0 (object position) must be in [-1,1] range, the given value is", y0);
        return -1;
    }
    if ((a < -1) || (a > 1)) {
        printf("%s %f\n", "a (object size) must be in [-1,1] range, the given value is", a);
        return -1;
    }
    if ((b < -1) || (b > 1)) {
        printf("%s %f\n", "b (object position) must be in [-1,1] range, the given value is", b);
        return -1;
    }
    return 0;
}
float parameters_check3D(float C0, float x0, float y0, float z0, float a, float b, float c)
{
    if (C0 <= 0) {
        printf("%s %f\n", "C0 (intensity) cannot be negative or equal to zero, the given value is", C0);
    }
    if ((x0 < -1) || (x0 > 1)) {
        printf("%s %f\n", "x0 (object position) must be in [-1,1] range, the given value is", x0);
        return -1;
    }
    if ((y0 < -1) || (y0 > 1)) {
        printf("%s %f\n", "y0 (object position) must be in [-1,1] range, the given value is", y0);
        return -1;
    }
    if ((z0 < -1) || (z0 > 1)) {
        printf("%s %f\n", "z0 (object position) must be in [-1,1] range, the given value is", z0);
        return -1;
    }
    if ((a < -1) || (a > 1)) {
        printf("%s %f\n", "a (object size) must be in [-1,1] range, the given value is", a);
        return -1;
    }
    if ((b < -1) || (b > 1)) {
        printf("%s %f\n", "b (object position) must be in [-1,1] range, the given value is", b);
        return -1;
    }
    if ((c < -1) || (c > 1)) {
        printf("%s %f\n", "c (object position) must be in [-1,1] range, the given value is", c);
        return -1;
    }
    return 0;
}

/* rotation matrix routine */
float matrot3(float *A, float psi1, float psi2, float psi3)
{
    A[0] = cosf(psi1)*cosf(psi2)*cosf(psi3)-sinf(psi1)*sinf(psi3);
    A[1] = sinf(psi1)*cosf(psi2)*cosf(psi3)+cosf(psi1)*sinf(psi3);
    A[2] = -sinf(psi2)*cosf(psi3);
    A[3] = -cosf(psi1)*cosf(psi2)*sinf(psi3)-sinf(psi1)*cosf(psi3);
    A[4] = -sinf(psi1)*cosf(psi2)*sinf(psi3)+cosf(psi1)*cosf(psi3);
    A[5] = sinf(psi2)*sinf(psi3);
    A[6] = cosf(psi1)*sinf(psi2);
    A[7] = sinf(psi1)*sinf(psi2);
    A[8] = cosf(psi2);    
    return *A;
}

/*matrix-vector multiplication*/
float matvet3(float *A, float *V1, float *V2)
{
    int i, j, counter;
    
    counter = 0;
    for(i=0; i<3; i++) {
        V2[i] = 0.0f;
        for(j=0; j<3; j++) {
            V2[i] += A[counter]*V1[j];
            counter++;
        }}
    return *V2;
}

/*matrix-matrix multiplication*/
float matmat3(float *A, float *B, float *C)
{
	float Am[3][3];
	float Bm[3][3];
	float Cm[3][3];
	
	 Am[0][0]=A[0];
     Am[0][1]=A[1];
     Am[0][2]=A[2];
     Am[1][0]=A[3];
     Am[1][1]=A[4];
     Am[1][2]=A[5];
     Am[2][0]=A[6];
     Am[2][1]=A[7];
     Am[2][2]=A[8];
     
     Bm[0][0]=B[0];
     Bm[0][1]=B[1];
     Bm[0][2]=B[2];
     Bm[1][0]=B[3];
     Bm[1][1]=B[4];
     Bm[1][2]=B[5];
     Bm[2][0]=B[6];
     Bm[2][1]=B[7];
     Bm[2][2]=B[8];
     
     int i, j, k;
     for(i=0; i<3; i++) {        
        for(j=0; j<3; j++) {
		Cm[i][j]=0.0f;
			for(k=0; k<3; k++) {
			Cm[i][j] += Am[i][k]*Bm[k][j];
			}
	}}
	
	 C[0] = Cm[0][0];
     C[1] = Cm[0][1];
     C[2] = Cm[0][2];
     C[3] = Cm[1][0];
	 C[4] = Cm[1][1];
     C[5] = Cm[1][2];
     C[6] = Cm[2][0];
     C[7] = Cm[2][1];
     C[8] = Cm[2][2];
	
	
    return *C;
}
