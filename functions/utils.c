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

#define MAXCHAR 1000

/***********************************************************************************************/
float checkParams2D(int *params_switch, int ModelSelected, char *ModelParametersFilename)
{
    /* parameters check which returns 'params_switch' vector with 0/1 values */
    params_switch[0] = 1; /* parameters file */
    params_switch[1] = 0; /* model */
    params_switch[2] = 1; /* components */
    params_switch[3] = 1; /* time steps (actual value if > 1)*/
    params_switch[4] = 1; /* model types */
    params_switch[5] = 1; /* C0 */
    params_switch[6] = 1; /* x0 */
    params_switch[7] = 1; /* y0 */
    params_switch[8] = 1; /* a */
    params_switch[9] = 1; /* b */
    
    
    FILE *fp = fopen(ModelParametersFilename, "r"); // read parameters file    
    int Model=0, Components=0, steps = 0, counter=0, ii;
    float C0 = 0.0f, x0 = 0.0f, y0 = 0.0f, a = 0.0f, b = 0.0f;
    
    if( fp == NULL ) {
        printf("%s %s\n", "Parameters file does not exist or cannot be read!", ModelParametersFilename);
        params_switch[0] = 0;
        return 0;
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
                            printf("%s \n", "Unexpected the end of the line (Components) in parameters file");
                            params_switch[0] = 0;
                            break; }
                        if  (strcmp(tmpstr1,"Components") == 0) Components = atoi(tmpstr2);
                        printf("%s %i\n", "Components:", Components);
                        if (Components <= 0) {
                            printf("%s %i\n", "Components cannot be negative, the given value is", Components);
                            params_switch[2] = 0;                           
                            break; }
                        if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %21[^;];", tmpstr1, tmpstr2);
                        else {
                            printf("%s \n", "Unexpected the end of the line (TimeSteps) in parameters file");
                            params_switch[0] = 0;
                            break; }
                        if  (strcmp(tmpstr1,"TimeSteps") == 0) steps = atoi(tmpstr2);
                        if (steps <= 0) {
                            printf("%s %i\n", "TimeSteps cannot be negative, the given value is", steps);
                            params_switch[3] = 0;
                            break; }
                        else {params_switch[3] = steps;}
                        printf("%s %i\n", "TimeSteps:", steps);
                        if (steps == 1) {
                            /**************************************************/
                            // printf("\n %s %i %s \n", "Stationary 2D model", ModelSelected, " is selected");
                            
                            /* loop over all components */
                            for(ii=0; ii<Components; ii++) {
                                if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %21s %15s %15s %15s %15s %15s %15[^;];", tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5, tmpstr6, tmpstr7, tmpstr8);
                                else {
                                    printf("%s \n", "Unexpected the end of the line (objects loop) in parameters file");
                                    params_switch[0] = 0;
                                    break; }
                                
                                if  (strcmp(tmpstr1,"Object") == 0) {
                                    C0 = (float)atof(tmpstr3); /* intensity */
                                    x0 = (float)atof(tmpstr4); /* x0 position */
                                    y0 = (float)atof(tmpstr5); /* y0 position */
                                    a = (float)atof(tmpstr6); /* a - size object */
                                    b = (float)atof(tmpstr7); /* b - size object */
                                }
                                else {
                                    printf("%s \n", "Cannot find 'Object' string in parameters file");
                                    params_switch[0] = 0;
                                    break; }
                                
                                if ((strcmp("gaussian",tmpstr2) != 0) && (strcmp("parabola",tmpstr2) != 0) && (strcmp("ellipse",tmpstr2) != 0) && (strcmp("parabola1",tmpstr2) != 0) && (strcmp("cone",tmpstr2) != 0) && (strcmp("rectangle",tmpstr2) != 0) ) {
                                    printf("%s %s\n", "Unknown name of the object, the given name is", tmpstr2);
                                    params_switch[4] = 0;
                                    break; }
                                if (C0 == 0) {
                                    printf("%s %f\n", "C0 should not be equal to zero, the given value is", C0);
                                    params_switch[5] = 0;
                                    break; }
                                if ((x0 < -1) || (x0 > 1)) {
                                    printf("%s %f\n", "x0 (object position) must be in [-1,1] range, the given value is", x0);
                                    params_switch[6] = 0;
                                    break; }
                                if ((y0 < -1) || (y0 > 1)) {
                                    printf("%s %f\n", "y0 (object position) must be in [-1,1] range, the given value is", y0);
                                    params_switch[7] = 0;
                                    break; }
                                if ((a <= 0) || (a > 2)) {
                                    printf("%s %f\n", "a (object size) must be positive in [0,2] range, the given value is", a);
                                    params_switch[8] = 0;
                                    break; }
                                if ((b <= 0) || (b > 2)) {
                                    printf("%s %f\n", "b (object size) must be positive in [0,2] range, the given value is", b);
                                    params_switch[9] = 0;
                                    break; }
                                printf("\nObject : %s \nC0 : %f \nx0 : %f \ny0 : %f \na : %f \nb : %f \n", tmpstr2, C0, x0, y0,  a, b);                           
                            }
                        }
                        else {
                            /**************************************************/
                            //printf("\n %s %i %s \n", "Temporal 2D + time model", ModelSelected, " is selected");
                            /* temporal phantom 2D + time (3D) */
                            
                            float C1 = 0.0f, x1 = 0.0f, y1 = 0.0f, a1 = 0.0f, b1 = 0.0f;
                            /* loop over all components */
                            for(ii=0; ii<Components; ii++) {
                                
                                if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %15s %15s %15s %15s %15s %15s %15[^;];", tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5, tmpstr6, tmpstr7, tmpstr8);
                                else {                                   
                                   printf("%s \n", "Unexpected the end of the line (objects loop) in parameters file");
                                   params_switch[0] = 0;
                                   break; }
                                
                                if  (strcmp(tmpstr1,"Object") == 0) {
                                    C0 = (float)atof(tmpstr3); /* intensity */
                                    x0 = (float)atof(tmpstr4); /* x0 position */
                                    y0 = (float)atof(tmpstr5); /* y0 position */
                                    a = (float)atof(tmpstr6); /* a - size object */
                                    b = (float)atof(tmpstr7); /* b - size object */
                                }
                                else {
                                    printf("%s \n", "Cannot find 'Object' string in parameters file");
                                    params_switch[0] = 0;
                                    break; }
                                
                                if ((strcmp("gaussian",tmpstr2) != 0) && (strcmp("parabola",tmpstr2) != 0) && (strcmp("ellipse",tmpstr2) != 0) && (strcmp("parabola1",tmpstr2) != 0) && (strcmp("cone",tmpstr2) != 0) && (strcmp("rectangle",tmpstr2) != 0) ) {
                                    printf("%s %s\n", "Unknown name of the object, the given name is", tmpstr2);
                                    params_switch[4] = 0;
                                    break; }
                                if (C0 == 0) {
                                    printf("%s %f\n", "C0 should not be equal to zero, the given value is", C0);
                                    params_switch[5] = 0;
                                    break; }
                                if ((x0 < -1) || (x0 > 1)) {
                                    printf("%s %f\n", "x0 (object position) must be in [-1,1] range, the given value is", x0);
                                    params_switch[6] = 0;
                                    break; }
                                if ((y0 < -1) || (y0 > 1)) {
                                    printf("%s %f\n", "y0 (object position) must be in [-1,1] range, the given value is", y0);
                                    params_switch[7] = 0;
                                    break; }
                                if ((a <= 0) || (a > 2)) {
                                    printf("%s %f\n", "a (object size) must be positive in [0,2] range, the given value is", a);
                                    params_switch[8] = 0;
                                    break; }
                                if ((b <= 0) || (b > 2)) {
                                    printf("%s %f\n", "b (object size) must be positive in [0,2] range, the given value is", b);
                                    params_switch[9] = 0;
                                    break; }
                                printf("\nObject : %s \nC0 : %f \nx0 : %f \ny0 : %f \na : %f \nb : %f \n", tmpstr2, C0, x0, y0, a, b);
                                
                                /* check Endvar relatedparameters */
                                if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %15s %15s %15s %15s %15s %15[^;];", tmpstr1, tmpstr3, tmpstr4, tmpstr5, tmpstr6, tmpstr7, tmpstr8);
                                else {
                                    printf("%s \n", "Unexpected the end of the line (Endvar loop) in parameters file");
                                    params_switch[0] = 0;
                                    break; }
                                
                                if  (strcmp(tmpstr1,"Endvar") == 0) {
                                    C1 = (float)atof(tmpstr3); /* intensity */
                                    x1 = (float)atof(tmpstr4); /* x0 position */
                                    y1 = (float)atof(tmpstr5); /* y0 position */
                                    a1 = (float)atof(tmpstr6); /* a - size object */
                                    b1 = (float)atof(tmpstr7); /* b - size object */
                                }
                                else {
                                    printf("%s \n", "Cannot find 'Endvar' string in parameters file");
                                    params_switch[0] = 0;
                                    break; }
                                
                                if (C1 == 0) {
                                    printf("%s %f\n", "Endvar C1 should not be equal to zero, the given value is", C1);
                                    params_switch[5] = 0;
                                    break; }
                                if ((x1 < -1) || (x1 > 1)) {
                                    printf("%s %f\n", "Endvar x1 (object position) must be in [-1,1] range, the given value is", x1);
                                    params_switch[6] = 0;
                                    break; }
                                if ((y1 < -1) || (y1 > 1)) {
                                    printf("%s %f\n", "Endvar y1 (object position) must be in [-1,1] range, the given value is", y1);
                                    params_switch[7] = 0;
                                    break; }
                                if ((a1 <= 0) || (a1 > 2)) {
                                    printf("%s %f\n", "Endvar a1 (object size) must be positive in [0,2] range, the given value is", a1);
                                    params_switch[8] = 0;
                                    break; }
                                if ((b1 <= 0) || (b1 > 2)) {
                                    printf("%s %f\n", "Endvar b1 (object size) must be positive in [0,2] range, the given value is", b1);
                                    params_switch[9] = 0;
                                    break; }
                             } /*components loop*/
                        }
                        params_switch[1] = 1;
                        counter++;
                    } /*model if*/
                }
            }
        }
    }    
    if  (params_switch[1] == 0) printf("%s %i\n", "No such model available, the given value is", Model);
    fclose(fp);
    return *params_switch;
}                
                            

/***********************************************************************************************/
float checkParams3D(int *params_switch, int ModelSelected, char *ModelParametersFilename)
{
    FILE *fp;
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
    
    params_switch[0] = 1; /* parameters file */
    params_switch[1] = 0; /* model */
    params_switch[2] = 1; /* components */
    params_switch[3] = 1; /* steps (actual value if > 1)*/
    params_switch[4] = 1; /* model types */
    params_switch[5] = 1; /* C0 */
    params_switch[6] = 1; /* x0 */
    params_switch[7] = 1; /* y0 */
    params_switch[8] = 1; /* z0 */
    params_switch[9] = 1; /* a */
    params_switch[10] = 1; /* b */
    params_switch[11] = 1; /* c */
    
    int Model=0, Components=0, steps = 0, counter=0, ii;
    float C0 = 0.0f, x0 = 0.0f, y0 = 0.0f, z0 = 0.0f, a = 0.0f, b = 0.0f, c = 0.0f;
    
    fp = fopen(ModelParametersFilename, "r");
    if (fp == NULL){
        printf("Could not open file %s",ModelParametersFilename);
        params_switch[0] = 0;
        return 0;
    }
    else {
        while (fgets(str, MAXCHAR, fp) != NULL)
        {
            /* work with non-# commented lines */
            if(str[0] != '#') {
                sscanf(str, "%15s : %21[^;];", tmpstr1, tmpstr2);
                if (strcmp(tmpstr1,"Model")==0) Model = atoi(tmpstr2);   {
                    if ((ModelSelected == Model) && (counter == 0)) {
                        /* check if we have a right model */
                        if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %21[^;];", tmpstr1, tmpstr2);
                        else {
                            params_switch[0] = 0;
                            printf("%s\n", "Unexpected the end of the line (Components) in parameters file");
                            break; }
                        if  (strcmp(tmpstr1,"Components") == 0) Components = atoi(tmpstr2);
                        printf("%s %i\n", "Components:", Components);
                        if (Components <= 0) {
                            params_switch[2] = 0; /*checking components*/
                            printf("%s %i\n", "Components cannot be negative, the given value is", Components);
                            break; }
                        if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %21[^;];", tmpstr1, tmpstr2);
                        else {
                            params_switch[0] = 0;
                            printf("%s\n", "Unexpected the end of the line (TimeSteps) in parameters file");
                            break; }
                        if  (strcmp(tmpstr1,"TimeSteps") == 0) steps = atoi(tmpstr2);
                        if (steps <= 0) {
                            params_switch[3] = 0;
                            printf("%s %i\n", "TimeSteps cannot be negative, the given value is", steps);
                            break; }
                        else {params_switch[3] = steps;}
                        printf("%s %i\n", "TimeSteps:", steps);
                        if (steps == 1) {
                            printf("\n %s %i %s \n", "Stationary 3D model", ModelSelected, " is selected");
                            /* loop over all components */
                            for(ii=0; ii<Components; ii++) {
                                
                                if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %21s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15[^;];", tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5, tmpstr6, tmpstr7, tmpstr8, tmpstr9, tmpstr10, tmpstr11, tmpstr12);
                                else {
                                    params_switch[0] = 0;
                                    printf("%s\n", "Unexpected the end of the line (objecta loop) in parameters file");
                                    break; }
                                
                                if  (strcmp(tmpstr1,"Object") == 0) {
                                    C0 = (float)atof(tmpstr3); /* intensity */
                                    x0 = (float)atof(tmpstr4); /* x0 position */
                                    y0 = (float)atof(tmpstr5); /* y0 position */
                                    z0 = (float)atof(tmpstr6); /* y0 position */
                                    a = (float)atof(tmpstr7); /* a - size object */
                                    b = (float)atof(tmpstr8); /* b - size object */
                                    c = (float)atof(tmpstr9); /* b - size object */
                                }
                                else {
                                    params_switch[0] = 0;
                                    printf("%s\n", "Cannot find 'Object' string in parameters file");
                                    break; }
                                
                                if ((strcmp("gaussian",tmpstr2) != 0) && (strcmp("paraboloid",tmpstr2) != 0) && (strcmp("ellipsoid",tmpstr2) != 0) && (strcmp("cone",tmpstr2) != 0) && (strcmp("cuboid",tmpstr2) != 0) && (strcmp("elliptical_cylinder",tmpstr2) != 0) ) {
                                    printf("%s %s\n", "Unknown name of the object, the given name is", tmpstr2);
                                    params_switch[4] = 0;
                                    break; }
                                if (C0 == 0) {
                                    params_switch[5] = 0;
                                    printf("%s %f\n", "C0 should not be equal to zero, the given value is", C0);
                                    break; }
                                if ((x0 < -1) || (x0 > 1)) {
                                    params_switch[6] = 0;
                                    printf("%s %f\n", "x0 (object position) must be in [-1,1] range, the given value is", x0);
                                    break; }
                                if ((y0 < -1) || (y0 > 1)) {
                                    params_switch[7] = 0;
                                    printf("%s %f\n", "y0 (object position) must be in [-1,1] range, the given value is", y0);
                                    break; }
                                if ((z0 < -1) || (z0 > 1)) {
                                    params_switch[8] = 0;
                                    printf("%s %f\n", "z0 (object position) must be in [-1,1] range, the given value is", z0);
                                    break; }
                                if ((a <= 0) || (a > 2)) {
                                    params_switch[9] = 0;
                                    printf("%s %f\n", "a (object size) must be positive in [0,2] range, the given value is", a);
                                    break; }
                                if ((b <= 0) || (b > 2)) {
                                    params_switch[10] = 0;
                                    printf("%s %f\n", "b (object size) must be positive in [0,2] range, the given value is", b);
                                    break; }
                                if ((c <= 0) || (c > 2)) {
                                    params_switch[11] = 0;
                                    printf("%s %f\n", "c (object size) must be positive in [0,2] range, the given value is", c);
                                    break; }
                                printf("\nObject : %s \nC0 : %f \nx0 : %f \ny0 : %f \nz0 : %f \na : %f \nb : %f \nc : %f \n", tmpstr2, C0, x0, y0, z0, a, b, c);
                            }
                        }
                        else {                            
                            printf("\n %s %i %s \n", "Temporal 3D+time model", ModelSelected, " is selected");
                            /* loop over all components */
                            for(ii=0; ii<Components; ii++) {
                                
                                if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %21s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15[^;];", tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5, tmpstr6, tmpstr7, tmpstr8, tmpstr9, tmpstr10, tmpstr11, tmpstr12);
                                else {
                                    params_switch[0] = 0;
                                    printf("%s\n", "Unexpected the end of the line (objecta loop) in parameters file");
                                    break; }
                                
                                if  (strcmp(tmpstr1,"Object") == 0) {
                                    C0 = (float)atof(tmpstr3); /* intensity */
                                    x0 = (float)atof(tmpstr4); /* x0 position */
                                    y0 = (float)atof(tmpstr5); /* y0 position */
                                    z0 = (float)atof(tmpstr6); /* y0 position */
                                    a = (float)atof(tmpstr7); /* a - size object */
                                    b = (float)atof(tmpstr8); /* b - size object */
                                    c = (float)atof(tmpstr9); /* b - size object */
                                }
                                else {
                                    params_switch[0] = 0;
                                    printf("%s\n", "Cannot find 'Object' string in parameters file");
                                    break; }
                                
                                if ((strcmp("gaussian",tmpstr2) != 0) && (strcmp("paraboloid",tmpstr2) != 0) && (strcmp("ellipsoid",tmpstr2) != 0) && (strcmp("cone",tmpstr2) != 0) && (strcmp("cuboid",tmpstr2) != 0) && (strcmp("elliptical_cylinder",tmpstr2) != 0) ) {
                                    printf("%s %s\n", "Unknown name of the object, the given name is", tmpstr2);
                                    params_switch[4] = 0;
                                    break; }
                                if (C0 == 0) {
                                    params_switch[5] = 0;
                                    printf("%s %f\n", "C0 should not be equal to zero, the given value is", C0);
                                    break; }
                                if ((x0 < -1) || (x0 > 1)) {
                                    params_switch[6] = 0;
                                    printf("%s %f\n", "x0 (object position) must be in [-1,1] range, the given value is", x0);
                                    break; }
                                if ((y0 < -1) || (y0 > 1)) {
                                    params_switch[7] = 0;
                                    printf("%s %f\n", "y0 (object position) must be in [-1,1] range, the given value is", y0);
                                    break; }
                                if ((z0 < -1) || (z0 > 1)) {
                                    params_switch[8] = 0;
                                    printf("%s %f\n", "z0 (object position) must be in [-1,1] range, the given value is", z0);
                                    break; }
                                if ((a <= 0) || (a > 2)) {
                                    params_switch[9] = 0;
                                    printf("%s %f\n", "a (object size) must be positive in [0,2] range, the given value is", a);
                                    break; }
                                if ((b <= 0) || (b > 2)) {
                                    params_switch[10] = 0;
                                    printf("%s %f\n", "b (object size) must be positive in [0,2] range, the given value is", b);
                                    break; }
                                if ((c <= 0) || (c > 2)) {
                                    params_switch[11] = 0;
                                    printf("%s %f\n", "c (object size) must be positive in [0,2] range, the given value is", c);
                                    break; }
                                // printf("\nObject : %s \nC0 : %f \nx0 : %f \ny0 : %f \nz0 : %f \na : %f \nb : %f \n", tmpstr2, C0, x0, y0, z0, a, b, c);
                                
                                /* check Endvar relatedparameters */
                                if (fgets(str, MAXCHAR, fp) != NULL) sscanf(str, "%15s : %15s %15s %15s %15s %15s %15s %15s %15s %15s %15[^;];", tmpstr1, tmpstr3, tmpstr4, tmpstr5, tmpstr6, tmpstr7, tmpstr8, tmpstr9, tmpstr10, tmpstr11, tmpstr12);
                                else {
                                    params_switch[0] = 0;
                                    printf("%s\n", "Unexpected the end of the line (Endvar loop) in parameters file");
                                    break; }
                                
                                if  (strcmp(tmpstr1,"Endvar") == 0) {
                                    C0 = (float)atof(tmpstr3); /* intensity */
                                    x0 = (float)atof(tmpstr4); /* x0 position */
                                    y0 = (float)atof(tmpstr5); /* y0 position */
                                    z0 = (float)atof(tmpstr6); /* y0 position */
                                    a = (float)atof(tmpstr7); /* a - size object */
                                    b = (float)atof(tmpstr8); /* b - size object */
                                    c = (float)atof(tmpstr9); /* b - size object */
                                }
                                else {
                                    params_switch[0] = 0;
                                    printf("%s\n", "Cannot find 'Endvar' string in parameters file");
                                    break; }
                                
                                if (C0 == 0) {
                                    params_switch[5] = 0;
                                    printf("%s %f\n", "Endvar C0 should not be equal to zero, the given value is", C0);
                                    break; }
                                if ((x0 < -1) || (x0 > 1)) {
                                    params_switch[6] = 0;
                                    printf("%s %f\n", "Endvar x0 (object position) must be in [-1,1] range, the given value is", x0);
                                    break; }
                                if ((y0 < -1) || (y0 > 1)) {
                                    params_switch[7] = 0;
                                    printf("%s %f\n", "Endvar y0 (object position) must be in [-1,1] range, the given value is", y0);
                                    break; }
                                if ((z0 < -1) || (z0 > 1)) {
                                    params_switch[8] = 0;
                                    printf("%s %f\n", "Endvar z0 (object position) must be in [-1,1] range, the given value is", z0);
                                    break; }
                                if ((a <= 0) || (a > 2)) {
                                    params_switch[9] = 0;
                                    printf("%s %f\n", "Endvar a (object size) must be positive in [0,2] range, the given value is", a);
                                    break; }
                                if ((b <= 0) || (b > 2)) {
                                    params_switch[10] = 0;
                                    printf("%s %f\n", "Endvar b (object size) must be positive in [0,2] range, the given value is", b);
                                    break; }
                                if ((c <= 0) || (c > 2)) {
                                    params_switch[11] = 0;
                                    printf("%s %f\n", "Endvar c (object size) must be positive in [0,2] range, the given value is", c);
                                    break; }
                                printf("\nObject : %s \nC0 : %f \nx0 : %f \ny0 : %f \nz0 : %f \na : %f \nb : %f \nc : %f \n", tmpstr2, C0, x0, y0, z0, a, b, c);
                            } /*components loop*/
                        }
                        params_switch[1] = 1;
                        counter++;
                    }
                }
            }
        }       
    }
    if  (params_switch[1] == 0) printf("%s %i\n", "No such model available, the given value is", Model);
    fclose(fp);
    return *params_switch;
}

/***********************************************************************************************/
/* rotation matrix routine */
float matrot3(float *Ad, float psi1, float psi2, float psi3)
{
    Ad[0] = cosf(psi1)*cosf(psi2)*cosf(psi3)-sinf(psi1)*sinf(psi3);
    Ad[1] = sinf(psi1)*cosf(psi2)*cosf(psi3)+cosf(psi1)*sinf(psi3);
    Ad[2] = -sinf(psi2)*cosf(psi3);
    Ad[3] = -cosf(psi1)*cosf(psi2)*sinf(psi3)-sinf(psi1)*cosf(psi3);
    Ad[4] = -sinf(psi1)*cosf(psi2)*sinf(psi3)+cosf(psi1)*cosf(psi3);
    Ad[5] = sinf(psi2)*sinf(psi3);
    Ad[6] = cosf(psi1)*sinf(psi2);
    Ad[7] = sinf(psi1)*sinf(psi2);
    Ad[8] = cosf(psi2);
    return *Ad;
}

/*matrix-vector multiplication*/
float matvet3(float *Ad, float *V1, float *V2)
{
    int l, m, counterT;
    
    counterT = 0;
    for(l=0; l<3; l++) {
        V2[l] = 0.0f;
        for(m=0; m<3; m++) {
            V2[l] += Ad[counterT]*V1[m];
            counterT++;
        }}
    return *V2;
}

/*matrix-matrix multiplication*/
// float matmat3(float *A, float *B, float *C)
// {
//     float Am[3][3];
//     float Bm[3][3];
//     float Cm[3][3];
//     
//     Am[0][0]=A[0];
//     Am[0][1]=A[1];
//     Am[0][2]=A[2];
//     Am[1][0]=A[3];
//     Am[1][1]=A[4];
//     Am[1][2]=A[5];
//     Am[2][0]=A[6];
//     Am[2][1]=A[7];
//     Am[2][2]=A[8];
//     
//     Bm[0][0]=B[0];
//     Bm[0][1]=B[1];
//     Bm[0][2]=B[2];
//     Bm[1][0]=B[3];
//     Bm[1][1]=B[4];
//     Bm[1][2]=B[5];
//     Bm[2][0]=B[6];
//     Bm[2][1]=B[7];
//     Bm[2][2]=B[8];
//     
//     int i, j, k;
//     for(i=0; i<3; i++) {
//         for(j=0; j<3; j++) {
//             Cm[i][j]=0.0f;
//             for(k=0; k<3; k++) {
//                 Cm[i][j] += Am[i][k]*Bm[k][j];
//             }
//         }}
//     
//     C[0] = Cm[0][0];
//     C[1] = Cm[0][1];
//     C[2] = Cm[0][2];
//     C[3] = Cm[1][0];
//     C[4] = Cm[1][1];
//     C[5] = Cm[1][2];
//     C[6] = Cm[2][0];
//     C[7] = Cm[2][1];
//     C[8] = Cm[2][2];
//     
//     
//     return *C;
// }
