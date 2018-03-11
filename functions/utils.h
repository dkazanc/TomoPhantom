/*
Copyright 2017 Daniil Kazantsev

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
#include "omp.h"
#ifdef __cplusplus
extern "C" {
#endif
float checkParams2D(int *params_switch, int ModelSelected, char *ModelParametersFilename);
float checkParams3D(int *params_switch, int ModelSelected, char *ModelParametersFilename);
float matrot3(float *Ad, float psi1, float psi2, float psi3);
float matvet3(float *Ad, float *V1, float *V2);
// float matmat3(float *A, float *B, float *C);
#ifdef __cplusplus
}
#endif
