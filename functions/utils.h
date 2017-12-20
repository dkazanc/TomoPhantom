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
float parameters_check2D(float C0, float x0, float y0, float a, float b, float phi_rot);
float parameters_check3D(float C0, float x0, float y0, float z0, float a, float b, float c);
float su3(float *A, float psi1, float psi2, float psi3);
float mmtvc(float *A, float *V1, float *V2);
#ifdef __cplusplus
}
#endif
