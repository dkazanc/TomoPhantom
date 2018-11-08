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
#include "CCPiDefines.h"

#ifdef __cplusplus
extern "C" {
#endif
CCPI_EXPORT float checkParams2D(int *params_switch, int ModelSelected, char *ModelParametersFilename);
CCPI_EXPORT float checkParams3D(int *params_switch, int ModelSelected, char *ModelParametersFilename);
CCPI_EXPORT float matrot3(float Ad[3][3], float psi1, float psi2, float psi3);
CCPI_EXPORT float matvet3(float Ad[3][3], float V1[3], float V2[3]);
CCPI_EXPORT float matmat3(float Am[3][3], float Bm[3][3], float Cm[3][3]);
#ifdef __cplusplus
}
#endif
