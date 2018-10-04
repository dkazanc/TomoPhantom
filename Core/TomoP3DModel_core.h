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

//#include <matrix.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
#include "omp.h"
#include "CCPiDefines.h"
#include "utils.h"

CCPI_EXPORT float TomoP3DModel_core(float *A, int ModelSelected, long N1, long N2, long N3, long Z1, long Z2, char *ModelParametersFilename);
CCPI_EXPORT float TomoP3DObject_core(float *A, long N1, long N2, long N3, long Z1, long Z2, char *Object,
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
	long tt /*temporal index, 0 - for stationary */);
