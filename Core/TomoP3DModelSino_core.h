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

#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
#include "omp.h"
#include "CCPiDefines.h"
#include "utils.h"

CCPI_EXPORT float TomoP3DModelSino_core(float *A, int ModelSelected, long Horiz_det, long Vert_det, long Z1, long Z2, long N, float *Theta_proj, int AngTot, char* ModelParametersFilename);
CCPI_EXPORT float TomoP3DObjectSino_core(float *A, long Horiz_det, long Vert_det, long Z1, long Z2, long N, float *Angl_vector, int AngTot, char *Object,
	float C0, /* longensity */
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
