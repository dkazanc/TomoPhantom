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

float checkParams2D(nb::ndarray<int> params_switch, int ModelSelected, char *ModelParametersFilename);
float checkParams3D(nb::ndarray<int> params_switch, int ModelSelected, char *ModelParametersFilename);
float matrot3(nb::ndarray<float, nb::shape<3, 3>> Ad, float psi1, float psi2, float psi3);
float matvet3(nb::ndarray<float, nb::shape<3, 3>> Ad, nb::ndarray<float, nb::shape<3>> V1, nb::ndarray<float, nb::shape<3>> V2);
float matmat3(nb::ndarray<float, nb::shape<3, 3>> Am, nb::ndarray<float, nb::shape<3, 3>> Bm, nb::ndarray<float, nb::shape<3, 3>> Cm);
