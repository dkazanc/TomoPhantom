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

float TomoP2DSinoNum_core(nb::ndarray<float> Sinogram, nb::ndarray<float> Phantom, int dimX, int DetSize, nb::ndarray<float> Theta, int ThetaLength, int sys);
float BilinearInterpolation(nb::ndarray<float> Phantom_pad, nb::ndarray<float> B, int DetSize, float ct, float st);
float padding(nb::ndarray<float> Phantom, nb::ndarray<float> Phantom_pad, int DetSize, int PhantSize, int padXY, int sys);
