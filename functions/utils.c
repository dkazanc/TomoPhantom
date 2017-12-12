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

#include "utils.h"

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

float parameters_check3D(float C0, float x0, float y0, float z0, float a, float b, float c, float phi_rot)
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
