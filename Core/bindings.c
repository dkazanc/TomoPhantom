#include "TomoP2DModelSino_core.h"
#include "TomoP2DModel_core.h"
#include "TomoP2DSinoNum_core.h"
#include "TomoP3DModelSino_core.h"
#include "TomoP3DModel_core.h"
#include "utils.h"

#include <nanobind/nanobind.h>

NB_MODULE(nanobind, m) {
  m.def("TomoP2DModel_core", &TomoP2DModel_core);
  m.def("TomoP2DObject_core", &TomoP2DObject_core);

  m.def("TomoP2DModelSino_core", &TomoP2DModelSino_core);
  m.def("TomoP2DObjectSino_core", &TomoP2DObjectSino_core);

  m.def("TomoP2DSinoNum_core", &TomoP2DSinoNum_core);
  m.def("BilinearInterpolation", &BilinearInterpolation);
  m.def("padding", &padding);

  m.def("TomoP3DModel_core", &TomoP3DModel_core);
  m.def("TomoP3DObject_core", &TomoP3DObject_core);

  m.def("TomoP3DModelSino_core", &TomoP3DModelSino_core);
  m.def("TomoP3DObjectSino_core", &TomoP3DObjectSino_core);

  m.def("checkParams2D", &checkParams2D);
  m.def("checkParams3D", &checkParams3D);
  m.def("matrot3", &matrot3);
  m.def("matvet3", &matvet3);
  m.def("matmat3", &matmat3);
}
