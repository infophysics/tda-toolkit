#include <pybind11/pybind11.h>
#include "../algorithms/CubicalRipser_2dim/src/cubicalripser_2dim.h"
#include "../algorithms/CubicalRipser_3dim/src/cubicalripser_3dim.h"
//#include "../algorithms/dipha/src/dipha.h"
#include "../algorithms/Perseus/Perseus.h"
#include "../algorithms/Ripser/ripser.h"
#include "../algorithms/bottleneck/CBottleneckDistance.h"
#include "Filter.h"
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>


namespace py = pybind11;

PYBIND11_MODULE(tda, m) {
  //	Binding for CubicalRipser2D
  py::class_<CubicalRipser2D>(m, "CubicalRipser2D")
		.def(py::init<>())
	  	.def("ComputeBarcode", &CubicalRipser2D::ComputeBarcode)
		.def("getBarcode", &CubicalRipser2D::getBarcode)
		;
  py::class_<CubicalRipser3D>(m, "CubicalRipser3D")
		  .def(py::init<>())
		  .def("ComputeBarcode", &CubicalRipser3D::ComputeBarcode)
		  ;
  py::class_<Perseus>(m, "Perseus")
		  .def(py::init<>())
		  .def("ComputeBarcode", &Perseus::ComputeBarcode)
		  .def("getBarcode", &Perseus::getBarcode)
		  ;

  py::class_<Ripser>(m, "Ripser")
		  .def(py::init<>())
		  .def("ComputeBarcode", &Ripser::ComputeBarcode)
		  .def("getBarcode", &Ripser::getBarcode)
		  .def("saveBarcodeToFile", &Ripser::saveBarcodeToFile)
		  ;
  
  py::class_<Filter2D>(m, "Filter2D")
		  .def(py::init<>())
		  .def("loadBinaryFromFile", &Filter2D::loadBinaryFromFile)
		  //	Various filterings
		  //	Binary filterings
		  .def("filterBinaryL1", &Filter2D::filterBinaryL1)
		  .def("filterBinaryL2", &Filter2D::filterBinaryL2)
		  .def("filterBinaryLinf", &Filter2D::filterBinaryLinf)
		  //	Save filtration
		  .def("saveBinaryFiltration", &Filter2D::saveBinaryFiltration)
  	  	  .def("filter3StateAsBinary", &Filter2D::filter3StateAsBinary)  
		  ;
  
  py::class_<Generator>(m, "Generator")
		  .def(py::init<>())
		  .def_readwrite("birth", &Generator::birth)
		  .def_readwrite("death", &Generator::death)
		  ;
  
  py::class_<CBottleneckDistance>(m, "BottleneckDistance")
		  .def(py::init<>())
		  .def("Distance", (double (CBottleneckDistance::*)(std::vector<Generator>,std::vector<Generator>,double)) &CBottleneckDistance::Distance)
		  .def("Distance", (double (CBottleneckDistance::*)(const char*,const char*,double)) &CBottleneckDistance::Distance)
		  ;
		  
  //py::class_<Dipha>(m, "Dipha")
  //		.def(py::init<>())
  //		.def("compute", &Dipha::compute)
  //		.def("createPhatFiltration", &Dipha::createPhatFiltration)
  //		.def("createSparseRepresentation", &Dipha::createSparseRepresentation)
  //		;
}
