#include <pybind11/pybind11.h>
#include "../algorithms/CubicalRipser_2dim/src/cubicalripser_2dim.h"
#include "../algorithms/CubicalRipser_3dim/src/cubicalripser_3dim.h"
//#include "../algorithms/dipha/src/dipha.h"
#include "../algorithms/Perseus/Perseus.h"
#include "../algorithms/Ripser/ripser.h"
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>


namespace py = pybind11;

PYBIND11_MODULE(tda, m) {
  //	Binding for CubicalRipser2D
  py::class_<CubicalRipser2D>(m, "CubicalRipser2D")
		.def(py::init<>())
	  	.def("ComputeBarcode", &CubicalRipser2D::ComputeBarcode)
		;
  py::class_<CubicalRipser3D>(m, "CubicalRipser3D")
		  .def(py::init<>())
		  .def("ComputeBarcode", &CubicalRipser3D::ComputeBarcode)
		  ;
  py::class_<Perseus>(m, "Perseus")
		  .def(py::init<>())
		  .def("ComputeBarcode", &Perseus::ComputeBarcode)
		  ;

  py::class_<Ripser>(m, "Ripser")
		  .def(py::init<>())
		  .def("ComputeBarcode", &Ripser::ComputeBarcode)
		  ;
  //py::class_<Dipha>(m, "Dipha")
  //		.def(py::init<>())
  //		.def("compute", &Dipha::compute)
  //		.def("createPhatFiltration", &Dipha::createPhatFiltration)
  //		.def("createSparseRepresentation", &Dipha::createSparseRepresentation)
  //		;
}
