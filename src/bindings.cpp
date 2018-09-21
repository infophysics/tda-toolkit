#include <pybind11/pybind11.h>
#include "../algorithms/CubicalRipser_2dim/src/cubicalripser_2dim.h"
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
}
