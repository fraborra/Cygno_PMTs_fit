#include "PMT_singleEvent.hpp"
#include "helper_lib.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(fitEvent, m) {
  m.doc() = "Python bindings for single-event PMT reconstruction";

  // Bind of the C++ struct SinglePointResult
  py::class_<SinglePointResult>(m, "SinglePointResult")
      .def(py::init<>())
      .def_readwrite("L_mean", &SinglePointResult::L_mean)
      .def_readwrite("L_std", &SinglePointResult::L_std)
      .def_readwrite("x_mean", &SinglePointResult::x_mean)
      .def_readwrite("x_std", &SinglePointResult::x_std)
      .def_readwrite("y_mean", &SinglePointResult::y_mean)
      .def_readwrite("y_std", &SinglePointResult::y_std)
      .def("__repr__", [](const SinglePointResult &s) {
        return "L_mean = " + std::to_string(s.L_mean) + "\n" +
               "L_std = " + std::to_string(s.L_std) + "\n" +
               "x_mean = " + std::to_string(s.x_mean) + "\n" +
               "x_std = " + std::to_string(s.x_std) + "\n" +
               "y_mean = " + std::to_string(s.y_mean) + "\n" +
               "y_std = " + std::to_string(s.y_std);
      });

  // bind the runSingleEvent function
  m.def(
      "runSingleEvent",
      [](const std::vector<double> &q_values, const std::vector<double> &calib,
         int NIterPrerun, int NIter, bool second_round) {
        if (q_values.size() != 4 || calib.size() != 4)
          throw std::runtime_error(
              "runSingleEvent expects two arrays of 4 doubles");

        double q[4] = {q_values[0], q_values[1], q_values[2], q_values[3]};
        double c[4] = {calib[0], calib[1], calib[2], calib[3]};

        return runSingleEvent(q, c, NIterPrerun, NIter, second_round);
      },
      py::arg("q_values"), py::arg("calib"), py::arg("NIterPrerun") = 100000,
      py::arg("NIter") = 10000, py::arg("second_round") = false,
      "Run single-event PMT reconstruction with optional parameters");
}
