#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "PMTsingleEvent.hpp"
#include "helper_lib.hpp"

namespace py = pybind11;

PYBIND11_MODULE(fitEvent, m) {
    m.doc() = "Python bindings for single-event PMT reconstruction";

    // Bind della struct C++
    py::class_<SinglePointResult>(m, "SinglePointResult")
        .def(py::init<>())
        .def_readwrite("L_mean", &SinglePointResult::L_mean)
        .def_readwrite("L_std",  &SinglePointResult::L_std)
        .def_readwrite("x_mean", &SinglePointResult::x_mean);
        .def_readwrite("x_std", &SinglePointResult::x_std);
        .def_readwrite("y_mean", &SinglePointResult::y_mean);
        .def_readwrite("y_std", &SinglePointResult::y_std);

    // Funzione con parametri opzionali
    m.def("run_single_event",
        [](const std::vector<double>& q_values,
           const std::vector<double>& calib,
           int NIterPrerun,
           int NIter,
           bool second_round) 
        {
            if(q_values.size() != 4 || calib.size() != 4)
                throw std::runtime_error("run_single_event expects two arrays of 4 doubles");

            double q[4] = { q_values[0], q_values[1], q_values[2], q_values[3] };
            double c[4] = { calib[0], calib[1], calib[2], calib[3] };

            return runSingleEvent(q, c, NIterPrerun, NIter, second_round);
        },
        py::arg("q_values"),
        py::arg("calib"),
        py::arg("NIterPrerun") = 100000,
        py::arg("NIter") = 10000,
        py::arg("second_round") = false,
        "Run single-event PMT reconstruction with optional parameters");
}