import os
import subprocess
import sys

import pybind11
from setuptools import Extension, setup

# Path to local pybind11
local_pybind11 = os.path.join("external", "pybind11", "include")
root_cflags = subprocess.check_output(["root-config", "--cflags"], text=True).split()
root_libs = subprocess.check_output(["root-config", "--libs"], text=True).split()


def find_bat():
    # 1. BAT_PREFIX explicitly provided
    if "BAT_PREFIX" in os.environ:
        pref = os.environ["BAT_PREFIX"]
        return (
            os.path.join(pref, "include"),
            os.path.join(pref, "lib"),
        )

    # 2. Search in LD_LIBRARY_PATH
    ld_paths = os.environ.get("LD_LIBRARY_PATH", "").split(":")
    for path in ld_paths:
        if os.path.exists(os.path.join(path, "libBAT.so")):
            include_guess = os.path.abspath(os.path.join(path, "..", "include"))
            return include_guess, path

    # 3. Fallback scan in common prefixes
    common_prefixes = [
        "/usr/local",
        "/usr",
        "/opt",
    ]
    for base in common_prefixes:
        libpath = os.path.join(base, "lib")
        if os.path.exists(os.path.join(libpath, "libBAT.so")):
            inc = os.path.join(base, "include")
            return inc, libpath

    # 4. If not found, error out
    print(
        "ERROR: Cannot find libBAT.so. "
        "Please set BAT_PREFIX or update LD_LIBRARY_PATH.",
        file=sys.stderr,
    )
    sys.exit(1)


BAT_INCLUDE, BAT_LIB = find_bat()

include_dirs = [local_pybind11, "src/", pybind11.get_include(), BAT_INCLUDE]  # fallback

ext_modules = [
    Extension(
        "fitEvent",
        sources=[
            "python_bindings/bindings.cpp",
            "src/PMT_association.cpp",
            "src/PMT_calibration.cpp",
            "src/PMT_FindAlpha.cpp",
            "src/PMT_singleEvent.cpp",
            "src/helper_lib.cpp",
        ],
        include_dirs=include_dirs,
        library_dirs=[BAT_LIB],
        libraries=["BAT"],  # linker automatically find libBAT.so
        extra_compile_args=["-std=c++11", "-O3", "-fPIC", "-fopenmp"] + root_cflags,
        extra_link_args=["-fopenmp"] + root_libs,
    )
]

setup(
    name="fitEvent",
    version="0.1.0",
    ext_modules=ext_modules,
)
