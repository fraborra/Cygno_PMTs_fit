# CYGNO PMT fits
[![License](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
![Language](https://img.shields.io/badge/language-C%2B%2B%20%2F%20Python-brightgreen.svg)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://pre-commit.com/)

## Overview

**CYGNO PMT Fits** is a C++ library developed for Bayesian reconstruction of cluster position and deposited energy in the LIME detector, using charge integrals from the photomultiplier tubes (PMTs).
The project relies on **BAT (Bayesian Analysis Toolkit)** as the fitting engine and also exposes **Python bindings** offering a single-point fitting interface.

The main executable is configured via a **configuration file** (`.conf`), which selects the fitting mode and controls input/output and runtime options.

## Requirements

### Mandatory
- A working **C++ compiler** compatible with BAT (e.g. GCC ≥ 4.8.x, C++11)
- **BAT – Bayesian Analysis Toolkit** must be installed and discoverable by the compiler and linker.

Please refer to the official BAT website for installation and documentation: <https://bat.mpp.mpg.de>

### Optional
- **Python >= 3.8** (if you want to use the Python interface)
- **pre-commit** (recommended for contributors)

---

## Installation

1. **Clone the repository:**
    ```bash
    git clone https://github.com/fraborra/Cygno_PMTs_fit.git
    cd Cygno_PMTs_fit
    ```

2. **Ensure BAT is installed:**

    Make sure BAT is installed on your system and that the compiler and linker can find:
    - BAT headers
    - BAT libraries

    Depending on your environment, you may need to set environment variables like:
    ```bash
    export LD_LIBRARY_PATH=/path/to/bat/lib:$LD_LIBRARY_PATH
    ```
    or use the tools provided by BAT (`bat-config`, etc.) when configuring the Makefile.

3. **Configure the Makefile**
    In the provided `Makefile`, ensure that:
    - BAT include paths are added to the compiler (e.g. `-I/path/to/bat/include/`)
    - BAT libraries are linked (e.g. `-L/path/to/bat/lib -lBAT` plus any additional dependencies)

    Adjust any other paths and flags according to your environment.

4. **Compile the project:**
    ```bash
    make
    ```
    This should produce the main executable, called:
    ```bash
    ./bin/runfit
    ```
---

## Usage
The main ececutable is invoked as:

```bash
./bin/runfit config/configuration.conf
```
Configuration files (`.conf`) define:
- the operating mode (`mode`) `['association', 'PMTcalibration', 'PMTfindalpha']`
- the input file
- the output file (only in `association` mode)
- fitting and output options

Example configuration files are provided in the config/ directory (e.g. `association.conf`, `calibration.conf`).

### Key configuration parameters

- `mode`: to specify the mode of the program (`['association', 'PMTcalibration', 'PMTfindalpha']`).
- `input_file`: path to the input data file.
- `output_file`: path to the output file where results are saved (used in `association` mode).
- `start_ind`, `end_ind`: define the data range (by index) to be processed:
    - `start_ind`: starting line index (0-based; `0` = first line)
    - `end_ind`: ending line index (`-1` = process until end of file)
- `plot`: boolean option to save diagnostic plots for the MCMC parameters.
- `write_log`: boolean option to save the log files.
- `write_chains`: boolean option to save MCMC chains for the model parameters:
    - In `association` mode, typically used for detailed diagnostics (forces a max of 5 events)
    - In `PMTcalibration` and `PMTfindalpha` modes chains are always written.
- `print_summary`: boolean option to print the summary of the MCMC integration summary of the MCMC integration to stdout:
    - Set to `false` in `association` to reduce runtime overhead.
    - always `true` in `PMTcalibration` and `PMTfindalpha`.
- `c1`, `c2`, `c3`, `c4`: calibration parameters for PMT 1–4, fitted from `PMTcalibration` mode.

---

## Input File Format
1. `association` mode

The input file is a TAB-separated text file. Each line corresponds to a single PMT cluster and must contain the following fields:

- `run` - run number
- `event` - event number
- `trigger` - trigger number
- `peak index` - index of the peak in the waveform
- `L1` - integral of **PMT 1** in **nC**
- `L2` - integral of **PMT 2** in **nC**
- `L3` - integral of **PMT 3** in **nC**
- `L4` - integral of **PMT 4** in **nC**

All fields must be separated by a tab character (`\t`)

An example file is provided in `data/golden_input.txt`

1. `PMTcalibration` mode:
The input format is identical to association, with three additional fields:
- `x` - x position of the cluster on the GEM plane, in **cm**.
- `y` - y position of the cluster on the GEM plane, in **cm**.
- `sc_integral`- camera integral of the cluster.

Again, fields are TAB-separated and follow the same ordering.

---

## Output File Format

1. `association` mode:
The output file will contain a series of lines, each representing processed data with the following fields:

- `run` - run number
- `event` - event number
- `trigger` - trigger number
- `peak index` - index of the peak in the waveform
- `L` – fitted parameter related to the light
- `Lstd`– standard deviation (uncertainty) of `L`
- `x`- reconstructed x position (cm).
- `xstd`– standard deviation (uncertainty) of `x`
- `y`- reconstructed y position (cm).
- `ystd`– standard deviation (uncertainty) of `y`

If the fit fails to converge for a given input line, all fitted quantities (`L`, `Lstd`, `x`, `xstd`, `y`,`ystd`) are set to `-1`.

An example output is provided in `data/golden_out.txt`

2. **PMTcalibration:**

In `PMTcalibration` and `PMTfindAlpha` modes, the program does not produce a simple line-by-line text output.
Instead, it writes the MCMC chains for the calibration parameters to file(s). These chains contain the sampled values of the calibration parameters over the MCMC iterations.

An example of how to read and analyze these chains is provided in: [How_to_read_BAT_output.ipynb](./notebooks/How_read_BAT_output.ipynb)

---

## Code style and pre-commit
This repository uses pre-commit to enforce consistent code style:
- **C++**: formatted with `clang-format`
- **Python**: formatted with `black`

To enable pre-commit locally:

```bash
pip install pre-commit
pre-commit install
```
After this, every git commit will automatically run the configured checks and formatters.
You can also run all hooks manually:

```bash
pre-commit run --all-files
```

---

## Citing This Work

This software was used to produce the results published in:

[**F. D. Amaro _et al_., “Bayesian network 3D event reconstruction in the Cygno optical TPC for dark matter direct detection,” Eur. Phys. J. C 85 (2025) 1261.**](https://doi.org/10.1140/epjc/s10052-025-14965-6)

If you use this code in academic work, please cite the article above.

---

## Authors and Contact
- **Francesco Borra** <br> Email: francesco.borra@uniroma3.it

---

## Acknowledgments
Special thanks to:
- Stefano Piacentini
- Matteo Folcarelli

for their contributions.

## License

This project is licensed under the GNU General Public License v3 (GPLv3).
See the LICENSE file for details: <https://www.gnu.org/licenses/>.
