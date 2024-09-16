# CYGNO PMT fits

## Project Description

CYGNO PMT fits is a project that uses a Bayesian fit to reconstruct the position and energy released by a cluster using the charge integral of PMTs for the CYGNO Collaboration.

## System Requirements

To use this project, you need to have BAT (Bayesian Analysis Toolkit) installed. You can download and install it from the official website: [BAT](https://bat.mpp.mpg.de).

## Installation

1. **Clone the repository:**
    ```bash
    git clone https://github.com/fraborra/Cygno_PMTs_fit.git
    cd Cygno_PMTs_fit
    ```
    
2. **Install dependencies:**
    Make sure BAT is installed and properly configured.

3. **Compile the project:**
    ```bash
    make
    ```
## Usage Instructions

1. **Run the main script:**
    To run the program use the following use the following command in the terminal: 
    ```bash
    ./runfit configuration.conf
    ```
    Where the `.conf` is the configuration file (examples in `calibration.conf` and `association.conf`).

    The options are:
    - `mode`: to specify the mode of the program (between **association** and **PMTcalibration**).
    - `input_file`: name of the input file to use for the program.
    - `output_file`: indicates the output file where the results will be saved for **association**.
    - `start_ind`: defines the starting row of the input file (`0` first row).
    - `end_ind`: defines the ending row of the input file (`-1` until the end).
    - `plot`: option to save some plots for the MCMC of the parameters (use only for 1 integration at time).
    - `write_log`: option to save the log files (use only for 1 integration at time).
    - `write_chains`: option to save MCMC chains of the parameters (use only for 1 integration at time in **association**, always `true` in **PMTcalibration**).
    - `print_summary`: option to print the summary of the MCMC integration on the screen (`false` to reduce compute time, always `true` in **PMTcalibration**).
    - `nPoints`: option to select how many events to integrate at once for the **association** (usually `1`).
    - `c1`: option to set calibratrion parameter for **PMT 1**
    - `c2`: option to set calibratrion parameter for **PMT 2**
    - `c3`: option to set calibratrion parameter for **PMT 3**
    - `c4`: option to set calibratrion parameter for **PMT 4**


## Input File Format
1. **association:**
    The input file should contain a series of lines, each representing a set of data with the following fields:

    - **run**: The run number.
    - **event**: The event number.
    - **trigger**: The trigger number.
    - **peak index**: The index indicating the position of the peak in the waveform.
    - **L1**: The integral of the **PMT 1** must be in **nC**.
    - **L2**: The integral of the **PMT 2** must be in **nC**.
    - **L3**: The integral of the **PMT 3** must be in **nC**.
    - **L4**: The integral of the **PMT 4** must be in **nC**.

    Each line in the input file should have these fields separated by a tab.

    Ensure that each line follows this format to allow the runfit program to correctly parse and process the input data.

    An example can be found in `golden_input.txt`

1. **PMTcalibration:**
    The same as the association but with three more variables:
    - **x**: x position of the cluster in the GEM plane, in **cm**.
    - **y**: y position of the cluster in the GEM plane, in **cm**.
    - **sc_integral**: camera integral of the cluster.

## Output File Format

1. **association:**
The output file will contain a series of lines, each representing processed data with the following fields:

- **run**: The run number.
- **event**: The event number.
- **trigger**: The trigger number.
- **peak index**: The index indicating the position of the peak in the waveform.
- **L**: The value of parameter L.
- **Lstd**: The standard deviation of parameter L.
- **x**: The value of parameter x.
- **xstd**: The standard deviation of parameter x.
- **y**: The value of parameter y.
- **ystd**: The standard deviation of parameter y.

Each line in the output file will have these fields separated by a tab.

Each line in the output file will contain the calculated values of these parameters for the corresponding input data line. If the fit didn't converge for a row then all the parameters values will be set with a **-1**.

An example can be found in `golden_out.txt`

2. **PMTcalibration:**
   
The PMTcalibration output will give only the chains of the computed "calibration parameters". How to read the chains can be found in the example [read_chains.ipynb](./read_chains.ipynb)

=======

## Authors and Contact

- **Name:** Francesco Borra
- **Email:** francesco.borra@uniroma3.it
## Credits

Special thanks Stefano Piacentini and Matteo Folcarelli.
