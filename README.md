# Code to reproduce figures from "HYPER: Group testing via hypergraph factorization applied to COVID-19"

Preprint available at: https://doi.org/10.1101/2021.02.24.21252394

Online tool available at: http://hyper.covid19-analysis.org

This repo contains Python and Julia code for HYPER:
+ `hyper/` contains example code for generating HYPER designs
+ `simulation/` contains data files (and code) for simulation under a realistic COVID-19 model
+ `simulation,figs/` generates figures from the outputs in `simulation`
+ `simulation,resource,analysis/` analyzes the effectiveness of various methods under resource constraints (using the outputs from `simulation`)
+ `stat,model,figs/` generates the figures to accompany the analysis under a common (non-COVID) statistical model
+ `supp,*/` contain code for associated supplementary figures

See subdirectories for more information.

*Note:* The Julia scripts are in fact also `Pluto.jl` notebooks (https://github.com/fonsp/Pluto.jl).

## System requirements

### Hardware Requirements

All codes can be run on a standard computer and have been tested on a 2019 Macbook Pro with the following specifications:
+ 2.4 GHz Quad-Core Intel Core i5
+ 16 GB RAM

The simulations in `simulation/` can be run in parallel so benefit greatly from running on a cluster;
PBS scripts are provided in the `simulation/` directory.
Similarly for some of the supplementary simulations.
See subdirectories for more information.

### Software Requirements

The code has been run on:
+ MacOS 11.6.4
+ Red Hat Enterprise Linux Server release 7.8

The Python code uses the Anaconda environment described here:
https://github.com/cleary-lab/covid19-group-tests/blob/master/code/group_test_simulations/README.md

The Julia code was run on Julia 1.7.1, which can be installed from here:
https://julialang.org/downloads/

The MATLAB code was run using the Python-Matlab connection (https://www.mathworks.com/help/matlab/matlab-engine-for-python.html)
with MATLAB R2019a.

## Installation

No further installation is required beyond the Anaconda (Python) environment,
Julia and MATLAB (with the Python-Matlab connection).

The Julia codes use packages installed via the Julia package manager (https://julialang.github.io/Pkg.jl/v1/)
and specified by `JuliaProject.toml` and `JuliaManifest.toml` files.
We give instructions where needed.
These installations can typically complete within 5-10 minutes.

## Demo / Instructions for use

Varies from directory to directory; see `README.md` files within.
