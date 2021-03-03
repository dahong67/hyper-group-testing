# Code to reproduce figures from "HYPER: Group testing via hypergraph factorization applied to COVID-19"

Preprint available at: https://www.medrxiv.org/content/10.1101/2021.02.24.21252394v1

Online tool available at: http://hyper.covid19-analysis.org

This repo contains Python and Julia code for HYPER:
+ `hyper/` contains example code for generating HYPER designs
+ `simulation/` contains data files (and code) for simulation under a realistic COVID-19 model
+ `simulation,figs/` generates figures from the outputs in `simulation`
+ `simulation,resource,analysis/` analyzes the effectiveness of various methods under resource constraints (using the outputs from `simulation`)
+ `theoretical,figs/` generates figures to accompany analysis under a theoretical model
See subdirectories for more information.

*Note:* The Julia scripts are in fact `Pluto.jl` notebooks (https://github.com/fonsp/Pluto.jl).
