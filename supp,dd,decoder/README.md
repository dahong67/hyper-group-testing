# Comparison with HYPER using DD decoding

Code in this directory carries out simulations of HYPER variants using DD decoding
and compares with the default conservative decoder.

The `designs/` directory has the main data files;
the subdirectories of `designs/` correspond to the two variants.
For example, `designs/hypergraph,dd,decode` contains the following files:
+ `Eff_avg.n-96_m-20_q-2.npy`: output file containing the average efficiency gains
+ `Recall_combined.n-96_m-20_q-2.npy`: output file containing the average sensitivity

## Steps to reproduce the output files

1. Generate (simulated) population files: https://github.com/cleary-lab/covid19-group-tests/tree/master/code/viral_kinetics  
The simulation scripts here expect a directory `../simulation/Simulated_populations/`
with the following files:
```
Simulated_populations
├── seir_viral_loads_swab.peak_times.npy
├── seir_viral_loads_swab.timepoints.npy
├── seir_viral_loads_swab.viral_loads.npz
└── used_pars_swab.csv
```
2. Run `design_sweep_dd_decode.pbs` cluster script.
This runs the simulation in `test_on_simulated_population_dd_decode.py` for the designs listed in `design_list.txt`.
3. Run `design_sweep_dd_skip.pbs` cluster script.
This runs the simulation in `test_on_simulated_population_dd_skip.py` for the designs listed in `design_list.txt`.

## Steps to create the figures

Julia code in this directory carries out the resource analyses and creates the figures.

To make sure the needed Julia packages are installed, run (from this directory):
```bash
julia --project=@. -E 'using Pkg; Pkg.instantiate(); Pkg.status()'
```
This *instantiates* the Julia environment: https://julialang.github.io/Pkg.jl/v1/environments/#Using-someone-else's-project

To generate Supplementary Fig. 21, run `figscript.jl`:
```bash
julia figscript.jl
```
This may take a few minutes on a laptop.
The expected output figures are:
+ `fig-s21.pdf`
+ `fig-s21.png`
