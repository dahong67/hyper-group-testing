# Comparison with additional designs

Code in this directory carries out some of the simulations from `../simulation/`
for additional designs.

The `designs/` directory has the main data files;
subdirectories of `designs/` correspond to methods.
For example, `designs/array,holes` includes the following files:
+ `Design.n-96_m-20_q-2.npy`: file defining the pooling design (input to simulation)
+ `Eff_avg.n-96_m-20_q-2.npy`: output file containing the average efficiency gains
+ `Recall_combined.n-96_m-20_q-2.npy`: output file containing the average sensitivity

The `double,pooling/` and `rand,assign/` subdirectories
don't contain `Design.*` files since those designs
are randomly generated in the simulation code.

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
2. Run `array_holes.pbs` cluster script.
This runs the simulation in `../simulation/covid19-group-tests/code/group_test_simulations/test_on_simulated_population.py` for the methods listed in `array_holes.txt`.
3. Run `rand_assign.pbs` cluster script.
This runs the simulation in `rand_test_on_simulated_population.py`
for the methods listed in `rand_assign.txt`.
4. Run `double_pooling.pbs` cluster script.
This runs the simulation in `rand_test_on_simulated_population.py`
for the methods listed in `double_pooling.txt`.

## Steps to create the figure

Julia code in this directory creates the accompanying figure.

To make sure the needed Julia packages are installed, run (from this directory):
```bash
julia --project=@. -E 'using Pkg; Pkg.instantiate(); Pkg.status()'
```
This *instantiates* the Julia environment: https://julialang.github.io/Pkg.jl/v1/environments/#Using-someone-else's-project

To generate the figures, run `figscript.jl`:
```bash
julia figscript.jl
```
This may take a couple minutes on a laptop.
The expected output figures are:
+ `fig-s5.pdf`
+ `fig-s5.png`
