# Study of impact of imbalance

Code in this directory studies the impact of imbalance on variability
under the COVID-19 model.

The `designs/` directory has the main data files;
subdirectories of `designs/` correspond to methods.
For example, `designs/consecutive/` includes the following files:
+ `Design.n-96_m-16_q-2.npy`: file defining the pooling design (input to simulation)
+ `num_stage_two.n-96_m-16_q-2.t-80.numpos-1.numiters-100000.randseed-1.npy`: output file containing the average number of stage two tests used as a function of where positive individuals are placed
+ `sens_avg_perm.n-96_m-16_q-2.t-80.numpos-1.numiters-100000.randseed-1.npy`: output file containing the average sensitivity as a function of where positive individuals are placed

The `double,pooling.seed-*/` and `rand,assign.seed-*/` subdirectories
don't contain `Design.*` files since those designs
are randomly generated in the simulation code.
Likewise the `hypergraph/`, `balarray/` and `rskscode,f2/` subdirectories
don't contain `Design.*` files since those designs
are available in `../simulation/designs/`, `../supp,comp,balarray/designs/` and `../supp,comp,rskscode/designs/`, respectively.

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
2. Run `design_sweep.pbs` cluster script.
This runs the simulation in `sample_test_on_simulated_population.py`
for the methods listed in `design_list.txt`.
3. Run `rand_assign.pbs` cluster script.
This runs the simulation in `sample_test_on_simulated_population.py`
for the random assignment design with the settings listed in `rand_list.txt`.
4. Run `double_pooling.pbs` cluster script.
This runs the simulation in `sample_test_on_simulated_population.py`
for double-pooling with the settings listed in `rand_list.txt`.

## Steps to create the figure

Julia code in this directory creates the accompanying figure.

To make sure the needed Julia packages are installed, run (from this directory):
```bash
julia --project=@. -E 'using Pkg; Pkg.instantiate(); Pkg.status()'
```
This *instantiates* the Julia environment: https://julialang.github.io/Pkg.jl/v1/environments/#Using-someone-else's-project

To generate the figure, run `figscript.jl`:
```bash
julia figscript.jl
```
This may take a couple minutes on a laptop.
The expected output figures are:
+ `fig-s7.pdf`
+ `fig-s7.png`
