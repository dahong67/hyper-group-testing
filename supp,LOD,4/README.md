# Simulation under a realistic COVID-19 model with reduced limit of detection (LOD)

Code in this directory repeats some of the simulations from `../simulation/`
with a reduced LOD.

The `designs/` directory has the main output data files;
subdirectories of `designs/` correspond to methods.
For example, `designs/array,8x12` contains the following files:
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
2. Run `array_hyper.pbs` cluster script.
This runs the simulation in `../simulation/covid19-group-tests/code/group_test_simulations/test_on_simulated_population.py` for the methods listed in `array_hyper.txt`.
3. Run `pbest.pbs` cluster script.
This runs the simulation in `pbest_test_on_simulated_population.py` and saves outputs to `pbest,results/`.
See the instructions in `../simulation/pbest/README.md` for how to setup (e.g., the Matlab-Python connection).
4. Run `../simulation/pbest/combine_results.py` to form combined results (across days) and place them in the `designs/pbest/` directory, for example:
```bash
conda activate covid-group-tests
python -u ../simulation/pbest/combine_results.py --resultspath ./pbest,results --savepath ./designs/pbest/ --start-time 10 --end-time 120
```
5. Run `individual_testing.py` to run the simulation for individual testing:
```bash
conda activate covid-group-tests
python -u individual_testing.py --viral-load-matrix ../simulation/Simulated_populations/seir_viral_loads_swab.viral_loads.npz --savepath designs/individual --start-time 0 --end-time 356 --LOD 4 | tee designs/individual/log.txt
```

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
+ `fig-s8.pdf`
+ `fig-s8.png`
