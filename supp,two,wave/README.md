# Simulation under a realistic COVID-19 model with a two-wave epidemic

Code in this directory repeats some of the simulations from `../simulation/`
with a two-wave epidemic.

The `designs/` directory has the main output data files;
subdirectories of `designs/` correspond to methods.
For example, `designs/array,8x12` contains the following files:
+ `Eff_avg.n-96_m-20_q-2.npy`: output file containing the average efficiency gains
+ `Recall_combined.n-96_m-20_q-2.npy`: output file containing the average sensitivity

## Steps to reproduce the output files

1. Generate (simulated) population files: https://github.com/cleary-lab/covid19-group-tests/tree/master/code/viral_kinetics  
The simulation scripts here expect a directory `Simulated_populations_two_wave/`
with the following files:
```
Simulated_populations_two_wave
├── seir_viral_loads_swab_switch_SEIR_new_1.peak_times.npy
├── seir_viral_loads_swab_switch_SEIR_new_1.timepoints.npy
├── seir_viral_loads_swab_switch_SEIR_new_1.viral_loads.npz
└── used_pars_swab_switch_SEIR_new_1.csv
```
2. Run `array_hyper.pbs` cluster script.
This runs the simulation in `../simulation/covid19-group-tests/code/group_test_simulations/test_on_simulated_population.py` for the methods listed in `array_hyper.txt`.
3. Run `pbest.pbs` cluster script.
This runs the simulation in `pbest_test_on_simulated_population.py` and saves outputs to `pbest,results/`.
See the instructions in `../simulation/pbest/README.md` for how to setup (e.g., the Matlab-Python connection).
4. Run `../simulation/pbest/combine_results.py` to form combined results (across days) and place them in the `designs/pbest/` directory, for example:
```bash
conda activate covid-group-tests
python -u ../simulation/pbest/combine_results.py --resultspath ./pbest,results --savepath ./designs/pbest/ --start-time 10 --end-time 200
```
5. Run `individual_testing.py` to run the simulation for individual testing:
```bash
conda activate covid-group-tests
python -u individual_testing.py --viral-load-matrix Simulated_populations_two_wave/seir_viral_loads_swab_switch_SEIR_new_1.viral_loads.npz --savepath designs/individual --start-time 0 --end-time 200 | tee designs/individual/log.txt
```
6. Run `posviralloads/extractscript.py` to extract the distribution of positive viral loads for days 82, 104 and 173:
```bash
conda activate covid-group-tests
python -u extractscript.py --viral-load-matrix ../Simulated_populations_two_wave/seir_viral_loads_swab_switch_SEIR_new_1.viral_loads.npz --savepath . --at-time 82
python -u extractscript.py --viral-load-matrix ../Simulated_populations_two_wave/seir_viral_loads_swab_switch_SEIR_new_1.viral_loads.npz --savepath . --at-time 104
python -u extractscript.py --viral-load-matrix ../Simulated_populations_two_wave/seir_viral_loads_swab_switch_SEIR_new_1.viral_loads.npz --savepath . --at-time 173
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
+ `fig-s9.pdf`
+ `fig-s9.png`
