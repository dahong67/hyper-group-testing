# Simulation under a realistic COVID-19 model

Python code in this directory simulates different group testing methods under a realistic COVID-19 model.

The `designs/` directory has the main data files;
subdirectories of `designs/` correspond to methods.
For example, `designs/array,8x12` contains the following files:
+ `Design.n-96_m-20_q-2.npy`: file defining the pooling design (input to simulation)
+ `Eff_avg.n-96_m-20_q-2.npy`: output file containing the average efficiency gains
+ `Recall_combined.n-96_m-20_q-2.npy`: output file containing the average sensitivity

The `covid19-group-tests` directory is a git submodule
pointing to the simulation framework we use.

## Steps to reproduce the output files

1. Generate (simulated) population files: https://github.com/cleary-lab/covid19-group-tests/tree/master/code/viral_kinetics  
The simulation scripts here expect a directory `./Simulated_populations/`
with the following files:
```
Simulated_populations
├── seir_viral_loads_swab.peak_times.npy
├── seir_viral_loads_swab.timepoints.npy
├── seir_viral_loads_swab.viral_loads.npz
└── used_pars_swab.csv
```
2. Run `design_sweep.pbs` cluster script.
This runs the simulation in `covid19-group-tests/code/group_test_simulations/test_on_simulated_population.py` for the methods listed in `design_list.txt`.
3. Follow instructions in `pbest/` to run the simulation for P-BEST (simulating for P-BEST required a modification of the code in `covid19-group-tests`).
