# Simulation under a realistic COVID-19 model

Python code in this directory simulates different methods under a realistic COVID-19 model.

The `designs/` directory has the main data files;
subdirectories of `designs/` correspond to methods.
For example, `designs/array,8x12` contains the following files:
+ `Design.n-96_m-20_q-2.npy`: file defining the pooling design (input to simulation)
+ `Eff_avg.n-96_m-20_q-2.npy`: output file containing the average efficiency gains
+ `Recall_combined.n-96_m-20_q-2.npy`: output file containing the average sensitivity

The `covid19-group-tests` directory is a git submodule
pointing to the simulation framework we use.

## Steps to produce the output files

1. Run `design_sweep.pbs` cluster script.
This runs the simulation in `covid19-group-tests/code/group_test_simulations/test_on_simulated_population.py` for the methods listed in `design_list.txt`.
2. Follow instructions in `pbest/` to run the simulation for P-BEST (simulating for P-BEST required a modification of the code in `covid19-group-tests`).
