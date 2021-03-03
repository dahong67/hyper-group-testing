# Simulation for P-BEST

Code in this directory produces the outputs in `../designs/pbest/`.

The modified simulation for P-BEST is implemented
in `test_on_simulated_population.py`.

Notes:
+ It requires a working Matlab-Python connection:
  https://www.mathworks.com/help/matlab/matlab-engine-for-python.html
+ It uses the P-BEST implementation provided online at https://github.com/NoamShental/PBEST loaded here as the git submodule `PBEST`

## Steps to produce the output

1. Run the `pbest.pbs` cluster script
with each desired day passed as job array environment variable `SGE_TASK_ID` (each run simulates a single day).
2. Outputs are placed in `results/`.
3. Run `combine_results.py` to form combined results (across days) and place them in the `../designs/pbest` directory, for example:
```bash
conda activate covid-group-tests
python -u combine_results.py --resultspath ./results --savepath ../designs/pbest/ --start-time 10 --end-time 120
```
