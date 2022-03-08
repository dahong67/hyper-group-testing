# Comparison with balanced arrays

Code in this directory carries out simulations for balanced arrays
and compares with HYPER.

The `designs/balarray` directory has the main data files.
For example:
+ `Design.n-96_m-20_q-2.npy`: file defining the pooling design (input to simulation)
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
2. Run `balarray.pbs` cluster script.
This runs the simulation in `../simulation/covid19-group-tests/code/group_test_simulations/test_on_simulated_population.py` for the designs listed in `balarray.txt`.

## Steps to create the figures

Julia code in this directory carries out the resource analyses and creates the figures.

To make sure the needed Julia packages are installed, run (from this directory):
```bash
julia --project=@. -E 'using Pkg; Pkg.instantiate(); Pkg.status()'
```
This *instantiates* the Julia environment: https://julialang.github.io/Pkg.jl/v1/environments/#Using-someone-else's-project

To generate Supplementary Figs. 6a-b, run `figscript.jl`:
```bash
julia figscript.jl
```
This may take a couple minutes on a laptop.
The expected output figures are:
+ `fig-s6ab.pdf`
+ `fig-s6ab.png`

To carry out the resource analysis and generate Supplementary Figs. 17-18, run `figscript-res.jl`:
```bash
julia figscript-res.jl
```
after un/commenting the corresponding line
in the following block of `figscript-res.jl`:
```julia
# ╔═╡ b58f8921-cd4c-4740-b17e-fc669a999b84
DAYS, SELECTED, SAVENAME =
	IdentityRange(53:53), CartesianIndex.([(6,1); (8,2);;]), "fig-s17"
	# IdentityRange(80:80), CartesianIndex.([(4,1); (9,5);;]), "fig-s18"
```
Each run may take a few minutes on a laptop.
The expected output figures are:
+ `fig-s17.pdf`
+ `fig-s17.png`
+ `fig-s18.pdf`
+ `fig-s18.png`
