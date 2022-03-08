# Comparison with Reed-Solomon Kautz-Singleton (RS-KS) code-based designs

Code in this directory carries out simulations
for Reed-Solomon Kautz-Singleton (RS-KS) code-based designs
and compares with HYPER.

The `designs/` directory has the main data files;
subdirectories of `designs/` correspond to RS-KS designs
with `f=2` and `f=3`, respectively.
For example, `designs/rskscode,f2` contains the following files:
+ `Design.n-96_m-22_q-2.npy`: file defining the pooling design (input to simulation)
+ `Eff_avg.n-96_m-22_q-2.npy`: output file containing the average efficiency gains
+ `Recall_combined.n-96_m-22_q-2.npy`: output file containing the average sensitivity

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
2. Run `rskscode.pbs` cluster script.
This runs the simulation in `../simulation/covid19-group-tests/code/group_test_simulations/test_on_simulated_population.py` for the designs listed in `rskscode.txt`.

## Steps to create the figures

Julia code in this directory carries out the resource analyses and creates the figures.

To make sure the needed Julia packages are installed, run (from this directory):
```bash
julia --project=@. -E 'using Pkg; Pkg.instantiate(); Pkg.status()'
```
This *instantiates* the Julia environment: https://julialang.github.io/Pkg.jl/v1/environments/#Using-someone-else's-project

To generate Supplementary Figs. 6c-d, run `figscript.jl`:
```bash
julia figscript.jl
```
This may take a couple minutes on a laptop.
The expected output figures are:
+ `fig-s6cd.pdf`
+ `fig-s6cd.png`

To carry out the resource analysis and generate Supplementary Figs. 19-20, run `figscript-res.jl`:
```bash
julia figscript-res.jl
```
after un/commenting the corresponding line
in the following block of `figscript-res.jl`:
```julia
# ╔═╡ b58f8921-cd4c-4740-b17e-fc669a999b84
DAYS, SELECTED, SAVENAME =
	IdentityRange(53:53), CartesianIndex.([(6,1); (8,2);;]), "fig-s19"
	# IdentityRange(80:80), CartesianIndex.([(4,1); (6,2);;]), "fig-s20"
```
Each run may take a few minutes on a laptop.
The expected output figures are:
+ `fig-s19.pdf`
+ `fig-s19.png`
+ `fig-s20.pdf`
+ `fig-s20.png`
