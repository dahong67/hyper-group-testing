# Resource analysis

Julia code in this directory uses the outputs in `../simulation/`
to analyze the effectiveness of the various methods under resource constraints.

To make sure the needed Julia packages are installed, run (from this directory):
```bash
julia --project=@. -E 'using Pkg; Pkg.instantiate(); Pkg.status()'
```
This *instantiates* the Julia environment: https://julialang.github.io/Pkg.jl/v1/environments/#Using-someone-else's-project

To generate Fig. 3 and Supplementary Fig. 10, run `figscript.jl`:
```bash
julia figscript.jl
```
This may take a few minutes on a laptop.
The expected output figures are:
+ `fig-3.pdf`
+ `fig-3.png`
+ `fig-s10.pdf`
+ `fig-s10.png`

To generate Supplementary Figs. 11-16, run `figscriptcomb.jl`:
```bash
julia figscriptcomb.jl
```
after un/commenting the corresponding line
in the following block of `figscriptcomb.jl`:
```julia
# ╔═╡ b58f8921-cd4c-4740-b17e-fc669a999b84
DAYS, SAVENAME = IdentityRange(53:53), "fig-s11"
# DAYS, SAVENAME = IdentityRange(80:80), "fig-s12"
# DAYS, SAVENAME = IdentityRange(83:83), "fig-s13"
# DAYS, SAVENAME = IdentityRange(84:84), "fig-s14"
# DAYS, SAVENAME = IdentityRange(90:90), "fig-s15"
# DAYS, SAVENAME = IdentityRange(93:93), "fig-s16"
```
Each run may take a few minutes on a laptop.
The expected output figures are:
+ `fig-s11.pdf`
+ `fig-s11.png`
+ ...
+ `fig-s16.pdf`
+ `fig-s16.png`
