# Simulation figures

Julia code in this directory uses the outputs in `../simulation/`
to create figures of the simulation traces for the paper.

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
+ `fig-2.pdf`
+ `fig-2.png`
