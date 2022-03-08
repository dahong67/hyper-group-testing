# Figure comparing array designs with additional HYPER designs

Julia code in this directory creates a figure
that compares array designs with additional HYPER designs.

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
+ `fig-s3.pdf`
+ `fig-s3.png`
