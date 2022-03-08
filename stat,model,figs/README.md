# Analysis under the statistical model

Julia code in this directory runs simulations and creates figures
to accompany the analysis under the general statistical model.

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
+ `fig-s1.pdf`
+ `fig-s1.png`

The script uses a saved cache of the simulation results: `simcache.bson`.
To regenerate these results simply delete this cache file and rerun:
```bash
rm simcache.bson
julia figscript.jl
```
or to have a progress bar, run `figscript.jl` with an appropriate logger:
```bash
julia --project=. -e 'using Logging, TerminalLoggers; with_logger(()->include("figscript.jl"),TerminalLogger());'
```
This may take roughly 6 minutes on a laptop.

Note that `figscript.jl` uses [`StableRNGs.jl`](https://github.com/JuliaRandom/StableRNGs.jl) to make the simulation reproducible.
