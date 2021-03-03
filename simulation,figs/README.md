# Simulation figures

Julia code in this directory uses the outputs in `../simulation/`
to create figures of the simulation traces for the paper.

To make sure the needed Julia packages are installed, run (from this directory):

```bash
julia --project=@. -E 'using Pkg; Pkg.instantiate(); Pkg.status()'
```

To generate the figures, run `figs.jl` (from the terminal):
```bash
julia --project=@. -E 'include("figs.jl")'
```
