# Figure comparing bounds

Julia code in this directory creates a figure
that compares two design problems and their bounds.

To make sure the needed Julia packages are installed, run (from this directory):

```bash
julia --project=@. -E 'using Pkg; Pkg.instantiate(); Pkg.status()'
```

To generate the figure, run `figs.jl` (from the terminal):
```bash
julia --project=@. -E 'include("figs.jl")'
```
The run may take roughly 30 second on a laptop.
The expected output figure is `fig.png`.
