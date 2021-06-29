# Analysis under the statistical model

Julia code in this directory creates figures to accompany
the analysis under a common statistical model.

To make sure the needed Julia packages are installed, run (from this directory):

```bash
julia --project=@. -E 'using Pkg; Pkg.instantiate(); Pkg.status()'
```

To generate the figures, run `figs.jl` (from the terminal):
```bash
julia --project=@. -E 'include("figs.jl")'
```
To generate figs a-c vs d-f,
un/comment lines in the following block of `figs.jl`:
```
# ╔═╡ 75a6fdb4-3aa0-11eb-38d5-f3c1c4ade957
opchar, figlabels = (α=0.05,β=0.90), (eff="a",sens="b",spec="c")
# opchar, figlabels = (α=0.05,β=0.80), (eff="d",sens="e",spec="f")
```
Each run may take roughly 5 minutes on a laptop.
The expected output figures for α=0.05,β=0.90 (first set)  are
+ `fig-a.png`
+ `fig-b.png`
+ `fig-c.png`
+ `fig-g.png`
+ `fig-h.png`