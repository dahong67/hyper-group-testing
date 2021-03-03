# Theoretical analysis figures

Julia code in this directory creates figures to accompany
the analysis under a theoretical model.

To make sure the needed Julia packages are installed, run (from this directory):

```bash
julia --project=@. -E 'using Pkg; Pkg.instantiate(); Pkg.status()'
```

To generate the figures, run `figs.jl` (from the terminal):
```bash
julia --project=@. -E 'include("figs.jl")'
```
To generate figs c-e vs f-h,
un/comment lines in the following block of `figs.jl`:
```
# ╔═╡ 75a6fdb4-3aa0-11eb-38d5-f3c1c4ade957
opchar, figlabels = (α=0.05,β=0.90), (eff="c",sens="d",spec="e")
# opchar, figlabels = (α=0.05,β=0.80), (eff="f",sens="g",spec="h")
```
Each run may take roughly 5 minutes on a laptop.
The expected output figures for α=0.05,β=0.90 (first set)  are
+ `fig-a.png`
+ `fig-b.png`
+ `fig-c.png`
+ `fig-d.png`
+ `fig-e.png`