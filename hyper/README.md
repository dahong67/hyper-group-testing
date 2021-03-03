# Generation of HYPER design

Julia code in this directory implements construction of HYPER designs
(see also our online tool: http://hyper.covid19-analysis.org).

To make sure the needed Julia packages are installed, run (from this directory):

```bash
julia --project=@. -E 'using Pkg; Pkg.instantiate(); Pkg.status()'
```

This *instantiates* the Julia environment: https://julialang.github.io/Pkg.jl/v1/environments/#Using-someone-else's-project

To generate a design with `n=10` individuals, `m=6` pools and `q=3` splits, run (from the terminal):

```bash
julia --project=@. -E 'include("gendesign.jl"); A = graphdesign(10,6,3);'
```

This typically takes a few seconds to run.

The output can be saved to be read in by numpy using

```bash
julia --project=@. -E 'include("gendesign.jl"); A = graphdesign(10,6,3); using NPZ; npzwrite("outs.npy",A)'
```

## Description of files

+ `gendesign.jl`: Julia script (in `Pluto.jl` notebook form) implementing the HYPER constructions
+ `outs.npy`: output from example above
+ `JuliaManifest.toml` & `JuliaProject.toml`: Julia environment files specifying the packages used.
+ `PGField/`: utility package for Beth's combinatorial construction