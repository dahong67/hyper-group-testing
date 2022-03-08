# Generation of HYPER design

Julia code in this directory implements construction of HYPER designs
(see also our online tool: http://hyper.covid19-analysis.org).

To make sure the needed Julia packages are installed, run (from this directory):
```bash
julia --project=@. -E 'using Pkg; Pkg.instantiate(); Pkg.status()'
```
This *instantiates* the Julia environment: https://julialang.github.io/Pkg.jl/v1/environments/#Using-someone-else's-project

To generate a design with `n=10` individuals, `m=6` pools and `q=3` splits,
run the `gendesign.jl` script:
```bash
julia gendesign.jl
```
This typically takes a few seconds to run and saves the output to `outs.npy` (which can be read in by numpy).

To change the design parameters simply modify the following line of the script:
```julia
A = hyperdesign(10,6,3)
```

## Description of files

+ `gendesign.jl`: Julia scipt (in [`Pluto.jl`](https://github.com/fonsp/Pluto.jl) notebook form) that generates a HYPER design and saves the output to `outs.npy`
+ `outs.npy`: output from example above
+ `JuliaManifest.toml` & `JuliaProject.toml`: Julia environment files specifying the packages (and versions) used.
+ `HyperGen/`: utility package that generates HYPER designs