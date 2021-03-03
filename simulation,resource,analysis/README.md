# Resource analysis

Python/Julia code in this directory uses the outputs in `../simulation/`
to analyze the effectiveness of the various methods under resource constraints.

To make sure the needed Julia packages are installed, run (from this directory):

```bash
julia --project=@. -E 'using Pkg; Pkg.instantiate(); Pkg.status()'
```

## Steps to produce the output files

1. Run `resource_analysis_hypergraph.py` to carry out the resource analyses used in the paper:
```bash
python -u resource_analysis_hypergraph.py --viral-load-matrix ../simulation/Simulated_populations/seir_viral_loads_swab.viral_loads.npz --resultspath ../simulation/designs/ --savepath . --start-time 40 --end-time 90 --n-swabs-list 12 24 48 96 192 384 768 1536 3072 6144 --m-kits-list 12 24 48 96 192 384 768 1536 3072 6144
python -u resource_analysis_hypergraph.py --viral-load-matrix ../simulation/Simulated_populations/seir_viral_loads_swab.viral_loads.npz --resultspath ../simulation/designs/ --savepath . --start-time 53 --end-time 53 --n-swabs-list 12 24 48 96 192 384 768 1536 3072 6144 --m-kits-list 12 24 48 96 192 384 768 1536 3072 6144
python -u resource_analysis_hypergraph.py --viral-load-matrix ../simulation/Simulated_populations/seir_viral_loads_swab.viral_loads.npz --resultspath ../simulation/designs/ --savepath . --start-time 80 --end-time 80 --n-swabs-list 12 24 48 96 192 384 768 1536 3072 6144 --m-kits-list 12 24 48 96 192 384 768 1536 3072 6144
python -u resource_analysis_hypergraph.py --viral-load-matrix ../simulation/Simulated_populations/seir_viral_loads_swab.viral_loads.npz --resultspath ../simulation/designs/ --savepath . --start-time 83 --end-time 83 --n-swabs-list 12 24 48 96 192 384 768 1536 3072 6144 --m-kits-list 12 24 48 96 192 384 768 1536 3072 6144
python -u resource_analysis_hypergraph.py --viral-load-matrix ../simulation/Simulated_populations/seir_viral_loads_swab.viral_loads.npz --resultspath ../simulation/designs/ --savepath . --start-time 84 --end-time 84 --n-swabs-list 12 24 48 96 192 384 768 1536 3072 6144 --m-kits-list 12 24 48 96 192 384 768 1536 3072 6144
python -u resource_analysis_hypergraph.py --viral-load-matrix ../simulation/Simulated_populations/seir_viral_loads_swab.viral_loads.npz --resultspath ../simulation/designs/ --savepath . --start-time 90 --end-time 90 --n-swabs-list 12 24 48 96 192 384 768 1536 3072 6144 --m-kits-list 12 24 48 96 192 384 768 1536 3072 6144
python -u resource_analysis_hypergraph.py --viral-load-matrix ../simulation/Simulated_populations/seir_viral_loads_swab.viral_loads.npz --resultspath ../simulation/designs/ --savepath . --start-time 93 --end-time 93 --n-swabs-list 12 24 48 96 192 384 768 1536 3072 6144 --m-kits-list 12 24 48 96 192 384 768 1536 3072 6144
```
This produces the files `summary.resource_t0-40_t1-90.csv`, etc.

2. Run `figs.jl` to generate the figures:
```bash
julia --project=@. -E 'include("figs.jl")'
```
To generate the figures for different sets of days,
un/comment lines in the following block of `figs.jl`:
```
# ╔═╡ e88985de-2f4e-11eb-32cf-1f7a851792bb
# prevalence ≈ 0.033808% - 2.459128%
t0, t1, q1pos, q2pos, q3pos, pbestpos = 40, 90, 8.5, 7, 4, missing

# prevalence ≈ 0.103552%
# t0, t1, q1pos, q2pos, q3pos, pbestpos = 53, 53, 8, 5.5, 3, missing

# prevalence ≈ 1.055872%
# t0, t1, q1pos, q2pos, q3pos, pbestpos = 80, 80, 8.5, 7, 4, missing

# prevalence ≈ 1.362224%
# t0, t1, q1pos, q2pos, q3pos, pbestpos = 83, 83, 8.5, 1, 4, missing

# prevalence ≈ 1.484104%
# t0, t1, q1pos, q2pos, q3pos, pbestpos = 84, 84, 8.5, 1.5, missing, 5

# prevalence ≈ 2.459128%
# t0, t1, q1pos, q2pos, q3pos, pbestpos = 90, 90, 9, 1.5, missing, 5

# prevalence ≈ 3.149712%
# t0, t1, q1pos, q2pos, q3pos, pbestpos = 93, 93, 9, 1.5, 5, missing
```
