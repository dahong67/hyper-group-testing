### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ 95eaf04a-08e5-11eb-043f-27646f1bf207
begin
	import Pkg; Pkg.activate(@__DIR__)
	using CairoMakie, Dictionaries, NPZ, Statistics
end

# ╔═╡ 849adbca-08e5-11eb-1956-cd02f2a6be7a
md"""
# Figure: simulation traces
"""

# ╔═╡ 210d331a-09c3-11eb-06ed-e1f4ddc023fc
md"""
## Load data
"""

# ╔═╡ 66dd9e7e-f948-4457-8bd6-45345be6abb5
DESDIR = relpath(joinpath(@__DIR__,"designs"))

# ╔═╡ dcfbfb7b-3d7b-42f5-9635-cbdc92453889
DES = (;n=96,m=16,q=2)

# ╔═╡ 89f93310-7924-46a9-890f-c155f3dace67
CONFIG = (;t=80,numpos=1,numiters=100000,randseed=1)

# ╔═╡ 0b62fdf6-9ffd-48b6-891d-4aa1c27f7f43
eff(numstage2) = DES.n/(DES.m+numstage2)

# ╔═╡ fd643c2c-08e5-11eb-0d25-014798a505d1
md"""
### Utility functions
"""

# ╔═╡ fe7cf842-08e5-11eb-0d07-41ed2902921a
function loadresults(dir, pattern)
	regmatches = filter(!isnothing, match.(pattern, readdir(dir)))
	pairs = map(regmatches) do regmatch
		design = (;
			n = parse(Int,regmatch[:n]),
			m = parse(Int,regmatch[:m]),
			q = parse(Int,regmatch[:q]),
		)
		config = (;
			t = parse(Int,regmatch[:t]),
			numpos = parse(Int,regmatch[:numpos]),
			numiters = parse(Int,regmatch[:numiters]),
			randseed = parse(Int,regmatch[:randseed]),
		)
		result = npzread(joinpath(dir,regmatch.match))
		dictionary([design => dictionary([config=>result])])
	end
	return reduce((d1,d2)->mergewith(merge,d1,d2),pairs)
end

# ╔═╡ 3f1b753e-0901-11eb-0060-1f1269588bef
load_numstagetwo(dir; pattern=r"^num_stage_two.n-(?<n>\d+)_m-(?<m>\d+)_q-(?<q>\d+).t-(?<t>\d+).numpos-(?<numpos>\d+).numiters-(?<numiters>\d+).randseed-(?<randseed>\d+).npy$") = loadresults(dir, pattern)

# ╔═╡ 3f1c5288-0901-11eb-391f-dd7dd02d3d9b
load_sensavgperm(dir; pattern=r"^sens_avg_perm.n-(?<n>\d+)_m-(?<m>\d+)_q-(?<q>\d+).t-(?<t>\d+).numpos-(?<numpos>\d+).numiters-(?<numiters>\d+).randseed-(?<randseed>\d+).npy$") = loadresults(dir, pattern)

# ╔═╡ 122151e0-08e6-11eb-18f5-3b277145bc9d
md"""
### HYPER
"""

# ╔═╡ 47d9acce-2fa6-11eb-0260-7767c0252b8c
hyper_eff = eff.(load_numstagetwo(joinpath(DESDIR,"hypergraph"))[DES][CONFIG])

# ╔═╡ f8274a32-2fa6-11eb-3836-9bbe0aaf63cf
hyper_sens = load_sensavgperm(joinpath(DESDIR,"hypergraph"))[DES][CONFIG]

# ╔═╡ 4864b93b-cf34-4a2e-81c9-bf53c0aefee3
md"""
### Consecutive pooling
"""

# ╔═╡ c4c5582b-2512-4e22-9a12-d2af7b3bece3
cons_eff = eff.(load_numstagetwo(joinpath(DESDIR,"consecutive"))[DES][CONFIG])

# ╔═╡ 1cf355b5-f626-4528-8e98-ddc2994f2fa3
cons_sens = load_sensavgperm(joinpath(DESDIR,"consecutive"))[DES][CONFIG]

# ╔═╡ 233f8121-60ea-4d35-ab3d-b7a2ff7b7c32
md"""
### Lexicographic pooling
"""

# ╔═╡ a326b7fc-0cf7-4712-9fa0-795e56df378b
lex_eff = eff.(load_numstagetwo(joinpath(DESDIR,"lexicographic"))[DES][CONFIG])

# ╔═╡ bdbc801d-0f0f-4d5a-ab6a-140f67b403fa
lex_sens = load_sensavgperm(joinpath(DESDIR,"lexicographic"))[DES][CONFIG]

# ╔═╡ 43da055c-58c9-42bb-9e0b-eb865ef15c92
md"""
### Random assignment
"""

# ╔═╡ 91954667-8373-4650-aa05-8652f94bbbfb
RANDASSIGN_SEEDS = 0:2

# ╔═╡ b32f3819-e32a-4867-ab71-ed8d327ab872
randassign_eff = map(dictionary(RANDASSIGN_SEEDS)) do seed
	eff.(load_numstagetwo(joinpath(DESDIR,"rand,assign.seed-$seed"))[DES][CONFIG])
end

# ╔═╡ 287438e0-9c89-4f87-8351-f211ea30b02b
randassign_sens = map(dictionary(RANDASSIGN_SEEDS)) do seed
	load_sensavgperm(joinpath(DESDIR,"rand,assign.seed-$seed"))[DES][CONFIG]
end

# ╔═╡ ae9d54fc-9be4-4b54-a105-1ae44768b0bb
md"""
### Double-pooling
"""

# ╔═╡ 6b84f456-756f-43d6-933b-5f1fd16316e4
DOUBLEPOOL_SEEDS = 0:2

# ╔═╡ 1a3408e1-4287-420c-bff8-24a5987653f1
doublepool_eff = map(dictionary(DOUBLEPOOL_SEEDS)) do seed
	eff.(load_numstagetwo(joinpath(DESDIR,"double,pooling.seed-$seed"))[DES][CONFIG])
end

# ╔═╡ d165a552-7ee3-4099-8059-a283b2cb0627
doublepool_sens = map(dictionary(DOUBLEPOOL_SEEDS)) do seed
	load_sensavgperm(joinpath(DESDIR,"double,pooling.seed-$seed"))[DES][CONFIG]
end

# ╔═╡ 2eb0f323-5c84-4002-9286-eedaba0a3ee9
md"""
### Balanced arrays
"""

# ╔═╡ 39daa00b-7819-4acf-83aa-f63e5d71f5bc
balarray_eff = eff.(load_numstagetwo(joinpath(DESDIR,"balarray"))[DES][CONFIG])

# ╔═╡ 28a7b979-f2f2-4229-8cca-0b9b8e9d7dd3
balarray_sens = load_sensavgperm(joinpath(DESDIR,"balarray"))[DES][CONFIG]

# ╔═╡ 727fba47-eccc-40bd-b864-893ff7b5fdee
md"""
### RS-KS code-based
"""

# ╔═╡ 0f50881c-53ad-4e32-8505-b3daf661a72d
rskscode_eff = eff.(load_numstagetwo(joinpath(DESDIR,"rskscode,f2"))[DES][CONFIG])

# ╔═╡ 3144baf3-beaa-4eee-b77b-84deeb37411b
rskscode_sens = load_sensavgperm(joinpath(DESDIR,"rskscode,f2"))[DES][CONFIG]

# ╔═╡ 0fd98c50-2203-11eb-0f58-d946dcf484b4
md"""
## Figures
"""

# ╔═╡ 635fde7b-3977-4161-8776-db116d6949ad
METHODS = [
	["Random assignment\n(draw $(i+1))" =>
		(;eff=randassign_eff[i],sens=randassign_sens[i]) for i in 0:2]...,
	["Double-pooling\n(draw $(i+1))" =>
		(;eff=doublepool_eff[i],sens=doublepool_sens[i]) for i in 0:2]...,
	"Consecutive"                   => (;eff=cons_eff,sens=cons_sens),
	"Lexicographic"                 => (;eff=lex_eff,sens=lex_sens),
	"Balanced arrays"               => (;eff=balarray_eff,sens=balarray_sens),
	"RS-KS code"                    => (;eff=rskscode_eff,sens=rskscode_sens),
	"HYPER"                         => (;eff=hyper_eff,sens=hyper_sens),
]

# ╔═╡ b170270d-fd2e-40dc-a038-5e343523315d
FIGIDX = vec(permutedims(CartesianIndices((2,6))))

# ╔═╡ 62d35fcc-47d2-44b2-9084-6c030acf902e
begin
	mm_to_units(mm) = floor(Int,mm/25.4*72/0.75)
	mm_to_units(mm1,mm2,mmrest...) = mm_to_units.((mm1,mm2,mmrest...))
end

# ╔═╡ bc933543-3ffe-4858-a06f-c110f65314c7
THEME = Theme(;
	font = "Arial", fontsize = 7,
	linewidth = 1, markersize = 4,
	Axis    = (;xticksize = 2.5, yticksize = 2.5),
	Text    = (;textsize = 7),
	BarPlot = (;label_size = 7, label_offset = 1, strokewidth = 0.5, gap = 0),
	Hist    = (;strokewidth = 0),
	BoxPlot = (;markersize = 3, whiskerwidth = 0.5),
	Legend  = (;framevisible = false, titlefont = "Arial Bold", patchsize = (8,8)),
	Lines   = (;linewidth = 2),
	Scatter = (;markersize = 2),
	Arrows  = (;linewidth = 1.5, arrowsize = 8),
)

# ╔═╡ c2d66877-4dc1-404b-9715-71f2b874e335
with_theme(THEME; Lines=(;linewidth=1)) do
	fig = Figure(; resolution=mm_to_units(180,120))

	for ((methodstr,(eff,sens)),figidx,figlbl) in zip(METHODS,FIGIDX,'a':'z')
		fig[Tuple(figidx)...] = gl = GridLayout()
		Label(gl[1,1,TopLeft()], string(figlbl);
			halign=:left, valign=:top, font="Arial Bold")

		# Efficiency
		ax_eff_ind, _ = lines(gl[1,1], eff)
		ax_eff_box, _ = boxplot(gl[1,2], fill(1,length(eff)), eff;
			axis=(;xautolimitmargin=(0.25f0,0.25f0)))

		# Sensitivity
		ax_sens_ind, _ = lines(gl[2,1], sens)
		ax_sens_box, _ = boxplot(gl[2,2], fill(1,length(sens)), sens;
			axis=(;xautolimitmargin=(0.25f0,0.25f0)))

		# Title
		Label(gl[1,:,Top()], methodstr; padding=(0,0,5,0))

		# Layout, axis limits, etc.
		linkyaxes!(ax_eff_ind,  ax_eff_box)
		linkyaxes!(ax_sens_ind, ax_sens_box)
		hidedecorations!(ax_eff_box)
		hidedecorations!(ax_sens_box)
		ax_eff_box.leftspinevisible = ax_sens_box.leftspinevisible = false
		ax_eff_ind.rightspinecolor  = ax_sens_ind.rightspinecolor  = :grey50
		colgap!(gl, 0)
		colsize!(gl, 2, Relative(0.2))

		ax_eff_ind.limits  = (1,DES.n,3.7,5.8)
		ax_sens_ind.limits = (1,DES.n,0.72,0.76)
		ax_eff_ind.xticks  = ax_sens_ind.xticks = [1; 32:32:96]
		ax_sens_ind.yticks = 0.72:0.01:0.76
		ax_sens_ind.ytickformat = y->string.(convert.(Int,100*y),'%')
		linkxaxes!(ax_eff_ind, ax_sens_ind)
		hidexdecorations!(ax_eff_ind; ticks=false, grid=false)
		rowgap!(gl, 5)
		ax_eff_ind.ylabel  = "Efficiency\n(rel. to indiv. testing)"
		ax_sens_ind.ylabel = "Sensitivity"
		ax_sens_ind.xlabel = "Location of the\npositive individual"
	end
	colgap!(fig.layout, 12)
	rowgap!(fig.layout, 12)

	save("fig-s7.png", fig; px_per_unit=2)
	save("fig-s7.pdf", fig)
	fig
end

# ╔═╡ c7c610b1-8524-4c17-94c7-3e51c797f971
md"""
**Maximums:**
"""

# ╔═╡ 32da154d-c5b8-475d-93cf-be2794fef8eb
(dictionary∘map)(METHODS) do (methodstr,effsens)
	methodstr=>map(findmax,effsens)
end

# ╔═╡ a8f9d744-5834-4168-9b57-3a99b7c12255
md"""
**Minimums:**
"""

# ╔═╡ d9008d4c-70eb-433b-b7cb-e34d527615e3
(dictionary∘map)(METHODS) do (methodstr,effsens)
	methodstr=>map(findmin,effsens)
end

# ╔═╡ 74416809-e943-4ce8-8142-822b9d6aba33
md"""
**Medians:**
"""

# ╔═╡ ee18c5b2-08fc-49e7-981b-192b4dda43c4
(dictionary∘map)(METHODS) do (methodstr,effsens)
	methodstr=>map(median,effsens)
end

# ╔═╡ Cell order:
# ╟─849adbca-08e5-11eb-1956-cd02f2a6be7a
# ╠═95eaf04a-08e5-11eb-043f-27646f1bf207
# ╟─210d331a-09c3-11eb-06ed-e1f4ddc023fc
# ╟─66dd9e7e-f948-4457-8bd6-45345be6abb5
# ╟─dcfbfb7b-3d7b-42f5-9635-cbdc92453889
# ╟─89f93310-7924-46a9-890f-c155f3dace67
# ╟─0b62fdf6-9ffd-48b6-891d-4aa1c27f7f43
# ╟─fd643c2c-08e5-11eb-0d25-014798a505d1
# ╟─fe7cf842-08e5-11eb-0d07-41ed2902921a
# ╟─3f1b753e-0901-11eb-0060-1f1269588bef
# ╟─3f1c5288-0901-11eb-391f-dd7dd02d3d9b
# ╟─122151e0-08e6-11eb-18f5-3b277145bc9d
# ╟─47d9acce-2fa6-11eb-0260-7767c0252b8c
# ╟─f8274a32-2fa6-11eb-3836-9bbe0aaf63cf
# ╟─4864b93b-cf34-4a2e-81c9-bf53c0aefee3
# ╟─c4c5582b-2512-4e22-9a12-d2af7b3bece3
# ╟─1cf355b5-f626-4528-8e98-ddc2994f2fa3
# ╟─233f8121-60ea-4d35-ab3d-b7a2ff7b7c32
# ╟─a326b7fc-0cf7-4712-9fa0-795e56df378b
# ╟─bdbc801d-0f0f-4d5a-ab6a-140f67b403fa
# ╟─43da055c-58c9-42bb-9e0b-eb865ef15c92
# ╟─91954667-8373-4650-aa05-8652f94bbbfb
# ╟─b32f3819-e32a-4867-ab71-ed8d327ab872
# ╟─287438e0-9c89-4f87-8351-f211ea30b02b
# ╟─ae9d54fc-9be4-4b54-a105-1ae44768b0bb
# ╟─6b84f456-756f-43d6-933b-5f1fd16316e4
# ╟─1a3408e1-4287-420c-bff8-24a5987653f1
# ╟─d165a552-7ee3-4099-8059-a283b2cb0627
# ╟─2eb0f323-5c84-4002-9286-eedaba0a3ee9
# ╟─39daa00b-7819-4acf-83aa-f63e5d71f5bc
# ╟─28a7b979-f2f2-4229-8cca-0b9b8e9d7dd3
# ╟─727fba47-eccc-40bd-b864-893ff7b5fdee
# ╟─0f50881c-53ad-4e32-8505-b3daf661a72d
# ╟─3144baf3-beaa-4eee-b77b-84deeb37411b
# ╟─0fd98c50-2203-11eb-0f58-d946dcf484b4
# ╟─635fde7b-3977-4161-8776-db116d6949ad
# ╟─b170270d-fd2e-40dc-a038-5e343523315d
# ╟─62d35fcc-47d2-44b2-9084-6c030acf902e
# ╟─bc933543-3ffe-4858-a06f-c110f65314c7
# ╟─c2d66877-4dc1-404b-9715-71f2b874e335
# ╟─c7c610b1-8524-4c17-94c7-3e51c797f971
# ╟─32da154d-c5b8-475d-93cf-be2794fef8eb
# ╟─a8f9d744-5834-4168-9b57-3a99b7c12255
# ╟─d9008d4c-70eb-433b-b7cb-e34d527615e3
# ╟─74416809-e943-4ce8-8142-822b9d6aba33
# ╟─ee18c5b2-08fc-49e7-981b-192b4dda43c4
