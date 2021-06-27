### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ 95eaf04a-08e5-11eb-043f-27646f1bf207
using DelimitedFiles, Dictionaries, IdentityRanges, LsqFit, NPZ, OffsetArrays, Plots, Plots.Measures, SparseArrays, Statistics

# ╔═╡ bfbc2e08-0c3b-11eb-10c0-7bd1d3a723b9
# Code to load CSC matrix from NPZ file
# based on description here: https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html
# see also:
# + https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.save_npz.html
# + https://docs.julialang.org/en/v1/stdlib/SparseArrays/#man-csc

module NPZ_CSC
using PyCall: pyimport, pyrepr
using SparseArrays: SparseMatrixCSC

np = pyimport("numpy")

function load(filename)
	contents = np.load(filename, allow_pickle=true)
	format = get(contents, "format")
	pyrepr(format) == "array(b'csc', dtype='|S3')" ||
		throw("Array format seems wrong. Got $(pyrepr(format))")
	
	shape = get(contents, "shape")
	data = get(contents, "data")
	indptr = get(contents, "indptr") .+ 1
	indices = get(contents, "indices") .+ 1
	
	SparseMatrixCSC(shape[1], shape[2], indptr, indices, data)
end

end

# ╔═╡ 849adbca-08e5-11eb-1956-cd02f2a6be7a
md"""
# Plots for reduced limit of detection (LOD=4.0)
"""

# ╔═╡ 34f9cdda-0901-11eb-1117-099323e644ef
gr(linewidth=2,markerstrokecolor=:auto,markersize=1.5,legend=nothing,
	dpi=200,size=(900,400),titlefontsize=9,guidefontsize=9,tickfontsize=7)

# ╔═╡ 289fbed4-21ff-11eb-261a-25eace4cbf4e
# polynomial fit
function polyfit(xx,yy,order)
	axes(xx) == axes(yy) || throw("Arrays must have same axes. Got: $(axes(xx)) and $(axes(yy))")
	x, y = OffsetArrays.no_offset_view(xx), OffsetArrays.no_offset_view(yy)
	
	fit_func = (x,p) -> evalpoly.(x,Ref(p))
	popt = curve_fit(fit_func,x,y,ones(order)).param
	return fit_func(xx,popt)
end

# ╔═╡ 6539e014-543f-11eb-17c3-991de18e923c
# plot a rectangle
rectangle(xlims,ylims) = Shape([xlims; reverse(xlims)],repeat(ylims,inner=2))

# ╔═╡ eb7bec9e-5453-11eb-0876-c93037e57988
hypername((n,m,q)) = "HYPER (m=$m, q=$q)"

# ╔═╡ 210d331a-09c3-11eb-06ed-e1f4ddc023fc
md"""
## Load data
"""

# ╔═╡ fd643c2c-08e5-11eb-0d25-014798a505d1
md"""
### Utility functions
"""

# ╔═╡ fe7cf842-08e5-11eb-0d07-41ed2902921a
function loadresults(dir, pattern)
	regmatches = filter(!isnothing, match.(pattern, readdir(dir)))
	pairs = map(regmatches) do regmatch
		design = (n=parse(Int,regmatch[:n]), m=parse(Int,regmatch[:m]), q=parse(Int,regmatch[:q]))
		result = OffsetVector(npzread(joinpath(dir,regmatch.match)),-1)
		design => result
	end
	return dictionary(pairs)
end

# ╔═╡ 3f1b753e-0901-11eb-0060-1f1269588bef
loadeff(dir; pattern=r"^Eff_avg.n-(?<n>\d+)_m-(?<m>\d+)_q-(?<q>\d+).npy$") = loadresults(dir, pattern)

# ╔═╡ 3f1c5288-0901-11eb-391f-dd7dd02d3d9b
loadsens(dir; pattern=r"^Recall_combined.n-(?<n>\d+)_m-(?<m>\d+)_q-(?<q>\d+).npy$") = loadresults(dir, pattern)

# ╔═╡ 05c6923e-08e6-11eb-1d9b-21afc86a6b04
md"""
### Viral load distributions and average prevalence
"""

# ╔═╡ c6c3f3ca-0c3b-11eb-2593-439620bd9e78
ViralLoad = NPZ_CSC.load("../simulation/Simulated_populations/seir_viral_loads_swab.viral_loads.npz")

# ╔═╡ cf1b2b1a-0c3b-11eb-2159-6ffeb682d595
prevalence = OffsetArray(vec(mean(>(0), ViralLoad, dims=1)), -1)

# ╔═╡ 122151e0-08e6-11eb-18f5-3b277145bc9d
md"""
### HYPER
"""

# ╔═╡ 47d9acce-2fa6-11eb-0260-7767c0252b8c
hypereff = loadeff("designs/hypergraph/")

# ╔═╡ f8274a32-2fa6-11eb-3836-9bbe0aaf63cf
hypersens = loadsens("designs/hypergraph/")

# ╔═╡ 2e0cc554-08e6-11eb-0186-bba8c2aee7ec
md"""
### Array
"""

# ╔═╡ 339ca8ce-08e6-11eb-19c3-5f0608ee2794
arrayeff = (dictionary∘Dict)(
	(8,12) => only(loadeff("designs/array,8x12/")),
	(16,24) => only(loadeff("designs/array,16x24/"))
)

# ╔═╡ 33a20596-0902-11eb-2d50-cb33f008ffd9
arraysens = (dictionary∘Dict)(
	(8,12) => only(loadsens("designs/array,8x12/")),
	(16,24) => only(loadsens("designs/array,16x24/"))
)

# ╔═╡ 3706cb8e-08e6-11eb-3c27-f1142c9d9159
md"""
### P-BEST
"""

# ╔═╡ 6f20ef1a-2216-11eb-37d8-cb470acc99e4
pbesteff = only(loadeff("designs/pbest/"))

# ╔═╡ 744ec278-2216-11eb-2aee-6138f60f7722
pbestsens = only(loadsens("designs/pbest/"))

# ╔═╡ a92db3c4-5436-11eb-181b-3f1a71204b12
md"""
### Individual testing
"""

# ╔═╡ 936c39ce-368e-11eb-0cb2-938376a6c18f
LOD = 4.0

# ╔═╡ 72f1ae0a-544c-11eb-0a55-f18f192d64c2
indiveff = [1.0 for _ in eachindex(prevalence)]

# ╔═╡ 21073ca8-3684-11eb-0a25-8d16a9799cf1
indivsens = mapslices(ViralLoad, dims=1) do vl
	posloads = filter(>(0),vl)
	isempty(posloads) ? 1.0 : mean(>(LOD),posloads)
end |> x->OffsetArray(vec(x),-1)

# ╔═╡ 0fd98c50-2203-11eb-0f58-d946dcf484b4
md"""
## Figures
"""

# ╔═╡ 222d4b16-7b8c-11eb-0622-797f8f7bba13
POLYORDER = 8

# ╔═╡ 19411dcc-220c-11eb-1b45-4d393e53f5bf
PLTWINDOW = IdentityRange(20:110)

# ╔═╡ a6c3504c-543f-11eb-2f20-791273bc06b8
SHADELIMS = [40,90]

# ╔═╡ e3dcaa08-5036-11eb-2f56-6f1d48bb9188
shadecolor = colorant"gray87"

# ╔═╡ 71950c84-543b-11eb-3c6b-4b4114ade7bf
md"""
### Fig a
"""

# ╔═╡ 6c66ebe0-543b-11eb-1619-294095a5ae7a
let compfigs = [96=>(n=96,m=16,q=2),384=>(n=384,m=32,q=2)]
	# Plot colors
	dcolors = distinguishable_colors(7, [RGB(1,1,1)], dropseed=true)
	indivcolor = dcolors[1]
	pbestcolor = dcolors[2]
	arraycolor = dcolors[3]
	hypercolor = Dict(
		16=>dcolors[4],
		32=>dcolors[5],
	)
	
	# Efficiency gain plots
	effplts = map(compfigs) do (n,hyperdesign)
		plot(rectangle(SHADELIMS,[0,n]),opacity=0.4,linewidth=0,color=shadecolor)

		# Array
		plot!(arrayeff[n == 96 ? (8,12) : (16,24)],linestyle=:solid,color=arraycolor)

		# P-BEST
		if n == 384
			plot!(pbesteff,linestyle=:solid,color=pbestcolor)
		end

		# HYPER
		plot!(hypereff[hyperdesign],linestyle=:dash,color=hypercolor[hyperdesign.m])

		# Individual
		plot!(indiveff,linestyle=:dot,color=indivcolor)
	end

	# Sensitivity plots
	sensplts = map(compfigs) do (n,hyperdesign)
		plot(rectangle(SHADELIMS,[0,1]),opacity=0.4,linewidth=0,color=shadecolor)

		# Array
		scatter!(arraysens[n == 96 ? (8,12) : (16,24)],opacity=0.6,color=arraycolor)
		plot!(polyfit(log10.(prevalence[PLTWINDOW]),arraysens[n == 96 ? (8,12) : (16,24)][PLTWINDOW],POLYORDER),
			linestyle=:solid,color=arraycolor)

		# P-BEST
		if n == 384
			scatter!(pbestsens,opacity=0.6,color=pbestcolor)
			plot!(polyfit(log10.(prevalence[PLTWINDOW]),pbestsens[PLTWINDOW],POLYORDER),
				linestyle=:solid,color=pbestcolor)
		end

		# HYPER
		scatter!(hypersens[hyperdesign],opacity=0.6,color=hypercolor[hyperdesign.m])
		plot!(polyfit(log10.(prevalence[PLTWINDOW]),hypersens[hyperdesign][PLTWINDOW],POLYORDER),
			linestyle=:dash,color=hypercolor[hyperdesign.m])

		# Individual
		scatter!(indivsens,opacity=0.6,color=indivcolor)
		plot!(polyfit(log10.(prevalence[PLTWINDOW]),indivsens[PLTWINDOW],POLYORDER),
			linestyle=:dot,color=indivcolor)
	end
	
	# Figure (titles, axes, annotations)
	title!.(effplts,string.("n=",first.(compfigs)," individuals"))
	
	plot!.(sensplts,xlims=(20,110),xticks=(20:10:110,["";["$day\n($(round(100*prevalence[day],digits=2))%)" for day in 30:10:100];""]),xlabel="day (prevalence)")
	plot!.(effplts, xlims=(20,110),xticks=(20:10:110,""))
	
	plot!(effplts[1],ylims=(0, 6*1.14),yticks=[0,1,96/20,96/16])
	plot!(effplts[2],ylims=(0,12*1.14),yticks=[0,1,384/40,384/48,384/32])
	plot!.(sensplts,ylims=(0.7,1.0),yticks=(0.75:0.05:0.95,string.(75:5:95,"%")))
	plot!(effplts[1], ylabel="individuals/test")
	plot!(sensplts[1],ylabel="sensitivity")
	
	annotate!(effplts[1],[
		(21.75, 6.1, text(hypername(last(compfigs[1])), 8, :left, :bottom, hypercolor[last(compfigs[1]).m])),
		(21.75, 4.9, text("8x12 array",                 8, :left, :bottom, arraycolor)),
		(21.75, 1.1, text("Individual testing",         8, :left, :bottom, indivcolor))
	])
	
	annotate!(effplts[2],[
		(21.75, 12.2, text(hypername(last(compfigs[2])), 8, :left, :bottom, hypercolor[last(compfigs[2]).m])),
		(21.75,  9.8, text("16x24 array",                8, :left, :bottom, arraycolor)),
		(21.75,  7.8, text("P-BEST",                     8, :left, :top,    pbestcolor)),
		(21.75,  1.2, text("Individual testing",         8, :left, :bottom, indivcolor))
	])
	
	annotate!(effplts[1],10,7.5,text("a",14))

	plot!.(effplts, bottom_margin=0mm)
	plot!.(sensplts,bottom_margin=5mm,top_margin=-2.5mm)
	plot!.([effplts; sensplts],right_margin=-1mm)
	
	plot(effplts...,sensplts...,layout=(2,2),link=:x)
	savefig("fig-a.png"); plot!()
end

# ╔═╡ 157b5eb0-5474-11eb-2bac-4daf68438917
md"""
**Average sensitivity loss w.r.t. 8x12 array in percentage points:**
"""

# ╔═╡ 0159c02a-5474-11eb-137e-8b07b5079fb7
100*mean((arraysens[(8,12)] .- hypersens[(n=96,m=16,q=2)])[40:90])

# ╔═╡ c7c95bee-5474-11eb-10aa-d3d7e819accb
md"""
**Average sensitivity loss w.r.t. 16x24 array in percentage points:**
"""

# ╔═╡ cc428054-5474-11eb-02bf-af72fb04320d
100*mean((arraysens[(16,24)] .- hypersens[(n=384,m=32,q=2)])[40:90])

# ╔═╡ 5b53ff88-5458-11eb-1313-b74a1235ad08
md"""
### Figs b-c
"""

# ╔═╡ 7eea53ac-5457-11eb-39e2-19c11bcff57d
let compfigs = [384=>[(n=384,m=32,q=2),(n=384,m=16,q=2),(n=384,m=12,q=2)],384=>[(n=384,m=12,q=1),(n=384,m=12,q=2),(n=384,m=12,q=3)]]
	# Plot colors and linestyles
	dcolors = distinguishable_colors(7, [RGB(1,1,1)], dropseed=true)
	indivcolor = dcolors[1]
	hypercolor = Dict(
		16=>dcolors[4],
		32=>dcolors[5],
		12=>dcolors[7],
	)
	hyperstyle = Dict(1=>:dot, 2=>:dash, 3=>:solid)
	
	# Efficiency gain plots
	effplts = map(compfigs) do (n,hyperdesigns)
		plot(rectangle(SHADELIMS,[0,n]),opacity=0.4,linewidth=0,color=shadecolor)

		# HYPER
		for design in hyperdesigns
			plot!(hypereff[design],linestyle=hyperstyle[design.q],color=hypercolor[design.m])
		end

		# Individual
		plot!(indiveff,linestyle=:dot,color=indivcolor)
	end

	# Sensitivity plots
	sensplts = map(compfigs) do (n,hyperdesigns)
		plot(rectangle(SHADELIMS,[0,1]),opacity=0.4,linewidth=0,color=shadecolor)

		# HYPER
		for design in hyperdesigns
			scatter!(hypersens[design],opacity=0.6,color=hypercolor[design.m])
			plot!(polyfit(log10.(prevalence[PLTWINDOW]),hypersens[design][PLTWINDOW],POLYORDER),
				linestyle=hyperstyle[design.q],color=hypercolor[design.m])
		end

		# Individual
		scatter!(indivsens,opacity=0.6,color=indivcolor)
		plot!(polyfit(log10.(prevalence[PLTWINDOW]),indivsens[PLTWINDOW],POLYORDER),
			linestyle=:dot,color=indivcolor)
	end
	
	# Figure (titles, axes, annotations)
	title!.(effplts,string.("n=",first.(compfigs)," individuals"))
	
	plot!.(sensplts,xlims=(20,110),xticks=(20:10:110,["";["$day\n($(round(100*prevalence[day],digits=2))%)" for day in 30:10:100];""]),xlabel="day (prevalence)")
	plot!.(effplts, xlims=(20,110),xticks=(20:10:110,""))
	
	plot!.(effplts, ylims=(0,32*1.02),yticks=[1,384/32,384/16,384/12])
	plot!.(sensplts,ylims=(0.7,1.0),yticks=(0.75:0.05:0.95,string.(75:5:95,"%")))
	plot!(effplts[1], ylabel="individuals/test")
	plot!(sensplts[1],ylabel="sensitivity")
	
	annotate!(effplts[1],[
		(21.75, 11.2, text(hypername(last(compfigs[1])[1]), 8, :left, :top, hypercolor[last(compfigs[1])[1].m])),
		(21.75, 22.4, text(hypername(last(compfigs[1])[2]), 8, :left, :top, hypercolor[last(compfigs[1])[2].m])),
		(45.75, 32.0, text(hypername(last(compfigs[1])[3]), 8, :left, :top, hypercolor[last(compfigs[1])[3].m])),
		(21.75,  1.4, text("Individual testing",            8, :left, :bottom, indivcolor))
	])
	annotate!(effplts[1],10,35,text("b",14))

	annotate!(effplts[2],[
		(21.75, 18.2, text(hypername(last(compfigs[2])[1]), 8, :left, :top, hypercolor[last(compfigs[2])[1].m])),
		(21.75, 28.4, text(hypername(last(compfigs[2])[2]), 8, :left, :top, hypercolor[last(compfigs[2])[2].m])),
		(51.75, 32.0, text(hypername(last(compfigs[2])[3]), 8, :left, :top, hypercolor[last(compfigs[2])[3].m])),
		(21.75,  1.4, text("Individual testing",            8, :left, :bottom, indivcolor))
	])
	annotate!(effplts[2],12,35,text("c",14))

	plot!.(effplts, bottom_margin=0mm)
	plot!.(sensplts,bottom_margin=5mm,top_margin=-2.5mm)
	plot!.([effplts; sensplts],right_margin=-1mm)
	
	plot(effplts...,sensplts...,layout=(2,2),link=:x)
	savefig("fig-b,c.png"); plot!()
end

# ╔═╡ Cell order:
# ╟─849adbca-08e5-11eb-1956-cd02f2a6be7a
# ╠═95eaf04a-08e5-11eb-043f-27646f1bf207
# ╠═34f9cdda-0901-11eb-1117-099323e644ef
# ╟─289fbed4-21ff-11eb-261a-25eace4cbf4e
# ╟─6539e014-543f-11eb-17c3-991de18e923c
# ╟─eb7bec9e-5453-11eb-0876-c93037e57988
# ╟─210d331a-09c3-11eb-06ed-e1f4ddc023fc
# ╟─fd643c2c-08e5-11eb-0d25-014798a505d1
# ╟─fe7cf842-08e5-11eb-0d07-41ed2902921a
# ╟─3f1b753e-0901-11eb-0060-1f1269588bef
# ╟─3f1c5288-0901-11eb-391f-dd7dd02d3d9b
# ╟─05c6923e-08e6-11eb-1d9b-21afc86a6b04
# ╟─bfbc2e08-0c3b-11eb-10c0-7bd1d3a723b9
# ╟─c6c3f3ca-0c3b-11eb-2593-439620bd9e78
# ╟─cf1b2b1a-0c3b-11eb-2159-6ffeb682d595
# ╟─122151e0-08e6-11eb-18f5-3b277145bc9d
# ╟─47d9acce-2fa6-11eb-0260-7767c0252b8c
# ╟─f8274a32-2fa6-11eb-3836-9bbe0aaf63cf
# ╟─2e0cc554-08e6-11eb-0186-bba8c2aee7ec
# ╟─339ca8ce-08e6-11eb-19c3-5f0608ee2794
# ╟─33a20596-0902-11eb-2d50-cb33f008ffd9
# ╟─3706cb8e-08e6-11eb-3c27-f1142c9d9159
# ╟─6f20ef1a-2216-11eb-37d8-cb470acc99e4
# ╟─744ec278-2216-11eb-2aee-6138f60f7722
# ╟─a92db3c4-5436-11eb-181b-3f1a71204b12
# ╟─936c39ce-368e-11eb-0cb2-938376a6c18f
# ╟─72f1ae0a-544c-11eb-0a55-f18f192d64c2
# ╟─21073ca8-3684-11eb-0a25-8d16a9799cf1
# ╟─0fd98c50-2203-11eb-0f58-d946dcf484b4
# ╟─222d4b16-7b8c-11eb-0622-797f8f7bba13
# ╟─19411dcc-220c-11eb-1b45-4d393e53f5bf
# ╟─a6c3504c-543f-11eb-2f20-791273bc06b8
# ╟─e3dcaa08-5036-11eb-2f56-6f1d48bb9188
# ╟─71950c84-543b-11eb-3c6b-4b4114ade7bf
# ╟─6c66ebe0-543b-11eb-1619-294095a5ae7a
# ╟─157b5eb0-5474-11eb-2bac-4daf68438917
# ╠═0159c02a-5474-11eb-137e-8b07b5079fb7
# ╟─c7c95bee-5474-11eb-10aa-d3d7e819accb
# ╠═cc428054-5474-11eb-02bf-af72fb04320d
# ╟─5b53ff88-5458-11eb-1313-b74a1235ad08
# ╟─7eea53ac-5457-11eb-39e2-19c11bcff57d
