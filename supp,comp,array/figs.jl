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
# Figure: comparison of array designs with additional HYPER designs
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
hypereff = loadeff("../simulation/designs/hypergraph/")

# ╔═╡ f8274a32-2fa6-11eb-3836-9bbe0aaf63cf
hypersens = loadsens("../simulation/designs/hypergraph/")

# ╔═╡ 2e0cc554-08e6-11eb-0186-bba8c2aee7ec
md"""
### Array
"""

# ╔═╡ 339ca8ce-08e6-11eb-19c3-5f0608ee2794
arrayeff = (dictionary∘Dict)(
	(8,12) => only(loadeff("../simulation/designs/array,8x12/")),
	(16,24) => only(loadeff("../simulation/designs/array,16x24/"))
)

# ╔═╡ 33a20596-0902-11eb-2d50-cb33f008ffd9
arraysens = (dictionary∘Dict)(
	(8,12) => only(loadsens("../simulation/designs/array,8x12/")),
	(16,24) => only(loadsens("../simulation/designs/array,16x24/"))
)

# ╔═╡ a92db3c4-5436-11eb-181b-3f1a71204b12
md"""
### Individual testing
"""

# ╔═╡ 936c39ce-368e-11eb-0cb2-938376a6c18f
LOD = 100.0

# ╔═╡ 72f1ae0a-544c-11eb-0a55-f18f192d64c2
indiveff = [1.0 for _ in eachindex(prevalence)]

# ╔═╡ 21073ca8-3684-11eb-0a25-8d16a9799cf1
indivsens = mapslices(ViralLoad, dims=1) do vl
	posloads = filter(>(0),vl)
	isempty(posloads) ? 1.0 : mean(>(LOD),posloads)
end |> x->OffsetArray(vec(x),-1)

# ╔═╡ 0fd98c50-2203-11eb-0f58-d946dcf484b4
md"""
## Figure
"""

# ╔═╡ 222d4b16-7b8c-11eb-0622-797f8f7bba13
POLYORDER = 8

# ╔═╡ 19411dcc-220c-11eb-1b45-4d393e53f5bf
PLTWINDOW = IdentityRange(20:110)

# ╔═╡ a6c3504c-543f-11eb-2f20-791273bc06b8
SHADELIMS = [40,90]

# ╔═╡ e3dcaa08-5036-11eb-2f56-6f1d48bb9188
shadecolor = colorant"gray87"

# ╔═╡ 6c66ebe0-543b-11eb-1619-294095a5ae7a
let compfigs = [96=>[(n=96,m=16,q=2),(n=96,m=20,q=2)],384=>[(n=384,m=32,q=2),(n=384,m=40,q=2)]]
	# Plot colors
	dcolors = distinguishable_colors(7, [RGB(1,1,1)], dropseed=true)
	indivcolor = dcolors[1]
	pbestcolor = dcolors[2]
	arraycolor = dcolors[3]
	hypercolor = Dict(
		16=>dcolors[4],
		32=>dcolors[5],
		20=>dcolors[6],
		40=>dcolors[7],
	)
	
	# Efficiency gain plots
	effplts = map(compfigs) do (n,hyperdesigns)
		plot(rectangle(SHADELIMS,[0,n]),opacity=0.4,linewidth=0,color=shadecolor)

		# Array
		plot!(arrayeff[n == 96 ? (8,12) : (16,24)],linestyle=:solid,color=arraycolor)

		# HYPER
		for hyperdesign in hyperdesigns
			plot!(hypereff[hyperdesign],linestyle=:dash,color=hypercolor[hyperdesign.m])
		end

		# Individual
		plot!(indiveff,linestyle=:dot,color=indivcolor)
	end

	# Sensitivity plots
	sensplts = map(compfigs) do (n,hyperdesigns)
		plot(rectangle(SHADELIMS,[0,1]),opacity=0.4,linewidth=0,color=shadecolor)

		# Array
		scatter!(arraysens[n == 96 ? (8,12) : (16,24)],opacity=0.6,color=arraycolor)
		plot!(polyfit(log10.(prevalence[PLTWINDOW]),arraysens[n == 96 ? (8,12) : (16,24)][PLTWINDOW],POLYORDER),
			linestyle=:solid,color=arraycolor)

		# HYPER
		for hyperdesign in hyperdesigns
			scatter!(hypersens[hyperdesign],opacity=0.6,color=hypercolor[hyperdesign.m])
			plot!(polyfit(log10.(prevalence[PLTWINDOW]),hypersens[hyperdesign][PLTWINDOW],POLYORDER),
				linestyle=:dash,color=hypercolor[hyperdesign.m])
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
	
	plot!(effplts[1],ylims=(0, 6*1.14),yticks=[0,1,96/20,96/16])
	plot!(effplts[2],ylims=(0,12*1.14),yticks=[0,1,384/40,384/48,384/32])
	plot!.(sensplts,ylims=(0.6,0.9),yticks=(0.65:0.05:0.85,string.(65:5:85,"%")))
	plot!(effplts[1], ylabel="individuals/test")
	plot!(sensplts[1],ylabel="sensitivity")
	
	annotate!(effplts[1],[
		(21.75, 6.1, text(hypername(last(compfigs[1])[1]), 8, :left, :bottom, hypercolor[last(compfigs[1])[1].m])),
		(21.75, 4.7, text(hypername(last(compfigs[1])[2]), 8, :left, :top,    hypercolor[last(compfigs[1])[2].m])),
		(21.75, 4.9, text("8x12 array",                    8, :left, :bottom, arraycolor)),
		(21.75, 1.1, text("Individual testing",            8, :left, :bottom, indivcolor))
	])
	
	annotate!(effplts[2],[
		(21.75, 12.2, text(hypername(last(compfigs[2])[1]), 8, :left, :bottom, hypercolor[last(compfigs[2])[1].m])),
		(21.75,  9.4, text(hypername(last(compfigs[2])[2]), 8, :left, :top,    hypercolor[last(compfigs[2])[2].m])),
		(21.75,  9.8, text("16x24 array",                   8, :left, :bottom, arraycolor)),
		(21.75,  1.2, text("Individual testing",            8, :left, :bottom, indivcolor))
	])
	
	plot!.(effplts, bottom_margin=0mm)
	plot!.(sensplts,bottom_margin=5mm,top_margin=-2.5mm)
	plot!.([effplts; sensplts],right_margin=-1mm)
	
	plot(effplts...,sensplts...,layout=(2,2),link=:x)
	savefig("fig.png"); plot!()
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
# ╟─a92db3c4-5436-11eb-181b-3f1a71204b12
# ╟─936c39ce-368e-11eb-0cb2-938376a6c18f
# ╟─72f1ae0a-544c-11eb-0a55-f18f192d64c2
# ╟─21073ca8-3684-11eb-0a25-8d16a9799cf1
# ╟─0fd98c50-2203-11eb-0f58-d946dcf484b4
# ╟─222d4b16-7b8c-11eb-0622-797f8f7bba13
# ╟─19411dcc-220c-11eb-1b45-4d393e53f5bf
# ╟─a6c3504c-543f-11eb-2f20-791273bc06b8
# ╟─e3dcaa08-5036-11eb-2f56-6f1d48bb9188
# ╟─6c66ebe0-543b-11eb-1619-294095a5ae7a
