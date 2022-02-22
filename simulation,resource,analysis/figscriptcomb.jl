### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ 95eaf04a-08e5-11eb-043f-27646f1bf207
begin
	import Pkg; Pkg.activate(@__DIR__)
	using Dictionaries, IdentityRanges, NPZ, OffsetArrays, Statistics
	using CairoMakie, Colors, ColorSchemes
end

# ╔═╡ 849adbca-08e5-11eb-1956-cd02f2a6be7a
md"""
# Figure: resource analyses (combined)
"""

# ╔═╡ b58f8921-cd4c-4740-b17e-fc669a999b84
DAYS, SAVENAME = IdentityRange(53:53), "fig-s11.png"
# DAYS, SAVENAME = IdentityRange(80:80), "fig-s12.png"
# DAYS, SAVENAME = IdentityRange(83:83), "fig-s13.png"
# DAYS, SAVENAME = IdentityRange(84:84), "fig-s14.png"
# DAYS, SAVENAME = IdentityRange(90:90), "fig-s15.png"
# DAYS, SAVENAME = IdentityRange(93:93), "fig-s16.png"

# ╔═╡ 210d331a-09c3-11eb-06ed-e1f4ddc023fc
md"""
## Load data
"""

# ╔═╡ a61d6b1b-a795-4883-b0be-55184cc79e58
DESDIR = "../../../hyper-group-testing/simulation/designs/"

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
### Average prevalence
"""

# ╔═╡ cf1b2b1a-0c3b-11eb-2159-6ffeb682d595
prevalence =
	OffsetVector(npzread(joinpath(DESDIR,"individual","Prevalence.npy")),-1)[DAYS]

# ╔═╡ 122151e0-08e6-11eb-18f5-3b277145bc9d
md"""
### HYPER
"""

# ╔═╡ 47d9acce-2fa6-11eb-0260-7767c0252b8c
hypereff = map(x->x[DAYS], loadeff(joinpath(DESDIR,"hypergraph")))

# ╔═╡ f8274a32-2fa6-11eb-3836-9bbe0aaf63cf
hypersens = map(x->x[DAYS], loadsens(joinpath(DESDIR,"hypergraph")))

# ╔═╡ 2e0cc554-08e6-11eb-0186-bba8c2aee7ec
md"""
### Array
"""

# ╔═╡ 339ca8ce-08e6-11eb-19c3-5f0608ee2794
arrayeff = (dictionary∘Dict)(
	(8,12) => only(loadeff(joinpath(DESDIR,"array,8x12")))[DAYS],
	(16,24) => only(loadeff(joinpath(DESDIR,"array,16x24")))[DAYS],
)

# ╔═╡ 33a20596-0902-11eb-2d50-cb33f008ffd9
arraysens = (dictionary∘Dict)(
	(8,12) => only(loadsens(joinpath(DESDIR,"array,8x12")))[DAYS],
	(16,24) => only(loadsens(joinpath(DESDIR,"array,16x24")))[DAYS],
)

# ╔═╡ 3706cb8e-08e6-11eb-3c27-f1142c9d9159
md"""
### P-BEST
"""

# ╔═╡ 6f20ef1a-2216-11eb-37d8-cb470acc99e4
pbesteff = only(loadeff(joinpath(DESDIR,"pbest")))[DAYS]

# ╔═╡ 744ec278-2216-11eb-2aee-6138f60f7722
pbestsens = only(loadsens(joinpath(DESDIR,"pbest/")))[DAYS]

# ╔═╡ a92db3c4-5436-11eb-181b-3f1a71204b12
md"""
### Individual testing
"""

# ╔═╡ 72f1ae0a-544c-11eb-0a55-f18f192d64c2
indiveff = OffsetVector(npzread(joinpath(DESDIR,"individual","Eff_avg.npy")),-1)[DAYS]

# ╔═╡ 21073ca8-3684-11eb-0a25-8d16a9799cf1
indivsens = OffsetVector(npzread(joinpath(DESDIR,"individual","Recall_combined.npy")),-1)[DAYS]

# ╔═╡ 1eb03185-c9c9-4ec7-aa2f-bccbcb32fc9c
md"""
## Resource analysis
"""

# ╔═╡ d56d4a09-6498-40f3-ac6e-ad7be1e7dbf4
SWABS = [12*2^j for j in 0:9]

# ╔═╡ 794090d7-65fe-41d3-9705-2f32a1fd92b2
TESTS = SWABS

# ╔═╡ ba58813a-818e-42fa-a177-bda658451233
function effcapacity(sampbudget::Number, testbudget::Number, n::Number, eff, sens)
	ntests = n./eff                             # tests/batch (for both stages)
	b = min.(sampbudget/n, testbudget./ntests)  # batches/day
	bavg = mean(b)
	effcap = bavg > 0.9 ? mean(n.*b.*sens) : 0.0
	return (;effcap, bavg)
end

# ╔═╡ 088ae14a-efa0-4428-9e09-204253c671c3
capsgrid = map(Iterators.product(SWABS,TESTS)) do budget
	# Individual testing
	indiv = (effcapacity(budget..., 1, indiveff, indivsens)..., des=())

	# HYPER (optimized across designs)
	hypercap, hyperdes = findmax(keys(hypereff)) do des
		effcapacity(budget..., des.n, hypereff[des], hypersens[des])
	end
	hyper = (hypercap...,des=hyperdes)

	# Arrays (optimized across designs)
	arraycap, arraydes = findmax(keys(arrayeff)) do des
		effcapacity(budget..., prod(des), arrayeff[des], arraysens[des])
	end
	array = (arraycap...,des=arraydes)

	# P-BEST
	pbest = (effcapacity(budget..., 384, pbesteff, pbestsens)..., des=())
	
	return (;indiv, hyper, array, pbest)
end

# ╔═╡ f5aafd70-e89f-4802-9ec7-d0ac81a17cd2
md"""
## Figures
"""

# ╔═╡ 4308103a-bf41-48f3-bdcb-d2bcf9cf9984
LEGEND = (;
	indiv  = ColorSchemes.OrRd_4[1] => "Individual testing",
	array  = colorant"dodgerblue"   => "Array design*",
	pbest  = colorant"forestgreen"  => "P-BEST*",
	hyper1 = ColorSchemes.OrRd_4[2] => "HYPER (q=1)",
	hyper2 = ColorSchemes.OrRd_4[3] => "HYPER (q=2)",
	hyper3 = ColorSchemes.OrRd_4[4] => "HYPER (q=3)",
)

# ╔═╡ 5b786730-e30e-4d6e-9229-d9087a08e5c3
METHODCOLOR = (;
	indiv = _   -> LEGEND[:indiv][1],
	hyper = des -> LEGEND[Symbol("hyper$(des.q)")][1],
	array = _   -> LEGEND[:array][1],
	pbest = _   -> LEGEND[:pbest][1],
)

# ╔═╡ 6d6411bb-2910-4c1c-82eb-154093782d33
SELECTED = [(9,1) (9,2); (9,3) (9,7); (4,2) (6,3)] .|> CartesianIndex

# ╔═╡ 0c901f1d-7da4-4688-b9a6-f5ddaa85ffc7
DAYSTR = if length(DAYS) == 1
	d = only(DAYS)
	p = round(prevalence[d].*100; digits=2)
	"(day $d: prevalence of $p%)"
else
	d1, d2 = extrema(DAYS)
	p1, p2 = round.(prevalence[[d1,d2]].*100; digits=2)
	"(days $d1-$d2: prevalence grows from $p1% to $p2%)"
end

# ╔═╡ 26c3d755-a2f6-40ef-9f8d-6e8942d4bd3f
with_theme() do
	fig = Figure(;resolution=(1200,1440))

	## Row 1: Selected cells + Best methods
	fig[1,1] = row1 = GridLayout(;alignmode=Outside())

	## Selected cells
	row1[1,1] = selfig = GridLayout()
	selidx = permutedims(collect(pairs(IndexCartesian(), SELECTED)))
	for ((figidx, grididx), figlbl) in zip(selidx, 'a':'z')
		# Compute bar values, labels, etc.
		bars = map(collect(pairs(capsgrid[grididx]))) do (method,(effcap,bavg,des))
			methodstr =
				(;indiv="Indiv.",hyper="HYPER",array="Array",pbest="P-BEST")[method]
			desstr = (;
				indiv = _   -> "",
				hyper = des -> des.q == 1 ?
					"n:$(convert(Int,des.n/des.m))\nm:1\nq:1" :
					"n:$(des.n)\nm:$(des.m)\nq:$(des.q)",
				array = dim -> join(dim,'x'),
				pbest = _   -> "",
			)[method](des)
			intround = x -> isinteger(x) ? convert(Int, x) : round(x; digits=1)
			bstr = (;
				indiv = _   -> "",
				hyper = des -> des.q == 1 ?
					"b:$(intround(des.m*bavg))" : "b:$(intround(bavg))",
				array = dim -> "b:$(intround(bavg))",
				pbest = _   -> "b:$(intround(bavg))",
			)[method](des)
			(;methodstr, effcap, desstr, bstr)
		end

		# Prepare axis
		ax = Axis(selfig[Tuple(figidx)...],
			title = "$(SWABS[grididx[1]]) samples, $(TESTS[grididx[2]]) tests",
			ylabel = "Effective\nscreening cap.",
			yticksvisible = false, yticklabelsvisible = false,
			xticks = (axes(bars,1), getindex.(bars,:methodstr)),
			xgridvisible = false, ygridvisible = false,
			yautolimitmargin = (0.0f0,0.25f0),
		)
		tightlimits!(ax, Left(), Right(), Bottom())

		# Draw barplot
		barplot!(ax, getindex.(bars,:effcap);
			bar_labels=:y, label_formatter=y->iszero(y) ? "" : round(y;digits=1),
			color=:dodgerblue, strokecolor=:black, strokewidth=1, gap=0,
			label_size=14f0)

		# Add design annotations
		for (idx,bar) in enumerate(bars)
			iszero(bar.effcap) && continue
			text!(ax, bar.desstr; position=(idx,bar.effcap), align=(:center,:top),
				offset=(0,-4), textsize=14f0, color=:white)
			desstrlines = isempty(bar.desstr) ? 0 : 1 + count('\n', bar.desstr)
			desstrlines+1 < 6*bar.effcap/maximum(getindex.(bars,:effcap)) &&
				text!(ax, bar.bstr; position=(idx,0), align=(:center,:bottom),
					offset=(0,2), textsize=12f0)
		end

		# Subfig label
		Label(selfig[Tuple(figidx)...,TopLeft()], string(figlbl);
			textsize=28f0, halign=:left, padding=(0,16,8,0))
	end

	## Grid of best methods
	ax = Axis(row1[1,2];
		title = "Best methods across the grid of resource constraints\n$DAYSTR",
		xlabel = "Daily sample collection capacity (average)",
		xticks = (axes(SWABS,1), string.(SWABS)),
		ylabel = "Daily test capacity (average)",
		yticks = (axes(TESTS,1), string.(TESTS)),
		xgridvisible = false, ygridvisible = false,
		xticksvisible = false, yticksvisible = false,
		aspect = DataAspect(),
	)
	tightlimits!(ax, Left(), Right(), Top(), Bottom())

	# Calculate colors, design strings, etc.
	cells = map(collect(pairs(IndexCartesian(),capsgrid))) do (grididx,caps)
		effcap, method = findmax(cap -> cap.effcap, caps)
		des = caps[method].des
		rect = Rect(grididx[1]-0.5,grididx[2]-0.5,1,1)
		color = METHODCOLOR[method](des)
		desstr = (;
			indiv = _   -> "",
			hyper = des -> des.q == 1 ?
				"n:$(convert(Int,des.n/des.m))\nm:1" : "n:$(des.n)\nm:$(des.m)",
			array = dim -> join(dim,'x'),
			pbest = _   -> "P-BEST",
		)[method](des)
		(;rect, color, effcap, desstr)
	end

	# Draw cell boxes
	poly!(ax, vec(getindex.(cells,:rect)); color=vec(getindex.(cells,:color)),
		strokecolor=(:black,0.1), strokewidth=1)

	# Add cell annotations to testing-limited regime
	testlimited = filter(idx->idx[1]>=idx[2],CartesianIndices(cells))
	text!(ax, string.(round.(getindex.(cells[testlimited],:effcap); digits=1));
		position=Point.(Tuple.(testlimited)), align=(:center,:bottom), offset=(0,5),
		textsize=13f0)
	text!(ax, getindex.(cells[testlimited],:desstr);
		position=Point.(Tuple.(testlimited)), align=(:center,:top), offset=(0,7),
		textsize=12f0, color=:white)

	# Add arrows and annotation to testing-rich regime
	ndiag = minimum(size(cells))
	arrows!(ax, 1:ndiag-1, (1:ndiag-1) .+ 0.4, fill(0,ndiag), fill(0.8,ndiag);
		linewidth=2, arrowsize=15)
	text!(ax, "Testing-rich regime:\nindividual testing optimal";
		position=((1+ndiag)/2-2,(1+ndiag)/2+2), align=(:center,:center),
		textsize=24f0, rotation=pi/4)

	# Mark selected cells
	poly!(ax,
		vec([Rect(grididx[1]-0.5,grididx[2]-0.5,1,1) for grididx in SELECTED]);
		color=(:white,0), strokecolor=:black, strokewidth=2.5)

	# Subfig label
	Label(row1[1,2,TopLeft()], string(('a':'z')[length(SELECTED)+1]);
		textsize=28f0, halign=:left)

	## Grid of all methods
	fig[2,1] = fullfig = GridLayout()

	# Prepare axes
	axs = [Axis(fullfig[i,j]) for j in 1:length(SWABS), i in 1:length(TESTS)]
	axs = reverse(axs; dims=2)

	# Fill axes
	for (ax, caps) in zip(axs, capsgrid)
		hidedecorations!(ax)
		tightlimits!(ax, Left(), Right())

		# Compute bar values, labels, etc.
		bars = map(collect(pairs(caps))) do (method,(effcap,bavg,des))
			color = METHODCOLOR[method](des)
			desstr = (;
				indiv = _   -> "Ind.",
				hyper = des -> des.q == 1 ?
					"n:$(convert(Int,des.n/des.m))\nm:1" : "n:$(des.n)\nm:$(des.m)",
				array = dim -> join(dim,'x'),
				pbest = _   -> "P-BEST",
			)[method](des)
			desstrcolor = (;
				indiv = _   -> :black,
				hyper = des -> des.q == 1 ? :black : :white,
				array = _   -> :white,
				pbest = _   -> :white,
			)[method](des)
			(;effcap, color, desstr, desstrcolor)
		end

		# Find best method and mark corresponding effective capacity
		maxbar = argmax(bar->bar.effcap, bars)
		hlines!(ax, [maxbar.effcap]; color=maxbar.color, linewidth=2)
		hlines!(ax, [maxbar.effcap]; color=:lightgray, linewidth=2, linestyle=:dash)
		text!(ax, "$(round(maxbar.effcap;digits=1))"; textsize=14f0,
			position=(mean(1:length(bars)),maxbar.effcap), align=(:center,:bottom))

		# Draw barplot
		barplot!(ax, getindex.(bars,:effcap); color=getindex.(bars,:color),
			strokecolor=:black, strokewidth=1, gap=0)

		# Add design annotations
		for (idx,bar) in enumerate(bars)
			iszero(bar.effcap) && continue
			textsize = 11f0 * min(1,
				6/maximum(length.(split(bar.desstr,'\n')))*bar.effcap/maxbar.effcap)
			textsize > 2f0 &&
				text!(ax, bar.desstr; position=(idx,bar.effcap/2), rotation=pi/2,
					align=(:center,:center), textsize, color=bar.desstrcolor)
		end

		ylims!(ax, (0,1.35*maxbar.effcap))
	end
	colgap!(fullfig, 0)
	rowgap!(fullfig, 0)

	# Add labels
	for (swabidx,swabs) in enumerate(SWABS)
		Label(fullfig[end,swabidx,Bottom()], "$swabs";
			padding=(0,0,0,8), textsize=16f0)
	end
	for (testidx,tests) in enumerate(TESTS)
		Label(fullfig[end-testidx+1,1,Left()], "$tests";
			padding=(0,8,0,0), textsize=16f0)
	end
	Label(fullfig[0,:], "Effective screening capacity for all methods \
		across the grid of resource constraints $DAYSTR"; textsize=16f0)
	Label(fullfig[end+1,:], "Daily sample collection budget (average)";
		textsize=16f0)
	Label(fullfig[2:end-1,0], "Daily test budget (average)";
		textsize=16f0, rotation=pi/2)
	colgap!(fullfig, 1, 8)
	rowgap!(fullfig, 1, 8)
	rowgap!(fullfig, length(TESTS)+1, 8)

	# Subfig label
	Label(fullfig[1,1], string(('a':'z')[length(SELECTED)+2]);
		textsize=28f0, halign=:left, valign=:top, padding=(0,0,0,0))

	## Legend
	LAYOUT = [:indiv, :hyper1, :hyper2, :hyper3, :array, :pbest]
	LEG = [LEGEND[m][1]=>rstrip(LEGEND[m][2],'*') for m in LAYOUT]
	Legend(fig[3,1], [PolyElement(;color) for (color,_) in LEG], last.(vec(LEG)),
		"Methods"; titlefont=assetpath("fonts", "NotoSans-Bold.ttf"),
		orientation=:horizontal, titleposition=:left,
		tellheight=true, tellwidth=false, framevisible=false,
		titlegap=36f0, colgap=24f0)

	## Layout sizes/gaps
	colgap!(selfig, 16)
	rowgap!(selfig, 0)
	colsize!(row1, 2, Fixed(500))
	colgap!(row1, 40)
	rowsize!(row1, 1, Fixed(500))
	rowgap!(fig.layout, 12)

	save(SAVENAME, fig)
	fig
end

# ╔═╡ Cell order:
# ╟─849adbca-08e5-11eb-1956-cd02f2a6be7a
# ╠═95eaf04a-08e5-11eb-043f-27646f1bf207
# ╠═b58f8921-cd4c-4740-b17e-fc669a999b84
# ╟─210d331a-09c3-11eb-06ed-e1f4ddc023fc
# ╟─a61d6b1b-a795-4883-b0be-55184cc79e58
# ╟─fd643c2c-08e5-11eb-0d25-014798a505d1
# ╟─fe7cf842-08e5-11eb-0d07-41ed2902921a
# ╟─3f1b753e-0901-11eb-0060-1f1269588bef
# ╟─3f1c5288-0901-11eb-391f-dd7dd02d3d9b
# ╟─05c6923e-08e6-11eb-1d9b-21afc86a6b04
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
# ╟─72f1ae0a-544c-11eb-0a55-f18f192d64c2
# ╟─21073ca8-3684-11eb-0a25-8d16a9799cf1
# ╟─1eb03185-c9c9-4ec7-aa2f-bccbcb32fc9c
# ╟─d56d4a09-6498-40f3-ac6e-ad7be1e7dbf4
# ╟─794090d7-65fe-41d3-9705-2f32a1fd92b2
# ╟─ba58813a-818e-42fa-a177-bda658451233
# ╟─088ae14a-efa0-4428-9e09-204253c671c3
# ╟─f5aafd70-e89f-4802-9ec7-d0ac81a17cd2
# ╟─4308103a-bf41-48f3-bdcb-d2bcf9cf9984
# ╟─5b786730-e30e-4d6e-9229-d9087a08e5c3
# ╟─6d6411bb-2910-4c1c-82eb-154093782d33
# ╟─0c901f1d-7da4-4688-b9a6-f5ddaa85ffc7
# ╟─26c3d755-a2f6-40ef-9f8d-6e8942d4bd3f
