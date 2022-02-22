### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ 95eaf04a-08e5-11eb-043f-27646f1bf207
begin
	import Pkg; Pkg.activate(@__DIR__)
	using Dictionaries, IdentityRanges, NPZ, OffsetArrays
	using CairoMakie, Colors, ColorSchemes, LsqFit, Statistics

	Makie.convert_arguments(P::PointBased, a::OffsetVector) =
		convert_arguments(P,collect(axes(a,1)),collect(a))
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
DESDIR = "../../../hyper-group-testing/simulation/designs/"

# ╔═╡ 301febdd-5ae1-4b4d-96cb-d4f7ea880d84
ADDDIR = "../../../hyper-group-testing/supp,dd,decoder/designs/"

# ╔═╡ 6b39ad12-5378-47dc-9a12-e1d5838b0dee
DAYS = IdentityRange(20:110)

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

# ╔═╡ 60bae43b-1ecf-4206-ae1f-5d8ee11a0ba4
md"""
### HYPER-DD-Decode
"""

# ╔═╡ f0bbadb2-9a72-4395-bf52-4d3adb694e4e
hyperdddeff = map(x->x[DAYS], loadeff(joinpath(ADDDIR,"hypergraph,dd,decode")))

# ╔═╡ 23c277b1-e68d-42ca-adb7-1a0c7c223dfc
hyperdddsens = map(x->x[DAYS], loadsens(joinpath(ADDDIR,"hypergraph,dd,decode")))

# ╔═╡ 95f938e5-1264-4e65-bed9-51fc1d57037e
md"""
### HYPER-DD-Skip
"""

# ╔═╡ 760f038e-cc3c-4a67-9c74-1a9c2e65b23c
hyperddseff = map(x->x[DAYS], loadeff(joinpath(ADDDIR,"hypergraph,dd,skip")))

# ╔═╡ 83610001-2459-4dec-a40c-345fa1b36dce
hyperddssens = map(x->x[DAYS], loadsens(joinpath(ADDDIR,"hypergraph,dd,skip")))

# ╔═╡ a92db3c4-5436-11eb-181b-3f1a71204b12
md"""
### Individual testing
"""

# ╔═╡ 72f1ae0a-544c-11eb-0a55-f18f192d64c2
indiveff = OffsetVector(npzread(joinpath(DESDIR,"individual","Eff_avg.npy")),-1)[DAYS]

# ╔═╡ 21073ca8-3684-11eb-0a25-8d16a9799cf1
indivsens = OffsetVector(npzread(joinpath(DESDIR,"individual","Recall_combined.npy")),-1)[DAYS]

# ╔═╡ 18b8b0a7-d969-4e8f-873c-d70ed9441889
md"""
## Resource analysis
"""

# ╔═╡ 95b8e800-8548-42b6-923e-91cb48c741d9
RESDAYS = IdentityRange(40:90)

# ╔═╡ 3cac1e88-a00e-4617-afa2-4da569b4da9c
SWABS = [12*2^j for j in 0:9]

# ╔═╡ 4da0f571-e58e-456c-8fbe-649e72d2a296
TESTS = SWABS

# ╔═╡ f63a10b3-1320-46ae-ad92-9a53a5beb41b
function effcapacity(sampbudget::Number, testbudget::Number, n::Number, eff, sens)
	ntests = n./eff                             # tests/batch (for both stages)
	b = min.(sampbudget/n, testbudget./ntests)  # batches/day
	bavg = mean(b)
	effcap = bavg > 0.9 ? mean(n.*b.*sens) : 0.0
	return (;effcap, bavg)
end

# ╔═╡ ad1555db-96c9-412f-9ffb-4d2f7427bbcf
capsgrid = map(Iterators.product(SWABS,TESTS)) do budget
	# Individual testing
	indiv = (;
		effcapacity(budget...,1,indiveff[RESDAYS], indivsens[RESDAYS])...,
		des=()
	)

	# HYPER (optimized across designs)
	hypercap, hyperdes = findmax(keys(hypereff)) do des
		effcapacity(budget..., des.n,
			hypereff[des][RESDAYS], hypersens[des][RESDAYS])
	end
	hyper = (hypercap...,des=hyperdes)

	# HYPER-DD-Decode (optimized across designs)
	hyperdddcap, hyperddddes = findmax(keys(hyperdddeff)) do des
		effcapacity(budget..., des.n,
			hyperdddeff[des][RESDAYS], hyperdddsens[des][RESDAYS])
	end
	hyperddd = (hyperdddcap...,des=hyperddddes)
	
	# HYPER-DD-Skip (optimized across designs)
	hyperddscap, hyperddsdes = findmax(keys(hyperddseff)) do des
		effcapacity(budget..., des.n,
			hyperddseff[des][RESDAYS], hyperddssens[des][RESDAYS])
	end
	hyperdds = (hyperddscap...,des=hyperddsdes)
	
	return (;indiv, hyper, hyperdds, hyperddd)
end

# ╔═╡ 0fd98c50-2203-11eb-0f58-d946dcf484b4
md"""
## Figures
"""

# ╔═╡ 581cfcf4-16e5-4662-80b6-54e202d9d1a6
# polynomial fit
function polyfit(xx,yy,order)
	axes(xx) == axes(yy) || throw("Arrays must have same axes. Got: $(axes(xx)) and $(axes(yy))")
	x, y = OffsetArrays.no_offset_view(xx), OffsetArrays.no_offset_view(yy)
	
	fit_func = (x,p) -> evalpoly.(x,Ref(p))
	popt = curve_fit(fit_func,x,y,ones(order)).param
	return fit_func(xx,popt)
end

# ╔═╡ 222d4b16-7b8c-11eb-0622-797f8f7bba13
POLYORDER = 8

# ╔═╡ 39c78658-ad44-444c-9808-c5011392d3e2
function plotsens!(ax,sens; prev=prevalence, days=DAYS, order=POLYORDER, style...)
	scatter!(ax, sens; color=(get(style,:color,:black),0.6))
	lines!(ax, polyfit(log10.(prev[days]),sens[days],order); style...)
end

# ╔═╡ a6c3504c-543f-11eb-2f20-791273bc06b8
SHADELIMS = extrema(RESDAYS)

# ╔═╡ c8c90422-4f35-4e90-90b8-4b30f3738f44
SHADECOLOR = (:gray87, 0.25)

# ╔═╡ 640b0974-da62-4cf7-963b-121a5030af16
desstr(n,m,q) = "n=$n individual, m=$m pools, q=$q splits"

# ╔═╡ 765191ba-fa31-4062-a1aa-bfcf85298791
LEGEND = (;
	indiv     = ColorSchemes.OrRd_4[1]   => "Individual testing",
	hyper1    = ColorSchemes.OrRd_4[2]   => "HYPER (q=1)",
	hyper2    = ColorSchemes.OrRd_4[3]   => "HYPER (q=2)",
	hyper3    = ColorSchemes.OrRd_4[4]   => "HYPER (q=3)",
	hyperdds1 = colorant"plum2"          => "HYPER-DD-Skip (q=1)",
	hyperdds2 = colorant"darkorchid2"    => "HYPER-DD-Skip (q=2)",
	hyperdds3 = colorant"purple4"        => "HYPER-DD-Skip (q=3)",
	hyperddd1 = colorant"palegreen2"     => "HYPER-DD-Decode (q=1)",
	hyperddd2 = colorant"mediumseagreen" => "HYPER-DD-Decode (q=2)",
	hyperddd3 = colorant"darkgreen"      => "HYPER-DD-Decode (q=3)",
)

# ╔═╡ eb0744cf-bd7a-4fb0-91ba-0eef39dc014d
METHODCOLOR = (;
	indiv = _   -> LEGEND[:indiv][1],
	hyper = des -> LEGEND[Symbol("hyper$(des.q)")][1],
	hyperdds = des -> LEGEND[Symbol("hyperdds$(des.q)")][1],
	hyperddd = des -> LEGEND[Symbol("hyperddd$(des.q)")][1],
)

# ╔═╡ 7e7d9de4-c3fc-439b-9b1b-4bab84f57947
DAYSTR = if length(RESDAYS) == 1
	d = only(RESDAYS)
	p = round(prevalence[d].*100; digits=2)
	"(day $d: prevalence of $p%)"
else
	d1, d2 = extrema(RESDAYS)
	p1, p2 = round.(prevalence[[d1,d2]].*100; digits=2)
	"(days $d1-$d2: prevalence grows from $p1% to $p2%)"
end

# ╔═╡ 27882286-ba06-4d0c-9876-ab8afe757d33
with_theme(; linewidth=3, markersize=3,
	Axis=(;xtickalign=1, ytickalign=1, xticklabelsize=12f0, yticklabelsize=12f0),
	Text=(;textsize=14f0),
) do
	fig = Figure(; resolution=(1200,1200))

	## Figs a-b
	fig[1,1] = row1 = GridLayout(;alignmode=Outside())
	
	# Plot styles
	indivstyle = (;color=colorant"black", linestyle=:dot)
	hyperstyle = (;color=ColorSchemes.OrRd_4[3])
	hyperdddstyle = (;color=colorant"mediumseagreen", linestyle=:dash)
	hyperddsstyle = (;color=colorant"darkorchid2", linestyle=:dash)
	
	# Fig a
	DES = (;n=96,m=16,q=2)
	row1[1,1] = GridLayout()
	ax_eff = Axis(row1[1,1][1,1]; title=desstr(DES...),
		limits=(nothing,(0,DES.n/DES.m*1.25)), yticks=[0,1,DES.n/DES.m])
	vspan!(ax_eff, SHADELIMS...; color=SHADECOLOR)
	lines!(ax_eff, hypereff[DES];    hyperstyle...)
	lines!(ax_eff, hyperdddeff[DES]; hyperdddstyle...)
	lines!(ax_eff, hyperddseff[DES]; hyperddsstyle...)
	lines!(ax_eff, indiveff;         indivstyle...)

	ax_sens = Axis(row1[1,1][2,1])
	vspan!(ax_sens, SHADELIMS...; color=SHADECOLOR)
	plotsens!(ax_sens, hypersens[DES];    hyperstyle...)
	plotsens!(ax_sens, hyperdddsens[DES]; hyperdddstyle...)
	plotsens!(ax_sens, hyperddssens[DES]; hyperddsstyle...)
	plotsens!(ax_sens, indivsens;         indivstyle...)

	text!(ax_eff, "HYPER"; color=hyperstyle.color,
		position=(20,DES.n/DES.m), offset=(8,4))
	text!(ax_eff, "HYPER-DD-Decode"; color=hyperdddstyle.color,
		position=(80,DES.n/DES.m), offset=(8,4))
	text!(ax_eff, "HYPER-DD-Skip"; color=hyperddsstyle.color,
		position=(20,DES.n/DES.m), offset=(8,-4), align=(:left,:top))
	text!(ax_eff, "Individual testing"; color=indivstyle.color,
		position=(20,1), offset=(8,4))

	# Fig b
	DES = (;n=384,m=32,q=2)
	row1[1,2] = GridLayout()
	ax_eff = Axis(row1[1,2][1,1]; title=desstr(DES...),
		limits=(nothing,(0,DES.n/DES.m*1.25)), yticks=[1,DES.n/DES.m])
	vspan!(ax_eff, SHADELIMS...; color=SHADECOLOR)
	lines!(ax_eff, hypereff[DES];    hyperstyle...)
	lines!(ax_eff, hyperdddeff[DES]; hyperdddstyle...)
	lines!(ax_eff, hyperddseff[DES]; hyperddsstyle...)
	lines!(ax_eff, indiveff;         indivstyle...)

	ax_sens = Axis(row1[1,2][2,1])
	vspan!(ax_sens, SHADELIMS...; color=SHADECOLOR)
	plotsens!(ax_sens, hypersens[DES];    hyperstyle...)
	plotsens!(ax_sens, hyperdddsens[DES]; hyperdddstyle...)
	plotsens!(ax_sens, hyperddssens[DES]; hyperddsstyle...)
	plotsens!(ax_sens, indivsens;         indivstyle...)

	text!(ax_eff, "HYPER"; color=hyperstyle.color,
		position=(20,DES.n/DES.m), offset=(8,4))
	text!(ax_eff, "HYPER-DD-Decode"; color=hyperdddstyle.color,
		position=(80,DES.n/DES.m), offset=(8,4))
	text!(ax_eff, "HYPER-DD-Skip"; color=hyperddsstyle.color,
		position=(20,DES.n/DES.m), offset=(8,-4), align=(:left,:top))
	text!(ax_eff, "Individual testing"; color=indivstyle.color,
		position=(20,1), offset=(8,4))
	
	# Common axis limits, ticks, etc.
	for (gl,sub) in zip(contents(row1),'a':'z')
		ax_eff, ax_sens = contents(gl)
		rowgap!(gl, 16)

		# x-axis
		xlims!(ax_eff, extrema(DAYS))
		linkxaxes!(ax_eff, ax_sens)
		hidexdecorations!(ax_eff; ticks=false, grid=false)
		ax_eff.xticks = ax_sens.xticks = 30:10:100
		daystr = day -> "$day\n($(round(100*prevalence[day];digits=2))%)"
		ax_sens.xtickformat = x->daystr.(convert.(Int,x))
		ax_sens.xlabel = "day (prevalence)"

		# y-axis
		ylims!(ax_sens, (0.5,0.9))
		ax_sens.yticks = 0.55:0.05:0.85
		ax_sens.ytickformat = y->string.(convert.(Int,round.(100*y;digits=10)),'%')
		ax_eff.ylabel = "individuals/test"
		ax_sens.ylabel = "sensitivity"

		# label
		Label(gl[1,1,TopLeft()], string(sub); textsize=28f0, halign=:left)
	end

	## Fig c: resource analysis
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
				hyperddd = des -> des.q == 1 ?
					"n:$(convert(Int,des.n/des.m))\nm:1" : "n:$(des.n)\nm:$(des.m)",
				hyperdds = des -> des.q == 1 ?
					"n:$(convert(Int,des.n/des.m))\nm:1" : "n:$(des.n)\nm:$(des.m)",
			)[method](des)
			desstrcolor = (;
				indiv = _   -> :black,
				hyper = des -> des.q == 1 ? :black : :white,
				hyperddd = des -> des.q == 1 ? :black : :white,
				hyperdds = des -> des.q == 1 ? :black : :white,
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
				5/maximum(length.(split(bar.desstr,'\n')))*bar.effcap/maxbar.effcap)
			textsize > 2f0 &&
				text!(ax, bar.desstr; position=(idx,bar.effcap/2), rotation=pi/2,
					align=(:center,:center), textsize, color=bar.desstrcolor)
		end

		ylims!(ax, (0,1.45*maxbar.effcap))
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
	Label(fullfig[1,1], "c";
		textsize=28f0, halign=:left, valign=:top, padding=(0,0,0,0))

	## Legend
	LAYOUT = [:indiv; nothing; nothing;; :hyper1; :hyper2; :hyper3;;
		:hyperdds1; :hyperdds2; :hyperdds3;; :hyperddd1; :hyperddd2; :hyperddd3]
	LEG = [isnothing(m) ? (:white,0.0)=>"" : LEGEND[m][1]=>rstrip(LEGEND[m][2],'*')
			for m in LAYOUT]
	Legend(fig[3,1], [PolyElement(;color) for (color,_) in LEG], last.(vec(LEG)),
		"Bars colored\nby method";
		titlefont=assetpath("fonts", "NotoSans-Bold.ttf"),
		orientation=:horizontal, nbanks=3, titleposition=:left,
		tellheight=true, tellwidth=false, framevisible=false,
		titlegap=36f0, colgap=24f0)

	## Layout
	rowsize!(row1, 1, Fixed(300))
	rowgap!(fig.layout, 2, 10)
	
	save("fig-s21.png", fig)
	fig
end

# ╔═╡ Cell order:
# ╟─849adbca-08e5-11eb-1956-cd02f2a6be7a
# ╠═95eaf04a-08e5-11eb-043f-27646f1bf207
# ╟─210d331a-09c3-11eb-06ed-e1f4ddc023fc
# ╟─66dd9e7e-f948-4457-8bd6-45345be6abb5
# ╟─301febdd-5ae1-4b4d-96cb-d4f7ea880d84
# ╟─6b39ad12-5378-47dc-9a12-e1d5838b0dee
# ╟─fd643c2c-08e5-11eb-0d25-014798a505d1
# ╟─fe7cf842-08e5-11eb-0d07-41ed2902921a
# ╟─3f1b753e-0901-11eb-0060-1f1269588bef
# ╟─3f1c5288-0901-11eb-391f-dd7dd02d3d9b
# ╟─05c6923e-08e6-11eb-1d9b-21afc86a6b04
# ╟─cf1b2b1a-0c3b-11eb-2159-6ffeb682d595
# ╟─122151e0-08e6-11eb-18f5-3b277145bc9d
# ╟─47d9acce-2fa6-11eb-0260-7767c0252b8c
# ╟─f8274a32-2fa6-11eb-3836-9bbe0aaf63cf
# ╟─60bae43b-1ecf-4206-ae1f-5d8ee11a0ba4
# ╟─f0bbadb2-9a72-4395-bf52-4d3adb694e4e
# ╟─23c277b1-e68d-42ca-adb7-1a0c7c223dfc
# ╟─95f938e5-1264-4e65-bed9-51fc1d57037e
# ╟─760f038e-cc3c-4a67-9c74-1a9c2e65b23c
# ╟─83610001-2459-4dec-a40c-345fa1b36dce
# ╟─a92db3c4-5436-11eb-181b-3f1a71204b12
# ╟─72f1ae0a-544c-11eb-0a55-f18f192d64c2
# ╟─21073ca8-3684-11eb-0a25-8d16a9799cf1
# ╟─18b8b0a7-d969-4e8f-873c-d70ed9441889
# ╟─95b8e800-8548-42b6-923e-91cb48c741d9
# ╟─3cac1e88-a00e-4617-afa2-4da569b4da9c
# ╟─4da0f571-e58e-456c-8fbe-649e72d2a296
# ╟─f63a10b3-1320-46ae-ad92-9a53a5beb41b
# ╟─ad1555db-96c9-412f-9ffb-4d2f7427bbcf
# ╟─0fd98c50-2203-11eb-0f58-d946dcf484b4
# ╟─581cfcf4-16e5-4662-80b6-54e202d9d1a6
# ╟─222d4b16-7b8c-11eb-0622-797f8f7bba13
# ╟─39c78658-ad44-444c-9808-c5011392d3e2
# ╟─a6c3504c-543f-11eb-2f20-791273bc06b8
# ╟─c8c90422-4f35-4e90-90b8-4b30f3738f44
# ╟─640b0974-da62-4cf7-963b-121a5030af16
# ╟─765191ba-fa31-4062-a1aa-bfcf85298791
# ╟─eb0744cf-bd7a-4fb0-91ba-0eef39dc014d
# ╟─7e7d9de4-c3fc-439b-9b1b-4bab84f57947
# ╟─27882286-ba06-4d0c-9876-ab8afe757d33
