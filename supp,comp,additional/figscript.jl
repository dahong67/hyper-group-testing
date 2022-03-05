### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ 95eaf04a-08e5-11eb-043f-27646f1bf207
begin
	import Pkg; Pkg.activate(@__DIR__)
	using Dictionaries, IdentityRanges, NPZ, OffsetArrays
	using CairoMakie, Colors, LsqFit, Statistics

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

# ╔═╡ 1c2f2782-f8cf-4f73-9881-6131b5e29c96
ADDDIR = "../../../hyper-group-testing/supp,comp,additional/designs/"

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

# ╔═╡ 7168f7dc-0377-4db0-8451-60645997827c
md"""
### Random assignment
"""

# ╔═╡ 2ac27e64-3c3a-46b5-af55-a177d5a37347
randassigneff = map(x->x[DAYS], loadeff(joinpath(ADDDIR,"rand,assign")))

# ╔═╡ 805bcefe-0e95-48d0-b902-bad29b7ce5a3
randassignsens = map(x->x[DAYS], loadsens(joinpath(ADDDIR,"rand,assign")))

# ╔═╡ 12f80cab-5feb-40e0-8653-037abd7c99a0
md"""
### Double-pooling
"""

# ╔═╡ 4377e981-8748-49ea-af76-34025e3030a2
doublepooleff = map(x->x[DAYS], loadeff(joinpath(ADDDIR,"double,pooling")))

# ╔═╡ e1c53a04-1625-4fa8-9372-de0bb3cba81c
doublepoolsens = map(x->x[DAYS], loadsens(joinpath(ADDDIR,"double,pooling")))

# ╔═╡ 2d9ac51e-50b6-445d-a2ff-4993d2d987e1
md"""
### Balanced square arrays with holes
"""

# ╔═╡ 9ef725c0-370b-4781-bbd9-c7547576cd51
balsqarrayeff = map(x->x[DAYS], loadeff(joinpath(ADDDIR,"array,holes")))

# ╔═╡ 8a7287e8-941f-448c-8e71-d9e7c01305f0
balsqarraysens = map(x->x[DAYS], loadsens(joinpath(ADDDIR,"array,holes")))

# ╔═╡ a92db3c4-5436-11eb-181b-3f1a71204b12
md"""
### Individual testing
"""

# ╔═╡ 72f1ae0a-544c-11eb-0a55-f18f192d64c2
indiveff = OffsetVector(npzread(joinpath(DESDIR,"individual","Eff_avg.npy")),-1)[DAYS]

# ╔═╡ 21073ca8-3684-11eb-0a25-8d16a9799cf1
indivsens = OffsetVector(npzread(joinpath(DESDIR,"individual","Recall_combined.npy")),-1)[DAYS]

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
SHADELIMS = [40,90]

# ╔═╡ c8c90422-4f35-4e90-90b8-4b30f3738f44
SHADECOLOR = (:gray87, 0.25)

# ╔═╡ cf761484-3d31-4294-bfc6-a1f897009a02
randassignstr(n,m,q) = "Random assignment (m=$m, q=$q)"

# ╔═╡ 4e5c92e8-059a-483b-97eb-038016cb39c9
doublepoolstr(n,m,q) = q == 2 ? "Double-pooling (m=$m)" : throw("q == $q != 2")

# ╔═╡ 1b91e41f-e346-4ea5-a6ef-f49d51634d31
function balsqarraystr(n,m,q)
	q == 2 || throw("q = $q != 2")
	iseven(m) || throw("m = $m is not even")
	platearraystr = Dict((96,20,2)=>"8x12",(384,40,2)=>"16x24")[(n,m,q)]
	return "$(m÷2)x$(m÷2) array with $((m÷2)^2-n) holes\n\
		(balanced variant of $platearraystr array)"
end

# ╔═╡ 640b0974-da62-4cf7-963b-121a5030af16
hyperstr(n,m,q) = "HYPER (m=$m, q=$q)"

# ╔═╡ fbab58ba-b9bf-4f82-bdd9-7b5a2a4ef5c2
begin
	mm_to_units(mm) = floor(Int,mm/25.4*72/0.75)
	mm_to_units(mm1,mm2,mmrest...) = mm_to_units.((mm1,mm2,mmrest...))
end

# ╔═╡ 2dff4440-1131-4e39-b55b-3dc0d46967e8
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

# ╔═╡ 27882286-ba06-4d0c-9876-ab8afe757d33
with_theme(THEME) do
	fig = Figure(; resolution=mm_to_units(180,75))

	# Plot styles
	dcolors = distinguishable_colors(7, [RGB(1,1,1)], dropseed=true)
	indivstyle = (;color=dcolors[1], linestyle=:dot)
	randassignstyle = (;color=dcolors[2], linewidth=3)
	doublepoolstyle = (;color=dcolors[6])
	balsqarraystyle = (;color=dcolors[3])
	hyperstyle = (n,m,q) -> (;
		color=get(Dict(16=>dcolors[4],32=>dcolors[5]),m,:black),
		linestyle=get(Dict(1=>:dot,2=>:dash,3=>:solid),q,:solid),
	)
	
	# Fig a
	fig[1,1] = GridLayout()
	ax_eff = Axis(fig[1,1][1,1]; title="n=96 individuals",
		limits=(nothing,(0,6*1.14)), yticks=[0,1,96/20,96/16])
	vspan!(ax_eff, SHADELIMS...; color=SHADECOLOR)
	lines!(ax_eff, randassigneff[(n=96,m=16,q=2)]; randassignstyle...)
	lines!(ax_eff, doublepooleff[(n=96,m=16,q=2)]; doublepoolstyle...)
	lines!(ax_eff, balsqarrayeff[(n=96,m=20,q=2)]; balsqarraystyle...)
	lines!(ax_eff, hypereff[(n=96,m=16,q=2)];      hyperstyle(96,16,2)...)
	lines!(ax_eff, indiveff;                       indivstyle...)

	ax_sens = Axis(fig[1,1][2,1])
	vspan!(ax_sens, SHADELIMS...; color=SHADECOLOR)
	plotsens!(ax_sens, randassignsens[(n=96,m=16,q=2)]; randassignstyle...)
	plotsens!(ax_sens, doublepoolsens[(n=96,m=16,q=2)]; doublepoolstyle...)
	plotsens!(ax_sens, balsqarraysens[(n=96,m=20,q=2)]; balsqarraystyle...)
	plotsens!(ax_sens, hypersens[(n=96,m=16,q=2)];      hyperstyle(96,16,2)...)
	plotsens!(ax_sens, indivsens;                       indivstyle...)

	text!(ax_eff, randassignstr(96,16,2); color=randassignstyle.color,
		position=(20,96/16), offset=(4,2))
	text!(ax_eff, doublepoolstr(96,16,2); color=doublepoolstyle.color,
		position=(20,96/16), offset=(4,-2), align=(:left,:top))
	text!(ax_eff, balsqarraystr(96,20,2); color=balsqarraystyle.color,
		position=(20,96/20), offset=(4,-2), align=(:left,:top))
	text!(ax_eff, hyperstr(96,16,2); color=hyperstyle(96,16,2).color,
		position=(83.0,5.1), offset=(4,2))
	text!(ax_eff, "Individual testing"; color=indivstyle.color,
		position=(20,1), offset=(4,2))

	# Fig b
	fig[1,2] = GridLayout()
	ax_eff = Axis(fig[1,2][1,1]; title="n=384 individuals",
		limits=(nothing,(0,12*1.14)), yticks=[0,1,384/40,384/32])
	vspan!(ax_eff, SHADELIMS...; color=SHADECOLOR)
	lines!(ax_eff, randassigneff[(n=384,m=32,q=2)]; randassignstyle...)
	lines!(ax_eff, doublepooleff[(n=384,m=32,q=2)]; doublepoolstyle...)
	lines!(ax_eff, balsqarrayeff[(n=384,m=40,q=2)]; balsqarraystyle...)
	lines!(ax_eff, hypereff[(n=384,m=32,q=2)];      hyperstyle(384,32,2)...)
	lines!(ax_eff, indiveff;                        indivstyle...)

	ax_sens = Axis(fig[1,2][2,1])
	vspan!(ax_sens, SHADELIMS...; color=SHADECOLOR)
	plotsens!(ax_sens, randassignsens[(n=384,m=32,q=2)]; randassignstyle...)
	plotsens!(ax_sens, doublepoolsens[(n=384,m=32,q=2)]; doublepoolstyle...)
	plotsens!(ax_sens, balsqarraysens[(n=384,m=40,q=2)]; balsqarraystyle...)
	plotsens!(ax_sens, hypersens[(n=384,m=32,q=2)];      hyperstyle(384,32,2)...)
	plotsens!(ax_sens, indivsens;                        indivstyle...)

	text!(ax_eff, randassignstr(384,32,2); color=randassignstyle.color,
		position=(20,384/32), offset=(4,2))
	text!(ax_eff, doublepoolstr(384,32,2); color=doublepoolstyle.color,
		position=(20,384/32), offset=(4,-2), align=(:left,:top))
	text!(ax_eff, balsqarraystr(384,40,2); color=balsqarraystyle.color,
		position=(20,384/40), offset=(4,-2), align=(:left,:top))
	text!(ax_eff, hyperstr(384,32,2); color=hyperstyle(384,32,2).color,
		position=(72.5,10.2), offset=(4,2))
	text!(ax_eff, "Individual testing"; color=indivstyle.color,
		position=(20,1), offset=(4,2))

	# Common axis limits, ticks, etc.
	for (gl,sub) in zip(contents(fig.layout),'a':'z')
		ax_eff, ax_sens = contents(gl)
		rowgap!(gl, 8)

		# x-axis
		xlims!(ax_eff, extrema(DAYS))
		linkxaxes!(ax_eff,ax_sens)
		hidexdecorations!(ax_eff; ticks=false, grid=false)
		ax_eff.xticks = ax_sens.xticks = 30:10:100
		daystr = day -> "$day\n($(round(100*prevalence[day];digits=2))%)"
		ax_sens.xtickformat = x->daystr.(convert.(Int,x))
		ax_sens.xlabel = "Day (Prevalence)"

		# y-axis
		ylims!(ax_sens, (0.6,0.9))
		ax_sens.yticks = 0.65:0.05:0.85
		ax_sens.ytickformat = y->string.(convert.(Int,100*y),'%')
		ax_eff.ylabel = "Efficiency\n(relative to individual testing)"
		ax_sens.ylabel = "Sensitivity"

		# label
		Label(gl[1,1,TopLeft()], string(sub);
			halign=:left, valign=:top, font="Arial Bold")
	end

	colgap!(fig.layout, 16)
	save("fig-s5.png", fig; px_per_unit=2)
	save("fig-s5.pdf", fig)
	fig
end

# ╔═╡ Cell order:
# ╟─849adbca-08e5-11eb-1956-cd02f2a6be7a
# ╠═95eaf04a-08e5-11eb-043f-27646f1bf207
# ╟─210d331a-09c3-11eb-06ed-e1f4ddc023fc
# ╟─66dd9e7e-f948-4457-8bd6-45345be6abb5
# ╟─1c2f2782-f8cf-4f73-9881-6131b5e29c96
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
# ╟─7168f7dc-0377-4db0-8451-60645997827c
# ╟─2ac27e64-3c3a-46b5-af55-a177d5a37347
# ╟─805bcefe-0e95-48d0-b902-bad29b7ce5a3
# ╟─12f80cab-5feb-40e0-8653-037abd7c99a0
# ╟─4377e981-8748-49ea-af76-34025e3030a2
# ╟─e1c53a04-1625-4fa8-9372-de0bb3cba81c
# ╟─2d9ac51e-50b6-445d-a2ff-4993d2d987e1
# ╟─9ef725c0-370b-4781-bbd9-c7547576cd51
# ╟─8a7287e8-941f-448c-8e71-d9e7c01305f0
# ╟─a92db3c4-5436-11eb-181b-3f1a71204b12
# ╟─72f1ae0a-544c-11eb-0a55-f18f192d64c2
# ╟─21073ca8-3684-11eb-0a25-8d16a9799cf1
# ╟─0fd98c50-2203-11eb-0f58-d946dcf484b4
# ╟─581cfcf4-16e5-4662-80b6-54e202d9d1a6
# ╟─222d4b16-7b8c-11eb-0622-797f8f7bba13
# ╟─39c78658-ad44-444c-9808-c5011392d3e2
# ╟─a6c3504c-543f-11eb-2f20-791273bc06b8
# ╟─c8c90422-4f35-4e90-90b8-4b30f3738f44
# ╟─cf761484-3d31-4294-bfc6-a1f897009a02
# ╟─4e5c92e8-059a-483b-97eb-038016cb39c9
# ╟─1b91e41f-e346-4ea5-a6ef-f49d51634d31
# ╟─640b0974-da62-4cf7-963b-121a5030af16
# ╟─fbab58ba-b9bf-4f82-bdd9-7b5a2a4ef5c2
# ╟─2dff4440-1131-4e39-b55b-3dc0d46967e8
# ╟─27882286-ba06-4d0c-9876-ab8afe757d33
