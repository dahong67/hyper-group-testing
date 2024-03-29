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
DESDIR = relpath(joinpath(@__DIR__,"..","simulation","designs"))

# ╔═╡ 301febdd-5ae1-4b4d-96cb-d4f7ea880d84
ADDDIR = relpath(joinpath(@__DIR__,"designs"))

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

# ╔═╡ 8c310c4a-9d6f-4dcb-9803-621e9f6c9697
md"""
### Balanced arrays
"""

# ╔═╡ 9d87fffc-0373-4d14-8c53-0d34c3236009
balarrayeff = map(x->x[DAYS], loadeff(joinpath(ADDDIR,"balarray")))

# ╔═╡ 2e0cc554-08e6-11eb-0186-bba8c2aee7ec
balarraysens = map(x->x[DAYS], loadsens(joinpath(ADDDIR,"balarray")))

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

# ╔═╡ 640b0974-da62-4cf7-963b-121a5030af16
hyperstr(n,m,q) = "HYPER (m=$m, q=$q)"

# ╔═╡ fdac5fcd-83f7-440b-af81-ac8658911048
balarraystr(n,m,q) = "Balanced array (m=$m, q=$q)"

# ╔═╡ 8939b5c4-cdf7-456e-b213-072d9e16b9b7
begin
	mm_to_units(mm) = floor(Int,mm/25.4*72/0.75)
	mm_to_units(mm1,mm2,mmrest...) = mm_to_units.((mm1,mm2,mmrest...))
end

# ╔═╡ 414b2041-d6ff-4c69-9083-d4cc4bb3373b
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
	balarraystyle = (n,m,q) -> (;
		color=get(Dict(8=>dcolors[6],16=>dcolors[3],32=>dcolors[7]),m,:black),
	)
	hyperstyle = (n,m,q) -> (;
		color=get(Dict(8=>dcolors[2],16=>dcolors[4],32=>dcolors[5]),m,:black),
		linestyle=get(Dict(1=>:dot,2=>:dash,3=>:solid),q,:solid),
	)
	
	# Fig a
	fig[1,1] = GridLayout()
	ax_eff = Axis(fig[1,1][1,1]; title="n=96 individuals",
		limits=(nothing,(0,12*1.05)), yticks=[0,1,96/16,96/8])
	vspan!(ax_eff, SHADELIMS...; color=SHADECOLOR)
	lines!(ax_eff, balarrayeff[(n=96,m=16,q=2)]; balarraystyle(96,16,2)...)
	lines!(ax_eff, balarrayeff[(n=96,m=8,q=2)];  balarraystyle(96,8,2)...)
	lines!(ax_eff, hypereff[(n=96,m=16,q=2)];    hyperstyle(96,16,2)...)
	lines!(ax_eff, hypereff[(n=96,m=8,q=2)];     hyperstyle(96,8,2)...)
	lines!(ax_eff, indiveff;                     indivstyle...)

	ax_sens = Axis(fig[1,1][2,1])
	vspan!(ax_sens, SHADELIMS...; color=SHADECOLOR)
	plotsens!(ax_sens, balarraysens[(n=96,m=16,q=2)]; balarraystyle(96,16,2)...)
	plotsens!(ax_sens, balarraysens[(n=96,m=8,q=2)];  balarraystyle(96,8,2)...)
	plotsens!(ax_sens, hypersens[(n=96,m=16,q=2)];    hyperstyle(96,16,2)...)
	plotsens!(ax_sens, hypersens[(n=96,m=8,q=2)];     hyperstyle(96,8,2)...)
	plotsens!(ax_sens, indivsens;                     indivstyle...)

	text!(ax_eff, balarraystr(96,16,2); color=balarraystyle(96,16,2).color,
		position=(20,96/16), offset=(4,-2), align=(:left,:top))
	text!(ax_eff, balarraystr(96,8,2); color=balarraystyle(96,8,2).color,
		position=(20,96/8), offset=(4,-4), align=(:left,:top))
	text!(ax_eff, hyperstr(96,16,2); color=hyperstyle(96,16,2).color,
		position=(20,96/16), offset=(4,2))
	text!(ax_eff, hyperstr(96,8,2); color=hyperstyle(96,8,2).color,
		position=(71.5,9.5), offset=(4,2))
	text!(ax_eff, "Individual testing"; color=indivstyle.color,
		position=(20,1), offset=(4,2))

	# Fig b
	fig[1,2] = GridLayout()
	ax_eff = Axis(fig[1,2][1,1]; title="n=384 individuals",
		limits=(nothing,(0,48*1.05)), yticks=[1,384/32,384/8])
	vspan!(ax_eff, SHADELIMS...; color=SHADECOLOR)
	lines!(ax_eff, balarrayeff[(n=384,m=32,q=2)]; balarraystyle(384,32,2)...)
	lines!(ax_eff, balarrayeff[(n=384,m=8,q=2)];  balarraystyle(384,8,2)...)
	lines!(ax_eff, hypereff[(n=384,m=32,q=2)];    hyperstyle(384,32,2)...)
	lines!(ax_eff, hypereff[(n=384,m=8,q=2)];     hyperstyle(384,8,2)...)
	lines!(ax_eff, indiveff;                      indivstyle...)

	ax_sens = Axis(fig[1,2][2,1])
	vspan!(ax_sens, SHADELIMS...; color=SHADECOLOR)
	plotsens!(ax_sens, balarraysens[(n=384,m=32,q=2)]; balarraystyle(384,32,2)...)
	plotsens!(ax_sens, balarraysens[(n=384,m=8,q=2)];  balarraystyle(384,8,2)...)
	plotsens!(ax_sens, hypersens[(n=384,m=32,q=2)];    hyperstyle(384,32,2)...)
	plotsens!(ax_sens, hypersens[(n=384,m=8,q=2)];     hyperstyle(384,8,2)...)
	plotsens!(ax_sens, indivsens;                      indivstyle...)

	text!(ax_eff, balarraystr(384,32,2); color=balarraystyle(384,32,2).color,
		position=(20,384/32), offset=(4,-2), align=(:left,:top))
	text!(ax_eff, balarraystr(384,8,2); color=balarraystyle(384,8,2).color,
		position=(20,384/8), offset=(8,-8), align=(:left,:top), rotation=-0.34)
	text!(ax_eff, hyperstr(384,32,2); color=hyperstyle(384,32,2).color,
		position=(20,384/32), offset=(4,2))
	text!(ax_eff, hyperstr(384,8,2); color=hyperstyle(384,8,2).color,
		position=(43.5,38.0), offset=(4,2))
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
	save("fig-s6ab.png", fig; px_per_unit=2)
	save("fig-s6ab.pdf", fig)
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
# ╟─8c310c4a-9d6f-4dcb-9803-621e9f6c9697
# ╟─9d87fffc-0373-4d14-8c53-0d34c3236009
# ╟─2e0cc554-08e6-11eb-0186-bba8c2aee7ec
# ╟─a92db3c4-5436-11eb-181b-3f1a71204b12
# ╟─72f1ae0a-544c-11eb-0a55-f18f192d64c2
# ╟─21073ca8-3684-11eb-0a25-8d16a9799cf1
# ╟─0fd98c50-2203-11eb-0f58-d946dcf484b4
# ╟─581cfcf4-16e5-4662-80b6-54e202d9d1a6
# ╟─222d4b16-7b8c-11eb-0622-797f8f7bba13
# ╟─39c78658-ad44-444c-9808-c5011392d3e2
# ╟─a6c3504c-543f-11eb-2f20-791273bc06b8
# ╟─c8c90422-4f35-4e90-90b8-4b30f3738f44
# ╟─640b0974-da62-4cf7-963b-121a5030af16
# ╟─fdac5fcd-83f7-440b-af81-ac8658911048
# ╟─8939b5c4-cdf7-456e-b213-072d9e16b9b7
# ╟─414b2041-d6ff-4c69-9083-d4cc4bb3373b
# ╟─27882286-ba06-4d0c-9876-ab8afe757d33
