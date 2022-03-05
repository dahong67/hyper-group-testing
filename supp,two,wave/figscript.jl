### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ 95eaf04a-08e5-11eb-043f-27646f1bf207
begin
	import Pkg; Pkg.activate(@__DIR__)
	using Dictionaries, IdentityRanges, NPZ, OffsetArrays
	using CairoMakie, Colors, Statistics

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
DESDIR = relpath(joinpath(@__DIR__,"designs"))

# ╔═╡ 6b39ad12-5378-47dc-9a12-e1d5838b0dee
DAYS = IdentityRange(10:190)

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

# ╔═╡ dabb2530-e060-4b8c-b03c-1cbeac98ea4b
md"""
### Positive viral loads
"""

# ╔═╡ de3ef15a-9196-4214-b694-6f1448fccf4a
posviralloads = let
	dir = "posviralloads"
	pattern = r"^posviralloads.t-(?<t>\d+).npy$"
	regmatches = filter(!isnothing, match.(pattern, readdir(dir)))
	pairs = map(regmatches) do regmatch
		t = parse(Int,regmatch[:t])
		posloads = npzread(joinpath(dir,regmatch.match))
		t => posloads
	end
	dictionary(sort(pairs))
end

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

# ╔═╡ 0fd98c50-2203-11eb-0f58-d946dcf484b4
md"""
## Figures
"""

# ╔═╡ 6adbc303-83dd-4505-b50d-a94c9c1834bc
HISTBINS = -3:0.1:25

# ╔═╡ 09e96bdc-ebde-4899-a089-6c0418fa7ecc
HISTDAYS = collect(keys(posviralloads))

# ╔═╡ a6c3504c-543f-11eb-2f20-791273bc06b8
SHADELIMS = extrema(HISTDAYS)

# ╔═╡ c8c90422-4f35-4e90-90b8-4b30f3738f44
SHADECOLOR = (:gray87, 0.25)

# ╔═╡ 640b0974-da62-4cf7-963b-121a5030af16
hyperstr(n,m,q) = "HYPER (m=$m, q=$q)"

# ╔═╡ 485f245c-f7ba-4ccb-86f7-08de5e9fd7e9
dorfmanstr(n,m) = string(
	"Dorfman with pools of size $(convert(Int,n//m))",
	"\n[i.e., ", hyperstr(n,m,1), "]"
)

# ╔═╡ 48b55bad-a658-4b22-9927-83e15bf900d1
begin
	mm_to_units(mm) = floor(Int,mm/25.4*72/0.75)
	mm_to_units(mm1,mm2,mmrest...) = mm_to_units.((mm1,mm2,mmrest...))
end

# ╔═╡ 125cf732-35fd-48f9-a2f1-6261d92969fa
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
with_theme(THEME; Legend=(;padding=(0,0,0,0))) do
	fig = Figure(; resolution=mm_to_units(180,150))

	# Plot styles
	dcolors = distinguishable_colors(7, [RGB(1,1,1)], dropseed=true)
	indivstyle = (;color=dcolors[1], linestyle=:dot)
	pbeststyle = (;color=dcolors[2])
	arraystyle = (;color=dcolors[3])
	hyperstyle = (n,m,q) -> (;
		color=get(
			Dict(
				(16,2)=>dcolors[4], (32,2)=>dcolors[5], (12,2)=>dcolors[7],
				(12,3)=>Makie.wong_colors()[1], (12,1)=>dcolors[6],
			),
			(m,q),:black),
		linestyle=get(Dict(1=>:dot,2=>:dash,3=>:solid),q,:solid),
	)
	dorfmanstyle = (;color=dcolors[6])
	
	# Fig a
	ax = Axis(fig[1,1])
	vspan!(ax, SHADELIMS...; color=SHADECOLOR)
	lines!(ax, prevalence; color=Makie.wong_colors()[2])
	ax.limits = (extrema(DAYS),(0,0.02))
	ax.xticks = 10:20:190
	ax.yticks = 0:0.005:0.02
	ax.ytickformat = y->string.(100*y,'%')
	ax.xlabel = "Day"
	ax.ylabel = "Prevalence"
	Label(fig[1,1,TopLeft()], "a"; halign=:left, valign=:top, font="Arial Bold")
	
	# Fig b
	ax = Axis(fig[1,2]; xlabel="log10(viral load)", xticks=-2:2:24,
		limits=(extrema(HISTBINS),nothing))
	for (day,color) in zip(HISTDAYS, Makie.wong_colors(0.5))
		hist!(ax, log10.(posviralloads[day]); bins=HISTBINS, color,
			label="Day $day: prevalence of $(round(100*prevalence[day],digits=1))%")
	end
	axislegend(ax; position=:rt)
	hideydecorations!(ax)
	hidespines!(ax, :l, :t, :r)
	Label(fig[1,2,TopLeft()], "b"; halign=:left, valign=:top, font="Arial Bold")
	
	# Fig c
	fig[2,1] = GridLayout()
	ax_eff = Axis(fig[2,1][1,1]; title="n=96 individuals",
		limits=(nothing,(0,6*1.7)), yticks=[0,1,96/20,96/16])
	vspan!(ax_eff, SHADELIMS...; color=SHADECOLOR)
	lines!(ax_eff, hypereff[(n=96,m=8,q=1)];  dorfmanstyle...)
	lines!(ax_eff, arrayeff[(8,12)];          arraystyle...)
	lines!(ax_eff, hypereff[(n=96,m=16,q=2)]; hyperstyle(96,16,2)...)
	lines!(ax_eff, indiveff;                  indivstyle...)

	ax_sens = Axis(fig[2,1][2,1])
	vspan!(ax_sens, SHADELIMS...; color=SHADECOLOR)
	scatter!(ax_sens, hypersens[(n=96,m=8,q=1)];  dorfmanstyle...)
	scatter!(ax_sens, arraysens[(8,12)];          arraystyle...)
	scatter!(ax_sens, hypersens[(n=96,m=16,q=2)]; hyperstyle(96,16,2)...)
	scatter!(ax_sens, indivsens;                  indivstyle...)

	text!(ax_eff, dorfmanstr(96,8); color=dorfmanstyle.color,
		position=(134,7.5), offset=(0,0), align=(:right,:bottom))
	text!(ax_eff, "8x12 array"; color=arraystyle.color,
		position=(10,96/20), offset=(4,-2), align=(:left,:top))
	text!(ax_eff, hyperstr(96,16,2); color=hyperstyle(96,16,2).color,
		position=(10,96/16), offset=(4,2))
	text!(ax_eff, "Individual testing"; color=indivstyle.color,
		position=(10,1), offset=(4,2))

	# Fig d
	fig[2,2] = GridLayout()
	ax_eff = Axis(fig[2,2][1,1]; title="n=384 individuals",
		limits=(nothing,(0,12*1.2)), yticks=[0,1,384/40,384/48,384/32])
	vspan!(ax_eff, SHADELIMS...; color=SHADECOLOR)
	lines!(ax_eff, hypereff[(n=384,m=16,q=1)]; dorfmanstyle...)
	lines!(ax_eff, arrayeff[(16,24)];          arraystyle...)
	lines!(ax_eff, hypereff[(n=384,m=32,q=2)]; hyperstyle(384,32,2)...)
	lines!(ax_eff, pbesteff;                   pbeststyle...)
	lines!(ax_eff, indiveff;                   indivstyle...)

	ax_sens = Axis(fig[2,2][2,1])
	vspan!(ax_sens, SHADELIMS...; color=SHADECOLOR)
	scatter!(ax_sens, hypersens[(n=384,m=16,q=1)]; dorfmanstyle...)
	scatter!(ax_sens, arraysens[(16,24)];          arraystyle...)
	scatter!(ax_sens, hypersens[(n=384,m=32,q=2)]; hyperstyle(384,32,2)...)
	scatter!(ax_sens, pbestsens;                   pbeststyle...)
	scatter!(ax_sens, indivsens;                   indivstyle...)

	text!(ax_eff, dorfmanstr(384,16); color=dorfmanstyle.color,
		position=(100,5.5), offset=(0,0), align=(:left,:top))
	text!(ax_eff, "16x24 array"; color=arraystyle.color,
		position=(10,384/40), offset=(4,2))
	text!(ax_eff, hyperstr(384,32,2); color=hyperstyle(384,32,2).color,
		position=(10,384/32), offset=(4,2))
	text!(ax_eff, "P-BEST"; color=pbeststyle.color,
		position=(10,384/48), offset=(4,-2), align=(:left,:top))
	text!(ax_eff, "Individual testing"; color=indivstyle.color,
		position=(10,1), offset=(4,2))

	# Fig e
	fig[3,1] = GridLayout()
	ax_eff = Axis(fig[3,1][1,1]; title="n=384 individuals",
		limits=(nothing,(0,32*1.04)), yticks=[1,384/32,384/16,384/12])
	vspan!(ax_eff, SHADELIMS...; color=SHADECOLOR)
	lines!(ax_eff, hypereff[(n=384,m=32,q=2)]; hyperstyle(384,32,2)...)
	lines!(ax_eff, hypereff[(n=384,m=16,q=2)]; hyperstyle(384,16,2)...)
	lines!(ax_eff, hypereff[(n=384,m=12,q=2)]; hyperstyle(384,12,2)...)
	lines!(ax_eff, indiveff;                   indivstyle...)

	ax_sens = Axis(fig[3,1][2,1])
	vspan!(ax_sens, SHADELIMS...; color=SHADECOLOR)
	scatter!(ax_sens, hypersens[(n=384,m=32,q=2)]; hyperstyle(384,32,2)...)
	scatter!(ax_sens, hypersens[(n=384,m=16,q=2)]; hyperstyle(384,16,2)...)
	scatter!(ax_sens, hypersens[(n=384,m=12,q=2)]; hyperstyle(384,12,2)...)
	scatter!(ax_sens, indivsens;                   indivstyle...)
	
	text!(ax_eff, hyperstr(384,32,2); color=hyperstyle(384,32,2).color,
		position=(10,384/32), offset=(4,-2), align=(:left,:top))
	text!(ax_eff, hyperstr(384,16,2); color=hyperstyle(384,16,2).color,
		position=(10,384/16), offset=(4,-6), align=(:left,:top))
	text!(ax_eff, hyperstr(384,12,2); color=hyperstyle(384,12,2).color,
		position=(55,28), offset=(0,0))
	text!(ax_eff, "Individual testing"; color=indivstyle.color,
		position=(10,1), offset=(4,2))
	
	# Fig f
	fig[3,2] = GridLayout()
	ax_eff = Axis(fig[3,2][1,1]; title="n=384 individuals",
		limits=(nothing,(0,32*1.04)), yticks=[1,384/32,384/16,384/12])
	vspan!(ax_eff, SHADELIMS...; color=SHADECOLOR)
	lines!(ax_eff, hypereff[(n=384,m=12,q=3)]; hyperstyle(384,12,3)...)
	lines!(ax_eff, hypereff[(n=384,m=12,q=2)]; hyperstyle(384,12,2)...)
	lines!(ax_eff, hypereff[(n=384,m=12,q=1)]; hyperstyle(384,12,1)...)
	lines!(ax_eff, indiveff;                   indivstyle...)

	ax_sens = Axis(fig[3,2][2,1])
	vspan!(ax_sens, SHADELIMS...; color=SHADECOLOR)
	scatter!(ax_sens, hypersens[(n=384,m=12,q=3)]; hyperstyle(384,12,3)...)
	scatter!(ax_sens, hypersens[(n=384,m=12,q=2)]; hyperstyle(384,12,2)...)
	scatter!(ax_sens, hypersens[(n=384,m=12,q=1)]; hyperstyle(384,12,1)...)
	scatter!(ax_sens, indivsens;                   indivstyle...)
	
	text!(ax_eff, hyperstr(384,12,3); color=hyperstyle(384,12,3).color,
		position=(60,28), offset=(0,0))
	text!(ax_eff, hyperstr(384,12,2); color=hyperstyle(384,12,2).color,
		position=(10,384/12), offset=(4,-4), align=(:left,:top), rotation=-0.13)
	text!(ax_eff, hyperstr(384,12,1); color=hyperstyle(384,12,1).color,
		position=(10,24), offset=(4,-18), align=(:left,:top))
	text!(ax_eff, "Individual testing"; color=indivstyle.color,
		position=(10,1), offset=(4,2))

	# Common axis limits, ticks, etc.
	for (gl,sub) in zip(contents(fig[2:end,:]),'c':'z')
		ax_eff, ax_sens = contents(gl)
		rowgap!(gl, 8)

		# x-axis
		xlims!(ax_eff, extrema(DAYS))
		linkxaxes!(ax_eff, ax_sens)
		hidexdecorations!(ax_eff; ticks=false, grid=false)
		ax_eff.xticks = ax_sens.xticks = 10:20:190
		daystr = day -> "$day\n($(round(100*prevalence[day];digits=2))%)"
		ax_sens.xtickformat = x->daystr.(convert.(Int,x))
		ax_sens.xlabel = "Day (Prevalence)"

		# y-axis
		ylims!(ax_sens, (0.37,0.98))
		ax_sens.yticks = 0.4:0.1:0.9
		ax_sens.ytickformat = y->string.(convert.(Int,100*y),'%')
		ax_eff.ylabel = "Efficiency\n(relative to indiv. testing)"
		ax_sens.ylabel = "Sensitivity"

		# label
		Label(gl[1,1,TopLeft()], string(sub);
			halign=:left, valign=:top, font="Arial Bold")
	end

	rowgap!(fig.layout, 16)
	colgap!(fig.layout, 16)
	rowsize!(fig.layout, 1, Relative(1/7))
	save("fig-s9.png", fig; px_per_unit=2)
	save("fig-s9.pdf", fig)
	fig
end

# ╔═╡ Cell order:
# ╟─849adbca-08e5-11eb-1956-cd02f2a6be7a
# ╠═95eaf04a-08e5-11eb-043f-27646f1bf207
# ╟─210d331a-09c3-11eb-06ed-e1f4ddc023fc
# ╟─66dd9e7e-f948-4457-8bd6-45345be6abb5
# ╟─6b39ad12-5378-47dc-9a12-e1d5838b0dee
# ╟─fd643c2c-08e5-11eb-0d25-014798a505d1
# ╟─fe7cf842-08e5-11eb-0d07-41ed2902921a
# ╟─3f1b753e-0901-11eb-0060-1f1269588bef
# ╟─3f1c5288-0901-11eb-391f-dd7dd02d3d9b
# ╟─05c6923e-08e6-11eb-1d9b-21afc86a6b04
# ╟─cf1b2b1a-0c3b-11eb-2159-6ffeb682d595
# ╟─dabb2530-e060-4b8c-b03c-1cbeac98ea4b
# ╟─de3ef15a-9196-4214-b694-6f1448fccf4a
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
# ╟─0fd98c50-2203-11eb-0f58-d946dcf484b4
# ╟─6adbc303-83dd-4505-b50d-a94c9c1834bc
# ╟─09e96bdc-ebde-4899-a089-6c0418fa7ecc
# ╟─a6c3504c-543f-11eb-2f20-791273bc06b8
# ╟─c8c90422-4f35-4e90-90b8-4b30f3738f44
# ╟─640b0974-da62-4cf7-963b-121a5030af16
# ╟─485f245c-f7ba-4ccb-86f7-08de5e9fd7e9
# ╟─48b55bad-a658-4b22-9927-83e15bf900d1
# ╟─125cf732-35fd-48f9-a2f1-6261d92969fa
# ╟─27882286-ba06-4d0c-9876-ab8afe757d33
