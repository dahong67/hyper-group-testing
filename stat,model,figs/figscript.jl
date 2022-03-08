### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ 571c959e-1aee-11eb-0e1e-a1b4d8f7762d
begin
	import Pkg; Pkg.activate(@__DIR__)
	using CacheVariables, Distributions, HyperGen, ProgressLogging, StableRNGs
	using CairoMakie
end

# ╔═╡ f7c7716c-3593-11eb-00c3-677a7dbc5fec
md"""
# Plots for the statistical model
"""

# ╔═╡ 12dbc40b-2495-459b-8fef-e80ed2a6bbb4
md"""
## Simulation according to the model
"""

# ╔═╡ b604d599-510b-43ed-8302-78f05f66fe1b
SimCounts = @NamedTuple begin
	numtests::Int
	falseneg::Int
	truepos::Int
	trueneg::Int
	falsepos::Int
end

# ╔═╡ 2d465f6e-5bbc-44fd-9c58-38e66d6501bc
Base.:+(x::SimCounts,y::SimCounts) = SimCounts((
	x.numtests+y.numtests,
	x.falseneg+y.falseneg,
	x.truepos+y.truepos,
	x.trueneg+y.trueneg,
	x.falsepos+y.falsepos,
))

# ╔═╡ f9d61ede-6e49-43b3-b973-4455fc6bdd97
function sim(rng,p,n,m,q;α=0,β=1)
	A = hyperdesign(n,m,q)
	x = rand(rng,Bernoulli(p),n)
	y = [any(x[sidx] for sidx in findall(pool)) ?
			rand(rng,Bernoulli(β)) : rand(rng,Bernoulli(α)) for pool in eachrow(A)]
	xput = [all(y[pidx] for pidx in findall(sample)) for sample in eachcol(A)]
	xhat = fill(false,n)
	for i in findall(xput)
		xhat[i] = x[i] ? rand(rng,Bernoulli(β)) : rand(rng,Bernoulli(α))
	end
	return SimCounts((
		m + count(xput),
		count(  x .& .!xhat),
		count(  x .&   xhat),
		count(.!x .& .!xhat),
		count(.!x .&   xhat),
	))
end

# ╔═╡ d94f7aae-50e6-4b3a-95db-2628b4b4d4db
md"""
### Overall simulation efficiency / sensitivity / specificity
"""

# ╔═╡ 62c3456a-75f0-4116-8149-3622bf81d6a8
eff(cts::SimCounts) =
	(cts.falseneg+cts.truepos+cts.trueneg+cts.falsepos)/cts.numtests

# ╔═╡ 64bd6c15-2c35-4b87-b145-af28173b0c05
sens(cts::SimCounts) = cts.truepos/(cts.truepos+cts.falseneg)

# ╔═╡ 89b642ac-5778-4ba8-ac44-4503a481318d
spec(cts::SimCounts) = cts.trueneg/(cts.trueneg+cts.falsepos)

# ╔═╡ 1ac80b88-3594-11eb-2a77-7321f527f86a
md"""
## Theoretical analysis
"""

# ╔═╡ 9e393e90-35f0-11eb-25bc-af3fc684f85f
Base.binomial(m::Float64,q::Int) =
	q == 0 ? 1 : prod(m-i for i in 0:q-1)/factorial(q)

# ╔═╡ fd32d722-3595-11eb-27ea-7d153568cb6a
md"""
Expected number of tests (bound for $q>2$):
"""

# ╔═╡ 1db68090-3594-11eb-3c00-1f948011570b
function numtestbound(p,n,m,q;α=0,β=1)
	k = n*q/m
	r = 1-p
	u = binomial(m-2,q-2)*ceil(n/binomial(m,q))
	p1 = 1 - β + (β-α)*r^k
	p2 = p1^2 + (β-α)^2*r^(2*k-u)*(1-r^u)
	return m + n*(1 - q*p1 + binomial(q,2)*p2)
end

# ╔═╡ fccdb23e-61fd-4336-9f6c-06c750f0acbc
effbound(p,n,m,q;α=0,β=1) = n/numtestbound(p,n,m,q;α,β)

# ╔═╡ 530c81aa-7bca-11eb-1769-a7fe4d74bb1b
md"""
Sensitivity/specificity:
"""

# ╔═╡ 133cd066-3a9d-11eb-3fa7-19b51f9943ff
sens(p,n,m,q;α=0,β=1) = β^(q+1)

# ╔═╡ 1d15d678-3a9d-11eb-06ef-8fd57eba76be
function spec(p,n,m,q;α=0,β=1)
	k = n*q/m
	r = 1-p
	γ = (β + (α-β)*r^(k-1))^q
	return 1-α*γ
end

# ╔═╡ 215f9960-3595-11eb-313d-ef186d77265b
md"""
Numerical (non-asymptotic) optimization of bound:
"""

# ╔═╡ 2101c74c-3595-11eb-059c-25b8a1a76ff5
mopt(p,n,q;α=0,β=1,maxm=n) = argmin(1:maxm) do m
	bound = numtestbound(p,n,m,q;α,β)
	return !isnan(bound) ? bound : Inf
end

# ╔═╡ 5c9a1260-3595-11eb-0e60-3fe0c174ca5a
effopt(p,n,q;α=0,β=1,maxm=n) = effbound(p,n,mopt(p,n,q;α,β,maxm),q;α,β)

# ╔═╡ c3bb886c-397f-4189-8ac3-6af46bfd9261
md"""
## Run/load simulation
"""

# ╔═╡ 2c512a1e-5ee0-4592-ad45-d54f26fa3aac
PRANGE = 10 .^ range(log10(0.1/100),log10(1/100),length=100)

# ╔═╡ 82bfa7d1-4ef5-4494-a776-e63493c76e9d
DESIGNS = [(n=3072,m=96, q=2), (n=3072,m=192,q=2), (n=3072,m=384,q=2)]

# ╔═╡ 2bd1d337-0d9a-496c-8543-530fbc5f6f79
CONFIGS = [(α=0.05,β=0.90), (α=0.05,β=0.80)]

# ╔═╡ e4d82079-b5ba-4f4a-b3e9-d36b98829e70
NRUNS = 25

# ╔═╡ b1d7979d-cfe0-4399-b07e-b02f070f2101
results = cache("simcache.bson") do
	map(CONFIGS) do (α,β)
		map(DESIGNS) do (n,m,q)
			rng = StableRNG(0)
			label = "α=$α,β=$β | n=$n,m=$m,q=$q | "
			@withprogress name=label map(enumerate(PRANGE)) do (idx,p)
				@logprogress idx/length(PRANGE)
				sum(sim(rng,p,n,m,q;α,β) for _ in 1:NRUNS)
			end
		end
	end
end

# ╔═╡ d2623c58-4c6e-4617-88e6-89204155c885
md"""
## Figures
"""

# ╔═╡ 93d64dbb-cc6f-4e9a-a84d-bf173b21ed6b
begin
	mm_to_units(mm) = floor(Int,mm/25.4*72/0.75)
	mm_to_units(mm1,mm2,mmrest...) = mm_to_units.((mm1,mm2,mmrest...))
end

# ╔═╡ 23acc613-bcd5-4f3f-97b1-b56c573271e7
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

# ╔═╡ 18514663-7393-4ec5-9010-48d293ea91f5
with_theme(THEME; Legend=(;padding=(0,0,0,0))) do
	fig = Figure(; resolution=mm_to_units(180,180))

	## Simulation figures
	PTICKS = [0.001, 0.002, 0.003, 0.005, 0.01]
	EFFTICKS, EFFLIMS = 5:5:25, (3.8,27)
	SENSTICKS, SENSLIMS = 0.5:0.1:0.9, (0.43,0.9)
	SPECTICKS, SPECLIMS = 0.99:0.002:1, (0.989,1)
	PALETTE = Makie.wong_colors(0.6)
	
	resultsitr = zip(Iterators.partition('a':'z',3),CONFIGS,results)
	for (row,(lbl,(α,β),res)) in enumerate(resultsitr)
		# Efficiency
		ax_eff = Axis(fig[row,1]; limits=(extrema(PTICKS), EFFLIMS),
			xlabel="Prevalence",
			xscale=log10, xticks=PTICKS, xtickformat=x->string.(x.*100,'%'),
			ylabel="Efficiency (relative to individual testing)",
			yscale=log10, yticks=EFFTICKS,
			title="β=$β, 1-α=$(1-α)")
		for (desres,(n,m,q),color) in zip(res,DESIGNS,PALETTE)
			scatter!(ax_eff, PRANGE, eff.(desres); color)
			lines!(ax_eff, PRANGE, p->effbound(p,n,m,q;α,β); color, label="m=$m")
		end
		axislegend(ax_eff)
		Label(fig[row,1,TopLeft()], string(lbl[1]);
			halign=:left, valign=:top, font="Arial Bold")

		# Sensitivity
		ax_sens = Axis(fig[row,2]; limits=(extrema(PTICKS), SENSLIMS),
			xlabel="Prevalence",
			xscale=log10, xticks=PTICKS, xtickformat=x->string.(x.*100,'%'),
			ylabel="Sensitivity",
			yticks=SENSTICKS, ytickformat=y->string.(convert.(Int,y.*100),'%'),
			title="β=$β, 1-α=$(1-α)")
		for (desres,(n,m,q),color) in zip(res,DESIGNS,PALETTE)
			scatter!(ax_sens, PRANGE, sens.(desres); color)
			lines!(ax_sens, PRANGE, p->sens(p,n,m,q;α,β); color, label="m=$m")
		end
		axislegend(ax_sens; position=:ct, orientation=:horizontal)
		Label(fig[row,2,TopLeft()], string(lbl[2]);
			halign=:left, valign=:top, font="Arial Bold")

		# Specificity
		ax_spec = Axis(fig[row,3]; limits=(extrema(PTICKS), SPECLIMS),
			xlabel="Prevalence",
			xscale=log10, xticks=PTICKS, xtickformat=x->string.(x.*100,'%'),
			ylabel="Specificity",
			yticks=SPECTICKS, ytickformat=y->string.(y.*100,'%'),
			title="β=$β, 1-α=$(1-α)")
		for (desres,(n,m,q),color) in zip(res,DESIGNS,PALETTE)
			scatter!(ax_spec, PRANGE, spec.(desres); color)
			lines!(ax_spec, PRANGE, p->spec(p,n,m,q;α,β); color, label="m=$m")
		end
		axislegend(ax_spec; position=:lb)
		Label(fig[row,3,TopLeft()], string(lbl[3]);
			halign=:left, valign=:top, font="Arial Bold")
	end

	## Optimization
	fig[end+1,:] = row = GridLayout()
	PRANGE_FINE = 10 .^ range(log10.(extrema(PRANGE))...,length=1000)
	
	ax_mopt = Axis(row[1,1]; limits=(extrema(PTICKS), (0.019,0.095)),
		xlabel="Prevalence",
		xscale=log10, xticks=PTICKS, xtickformat=x->string.(x.*100,'%'),
		ylabel="Optimal choice of pools per individual", yticks=0.02:0.01:0.09)
	for n in [96,384,1536,6144]
		lines!(ax_mopt, PRANGE_FINE, p->mopt(p,n,2)/n; label="n=$n")
	end
	lines!(ax_mopt, PRANGE_FINE, p->2*p^(2/3)-p; label="2p²ᐟ³-p",
		color=:black, linestyle=:dash)
	axislegend(ax_mopt; position=:lt)
	Label(row[1,1,TopLeft()], string(('a':'z')[length(results)*3+1]);
		halign=:left, valign=:top, font="Arial Bold")

	EFFOPTTICKS = inv.(3*PTICKS.^(2/3))
	ax_effopt = Axis(row[1,2]; limits=(extrema(PTICKS), nothing),
		xlabel="Prevalence",
		xscale=log10, xticks=PTICKS, xtickformat=x->string.(x.*100,'%'),
		ylabel="Efficiency (relative to individual testing)",
		yscale=log10, yticks=EFFOPTTICKS,
		ytickformat=y->string.(round.(y;digits=1)))
	for n in [96,384,1536,6144]
		lines!(ax_effopt, PRANGE_FINE, p->effopt(p,n,2); label="n=$n")
	end
	lines!(ax_effopt, PRANGE_FINE, p->inv(3*p^(2/3)); label="p⁻²ᐟ³/3",
		color=:black, linestyle=:dash)
	axislegend(ax_effopt)
	Label(row[1,2,TopLeft()], string(('a':'z')[length(results)*3+2]);
		halign=:left, valign=:top, font="Arial Bold")

	rowgap!(fig.layout, 16)
	colgap!(fig.layout, 16)
	colgap!(row, 16)
	save("fig-s1.png", fig; px_per_unit=2)
	save("fig-s1.pdf", fig)
	fig
end

# ╔═╡ Cell order:
# ╟─f7c7716c-3593-11eb-00c3-677a7dbc5fec
# ╠═571c959e-1aee-11eb-0e1e-a1b4d8f7762d
# ╟─12dbc40b-2495-459b-8fef-e80ed2a6bbb4
# ╠═b604d599-510b-43ed-8302-78f05f66fe1b
# ╠═2d465f6e-5bbc-44fd-9c58-38e66d6501bc
# ╠═f9d61ede-6e49-43b3-b973-4455fc6bdd97
# ╟─d94f7aae-50e6-4b3a-95db-2628b4b4d4db
# ╠═62c3456a-75f0-4116-8149-3622bf81d6a8
# ╠═64bd6c15-2c35-4b87-b145-af28173b0c05
# ╠═89b642ac-5778-4ba8-ac44-4503a481318d
# ╟─1ac80b88-3594-11eb-2a77-7321f527f86a
# ╠═9e393e90-35f0-11eb-25bc-af3fc684f85f
# ╟─fd32d722-3595-11eb-27ea-7d153568cb6a
# ╠═1db68090-3594-11eb-3c00-1f948011570b
# ╠═fccdb23e-61fd-4336-9f6c-06c750f0acbc
# ╟─530c81aa-7bca-11eb-1769-a7fe4d74bb1b
# ╠═133cd066-3a9d-11eb-3fa7-19b51f9943ff
# ╠═1d15d678-3a9d-11eb-06ef-8fd57eba76be
# ╟─215f9960-3595-11eb-313d-ef186d77265b
# ╠═2101c74c-3595-11eb-059c-25b8a1a76ff5
# ╠═5c9a1260-3595-11eb-0e60-3fe0c174ca5a
# ╟─c3bb886c-397f-4189-8ac3-6af46bfd9261
# ╟─2c512a1e-5ee0-4592-ad45-d54f26fa3aac
# ╟─82bfa7d1-4ef5-4494-a776-e63493c76e9d
# ╟─2bd1d337-0d9a-496c-8543-530fbc5f6f79
# ╟─e4d82079-b5ba-4f4a-b3e9-d36b98829e70
# ╠═b1d7979d-cfe0-4399-b07e-b02f070f2101
# ╟─d2623c58-4c6e-4617-88e6-89204155c885
# ╟─93d64dbb-cc6f-4e9a-a84d-bf173b21ed6b
# ╟─23acc613-bcd5-4f3f-97b1-b56c573271e7
# ╟─18514663-7393-4ec5-9010-48d293ea91f5
