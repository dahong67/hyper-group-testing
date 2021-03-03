### A Pluto.jl notebook ###
# v0.12.16

using Markdown
using InteractiveUtils

# ╔═╡ 571c959e-1aee-11eb-0e1e-a1b4d8f7762d
using Distributions, IterTools, LaTeXStrings, PGField, Plots, Plots.Measures, Primes, Random, SparseArrays

# ╔═╡ f7c7716c-3593-11eb-00c3-677a7dbc5fec
md"""
# Theoretical model figures
"""

# ╔═╡ 1ac80b88-3594-11eb-2a77-7321f527f86a
md"""
## Theoretical analysis
"""

# ╔═╡ 9e393e90-35f0-11eb-25bc-af3fc684f85f
Base.binomial(m::Float64,q::Int) = q == 0 ? 1 : prod(m-i for i in 0:q-1)/factorial(q)

# ╔═╡ fd32d722-3595-11eb-27ea-7d153568cb6a
md"""
Expected number of tests (bound for $q>2$):
"""

# ╔═╡ 1db68090-3594-11eb-3c00-1f948011570b
function effbound(p,n,m,q;α=0,β=1)
	k = n*q/m
	r = 1-p
	u = binomial(m-2,q-2)*ceil(n/binomial(m,q))
	p1 = 1 - β + (β-α)*r^k
	p2 = p1^2 + (β-α)^2*r^(2*k-u)*(1-r^u)
	return (m + n*(1 - q*p1 + binomial(q,2)*p2))/n
end

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
mopt(p,n,q;α=0,β=1,maxm=n) = argmin([effbound(p,n,m,q;α=α,β=β) for m in 1:maxm] .|> x->isnan(x) ? Inf : x)

# ╔═╡ 5c9a1260-3595-11eb-0e60-3fe0c174ca5a
effopt(p,n,q;α=0,β=1,maxm=n) = effbound(p,n,mopt(p,n,q;α=α,β=β,maxm=maxm),q;α=α,β=β)

# ╔═╡ 055a9656-35a5-11eb-2094-1d0f43ac2d1d
md"""
## Simulating according to the theoretical model
"""

# ╔═╡ 68679664-7bcc-11eb-2b20-c384ea767ada
md"""
### Generation of HYPER designs
"""

# ╔═╡ fb1e305a-1aef-11eb-0c09-dfdb058d559c
function gen2factors(m)
	(m % 2 == 0) || throw("m must be a multiple of 2, got m=$m.")
	
	r = m-1
	starter = vcat(
		[[GFElem{r}(0),PGInf{r}()]],
		[[GFElem{r}(i),GFElem{r}(-i)] for i in 1:(r-1)÷2]
	)
	
	return ([Ai.+GFElem{r}(g) for Ai in starter] for g in 0:r-1)
end

# ╔═╡ fed96cbe-1aef-11eb-18b7-bbf7c302625e
function gen3factors(m)
	(m % 6 == 0) || throw("m must be a multiple of 6, got m=$m.")
	isprime(m-1) || throw("m-1 must be prime, got m=$m.")
	
	r = m-1
	F = [GFElem{r}(i) for i in 0:r-1]
	PG = [F; PGInf{r}()]
	
	O = (Set∘map)(PG) do i
		orbit = iterated(x -> -(inv(x)*(one(x)+x)), i)
		Set(collect(Iterators.take(orbit,r+1)))
	end
	ω = F[findfirst(i -> order(i) == r-1, F)]
	
	return ([(λ.*Ai).+g for Ai in O] for g in F, λ in ω.^(1:(r-1)÷2))
end

# ╔═╡ 0298df42-1af0-11eb-151b-ad0cae6f5363
function graphdesign(n,m,q; reorder=true)
	# Generate factors
	if q == 1
		r = m-1
		factors = vcat([[GFElem{r}(i)] for i in 0:r-1], [[PGInf{r}()]])
	elseif q == 2
		factors = Iterators.flatten(gen2factors(m))
	elseif q == 3
		factors = Iterators.flatten(gen3factors(m))
	else
		throw("q must be 1, 2 or 3, got q=$q.")
	end
	
	# Form pools and (optionally) reorder to group individuals in identical pools together
	pools = collect(Iterators.take(Iterators.cycle(factors),n))
	if reorder
		nfactors = binomial(m,q)
		pools = collect(Iterators.flatten(pools[i:nfactors:end] for i in 1:nfactors))
	end
	
	# Convert to pool indices and create sparse array
	poolidx = map(pool -> [num isa GFElem ? num.value+1 : m for num in pool], pools)
	return sparse(collect(Iterators.flatten(poolidx)),repeat(1:n,inner=q),1,m,n)
end

# ╔═╡ c085ecb2-1af4-11eb-1f1b-45f2ced5ef31
graphmatrix(n,m,q; reorder=true) = convert(Array{Bool},graphdesign(n,m,q,reorder=reorder))

# ╔═╡ 7b4b8da8-7bcc-11eb-2976-4b2155781857
md"""
### Simulation following the model
"""

# ╔═╡ b7046b98-3aa5-11eb-31c4-95119c9a90fa
struct SimCounts{T<:Integer}
	numtests::T
	falseneg::T
	truepos::T
	trueneg::T
	falsepos::T
end

# ╔═╡ 1038eba8-3aa6-11eb-0bf1-37650421e05e
Base.:+(x::SimCounts,y::SimCounts) = SimCounts(
	x.numtests+y.numtests,
	x.falseneg+y.falseneg,
	x.truepos+y.truepos,
	x.trueneg+y.trueneg,
	x.falsepos+y.falsepos,
)

# ╔═╡ 99ccff42-35a3-11eb-1ef1-bb3d0a6f89c1
function sim(p,n,m,q;α=0,β=1)
	A = graphmatrix(n,m,q)
	x = rand(Bernoulli(p),n)
	y = [any(x[sidx] for sidx in findall(pool)) ? rand(Bernoulli(β)) : rand(Bernoulli(α)) for pool in eachrow(A)]
	xput = [all(y[pidx] for pidx in findall(sample)) for sample in eachcol(A)]
	xhat = fill(false,n)
	for i in findall(xput)
		xhat[i] = x[i] ? rand(Bernoulli(β)) : rand(Bernoulli(α))
	end
	return SimCounts(
		m + count(xput),
		count(  x .& .!xhat),
		count(  x .&   xhat),
		count(.!x .& .!xhat),
		count(.!x .&   xhat),
	)
end

# ╔═╡ c5735898-7bcc-11eb-0e6a-ed0c661835ad
md"""
### Overall simulation efficiency / sensitivity / specificity
"""

# ╔═╡ 2eb94ae2-3a73-11eb-0a46-5585fa02f634
eff(cts::SimCounts) = cts.numtests/(cts.falseneg+cts.truepos+cts.trueneg+cts.falsepos)

# ╔═╡ 62d13952-3a73-11eb-2e22-4db2910fdb48
sens(cts::SimCounts) = cts.truepos/(cts.truepos+cts.falseneg)

# ╔═╡ ddb734be-3a73-11eb-05f6-b30cd2f9e57b
spec(cts::SimCounts) = cts.trueneg/(cts.trueneg+cts.falsepos)

# ╔═╡ 2383c37a-35a8-11eb-3f37-ed1672ed16dc
md"""
## Figures a-b
"""

# ╔═╡ 9a21144c-7bcb-11eb-3757-a5e3b58e0fde
md"""
### Fig a
"""

# ╔═╡ 249a55a6-3677-11eb-0b55-091be1cbf0a2
let pticks=[0.1,0.2,0.3,0.4,0.5,1]./100
	mticks = 0.02:0.01:0.09
	plot(size=(400,300),dpi=200,framestyle=:box,
		xscale=:log10,legend=:bottomright,
		xlabel=L"\textrm{Prevalence } (\%)",ylabel=L"m/n",
		xticks=(pticks,string.(pticks*100)),xlim=extrema(pticks),
		yticks=(mticks,string.(mticks)),ylim=(0.019,0.095))
	prange = 10 .^ range(log10.(extrema(pticks))...,length=1000)
	
	for n in [96,384,1536,6144]
		plot!(p->mopt(p,n,2)/n,prange,label=latexstring("n = $n"),linewidth=3)
	end
	plot!(p->2*p^(2/3)-p,prange,label=L"2p^{2/3}-p",color=:black,linewidth=3,linestyle=:dash)
	annotate!(0.075/100,0.095,text("a"),top_margin=1mm)
	savefig("fig-a.png"); plot!()
end

# ╔═╡ a5f9a004-7bcb-11eb-2d97-b3b7a54a54e5
md"""
### Fig b
"""

# ╔═╡ adca3628-506e-11eb-3c52-5d8690c0a184
let pticks=[0.1,0.2,0.3,0.4,0.5,1]./100
	effticks = round.(inv.(3*pticks.^(2/3)),digits=1)
	plot(size=(400,300),dpi=200,framestyle=:box,
		xscale=:log10,yscale=:log10,legend=:topright,
		xlabel=L"\textrm{Prevalence } (\%)",ylabel=L"\textrm{Efficiency gain}: n/\mathrm{E}\;T",
		xticks=(pticks,string.(pticks*100)),xlim=extrema(pticks),
		yticks=(effticks,string.(effticks)))
	prange = 10 .^ range(log10.(extrema(pticks))...,length=1000)
	
	for n in [96,384,1536,6144]
		plot!(p->1/effopt(p,n,2),prange,label=latexstring("n = $n"),linewidth=3)
	end
	plot!(p->1/(3*p^(2/3)),prange,label=L"p^{-2/3}/3",color=:black,linewidth=3,linestyle=:dash)
	annotate!(0.075/100,36,text("b"))
	savefig("fig-b.png"); plot!()
end

# ╔═╡ 0f4e5394-7bcc-11eb-3521-77679d12146e
md"""
## Figs c-e / f-h

Figs c-e are generated by same code as Figs f-h,
with only the setup changed
(un/comment the relevant line to choose which set).
"""

# ╔═╡ e2a7aad6-7bcc-11eb-3383-65908f153f39
md"""
### Setup
"""

# ╔═╡ 75a6fdb4-3aa0-11eb-38d5-f3c1c4ade957
# opchar, figlabels = (α=0.05,β=0.90), (eff="c",sens="d",spec="e")
opchar, figlabels = (α=0.05,β=0.80), (eff="f",sens="g",spec="h")

# ╔═╡ f576fd34-3a76-11eb-164f-e545ad1f1ba4
prange = 10 .^ range(log10(0.1/100),log10(1/100),length=100)

# ╔═╡ 633999fe-3abc-11eb-113a-ab1c5b66d0fd
designs = [
	(n=3072,m=96, q=2),
	(n=3072,m=192,q=2),
	(n=3072,m=384,q=2),
]

# ╔═╡ 0ce1d1dc-3aaf-11eb-2512-355867435800
Random.seed!(0); results = [[sum(sim(p,des.n,des.m,des.q;α=opchar.α,β=opchar.β) for _ in 1:25) for p in prange] for des in designs]

# ╔═╡ 1a8fd33e-7bcc-11eb-2e76-eb8cb430ec07
md"""
### Fig c / f
"""

# ╔═╡ bcf081c2-506e-11eb-31f2-2743dbbf1bbf
let pticks=[0.1,0.2,0.3,0.4,0.5,1]./100
	effticks = 5:5:25
	plot(size=(300,300),dpi=200,framestyle=:box,
		xscale=:log10,yscale=:log10,legend=:topright,
		xlabel=L"\textrm{Prevalence } (\%)",ylabel=L"\textrm{Efficiency gain}: n/\mathrm{E}\;T",
		xticks=(pticks,string.(pticks*100)),xlim=extrema(pticks),
		yticks=(effticks,string.(effticks)),ylim=(3.8,27))
	
	for (des,res,lcolor) in zip(designs,results,palette(:default))
		scatter!(prange,inv.(eff.(res)),color=lcolor,markerstrokecolor=:auto,markersize=2,opacity=0.6,label="")
		plot!(prange,p->1/effbound(p,des.n,des.m,des.q;α=opchar.α,β=opchar.β),linewidth=2,color=lcolor,
			label=latexstring("m=$(des.m)"))
	end
	annotate!(0.07/100,32,text(figlabels.eff),top_margin=1mm)
	title!(latexstring("\\alpha = $(opchar.α), \\beta = $(opchar.β)"),titlefontsize=9)
	savefig("fig-$(figlabels.eff).png"); plot!()
end

# ╔═╡ 2c59b38c-7bcc-11eb-2a83-458f3d104481
md"""
### Fig d/g
"""

# ╔═╡ d50cdc90-3b1f-11eb-26d5-b9e3c3869715
let pticks=[0.1,0.2,0.3,0.4,0.5,1]./100
	sensticks = 0.5:0.1:0.9
	plot(size=(300,300),dpi=200,framestyle=:box,
		xscale=:log10,legend=:bottomright,
		xlabel=L"\textrm{Prevalence } (\%)",ylabel=L"\textrm{Sensitivity}",
		xticks=(pticks,string.(pticks*100)),xlim=extrema(pticks),
		yticks=(sensticks,string.(sensticks)),ylim=(0.43,0.9))
	
	for (des,res,lcolor) in zip(designs,results,palette(:default))
		scatter!(prange,sens.(res),color=lcolor,markerstrokecolor=:auto,markersize=2,opacity=0.6,label="")
		plot!(prange,p->sens(p,des.n,des.m,des.q;α=opchar.α,β=opchar.β),linewidth=2,color=lcolor,
			label=latexstring("m=$(des.m)"))
	end
	annotate!(0.065/100,0.94,text(figlabels.sens),top_margin=1mm)
	title!(latexstring("\\alpha = $(opchar.α), \\beta = $(opchar.β)"),titlefontsize=9)
	savefig("fig-$(figlabels.sens).png"); plot!()
end

# ╔═╡ 332fda94-7bcc-11eb-1694-7baf75b3a709
md"""
### Fig e/h
"""

# ╔═╡ c588b5de-3b20-11eb-012a-6bfd63359b4d
let pticks=[0.1,0.2,0.3,0.4,0.5,1]./100
	specticks = 0.99:0.002:1
	plot(size=(300,300),dpi=200,framestyle=:box,
		xscale=:log10,legend=:bottomleft,
		xlabel=L"\textrm{Prevalence } (\%)",ylabel=L"\textrm{Specificity}",
		xticks=(pticks,string.(pticks*100)),xlim=extrema(pticks),
		yticks=(specticks,string.(specticks)),ylim=(0.989,1),
	)
	
	for (des,res,lcolor) in zip(designs,results,palette(:default))
		scatter!(prange,spec.(res),color=lcolor,markerstrokecolor=:auto,markersize=2,opacity=0.6,label="")
		plot!(prange,p->spec(p,des.n,des.m,des.q;α=opchar.α,β=opchar.β),linewidth=2,color=lcolor,
			label=latexstring("m=$(des.m)"))
	end
	annotate!(0.055/100,1.001,text(figlabels.spec),top_margin=1mm)
	title!(latexstring("\\alpha = $(opchar.α), \\beta = $(opchar.β)"),titlefontsize=9)
	savefig("fig-$(figlabels.spec).png"); plot!()
end

# ╔═╡ Cell order:
# ╟─f7c7716c-3593-11eb-00c3-677a7dbc5fec
# ╠═571c959e-1aee-11eb-0e1e-a1b4d8f7762d
# ╟─1ac80b88-3594-11eb-2a77-7321f527f86a
# ╠═9e393e90-35f0-11eb-25bc-af3fc684f85f
# ╟─fd32d722-3595-11eb-27ea-7d153568cb6a
# ╠═1db68090-3594-11eb-3c00-1f948011570b
# ╟─530c81aa-7bca-11eb-1769-a7fe4d74bb1b
# ╠═133cd066-3a9d-11eb-3fa7-19b51f9943ff
# ╠═1d15d678-3a9d-11eb-06ef-8fd57eba76be
# ╟─215f9960-3595-11eb-313d-ef186d77265b
# ╠═2101c74c-3595-11eb-059c-25b8a1a76ff5
# ╠═5c9a1260-3595-11eb-0e60-3fe0c174ca5a
# ╟─055a9656-35a5-11eb-2094-1d0f43ac2d1d
# ╟─68679664-7bcc-11eb-2b20-c384ea767ada
# ╟─fb1e305a-1aef-11eb-0c09-dfdb058d559c
# ╟─fed96cbe-1aef-11eb-18b7-bbf7c302625e
# ╟─0298df42-1af0-11eb-151b-ad0cae6f5363
# ╟─c085ecb2-1af4-11eb-1f1b-45f2ced5ef31
# ╟─7b4b8da8-7bcc-11eb-2976-4b2155781857
# ╠═b7046b98-3aa5-11eb-31c4-95119c9a90fa
# ╠═1038eba8-3aa6-11eb-0bf1-37650421e05e
# ╠═99ccff42-35a3-11eb-1ef1-bb3d0a6f89c1
# ╟─c5735898-7bcc-11eb-0e6a-ed0c661835ad
# ╠═2eb94ae2-3a73-11eb-0a46-5585fa02f634
# ╠═62d13952-3a73-11eb-2e22-4db2910fdb48
# ╠═ddb734be-3a73-11eb-05f6-b30cd2f9e57b
# ╟─2383c37a-35a8-11eb-3f37-ed1672ed16dc
# ╟─9a21144c-7bcb-11eb-3757-a5e3b58e0fde
# ╟─249a55a6-3677-11eb-0b55-091be1cbf0a2
# ╟─a5f9a004-7bcb-11eb-2d97-b3b7a54a54e5
# ╟─adca3628-506e-11eb-3c52-5d8690c0a184
# ╟─0f4e5394-7bcc-11eb-3521-77679d12146e
# ╟─e2a7aad6-7bcc-11eb-3383-65908f153f39
# ╠═75a6fdb4-3aa0-11eb-38d5-f3c1c4ade957
# ╟─f576fd34-3a76-11eb-164f-e545ad1f1ba4
# ╟─633999fe-3abc-11eb-113a-ab1c5b66d0fd
# ╠═0ce1d1dc-3aaf-11eb-2512-355867435800
# ╟─1a8fd33e-7bcc-11eb-2e76-eb8cb430ec07
# ╟─bcf081c2-506e-11eb-31f2-2743dbbf1bbf
# ╟─2c59b38c-7bcc-11eb-2a83-458f3d104481
# ╟─d50cdc90-3b1f-11eb-26d5-b9e3c3869715
# ╟─332fda94-7bcc-11eb-1694-7baf75b3a709
# ╟─c588b5de-3b20-11eb-012a-6bfd63359b4d
