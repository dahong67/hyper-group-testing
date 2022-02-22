### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ 571c959e-1aee-11eb-0e1e-a1b4d8f7762d
begin
	import Pkg; Pkg.activate(@__DIR__)
	using CairoMakie, Distributions, LaTeXStrings
end

# ╔═╡ f7c7716c-3593-11eb-00c3-677a7dbc5fec
md"""
# Plots comparing bounds
"""

# ╔═╡ 6e5aab50-b74b-11eb-33c3-5d53afe92e95
md"""
## Problem 1: any number of stages, any decoder, etc.
"""

# ╔═╡ 7fdcfedc-b74b-11eb-399f-95f6fdf61e67
md"""
### Lower bound (binary entropy)
"""

# ╔═╡ 995ee74e-b74b-11eb-20ef-0fac0d588eac
prob1_lb(p) = entropy(Bernoulli(p),2)

# ╔═╡ fb810c18-b74b-11eb-39a6-5b769158670f
md"""
### Upper bound for fully adaptive

Theorem 3 from:
Aldridge, Matthew. "Rates of adaptive group testing in the linear regime." 2019 IEEE International Symposium on Information Theory (ISIT). IEEE, 2019.
"""

# ╔═╡ bd23cdb0-b74c-11eb-1d4e-2bd846b6d795
function adap_ub(p)
	q = 1-p
	
	m = ceil(-log(2-p)/log(1-p))
	a = floor(log(2,m))
	b = m - 2^a
	
	F = q^m + (1-q^(m-2*b))*(a+1) + (q^(m-2*b)-q^m)*(a+2)
	G = m*q^m + 1/p*(1 + m*q^(m+1) - (m+1)*q^m)
	return F/G
end

# ╔═╡ aac6b044-b74a-11eb-2199-1b6613ccaa1b
md"""
## Problem 2: two-stage conservative testing
"""

# ╔═╡ 1ac80b88-3594-11eb-2a77-7321f527f86a
md"""
### Lower bound

Theorem 5 from:
Aldridge, Matthew. "Conservative two-stage group testing." arXiv preprint arXiv:2005.06617 (2020).
"""

# ╔═╡ 0919ebde-b74b-11eb-3947-878a1f65f615
function prob2_lb(p;wrange=2:100000)
	if p >= (3-sqrt(5))/2
		return 1.0
	else
		fp = maximum(w->-w*log(1-(1-p)^(w-1)),wrange)
		gp = maximum(w->-w*log(1-(1-p)^w    ),wrange)
		bound2 = 1/gp*(log(gp)+1)
		bound3 = p + 1/fp*(log((1-p)*fp)+1)
		return max(bound2,bound3)
	end
end

# ╔═╡ 95bfab7c-b74a-11eb-04d7-ab9d03e13c5b
md"""
### Upper bound for HYPER
"""

# ╔═╡ 9e393e90-35f0-11eb-25bc-af3fc684f85f
Base.binomial(m::Float64,q::Int) =
	q == 0 ? 1 : prod(m-i for i in 0:q-1)/factorial(q)

# ╔═╡ 1db68090-3594-11eb-3c00-1f948011570b
function hyper_ub(p,n,m,q)
	α = 0
	β = 1
	k = n*q/m
	r = 1-p
	u = binomial(m-2,q-2)*ceil(n/binomial(m,q))
	p1 = 1 - β + (β-α)*r^k
	p2 = p1^2 + (β-α)^2*r^(2*k-u)*(1-r^u)
	return (m + n*(1 - q*p1 + binomial(q,2)*p2))/n
end

# ╔═╡ 215f9960-3595-11eb-313d-ef186d77265b
md"""
Numerical optimization of bound:
"""

# ╔═╡ 5c9a1260-3595-11eb-0e60-3fe0c174ca5a
hyper_ub(p,n,q;maxm=n) = minimum(hyper_ub(p,n,m,q) for m in q:maxm)

# ╔═╡ 2e3ff118-b74b-11eb-2204-034129ec8abe
hyper_ub(p,n;qrange=1:3) = minimum(hyper_ub(p,n,q) for q in qrange)

# ╔═╡ 2383c37a-35a8-11eb-3f37-ed1672ed16dc
md"""
## Plot
"""

# ╔═╡ 44682866-b1ee-11eb-15dc-69a15ac82723
with_theme(; linewidth=5, Legend=(;framevisible=false)) do
	# Config
	PRANGE = 0.001:0.001:0.2
	PTICKS = 0.01:0.01:0.2
	EFFTICKS = 0:0.2:1
	LSIZE = 20f0 # size for LaTeXStrings

	# Figure
	fig = Figure(;resolution=(1200,600))
	ax = Axis(fig[1,1]; limits=((0,maximum(PRANGE)), (0,1)),
		xlabel=LaTeXString("Prevalence"), xlabelsize=LSIZE,
		xticks=PTICKS, yticks=EFFTICKS,
		xtickformat=x->string.(convert.(Int,round.(x.*100;digits=10)),'%'),
		ylabel=L"tests/individual: $\mathbf{\mathrm{E}}T/n$", ylabelsize=LSIZE)

	# Problem 2
	n = 6144
	lines!(ax, PRANGE, p->hyper_ub(p,n); label="Upper bound for HYPER method")
	lines!(ax, PRANGE, prob2_lb; linestyle=:dash,
		label="Lower bound for problem 2: \
		two-stage testing with conservative decoder")
	
	# Problem 1
	lines!(ax, PRANGE, adap_ub; label="Upper bound for fully adaptive method")
	lines!(ax, PRANGE, prob1_lb; linestyle=:dash,
		label="Lower bound for problem 1: \
		any number of stages, any decoder, etc.")

	axislegend(ax; position=:rb)
	save("fig-s2.png", fig)
	fig
end

# ╔═╡ Cell order:
# ╟─f7c7716c-3593-11eb-00c3-677a7dbc5fec
# ╠═571c959e-1aee-11eb-0e1e-a1b4d8f7762d
# ╟─6e5aab50-b74b-11eb-33c3-5d53afe92e95
# ╟─7fdcfedc-b74b-11eb-399f-95f6fdf61e67
# ╠═995ee74e-b74b-11eb-20ef-0fac0d588eac
# ╟─fb810c18-b74b-11eb-39a6-5b769158670f
# ╠═bd23cdb0-b74c-11eb-1d4e-2bd846b6d795
# ╟─aac6b044-b74a-11eb-2199-1b6613ccaa1b
# ╟─1ac80b88-3594-11eb-2a77-7321f527f86a
# ╠═0919ebde-b74b-11eb-3947-878a1f65f615
# ╟─95bfab7c-b74a-11eb-04d7-ab9d03e13c5b
# ╠═9e393e90-35f0-11eb-25bc-af3fc684f85f
# ╠═1db68090-3594-11eb-3c00-1f948011570b
# ╟─215f9960-3595-11eb-313d-ef186d77265b
# ╠═5c9a1260-3595-11eb-0e60-3fe0c174ca5a
# ╠═2e3ff118-b74b-11eb-2204-034129ec8abe
# ╟─2383c37a-35a8-11eb-3f37-ed1672ed16dc
# ╟─44682866-b1ee-11eb-15dc-69a15ac82723
