### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ be6b935c-0918-11eb-227a-f7ffbb70c325
using IterTools, PGField, Primes, SparseArrays

# ╔═╡ 2f42420e-7bf8-11eb-2529-e3a58c1ee4b4
md"""
# HYPER
Functions to generate design.
"""

# ╔═╡ 9268d880-09a9-11eb-2d79-af4d3eaa85f8
function gen2factors(m)
	(m % 2 == 0) || throw("m must be a multiple of 2, got m=$m.")
	
	r = m-1
	starter = vcat(
		[[GFElem{r}(0),PGInf{r}()]],
		[[GFElem{r}(i),GFElem{r}(-i)] for i in 1:(r-1)÷2]
	)
	
	return ([Ai.+GFElem{r}(g) for Ai in starter] for g in 0:r-1)
end

# ╔═╡ c2fde6ba-0993-11eb-066b-6b54080e7bc1
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

# ╔═╡ 770a9a5c-09a0-11eb-1043-852c0087f749
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
	return sparse(collect(Iterators.flatten(poolidx)),repeat(1:n,inner=q),1,m,n) |> collect
end

# ╔═╡ Cell order:
# ╟─2f42420e-7bf8-11eb-2529-e3a58c1ee4b4
# ╠═be6b935c-0918-11eb-227a-f7ffbb70c325
# ╠═9268d880-09a9-11eb-2d79-af4d3eaa85f8
# ╠═c2fde6ba-0993-11eb-066b-6b54080e7bc1
# ╠═770a9a5c-09a0-11eb-1043-852c0087f749
