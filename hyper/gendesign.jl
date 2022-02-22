### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ be6b935c-0918-11eb-227a-f7ffbb70c325
begin
	import Pkg; Pkg.activate(@__DIR__)
	using HyperGen, NPZ
end

# ╔═╡ 2f42420e-7bf8-11eb-2529-e3a58c1ee4b4
md"""
# Generate HYPER design and save to file
"""

# ╔═╡ 3912dc07-f4a1-4385-992f-f33536e167b9
md"""
## Generate design
"""

# ╔═╡ 9268d880-09a9-11eb-2d79-af4d3eaa85f8
A = hyperdesign(10,6,3)

# ╔═╡ 4e501d9d-4cef-41d4-a110-3e5e97c185a5
md"""
## Save to `*.npy` file
"""

# ╔═╡ c2fde6ba-0993-11eb-066b-6b54080e7bc1
npzwrite("outs.npy", convert(Matrix{Int}, A))

# ╔═╡ 4754db29-a42d-4dab-a499-5e68e66665be
npzread("outs.npy")

# ╔═╡ Cell order:
# ╟─2f42420e-7bf8-11eb-2529-e3a58c1ee4b4
# ╠═be6b935c-0918-11eb-227a-f7ffbb70c325
# ╟─3912dc07-f4a1-4385-992f-f33536e167b9
# ╠═9268d880-09a9-11eb-2d79-af4d3eaa85f8
# ╟─4e501d9d-4cef-41d4-a110-3e5e97c185a5
# ╠═c2fde6ba-0993-11eb-066b-6b54080e7bc1
# ╠═4754db29-a42d-4dab-a499-5e68e66665be
