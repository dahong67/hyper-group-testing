### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ b53dd37e-0dc5-11eb-37ee-555eed0b2c9f
using CSV, OffsetArrays, Plots, Plots.Measures, PyCall, SplitApplyCombine, Statistics

# ╔═╡ 4592edbe-2edd-11eb-2ed0-679711b2d714
# Code to load CSC matrix from NPZ file
# based on description here: https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html
# see also:
# + https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.save_npz.html
# + https://docs.julialang.org/en/v1/stdlib/SparseArrays/#man-csc

module NPZ_CSC
using PyCall: pyimport, pyrepr
using SparseArrays: SparseMatrixCSC

np = pyimport("numpy")

function load(filename)
	contents = np.load(filename, allow_pickle=true)
	format = get(contents, "format")
	pyrepr(format) == "array(b'csc', dtype='|S3')" ||
		throw("Array format seems wrong. Got $(pyrepr(format))")
	
	shape = get(contents, "shape")
	data = get(contents, "data")
	indptr = get(contents, "indptr") .+ 1
	indices = get(contents, "indices") .+ 1
	
	SparseMatrixCSC(shape[1], shape[2], indptr, indices, data)
end

end

# ╔═╡ eaffbf72-7bed-11eb-2b72-b759220943ef
md"""
# Resource analysis figure

Un/comment the relevant line for each choice of days.
"""

# ╔═╡ e88985de-2f4e-11eb-32cf-1f7a851792bb
# prevalence ≈ 0.033808% - 2.459128%
t0, t1, q1pos, q2pos, q3pos, pbestpos = 40, 90, 8.5, 7, 4, missing

# prevalence ≈ 0.103552%
# t0, t1, q1pos, q2pos, q3pos, pbestpos = 53, 53, 8, 5.5, 3, missing

# prevalence ≈ 1.055872%
# t0, t1, q1pos, q2pos, q3pos, pbestpos = 80, 80, 8.5, 7, 4, missing

# prevalence ≈ 1.362224%
# t0, t1, q1pos, q2pos, q3pos, pbestpos = 83, 83, 8.5, 1, 4, missing

# prevalence ≈ 1.484104%
# t0, t1, q1pos, q2pos, q3pos, pbestpos = 84, 84, 8.5, 1.5, missing, 5

# prevalence ≈ 2.459128%
# t0, t1, q1pos, q2pos, q3pos, pbestpos = 90, 90, 9, 1.5, missing, 5

# prevalence ≈ 3.149712%
# t0, t1, q1pos, q2pos, q3pos, pbestpos = 93, 93, 9, 1.5, 5, missing

# ╔═╡ a17ac18e-7bee-11eb-1624-fde146fd6293
swabsweep = [12,24,48,96,192,384,768,1536,3072,6144]

# ╔═╡ a45a36b4-7bee-11eb-18d7-9f83f59c66fa
testsweep = [12,24,48,96,192,384,768,1536,3072,6144]

# ╔═╡ 6e2887ba-0e95-11eb-1475-e7ea7e5cf740
budgetsweep = collect(Iterators.product(swabsweep, testsweep));

# ╔═╡ 08baf128-2958-11eb-294d-35b2a6eeef5d
selectedcells = [(9,1),(9,2),(9,3),(9,7),(4,2),(6,3)]

# ╔═╡ 9b6802e0-1558-11eb-1f6e-b5f67e4cafd3
methodlist = ["individual","hypergraph","array_comb","pbest"]

# ╔═╡ f7be67a4-7bed-11eb-032c-e9f443225026
md"""
## Utility functions
"""

# ╔═╡ 213f310c-155e-11eb-0720-b748da5a582d
module FindMin
import Base.findmin, Base.findmax, Base.argmin, Base.argmax
isgreater(x, y) = _is_reflexive(x) && _is_reflexive(y) ? isless(y, x) : isless(x, y)
_is_reflexive(x) = let eq = x == x
	isa(eq, Bool) && eq
end
findmax(f, domain) = mapfoldl(x -> (f(x), x), _rf_findmax, domain)
_rf_findmax((fm, m), (fx, x)) = isless(fm, fx) ? (fx, x) : (fm, m)
findmin(f, domain) = mapfoldl(x -> (f(x), x), _rf_findmin, domain)
_rf_findmin((fm, m), (fx, x)) = isgreater(fm, fx) ? (fx, x) : (fm, m)
argmax(f, domain) = findmax(f, domain)[2]
argmin(f, domain) = findmin(f, domain)[2]
end

# ╔═╡ 2e6ae32c-7bee-11eb-12b2-a146668beb80
md"""
## Load results
"""

# ╔═╡ 685aa44c-2edd-11eb-237f-6715b380cb39
ViralLoad = NPZ_CSC.load("../simulation/Simulated_populations/seir_viral_loads_swab.viral_loads.npz")

# ╔═╡ 6ec783d4-2edd-11eb-2be6-7f0631417b11
prevalence = OffsetArray(vec(mean(>(0), ViralLoad, dims=1)), -1)

# ╔═╡ b37b43e4-2edd-11eb-23a1-e1c242eb0029
prev0 = round(prevalence[t0]*100,digits=2)

# ╔═╡ bf7652c6-2edd-11eb-05df-33027bb8ab4e
prev1 = round(prevalence[t1]*100,digits=2)

# ╔═╡ 0faadc44-0dc6-11eb-0c3f-db78ba861a01
results = CSV.File("summary.resource_t0-$(t0)_t1-$(t1).csv")

# ╔═╡ d33a506c-0e2b-11eb-24a6-41ddae38afa5
results_methods = group(res -> res[Symbol("Design method")], results)

# ╔═╡ d1b9df5a-0e2b-11eb-1767-e7762414041c
begin
	results_methods_budget = map(results_methods) do results_method
		map(budgetsweep) do (swabs, tests)
			filter(results_method) do res
				res[Symbol("Sample budget")] == swabs && res[Symbol("Test budget")] == tests
			end
		end
	end
	array_comb = map(results_methods_budget["array_8x12"],results_methods_budget["array_16x24"]) do a1,a2
		a = [a1;a2]
		length(a) <= 1 ? a : [FindMin.argmax(r->r[Symbol("Design effectiveness")],a)]
	end
	insert!(results_methods_budget,"array_comb",array_comb)
end

# ╔═╡ ccdb2506-0e95-11eb-34b4-8db7c3d43ead
efftests = map(results_methods_budget) do res
	map(r -> length(r) == 0 ? 0.0 : only(r)[Symbol("Design effectiveness")], res)
end

# ╔═╡ c3cc7816-2943-11eb-2807-c742e3307662
bestresults = map(Tuple.(CartesianIndices(budgetsweep))) do (swabidx,testidx)
	bestmethod = FindMin.argmax(method->efftests[method][swabidx,testidx], methodlist)
	only(results_methods_budget[bestmethod][swabidx,testidx])
end

# ╔═╡ e3aa4e4a-7be8-11eb-1a29-239b9fbce803
md"""
## Figs a-f
"""

# ╔═╡ 3d2aa366-295f-11eb-24be-f93fa0a0a9de
methodlabels = Dict(
	"individual"=>"Indiv.",
	"hypergraph"=>"HYPER",
	"array_comb"=>"Array",
	"pbest"=>"P-BEST"
)

# ╔═╡ eb95975e-2959-11eb-22af-8dc3c3bd4a3c
selectedplts = map(enumerate(selectedcells)) do (t,(swabidx,testidx))
	result = bestresults[swabidx,testidx]
	method = result[Symbol("Design method")]
	if method == "individual"
		bestq = 0
	elseif method == "hypergraph"
		bestq = result[Symbol("Design sample split")]
	else
		bestq = missing
	end
	ymax = 1.2*result[Symbol("Design effectiveness")]
	
	methodeffs = [efftests[method][swabidx,testidx] for method in methodlist]
	effticks = [0,round(result[Symbol("Design effectiveness")],digits=1)]
	effticklabels = effticks
	bar([methodlabels[method] for method in methodlist],methodeffs,
		bordercolor=:black,
		framestyle=:box,xgrid=nothing,label="",size=(325,250),
		bar_width=0.99,
		ylims=(0,ymax),xlims=(-0.01,4.01),
		ylabel="Effective screening capacity",
		yticks=(effticks,string.(" ",effticklabels)),tick_direction=:out)
	
	swabs, tests = result[Symbol("Sample budget")], result[Symbol("Test budget")]
	title!("$swabs samples, $tests tests")
	
	for (idx,method) in enumerate(methodlist)
		iszero(methodeffs[idx]) && continue
		annotate!([(idx-0.5,methodeffs[idx]+0.05*ymax,Plots.text("$(round(methodeffs[idx],digits=1))",6,:black,:center))])
		
		method == "individual" && continue
		res = only(results_methods_budget[method][swabidx,testidx])
		n = res[Symbol("Design samples")]
		m = res[Symbol("Design pools")]
		q = res[Symbol("Design sample split")]
		b = res[Symbol("Design runs per day")]
		
		if method == "hypergraph"
			if q == 1
				n = convert(Int,n / m)
				b = b*m
				m = 1
			end
			annotate!([(idx-0.5,methodeffs[idx]-0.07ymax,Plots.text("n:$n",5,:white,:center))])
			annotate!([(idx-0.5,methodeffs[idx]-0.15ymax,Plots.text("m:$m",5,:white,:center))])
			annotate!([(idx-0.5,methodeffs[idx]-0.23ymax,Plots.text("q:$q",5,:white,:center))])
		elseif method == "array_comb"
			arraylabel = split(only(results_methods_budget["array_comb"][swabidx,testidx])[Symbol("Design method")],'_')[end]
			annotate!([(idx-0.5,methodeffs[idx]-0.07*ymax,Plots.text("$arraylabel",5,:white,:center))])
		end
		annotate!([(idx-0.5,0.07*ymax,Plots.text("b:$(isinteger(b) ? convert(Int,b) : round(b,digits=1))",5,:black,:center))])
	end
	
	annotate!(-1.05,1.07*ymax,text("ABCDEFGHIJKLMNOPQRSTUVWXYZ"[t] |> lowercase,12))

	plot!()
end |> plts -> plot(plts...,layout=(3,2),link=:x,xtickfontsize=6,ytickfontsize=6,titlefontsize=7,guidefontsize=6,dpi=150,size=(375,500),top_margin=-2mm,bottom_margin=0.5mm,left_margin=0.0mm)

# ╔═╡ cc143562-2967-11eb-06fe-d7f4cd31664d
savefig(selectedplts,"res,analysis,selected,t0-$(t0),t1-$(t1).png")

# ╔═╡ e2fff794-7be8-11eb-0bfb-a17b3f02641e
md"""
## Fig g
"""

# ╔═╡ 8b5a9e38-7be9-11eb-1d71-f5529db9d4e8
HEATMAP_PALETTE = palette(:OrRd_4,4)

# ╔═╡ 422424ee-3032-11eb-2454-972aa1343745
pbestcolor=colorant"forestgreen"

# ╔═╡ fc03d208-5f2e-11eb-366e-45ca469d44ed
arraycolor=colorant"dodgerblue"

# ╔═╡ 17caf5de-2951-11eb-10b4-89609425f75a
bestfig = let cpal = palette([collect(HEATMAP_PALETTE); pbestcolor])
	bestq = map(bestresults) do result
		method = result[Symbol("Design method")]
		if method == "individual"
			return 0
		elseif method == "hypergraph" || method == "dorfman"
			return result[Symbol("Design sample split")]
		elseif method == "pbestdecoder"
			return 4
		else
			return missing
		end
	end
	
	heatmap(permutedims(bestq),colorbar=nothing,color=cpal,clims=(0,4),aspect_ratio=1,
		size=(550,500),dpi=150,right_margin=15mm,grid=nothing,
		xlims=(0.5,length(swabsweep)+0.5),ylims=(0.5,length(testsweep)+0.5),
		xlabel="Daily sample collection capacity (average)",ylabel="Daily test capacity (average)",guidefontsize=10,
		xticks=(axes(swabsweep,1),swabsweep),yticks=(axes(testsweep,1),testsweep),
		framestyle=:box)
	
	if !ismissing(pbestpos)
		annotate!([(10.75,pbestpos, Plots.text("P-BEST",10,:left,cpal[5]))])
	end

	if !ismissing(q1pos)
		annotate!([(10.75,q1pos+0.2, Plots.text("HYPER",10,:left,cpal[2]))])
		annotate!([(10.75,q1pos-0.2, Plots.text("(q=1)",10,:left,cpal[2]))])
	end

	if !ismissing(q2pos)
		annotate!([(10.75,q2pos+0.2, Plots.text("HYPER",10,:left,cpal[3]))])
		annotate!([(10.75,q2pos-0.2, Plots.text("(q=2)",10,:left,cpal[3]))])
	end

	if !ismissing(q3pos)
		annotate!([(10.75,q3pos+0.2, Plots.text("HYPER",10,:left,cpal[4]))])
		annotate!([(10.75,q3pos-0.2, Plots.text("(q=3)",10,:left,cpal[4]))])
	end
	
	bestefftests = getindex.(bestresults,Symbol("Design effectiveness"))
	for idx in CartesianIndices(bestefftests)
		plot!(Shape([(idx[1]-0.5,idx[2]-0.5),(idx[1]+0.5,idx[2]-0.5),(idx[1]+0.5,idx[2]+0.5),(idx[1]-0.5,idx[2]+0.5)]),
			label="",fillcolor=nothing,linewidth=0.5,linealpha=0.2)
	end

	for (t,idx) in enumerate(selectedcells)
		plot!(Shape([(idx[1]-0.5,idx[2]-0.5),(idx[1]+0.5,idx[2]-0.5),(idx[1]+0.5,idx[2]+0.5),(idx[1]-0.5,idx[2]+0.5)]),label="",
			fillcolor=nothing,linecolor=:black,linewidth=2)
	end
	
	for idx in CartesianIndices(bestefftests)
		idx[1] >= idx[2] && annotate!([(idx[1],idx[2]+0.2,Plots.text(round(bestefftests[idx],digits=1),7,:black,:center))])
		if idx[1] > idx[2]
			n, m, q = bestresults[idx][Symbol("Design samples")], bestresults[idx][Symbol("Design pools")], bestresults[idx][Symbol("Design sample split")]
			method = bestresults[idx][Symbol("Design method")]
			if method == "pbestdecoder"
				annotate!([(idx[1],idx[2]-0.20,Plots.text("P-BEST",5,:white,:center))])
				continue
			end
			if q == 1
				n = convert(Int,n / m)
				m = 1
			end
			annotate!([(idx[1],idx[2]-0.11,Plots.text("n=$n",5,ismissing(bestq[idx]) ? :black : :white,:center))])
			annotate!([(idx[1],idx[2]-0.29,Plots.text("m=$m",5,ismissing(bestq[idx]) ? :black : :white,:center))])
		end
	end
	for idx in axes(bestresults,1)[begin:end-1]
		plot!([(idx,idx+0.35),(idx,idx+1.25)],arrow=1,color=:black,label="")
	end
	
	annotate!([(3.5    ,7.5    ,Plots.text("Testing-rich regime:",     12,:black,:center,rotation=45))])
	annotate!([(3.5+0.4,7.5-0.4,Plots.text("individual testing optimal",12,:black,:center,rotation=45))])
	
	annotate!(-0.5,10.95,text("ABCDEFGHIJKLMNOPQRSTUVWXYZ"[length(selectedcells)+1] |> lowercase,12))

	if t0 == t1
		plot!(title="Day $t0: prevalence of $prev0%",titlefontsize=8)
	else
		plot!(title="Days $t0-$t1: prevalence grows exponentially from $prev0% to $prev1%",titlefontsize=8)
	end
	plot!(size=(550,500),bottom_margin=-8mm,left_margin=-2mm,top_margin=-8mm)
end

# ╔═╡ 20d151f0-2951-11eb-124b-b124fe93a1c8
savefig(bestfig,"res,analysis,bestmethods,t0-$(t0),t1-$(t1).png")

# ╔═╡ 4526a4f6-7bef-11eb-17e6-61b5bd0a2f00
md"""
## Fig h
"""

# ╔═╡ 85b7b0d6-1962-11eb-33c2-b759b672eac8
resfig = map(Tuple.(CartesianIndices(budgetsweep))) do (swabidx,testidx)
	swabs, tests = budgetsweep[swabidx,testidx]
	bestmethod = FindMin.argmax(method->efftests[method][swabidx,testidx], methodlist)
	bestresult = only(results_methods_budget[bestmethod][swabidx,testidx])
	n, m, q = bestresult[Symbol("Design samples")], bestresult[Symbol("Design pools")], bestresult[Symbol("Design sample split")]
	if q == 1
		n = convert(Int,n / m)
		m = 1
	end
	besteff = bestresult[Symbol("Design effectiveness")]
	hyperq = only(results_methods_budget["hypergraph"][swabidx,testidx])[Symbol("Design sample split")]
	bar([efftests[method][swabidx,testidx] for method in methodlist],
		color=[HEATMAP_PALETTE[1],HEATMAP_PALETTE[hyperq+1],arraycolor,pbestcolor],
		bar_width=0.99,xlims=(0.47,4.53),
		size=(100,100), ylims=(0,1.3*besteff), label="", framestyle=:box,grid=nothing,
		xticks=(testidx == 1 ? (2.5,swabs) : nothing),yticks=(swabidx == 1 ? (0.65*besteff,tests) : nothing),
	)
	hline!([besteff],color=:black,opacity=0.2,label="")
	annotate!([(2.5,1.16*besteff,Plots.text(round(besteff,digits=1),7,:black,:center))])
	annotate!([(2,efftests["hypergraph"][swabidx,testidx]/2, text("n:$n, m:$m", 5, :white, :center, rotation=90))])
	efftests["individual"][swabidx,testidx] > 0.25*besteff && annotate!([(1,efftests["individual"][swabidx,testidx]/2, text("Ind.", 5, :black, :center, rotation=90))])
	efftests["pbest"][swabidx,testidx] > 0 && annotate!([(4,efftests["pbest"][swabidx,testidx]/2, text("P-BEST", 5, :white, :center, rotation=90))])
	if efftests["array_comb"][swabidx,testidx] > 0
		arraylabel = split(only(results_methods_budget["array_comb"][swabidx,testidx])[Symbol("Design method")],'_')[end]
		annotate!([(3,efftests["array_comb"][swabidx,testidx]/2, text(arraylabel, 5, :white, :center, rotation=90))])
	end
	xlabel!(swabidx == 5 && testidx == 1 ? "Daily sample collection capacity (average)" : "",labelfontsize=10)
	ylabel!(testidx == 5 && swabidx == 1 ? "Daily test capacity (average)" : "",labelfontsize=10)
	
	swabidx == 1 && testidx == 10 && !(t0 == 40 && t1 == 90) && annotate!(-2.4,1.6*besteff,text("ABCDEFGHIJKLMNOPQRSTUVWXYZ"[length(selectedcells)+2] |> lowercase,12))

	if swabidx == 5 && testidx == 10
		if t0 == t1
			plot!(title="Day $t0: prevalence of $prev0%",titlefontsize=10)
		else
			plot!(title="Days $t0-$t1: prevalence grows exponentially from $prev0% to $prev1%",titlefontsize=10)
		end
	end

	plot!(leftmargin=(swabidx == 1 ? 0mm : -3.25mm),rightmargin=0mm,topmargin=-1mm,bottommargin=(testidx == 1 ? 0mm : -2.25mm),guidefontsize=10)
end |> plts->plot(plts[:,end:-1:begin]...,layout=size(budgetsweep),size=(800,700),legend=nothing,dpi=150)

# ╔═╡ e425af4e-1ad2-11eb-3317-e7dcb97dc924
savefig(resfig,"res,analysis,detail,t0-$(t0),t1-$(t1).png")

# ╔═╡ Cell order:
# ╟─eaffbf72-7bed-11eb-2b72-b759220943ef
# ╠═b53dd37e-0dc5-11eb-37ee-555eed0b2c9f
# ╠═e88985de-2f4e-11eb-32cf-1f7a851792bb
# ╠═a17ac18e-7bee-11eb-1624-fde146fd6293
# ╠═a45a36b4-7bee-11eb-18d7-9f83f59c66fa
# ╠═6e2887ba-0e95-11eb-1475-e7ea7e5cf740
# ╠═08baf128-2958-11eb-294d-35b2a6eeef5d
# ╠═9b6802e0-1558-11eb-1f6e-b5f67e4cafd3
# ╟─f7be67a4-7bed-11eb-032c-e9f443225026
# ╠═213f310c-155e-11eb-0720-b748da5a582d
# ╠═4592edbe-2edd-11eb-2ed0-679711b2d714
# ╟─2e6ae32c-7bee-11eb-12b2-a146668beb80
# ╠═685aa44c-2edd-11eb-237f-6715b380cb39
# ╠═6ec783d4-2edd-11eb-2be6-7f0631417b11
# ╠═b37b43e4-2edd-11eb-23a1-e1c242eb0029
# ╠═bf7652c6-2edd-11eb-05df-33027bb8ab4e
# ╠═0faadc44-0dc6-11eb-0c3f-db78ba861a01
# ╠═d33a506c-0e2b-11eb-24a6-41ddae38afa5
# ╠═d1b9df5a-0e2b-11eb-1767-e7762414041c
# ╠═ccdb2506-0e95-11eb-34b4-8db7c3d43ead
# ╠═c3cc7816-2943-11eb-2807-c742e3307662
# ╟─e3aa4e4a-7be8-11eb-1a29-239b9fbce803
# ╟─3d2aa366-295f-11eb-24be-f93fa0a0a9de
# ╠═eb95975e-2959-11eb-22af-8dc3c3bd4a3c
# ╠═cc143562-2967-11eb-06fe-d7f4cd31664d
# ╟─e2fff794-7be8-11eb-0bfb-a17b3f02641e
# ╠═8b5a9e38-7be9-11eb-1d71-f5529db9d4e8
# ╠═422424ee-3032-11eb-2454-972aa1343745
# ╠═fc03d208-5f2e-11eb-366e-45ca469d44ed
# ╠═17caf5de-2951-11eb-10b4-89609425f75a
# ╠═20d151f0-2951-11eb-124b-b124fe93a1c8
# ╟─4526a4f6-7bef-11eb-17e6-61b5bd0a2f00
# ╠═85b7b0d6-1962-11eb-33c2-b759b672eac8
# ╠═e425af4e-1ad2-11eb-3317-e7dcb97dc924
