### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ 95eaf04a-08e5-11eb-043f-27646f1bf207
using Dictionaries, NPZ, StatsPlots, StatsPlots.PlotMeasures, Statistics

# ╔═╡ 849adbca-08e5-11eb-1956-cd02f2a6be7a
md"""
# Figure: impact of imbalance
"""

# ╔═╡ eb7bec9e-5453-11eb-0876-c93037e57988
hypername((n,m,q)) = "HYPER (m=$m, q=$q)"

# ╔═╡ 210d331a-09c3-11eb-06ed-e1f4ddc023fc
md"""
## Load data
"""

# ╔═╡ fd643c2c-08e5-11eb-0d25-014798a505d1
md"""
### Utility functions
"""

# ╔═╡ fe7cf842-08e5-11eb-0d07-41ed2902921a
function loadresults(dir, pattern)
	regmatches = filter(!isnothing, match.(pattern, readdir(dir)))
	pairs = map(regmatches) do regmatch
		design = (n=parse(Int,regmatch[:n]), m=parse(Int,regmatch[:m]), q=parse(Int,regmatch[:q]))
		config = (t=parse(Int,regmatch[:t]), numpos=parse(Int,regmatch[:numpos]), numiters=parse(Int,regmatch[:numiters]), randseed=parse(Int,regmatch[:randseed]))
		result = npzread(joinpath(dir,regmatch.match))
		dictionary([design => dictionary([config=>result])])
	end
	return reduce((d1,d2)->mergewith(merge,d1,d2),pairs)
end

# ╔═╡ c4415a4a-ade7-11eb-16a4-1d071360686f
loadnumstage2(dir; pattern=r"^num_stage_two.n-(?<n>\d+)_m-(?<m>\d+)_q-(?<q>\d+).t-(?<t>\d+).numpos-(?<numpos>\d+).numiters-(?<numiters>\d+).randseed-(?<randseed>\d+).npy$") = loadresults(dir, pattern)

# ╔═╡ 981df61c-ade7-11eb-2524-c77bcf95ebd0
loadsens_perm(dir; pattern=r"^sens_avg_perm.n-(?<n>\d+)_m-(?<m>\d+)_q-(?<q>\d+).t-(?<t>\d+).numpos-(?<numpos>\d+).numiters-(?<numiters>\d+).randseed-(?<randseed>\d+).npy$") = loadresults(dir, pattern)

# ╔═╡ 122151e0-08e6-11eb-18f5-3b277145bc9d
md"""
### HYPER
"""

# ╔═╡ 4d02357e-ade9-11eb-3065-05222f3dfe21
hyp_numstage2 = loadnumstage2("designs/hypergraph")

# ╔═╡ 497b1ba0-ade9-11eb-39ac-3774907d79dd
hyp_sens_perm = loadsens_perm("designs/hypergraph")

# ╔═╡ 2e0cc554-08e6-11eb-0186-bba8c2aee7ec
md"""
### Lexicographic pooling
"""

# ╔═╡ dbc88970-ade9-11eb-2f5d-37beac0ae0de
lex_numstage2 = loadnumstage2("designs/lexicographic")

# ╔═╡ def24668-ade9-11eb-0e67-47bd5c7f64d7
lex_sens_perm = loadsens_perm("designs/lexicographic")

# ╔═╡ 0e4bfb8a-a21d-11eb-07d9-51acf99015b3
md"""
### Consecutive pooling
"""

# ╔═╡ fbb979f6-ade9-11eb-078e-cd247e41d9d0
con_numstage2 = loadnumstage2("designs/consecutive")

# ╔═╡ 03fae046-adea-11eb-1418-fbfb4abf42e0
con_sens_perm = loadsens_perm("designs/consecutive")

# ╔═╡ 1b3b908a-a21d-11eb-1ed1-318a1a5fd572
md"""
### Random assignment design
"""

# ╔═╡ b655856e-aeaa-11eb-2067-9140a3dd8094
rad_seeds = 0:2

# ╔═╡ 203f9d60-a21d-11eb-3d17-47daf91ce8db
rad_numstage2 = map(dictionary(rad_seeds)) do seed
	loadnumstage2("designs/rand,assign.seed-$seed")
end

# ╔═╡ 1ee354b0-adea-11eb-2226-8f78a0301e1c
rad_sens_perm = map(dictionary(rad_seeds)) do seed
	loadsens_perm("designs/rand,assign.seed-$seed")
end

# ╔═╡ 13866a4a-a21d-11eb-18f0-19efe13c34d1
md"""
### Double-pooling
"""

# ╔═╡ 0f4a798e-aead-11eb-05c9-db633b438e4b
dbl_seeds = 0:2

# ╔═╡ 2bb38bce-adea-11eb-3f9e-a1a101f25278
dbl_numstage2 = map(dictionary(dbl_seeds)) do seed
	loadnumstage2("designs/double,pooling.seed-$seed")
end

# ╔═╡ 2b5771e0-adea-11eb-32cd-99191a065b00
dbl_sens_perm = map(dictionary(dbl_seeds)) do seed
	loadsens_perm("designs/double,pooling.seed-$seed")
end

# ╔═╡ b25dd0a0-aeb0-11eb-07b6-2bcfce3c132d
md"""
## Collect results from designs
"""

# ╔═╡ 480a805a-aeaf-11eb-24fc-b7759b08b9a8
viz_numstage2 = dictionary([
	"Consecutive"    => con_numstage2,
	"Lexicographic"  => lex_numstage2,
	["Random assign. (draw $(i+1))" => rad_numstage2[i] for i in rad_seeds]...,
	["Double-pooling (draw $(i+1))" => dbl_numstage2[i] for i in dbl_seeds]...,
	"HYPER"          => hyp_numstage2,
])

# ╔═╡ 380150da-aeaf-11eb-0665-23d1cd1fede3
viz_sens_perm = dictionary([
	"Consecutive"    => con_sens_perm,
	"Lexicographic"  => lex_sens_perm,
	["Random assign. (draw $(i+1))" => rad_sens_perm[i] for i in rad_seeds]...,
	["Double-pooling (draw $(i+1))" => dbl_sens_perm[i] for i in dbl_seeds]...,
	"HYPER"          => hyp_sens_perm,
])

# ╔═╡ 0fd98c50-2203-11eb-0f58-d946dcf484b4
md"""
## Figures
"""

# ╔═╡ ca673ac2-adea-11eb-0289-07b30ddadda9
DES = (n=96,m=16,q=2)

# ╔═╡ 465ff718-aeb3-11eb-0a6b-dd01373ff488
CONFIG = (t=80,numpos=1,numiters=100000,randseed=1)

# ╔═╡ a2fd7cae-b7e9-11eb-32d9-21a05872a518
md"""
### Figs a-b
"""

# ╔═╡ f999303e-d7c0-11eb-034e-231a56f3e828
eff(numstage2) = DES.n/(DES.m+numstage2)

# ╔═╡ 99742f0c-b7e9-11eb-32a4-bd2e2de185e7
let
	plts = map(collect(pairs(viz_numstage2))) do (label,numstage2)
		plot(numstage2[DES][CONFIG] .|> eff;
			label="",title=replace(label," ("=>"\n("),titlefontsize=8,
			ylims=(3.7,5.8),yticks=4:0.5:5.5,
			xlims=(1,DES.n),xticks=[1; 32:32:96],
			xguidefontsize=7,
			xlabel="location of\npositive individual",ylabel="individuals/test",
			linewidth=1.5,
		)
	end
	plot!.(plts[2:end],ylabel="",yticks=(4:0.5:5.5,""),left_margin=-4mm)
	annotate!(plts[1],-50,6.21,text("a",14),left_margin=1.5mm)
	plot(plts...,layout=(1,9),size=(900,170),
		bottom_margin=5mm,top_margin=2mm,right_margin=1mm)
	savefig("fig-a.png"); plot!()
end

# ╔═╡ 9d432a66-b7e9-11eb-131f-1507c7391ef5
let
	numstage2list = map(viz_numstage2) do numstage2
		numstage2[DES][CONFIG]
	end
	boxplot(
		reduce(vcat,fill.(1:length(numstage2list),collect(length.(numstage2list)))),
		reduce(vcat,numstage2list) .|> eff;
		label=nothing,
		xlims=(0.5,9.5),
		xticks=(1:length(numstage2list),replace.(keys(numstage2list)," ("=>"\n(")),
		size=(900,170),right_margin=-4mm,bottom_margin=5mm,
		fillopacity=0.5,ylims=(3.7,5.8),
		ylabel="individuals/test",
		bar_width=0.5,
	)
	annotate!(0.060,6,text("b",14),top_margin=3mm,left_margin=1.5mm)
	savefig("fig-b.png"); plot!()
end

# ╔═╡ aaea1eba-b7e9-11eb-39d0-a5b9fd85c26b
md"""
### Figs c-d
"""

# ╔═╡ b0b8bb60-b7e9-11eb-3bb4-69c931bea667
let
	plts = map(collect(pairs(viz_sens_perm))) do (label,sens_perm)
		plot(sens_perm[DES][CONFIG];
			label="",title=replace(label," ("=>"\n("),titlefontsize=8,
			ylims=(0.72,0.76),yticks=0.72:0.01:0.76,
			xlims=(1,DES.n),xticks=[1; 32:32:96],
			xguidefontsize=7,
			xlabel="location of\npositive individual",ylabel="sensitivity",
			linewidth=1.5,
		)
	end
	plot!.(plts[2:end],ylabel="",yticks=(0.72:0.01:0.76,""),left_margin=-4mm)
	annotate!(plts[1],-58,0.77,text("c",14))
	plot(plts...,layout=(1,9),size=(900,170),
		bottom_margin=5mm,top_margin=3mm,right_margin=1mm)
	savefig("fig-c.png"); plot!()
end

# ╔═╡ eb6bae70-b7e9-11eb-22d7-af685cb504fb
let
	sens_permlist = map(viz_sens_perm) do sens_perm
		sens_perm[DES][CONFIG]
	end
	boxplot(
		reduce(vcat,fill.(1:length(sens_permlist),collect(length.(sens_permlist)))),
		reduce(vcat,sens_permlist);
		label=nothing,
		xlims=(0.5,9.5),
		xticks=(1:length(sens_permlist),replace.(keys(viz_sens_perm)," ("=>"\n(")),
		size=(900,170),right_margin=-4mm,bottom_margin=5mm,
		fillopacity=0.5,ylims=(0.72,0.76),
		ylabel="sensitivity",
		bar_width=0.5,
	)
	annotate!(-0.01,0.763,text("d",14),top_margin=3mm)
	savefig("fig-d.png"); plot!()
end

# ╔═╡ Cell order:
# ╟─849adbca-08e5-11eb-1956-cd02f2a6be7a
# ╠═95eaf04a-08e5-11eb-043f-27646f1bf207
# ╟─eb7bec9e-5453-11eb-0876-c93037e57988
# ╟─210d331a-09c3-11eb-06ed-e1f4ddc023fc
# ╟─fd643c2c-08e5-11eb-0d25-014798a505d1
# ╟─fe7cf842-08e5-11eb-0d07-41ed2902921a
# ╟─c4415a4a-ade7-11eb-16a4-1d071360686f
# ╟─981df61c-ade7-11eb-2524-c77bcf95ebd0
# ╟─122151e0-08e6-11eb-18f5-3b277145bc9d
# ╟─4d02357e-ade9-11eb-3065-05222f3dfe21
# ╟─497b1ba0-ade9-11eb-39ac-3774907d79dd
# ╟─2e0cc554-08e6-11eb-0186-bba8c2aee7ec
# ╟─dbc88970-ade9-11eb-2f5d-37beac0ae0de
# ╟─def24668-ade9-11eb-0e67-47bd5c7f64d7
# ╟─0e4bfb8a-a21d-11eb-07d9-51acf99015b3
# ╟─fbb979f6-ade9-11eb-078e-cd247e41d9d0
# ╟─03fae046-adea-11eb-1418-fbfb4abf42e0
# ╟─1b3b908a-a21d-11eb-1ed1-318a1a5fd572
# ╠═b655856e-aeaa-11eb-2067-9140a3dd8094
# ╟─203f9d60-a21d-11eb-3d17-47daf91ce8db
# ╟─1ee354b0-adea-11eb-2226-8f78a0301e1c
# ╟─13866a4a-a21d-11eb-18f0-19efe13c34d1
# ╠═0f4a798e-aead-11eb-05c9-db633b438e4b
# ╟─2bb38bce-adea-11eb-3f9e-a1a101f25278
# ╟─2b5771e0-adea-11eb-32cd-99191a065b00
# ╟─b25dd0a0-aeb0-11eb-07b6-2bcfce3c132d
# ╟─480a805a-aeaf-11eb-24fc-b7759b08b9a8
# ╟─380150da-aeaf-11eb-0665-23d1cd1fede3
# ╟─0fd98c50-2203-11eb-0f58-d946dcf484b4
# ╠═ca673ac2-adea-11eb-0289-07b30ddadda9
# ╠═465ff718-aeb3-11eb-0a6b-dd01373ff488
# ╟─a2fd7cae-b7e9-11eb-32d9-21a05872a518
# ╟─f999303e-d7c0-11eb-034e-231a56f3e828
# ╟─99742f0c-b7e9-11eb-32a4-bd2e2de185e7
# ╟─9d432a66-b7e9-11eb-131f-1507c7391ef5
# ╟─aaea1eba-b7e9-11eb-39d0-a5b9fd85c26b
# ╟─b0b8bb60-b7e9-11eb-3bb4-69c931bea667
# ╟─eb6bae70-b7e9-11eb-22d7-af685cb504fb
