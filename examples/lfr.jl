using Pkg
function useit(list::Array{Symbol})
        installed = [key for key in keys(Pkg.installed())]
        strpackages = @. string(list)
        uninstalled = setdiff(strpackages,installed)

        map(Pkg.add,uninstalled)
        for package ∈ list
            @eval using $package
        end
end

packages = [:OhMyREPL, :DelimitedFiles, :SNAPDatasets, :Clustering, :LightGraphs, :SimpleWeightedGraphs, :SparseArrays, :Random, :StatsBase, :GraphPlot, :Plots, :GraphRecipes, :Statistics, :LinearAlgebra, :Arpack]
useit(packages)

function lectura(red)
    Nodes = union(unique(red[:,1]),unique(red[:,2]))
    g = SimpleWeightedDiGraph()
    last_node = Int64(length(Nodes))
    add_vertices!(g,last_node)
    for n in 1:size(red)[1]
        add_edge!(g,red[n,1],red[n,2],red[n,3])
    end
    return g, Nodes
end

function lectura_uw(red)
    Nodes = union(unique(red[:,1]),unique(red[:,2]))
    g = SimpleGraph()
    last_node = Int64(length(Nodes))
    add_vertices!(g,last_node)
    for n in 1:size(red)[1]
        add_edge!(g,red[n,1],red[n,2])
        add_edge!(g,red[n,2],red[n,1])
    end
    return g, Nodes
end

function average_vmeasure(times,step)
	vihara = Array{Tuple{Float64,Float64}}(undef,0)
	vnbm = Array{Tuple{Float64,Float64}}(undef,0)
	vflux = Array{Tuple{Float64,Float64}}(undef,0)
	vreluct = Array{Tuple{Float64,Float64}}(undef,0)
	vnorm = Array{Tuple{Float64,Float64}}(undef,0)
	for muti in 0.01:step:1
		vm_ihara = Array{Float64}(undef,0)
		vm_nbm = Array{Float64}(undef,0)
		vm_flux = Array{Float64}(undef,0)
		vm_reluct = Array{Float64}(undef,0)
		vm_norm = Array{Float64}(undef,0)
		for i in 1:times
			run(`./lfr_benchmark/weighted_directed_nets/benchmark -N 100 -k 20 -mut $muti -muw 0.01 -maxk 30`)

			redecita = readdlm("network.dat")
			true_com = readdlm("community.dat")
			true_com = Int64.(true_com)

			g, Nodes = lectura(redecita)
            g1, Nodes = lectura_uw(redecita)
            ihara = communities(g)
            nbm = communities(g1;matrix=NB_matrix)
            flux = communities(g1;matrix=flux_matrix)
            reluctant = communities(g1;matrix=reluctant_matrix)
            norm_rel = communities(g1;matrix=normalized_reluctant)

			push!(vm_ihara,vmeasure(true_com[:,2],ihara))
			push!(vm_nbm,vmeasure(true_com[:,2],nbm))
			push!(vm_flux,vmeasure(true_com[:,2],flux))
			push!(vm_reluct,vmeasure(true_com[:,2],reluctant))
			push!(vm_norm,vmeasure(true_com[:,2],norm_rel))
		end
		push!(vihara,mean_and_std(vm_ihara))
		push!(vnbm,mean_and_std(vm_nbm))
		push!(vflux,mean_and_std(vm_flux))
		push!(vreluct,mean_and_std(vm_reluct))
		push!(vnorm,mean_and_std(vm_norm))
	end
    vihara, vnbm, vflux, vreluct, vnorm
end

vihara, vnbm, vflux, vreluct, vnorm = average_vmeasure(10,0.02)
the_time = 0.01:0.02:1 |> collect

theavg=map(x -> x[1],vihara) |> collect
σ=map(x->x[2],vihara)|> collect

theavg2=map(x -> x[1],vnbm) |> collect
σ2=map(x->x[2],vnbm)|> collect

theavg3=map(x -> x[1],vflux) |> collect
σ3=map(x->x[2],vflux)|> collect

theavg4=map(x -> x[1],vreluct) |> collect
σ4=map(x->x[2],vreluct)|> collect

theavg5=map(x -> x[1],vnorm) |> collect
σ5=map(x->x[2],vnorm)|> collect


plot(the_time,theavg,ribbon=σ,fillalpha=.1,lab="Ihara")
plot!(the_time,theavg2,ribbon=σ2,fillalpha=.1,lab="NBM")
plot!(the_time,theavg3,ribbon=σ3,fillalpha=.1,lab="Flux")
plot!(the_time,theavg4,ribbon=σ4,fillalpha=.1,lab="Reluctant")
plot!(the_time,theavg5,ribbon=σ5,fillalpha=.1,lab="Norm-Rel")
yaxis!("v-measure")
#  plt = twinx()
#  plot!(plt,thetime, dlws,lab="Diameter", color=:red)
title!("v-measure vs topologial mixture parameter")
xaxis!("mixture - p")

savefig("figs/vmeasure.png")







