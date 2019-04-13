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

packages = [:OhMyREPL, :SNAPDatasets, :Clustering, :LightGraphs, :SimpleWeightedGraphs, :Random, :GraphPlot, :Plots, :GraphRecipes, :Statistics, :LinearAlgebra, :Arpack]
useit(packages)

include("src/utils.jl")

g1 = SimpleDiGraph(128, 3000)

g = SimpleWeightedDiGraph(g1)

for i ∈ eachindex(g.weights.nzval)
    g.weights.nzval[i] = rand()
end


edgeidmap, m, aristas = mapa(g)
NBM, edgeidmap = ollin_matrix(g)
#  grad = map(x -> strength(g,x), vertices(g))
#  gradi = map(x -> strengthin(g,x), vertices(g))
#  grado = map(x -> strengthout(g,x), vertices(g))
#  threshold = maximum(NBM)

#  the_values, the_vectors = eigens(NBM,threshold)
the_values = eigvals(NBM)
the_real = real.(the_values)
the_imag = imag.(the_values)
the_norm = sqrt.((the_real.^2)+(the_imag.^2))
both = hcat(the_imag,the_norm)
sorted = sortslices(both, dims=1,by=x->x[2],rev=true)
threshold = sorted[findall(x->x≠0,sorted[:,1]),2][1]

scatter(the_real,the_imag,title="Eigenvalues of weighted NBM on
complex plane",aspect_ratio=1)

θ = 0:π/50:4π

r = threshold
#  r = 0.038
the_x = r*cos.(θ)
the_y = r*sin.(θ)
plot!(the_x,the_y,lab="Threshold")
xlabel!("Real")
ylabel!("Imaginary")

savefig("/figs/del_eigval_in_complex.png")


R = NBM
#  threshold=1.74
valores, vectores = eigens(R,threshold)
cuantos, index = num_com(valores,threshold)
cuantos, index = 2,[1,2]

index = index[2:length(index)]
matriz_embedded = real(vectores[:,index])
contraida = contraccion(g,index,matriz_embedded,edgeidmap)

grupos = kmeans(contraida',length(index)+1;init=:kmcen)
grupos = assignments(grupos)

membership = grupos
nodecolor = [colorant"red",colorant"yellow",colorant"blue",colorant"violet",colorant"orange",colorant"green"]
nodefillc =  nodecolor[membership]


graphplot(g)
graphplot(g,method=:stress,markercolor=nodefillc)

Nodes[findall(x->x==1,membership)]
Nodes[findall(x->x==2,membership)]
Nodes[findall(x->x==3,membership)]
Nodes[findall(x->x==4,membership)]
Nodes[findall(x->x==5,membership)]





