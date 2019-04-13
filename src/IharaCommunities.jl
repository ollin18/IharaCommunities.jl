module IharaCommunities

    function communities(g; matrix = ollin_matrix)
        edgeidmap, m, aristas = mapa(g)
        NBM, edgeidmap = matrix(g)

        the_values = eigvals(NBM)
        the_real = real.(the_values)
        the_imag = imag.(the_values)
        the_norm = sqrt.((the_real.^2)+(the_imag.^2))
        both = hcat(the_imag,the_norm)
        sorted = sortslices(both, dims=1,by=x->x[2],rev=true)
        threshold = sorted[findall(x->x≠0,sorted[:,1]),2][1]

        #  the_values, vectores = eigens(NBM,threshold)
        cuantos, index = num_com(the_values,threshold)
        if length(cuantos) > 1
            vectores = eigvecs(NBM)

            index = index[2:length(index)]
            matriz_embedded = real(vectores[:,index])
            contraida = contraccion(g,index,matriz_embedded,edgeidmap)

            grupos = kmeans(contraida',length(index)+1;init=:kmcen)
            grupos = assignments(grupos)
        else
            grupos = Int64.(ones(nv(g)))
        end
        grupos
    end

end # module
