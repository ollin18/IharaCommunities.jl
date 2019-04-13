#!/usr/bin/env julia

#  function NB_matrix(g)
#      A = Matrix(adjacency_matrix(g))
#      ceros = 0*(A)
#      D = 0*(A)
#      menos = -1 * one(A)
#      for n in 1:nv(g)
#          D[n,n] = degree(g,n)-1
#      end
#      sparse(hcat(vcat(ceros,menos),vcat(D,A)))
#  end

function NB_matrix(g)
    edgeidmap, m, aristas = mapa(g)
    B = zeros(Float64, 2*aristas, 2*aristas)
    for (e,u) in edgeidmap
        i, j = src(e), dst(e)
        eles = neighbors(g,j)
        k = j
        for l in eles
            B[edgeidmap[Edge(k,l)],u] = (kron_δ(j,k)*(1-kron_δ(i,l)))
        end
    end
    return B, edgeidmap
end


function mapa(g)
    edgeidmap = Dict{Edge, Int}()
    aristas = ne(g)
    m = 0
    if is_directed(g)
        for e in edges(g)
            m += 1
            edgeidmap[Edge(e.src,e.dst)] = m
            #  edgeidmap[Edge(e.dst,e.src)] = m + aristas
        end
    else
        for e in edges(g)
            m += 1
            edgeidmap[Edge(e.src,e.dst)] = m
            edgeidmap[Edge(e.dst,e.src)] = m + aristas
        end
    end
    edgeidmap, m, aristas
end

function flux_matrix(g)
    edgeidmap, m, aristas = mapa(g)
    B = zeros(Float64, 2*aristas, 2*aristas)
    for (e,u) in edgeidmap
        i, j = src(e), dst(e)
        eles = neighbors(g,j)
        k = j
        for l in eles
            B[edgeidmap[Edge(k,l)],u] = (kron_δ(j,k)*(1-kron_δ(i,l)))*(1/(degree(g,j)))
        end
    end
    return B, edgeidmap
end

function reluctant_matrix(g)
    edgeidmap, m, aristas = mapa(g)
    B = zeros(Float64, 2*aristas, 2*aristas)
    for (e,u) in edgeidmap
        i, j = src(e), dst(e)
        eles = neighbors(g,j)
        k = j
        for l in eles
            B[edgeidmap[Edge(k,l)],u] = (kron_δ(j,k)*(1-kron_δ(l,i)))+(kron_δ(j,k)*kron_δ(l,i)*inv_degree(g,j))
        end
    end
    return B, edgeidmap
end

function normalized_reluctant(g)
    edgeidmap, m, aristas = mapa(g)
    B = zeros(Float64, 2*aristas, 2*aristas)
    for (e,u) in edgeidmap
        i, j = src(e), dst(e)
        eles = neighbors(g,j)
        k = j
        for l in eles
            B[edgeidmap[Edge(k,l)],u] = ((kron_δ(j,k)*(1-kron_δ(i,l)))+(kron_δ(j,k)*kron_δ(l,i)*inv_degree(g,j)))*(1/((degree(g,i)-1)+inv_degree(g,j)))
        end
    end
    return B, edgeidmap
end

function ollin_matrix(g)
    v = g.weights
    edgeidmap, m, aristas = mapa(g)
    B = zeros(Float64, aristas, aristas)
    for (e,u) in edgeidmap
        i, j = src(e), dst(e)
        eles = outneighbors(g,j)
        k = j
        for l in eles
            if strengthin(g,l)>0
                #  B[edgeidmap[Edge(k,l)],u] = (kron_δ(j,k)*(1-kron_δ(i,l)))*(2/(1/v[i,j]+1/v[k,l]))*(1/strengthin(g,l))
                B[edgeidmap[Edge(k,l)],u] = (kron_δ(j,k)*(1-kron_δ(i,l)))*((v[i,j]+v[k,l])/2)*(1/strengthin(g,l))
                #  B[edgeidmap[Edge(k,l)],u] = (2/(1/v[i,j]+1/v[k,l]))*(1/strengthin(g,l))
            else
                B[edgeidmap[Edge(k,l)],u] = 0
            end
        end
    end
    return B, edgeidmap
end

function ollin_reluctant(g)
    v = g.weights
    edgeidmap, m, aristas = mapa(g)
    B = zeros(Float64, aristas, aristas)
    for (e,u) in edgeidmap
        i, j = src(e), dst(e)
        eles = outneighbors(g,j)
        k = j
        for l in eles
            if strengthin(g,l)>0
                B[edgeidmap[Edge(k,l)],u] = (kron_δ(j,k)*(1-kron_δ(i,l)))*((v[i,j]+v[k,l])/2)*(1/strengthin(g,l)) + rel_single(g,k)
            else
                B[edgeidmap[Edge(k,l)],u] = 0
            end
        end
    end
    return B, edgeidmap
end
