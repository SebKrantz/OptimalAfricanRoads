
function res_to_vec(Ijk, graph)
    n = size(graph, 1)
    rv = zeros(n)
    for i in 1:n
        rv[i] = (Ijk[graph.from[i], graph.to[i]] + Ijk[graph.to[i], graph.from[i]]) / 2
    end
    return rv
end

function vec_to_res(n, rv, graph)
    k = size(graph, 1)
    Ijk = zeros(n, n)
    for i in 1:k
        Ijk[graph.from[i], graph.to[i]] = Ijk[graph.to[i], graph.from[i]] = rv[i]
    end
    return Ijk
end
