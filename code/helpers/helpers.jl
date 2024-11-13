
function res_to_vec(Ijk, edges)
    n = size(edges, 1)
    rv = zeros(n)
    for i in 1:n
        rv[i] = (Ijk[edges.from[i], edges.to[i]] + Ijk[edges.to[i], edges.from[i]]) / 2
    end
    return rv
end

function vec_to_res(n, rv, edges)
    k = size(edges, 1)
    Ijk = zeros(n, n)
    for i in 1:k
        Ijk[edges.from[i], edges.to[i]] = Ijk[edges.to[i], edges.from[i]] = rv[i]
    end
    return Ijk
end
