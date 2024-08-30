function Ijk = vec_to_res(n, rv, graph)

    k = size(graph, 1);
    Ijk = zeros(n, n);
    
    for i = 1:k
         Ijk(graph.from(i), graph.to(i)) = rv(i);
         Ijk(graph.to(i), graph.from(i)) = rv(i);
    end
end