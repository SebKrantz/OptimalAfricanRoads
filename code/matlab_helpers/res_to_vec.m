function rv = res_to_vec(Ijk, graph)

    n = size(graph, 1);
    rv = zeros(n,1);
    
    for i = 1:n
        rv(i) = (Ijk(graph.from(i), graph.to(i)) + Ijk(graph.to(i), graph.from(i))) / 2;
    end
    
end