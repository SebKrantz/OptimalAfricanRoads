# Please also see performance notes at https://github.com/SebKrantz/OptimalTransportNetworks.jl
using DataFrames, CSV, LinearAlgebra, Statistics, Plots
using OptimalTransportNetworks
include("../helpers/helpers.jl")

# Read Undirected Graph
graph = CSV.read("data/transport_network/largest_pcities/fastest_routes_graph_edges.csv", DataFrame)
histogram(graph.ug_cost, bins=100)
# Adjusting 
graph.distance /= 1000    # Convert to km
graph.border_dist /= 1000 # Convert to km
graph.duration /= 60      # Convert to hours
graph.total_cost = graph.ug_cost / 1e6  # Convert to millions

n = maximum([maximum(graph.from), maximum(graph.to)])

# Create Adjacency Matrix
adj_matrix = falses(n, n)
for i in 1:size(graph, 1)
    adj_matrix[graph.from[i], graph.to[i]] = adj_matrix[graph.to[i], graph.from[i]] = true
end

# Read Nodes Data
nodes = CSV.read("data/transport_network/largest_pcities/fastest_routes_graph_nodes.csv", DataFrame)
nodes.population /= 1000 # Convert to thousands
nodes.outflows /= 1000   # Adjust in line with population
describe(nodes)

# Create Infrastructure Matrix: Following Graff (2024) = average speed in km/h: length of route is accounted for in cost function
infra_matrix = zeros(n, n)
for i in 1:size(graph, 1)
    speed = graph.distance[i] / graph.duration[i]
    infra_matrix[graph.from[i], graph.to[i]] = infra_matrix[graph.to[i], graph.from[i]] = speed
end
describe(graph.distance ./ graph.duration)

# Create Iceberg Trade Cost Matrix: Following Graff (2024)
# graph.distance += graph.border_dist # uncomment to enable frictions
iceberg_matrix = zeros(n, n)
for i in 1:size(graph, 1)
    iceberg_matrix[graph.from[i], graph.to[i]] = iceberg_matrix[graph.to[i], graph.from[i]] = 0.1158826 * log(graph.distance[i] / 1.609)
end
iceberg_matrix[iceberg_matrix .< 0] .= 0

# From Collier et al. (2016):
infra_building_matrix = zeros(n, n)
for i in 1:size(graph, 1)
    infra_building_matrix[graph.from[i], graph.to[i]] = infra_building_matrix[graph.to[i], graph.from[i]] = graph.total_cost[i]
end

# Basic characteristics of the economy
population = nodes.population;
population += (population .== 0) * 1e-6

# Productivity Matrix 
# Check largest pcities
sum(population .> 2000 .|| nodes.outflows .> 1000) == 47
productivity = zeros(n, maximum(nodes.product))
with_ports = true
for i in 1:n
    productivity[i, nodes.product[i]] = nodes.IWI[i]
    if with_ports && nodes.outflows[i] > 0
        productivity[i, nodes.product[i]] += (37 * nodes.outflows[i]) / population[i] # Productivity of ports
    end
end
extrema(productivity)
all(sum(productivity .> 0, dims = 2) .== 1) # Check armington assumption

J, N = size(productivity)

# Infrastructure Bounds
min_mask = infra_matrix                          
max_mask = max.(infra_matrix, adj_matrix .* 100)  

# Cost of all Road Activities (33 billion) 
K_inv = sum(infra_building_matrix) / 2 
# Total implied costs
infra_building_matrix ./= (max_mask - infra_matrix)
infra_building_matrix[isinf.(infra_building_matrix)] .= 0
infra_building_matrix[isnan.(infra_building_matrix)] .= 0
K_base = sum(infra_building_matrix .* infra_matrix) / 2
# Now setting budget
K = (K_base + 10e3) * 2 # 10 or 20 billion USD'15 investment volume (*2 because symmetric)

K / sum(infra_building_matrix .* max_mask)
K / sum(infra_building_matrix .* min_mask)
if K > sum(infra_building_matrix .* max_mask)
    error("Infrastructure budget is too large for the network")
end

# Parameters
alpha = 0.7 # Spending share on traded goods in utility (curvature parameter)
gamma = 0.946 # F&S: 0.1 ; TG: 0.946; Parameter governing intensity of congestion in transport
beta = 1.2446 * gamma # F&S: 0.13 ; TG: 1.2446 * gamma; Parameter governing returns to scale in infrastructure investment
# gamma = beta^2/gamma # IRS case
sigma = 3.8 #5; # Elasticity of substitution parameter (3.8 = Armington)
a = 1 # F&S: 1; TG: 0.7; Returns to scale to labor in production function Zn * Ln^a
rho = 0 # 2 # inequality aversion

# Initialise geography
param = init_parameters(annealing = false, labor_mobility = false, cross_good_congestion = true, 
                        a = a, sigma = sigma, N = N, alpha = alpha, beta = beta, gamma = gamma, rho = rho, 
                        K = K, tol = 1e-5, min_iter = 15, max_iter = 45, verbose = true)

param, g = create_graph(param, type = "custom", x = nodes.lon, y = nodes.lat, adjacency = adj_matrix, 
                        Lj = population, Zjn = productivity, Hj = population .* (1-alpha)) # TG: I normalise this because the general utility function has a (h_j/(1-alpha))^(1-alpha) thing with it

g[:delta_i] = infra_building_matrix;
g[:delta_tau] = iceberg_matrix;

# Naming conventions: if IRS -> "cgc_irs" or "irs_na" without annealing; if alpha = 0.1 -> add "_alpha01"; if with_ports = false -> add "_noport"; if frictions -> add "_bc" (border cost)
filename = "22g_10b_fixed_cgc_sigma38_rho0_julia" # adjust if sigma != 1.5
println("File extension: $filename")

# Recommended to use coin HSL linear solvers. See README of OptimalTransportNetworks.jl and Ipopt.jl
# param[:optimizer_attr] = Dict(:hsllib => "/usr/local/lib/libhsl.dylib", :linear_solver => "ma57") # Use ma86 for optimal performance on big machine 

# Solve allocation from existing infrastructure
@time res_stat = optimal_network(param, g, I0 = infra_matrix, solve_allocation = true, verbose = true)

# Solve Optimal Network
@time res_opt = optimal_network(param, g, I0 = infra_matrix, Il = min_mask, Iu = max_mask, verbose = false)
# # Run Annealing Separately
# @time res_opt, model, recover_allocation = optimal_network(param, g, I0 = infra_matrix, Il = min_mask, Iu = max_mask, verbose = false, return_model = 2)
# @time res_opt = annealing(param, g, res_opt[:Ijk], final_model = model, recover_allocation = recover_allocation, allocation = res_opt, verbose = true)

# Check: should be 1
K / sum(res_opt[:Ijk] .* g[:delta_i])

# Plot Network
plot_graph(g, res_opt[:Ijk], height = 800) # , node_sizes = res[:Cj])
plot_graph(g, res_opt[:Ijk] - infra_matrix, height = 800)

# Saving: Nodes
res_nodes = deepcopy(nodes)
res_nodes.uj_orig = res_stat[:uj]
res_nodes.Lj_orig = res_stat[:Lj]
res_nodes.Cj_orig = res_stat[:Cj]
res_nodes.Dj_orig = res_stat[:Dj]
res_nodes.PCj_orig = res_stat[:PCj]
# res_nodes.welfare = res_opt.welfare; # Same as: sum(res_opt.Lj .* res_opt.uj)
res_nodes.uj = res_opt[:uj]
res_nodes.Lj = res_opt[:Lj]
res_nodes.Cj = res_opt[:Cj]
res_nodes.Dj = res_opt[:Dj]
res_nodes.PCj = res_opt[:PCj]
for n in 1:N
   res_nodes[Symbol("Lj_$(n)")] = res_opt[:Ljn][:,n]
   res_nodes[Symbol("Dj_$(n)")] = res_opt[:Djn][:,n]
   res_nodes[Symbol("Yj_$(n)")] = res_opt[:Yjn][:,n]
   res_nodes[Symbol("Pj_$(n)")] = res_opt[:Pjn][:,n]
end
res_nodes |> CSV.write("results/transport_network/GE/largest_pcities/nodes_results_$(filename).csv")

# Saving: Graph / Edges
res_graph = deepcopy(graph)
res_graph.Ijk_orig = res_to_vec(infra_matrix, graph)
res_graph.Ijk = res_to_vec(res_opt[:Ijk], graph)
for n in 1:N
    res_graph[!, Symbol("Qjk_$(n)")] = res_to_vec(res_opt[:Qjkn][:,:,n], graph)
end
res_graph |> CSV.write("results/transport_network/GE/largest_pcities/edges_results_$(filename).csv")





