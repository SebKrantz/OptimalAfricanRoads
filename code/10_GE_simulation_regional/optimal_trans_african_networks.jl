# Please also see performance notes at https://github.com/SebKrantz/OptimalTransportNetworks.jl
using DataFrames, CSV, LinearAlgebra, Statistics, Plots
using OptimalTransportNetworks
include("../helpers/helpers.jl")

# Read Undirected Graph
edges_orig = CSV.read("data/transport_network/csv/graph_orig.csv", DataFrame)
edges_orig.add .= false
# Read Additional Links 
edges_add = CSV.read("data/transport_network/csv/graph_add.csv", DataFrame)
edges_add.add .= true
# Combining
edges = vcat(select(edges_orig, :add, :from, :to, :distance, :duration, :ug_cost_km => :cost_per_km, :border_dist), 
             select(edges_add, :add, :from, :to, :distance, :duration_100kmh => :duration, :cost_km => :cost_per_km, :border_dist))
histogram(edges.cost_per_km, bins=100)

# Adjusting 
edges.distance /= 1000    # Convert to km
edges.border_dist /= 1000 # Convert to km
edges.duration /= 60      # Convert to hours
edges.cost_per_km /= 1e6  # Convert to millions
edges.total_cost = edges.cost_per_km .* edges.distance

n = maximum([maximum(edges.from), maximum(edges.to)])

# Create Adjacency Matrix
adj_matrix = falses(n, n)
for i in 1:size(edges, 1)
    adj_matrix[edges.from[i], edges.to[i]] = adj_matrix[edges.to[i], edges.from[i]] = true
end

# Read Nodes Data
nodes = CSV.read("data/transport_network/csv/graph_nodes.csv", DataFrame)
nodes.population /= 1000 # Convert to thousands
nodes.outflows /= 1000   # Adjust in line with population
describe(nodes)

# Create Infrastructure Matrix: Following Graff (2024) = average speed in km/h: length of route is accounted for in cost function
infra_matrix = zeros(n, n)
for i in 1:size(edges, 1)
    if edges.add[i]
        infra_matrix[edges.from[i], edges.to[i]] = infra_matrix[edges.to[i], edges.from[i]] = 0.0
    else
        speed = edges.distance[i] / edges.duration[i]
        infra_matrix[edges.from[i], edges.to[i]] = infra_matrix[edges.to[i], edges.from[i]] = speed
    end
end

# Create Iceberg Trade Cost Matrix: Following Graff (2024)
# edges.distance += edges.border_dist # uncomment to enable frictions
iceberg_matrix = zeros(n, n)
for i in 1:size(edges, 1)
    iceberg_matrix[edges.from[i], edges.to[i]] = iceberg_matrix[edges.to[i], edges.from[i]] = 0.1158826 * log(edges.distance[i] / 1.609)
end
iceberg_matrix[iceberg_matrix .< 0] .= 0

# Create Infrastructure Building Cost Matrix: Following Graff (2024)
infra_building_matrix = zeros(n, n)
for i in 1:size(edges, 1)
    infra_building_matrix[edges.from[i], edges.to[i]] = infra_building_matrix[edges.to[i], edges.from[i]] = edges.total_cost[i]
end

# Basic characteristics of the economy
population = nodes.population
population += (population .== 0) * 1e-6 # this is needed as ipopt crashes otherwise (need to have population > 0 always), so I add an infinitisimal person

# Productivity Matrix 
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
all(sum(productivity .> 0, dims = 2) .== 1)

size(productivity)
extrema(productivity)
for i in 1:size(productivity, 2)
   println(extrema(productivity[:, i]))
end

J, N = size(productivity)

# Infrastructure Bounds
min_mask = infra_matrix                          
max_mask = max.(infra_matrix, adj_matrix .* 100)

# Cost of all Road Activities (162 billion)
K_inv = sum(infra_building_matrix) / 2 
# Total implied costs
infra_building_matrix ./= (max_mask - infra_matrix)
infra_building_matrix[isinf.(infra_building_matrix)] .= 0
infra_building_matrix[isnan.(infra_building_matrix)] .= 0
K_base = sum(infra_building_matrix .* infra_matrix) / 2 # Baseline budget = cost of existing infrastructure
# Now setting budget
K = (K_base + 50e3) * 2 # 50 billion USD'15 investment volume (*2 because symmetric)

K / sum(infra_building_matrix .* max_mask)
K / sum(infra_building_matrix .* min_mask)
if K > sum(infra_building_matrix .* max_mask)
    error("Infrastructure budget is too large for the network")
end

# Parameters
alpha = 0.7 # Spending share on traded goods in utility (curvature parameter)
gamma = 0.946 # F&S: 0.1 ; TG: 0.946; Parameter governing intensity of congestion in transport
beta_old = 1.2446 * gamma
beta = 1 # 1.2446 * gamma # F&S: 0.13 ; TG: 1.2446 * gamma; Parameter governing returns to scale in infrastructure investment
gamma = beta_old^2/gamma # IRS case
sigma = 3.8 #5; # Elasticity of substitution parameter
a = 1 # F&S: 1; TG: 0.7; Returns to scale to labor in production function Zn * Ln^a
rho = 0 # inequality aversion: not possible to solve model if enabled (= 2) -> set alpha = 0.1 instead and resolve allocation afterwards

# Initialise geography
param = init_parameters(annealing = false, labor_mobility = false, cross_good_congestion = true, duality = true,
                        a = a, sigma = sigma, N = N, alpha = alpha, beta = beta, gamma = gamma, rho = rho, 
                        K = K, tol = 1e-5, min_iter = 15, max_iter = 45, verbose = true)

graph = create_graph(param, type = "custom", x = nodes.lon, y = nodes.lat, adjacency = adj_matrix, 
                     Lj = population, Zjn = productivity, Hj = population .* (1-alpha)); # TG: I normalise this because the general utility function has a (h_j/(1-alpha))^(1-alpha) thing with it

graph[:delta_i] = infra_building_matrix
graph[:delta_tau] = iceberg_matrix

# Naming conventions: if IRS -> "cgc_irs"; if alpha = 0.1 -> add "_alpha01"; if with_ports = false -> add "_noport"; if frictions -> add "_bc" (border cost)
filename = "22g_50b_fixed_irs_sigma38_duality_julia" # adjust if sigma != 1.5
println("File extension: $filename")

# Recommended to use coin HSL linear solvers. See README of OptimalTransportNetworks.jl and Ipopt.jl
param[:optimizer_attr] = Dict(:hsllib => "/usr/local/lib/libhsl.dylib", :linear_solver => "ma57") # Use ma86 for optimal performance on big machine 

# Solve allocation from existing infrastructure
@time res_stat = optimal_network(param, graph, I0 = infra_matrix, solve_allocation = true, verbose = true)

# Solve Optimal Network (this can take long - up to 48h)
@time res_opt = optimal_network(param, graph, I0 = infra_matrix, Il = min_mask, Iu = max_mask, verbose = false)

# Check: should be 1
K / sum(res_opt[:Ijk] .* graph[:delta_i])

# Plot Network
plot_graph(graph, res_opt[:Ijk], height = 800) # , node_sizes = res_opt[:Cj])
plot_graph(graph, res_opt[:Ijk] - infra_matrix, height = 800)

# Saving: Nodes
res_nodes = deepcopy(nodes)
res_nodes.uj_orig = vec(res_stat[:uj])
res_nodes.Lj_orig = vec(res_stat[:Lj])
res_nodes.Cj_orig = vec(res_stat[:Cj])
res_nodes.Dj_orig = vec(res_stat[:Dj])
res_nodes.PCj_orig = vec(res_stat[:PCj])
# res_nodes.welfare = res_opt.welfare; # Same as: sum(res_opt.Lj .* res_opt.uj)
res_nodes.uj = vec(res_opt[:uj])
res_nodes.Lj = vec(res_opt[:Lj])
res_nodes.Cj = vec(res_opt[:Cj])
res_nodes.Dj = vec(res_opt[:Dj])
res_nodes.PCj = vec(res_opt[:PCj])
for n in 1:N
   res_nodes[!, Symbol("Lj_$(n)")] = res_opt[:Ljn][:,n]
   res_nodes[!, Symbol("Dj_$(n)")] = res_opt[:Djn][:,n]
   res_nodes[!, Symbol("Yj_$(n)")] = res_opt[:Yjn][:,n]
   res_nodes[!, Symbol("Pj_$(n)")] = res_opt[:Pjn][:,n]
end
res_nodes |> CSV.write("results/transport_network/GE/total/nodes_results_$(filename).csv")

# Saving: Graph / Edges
res_edges = deepcopy(edges)
res_edges.Ijk_orig = res_to_vec(infra_matrix, edges)
res_edges.Ijk = res_to_vec(res_opt[:Ijk], edges)
for n in 1:N
   res_edges[!, Symbol("Qjk_$(n)")] = res_to_vec(res_opt[:Qjkn][:,:,n], edges)
end
res_edges |> CSV.write("results/transport_network/GE/total/edges_results_$(filename).csv")

