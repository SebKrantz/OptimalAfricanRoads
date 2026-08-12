#####################################################################
# PIDA Static GE Allocation: Trans-African Evaluation
# -------------------------------------------------------------------
# GE counterpart to the PE analysis in
# `code/12_PIDA_evaluation/PIDA_PE_analysis_trans_african.R`.
#
# On the same 47-largest-port-cities fastest-routes graph used in
# `code/11_GE_simulation_trans_african/optimal_trans_african_networks_largest_pcities_dual_loop.jl`
# (330 condensed corridor edges, 212 nodes), we solve two Fajgelbaum-
# Schaal STATIC allocations per (gamma, rho, sigma) spec (no dual-loop
# optimisation of `Ijk`):
#
#   1. Baseline   `res_stat_base`: I0 = infra_matrix          (observed speeds)
#   2. PIDA       `res_stat_pida`: I0 = infra_matrix_pida     (PIDA links upgraded to >= 100 km/h)
#
# This mirrors the PE exercise in PIDA_PE_analysis_trans_african.R
# (Section 5): the PE code replaces `duration` with `duration_imp` on
# PIDA links and recomputes MA; here we replace observed speeds with
# max(speed, 100 km/h) on PIDA links and recompute the FS equilibrium.
# The welfare / consumption / real-output gain from res_base -> res_pida
# is the GE welfare gain "from just the PIDA links" - a static analogue
# of the optimally allocated PIDA-equivalent budget in the dual-loop
# file (Ki = 6.759e3).
#
# PIDA edge indices are identical to the R sibling files (58 links,
# coverage-based match at buffer=25km, threshold=50%, computed in
# `code/12_PIDA_evaluation/match_PIDA_to_trans_african.R`).
#
# Outputs (per spec):
#   results/transport_network/GE_dual/trans_african/pida/
#     nodes_results_pida_static_22g[_noport]_fixed_cgc[_irs_na]_sigma{s}_rho{r}_duality_julia.csv
#     edges_results_pida_static_22g[_noport]_fixed_cgc[_irs_na]_sigma{s}_rho{r}_duality_julia.csv
#####################################################################

using DataFrames, CSV, LinearAlgebra, Statistics, Plots
using OptimalTransportNetworks
import HSL_jll
include("../helpers/helpers.jl")

# -------------------------------------------------------------------
# 1. Read Undirected Graph and Nodes  (identical to dual-loop file)
# -------------------------------------------------------------------
edges = CSV.read("data/transport_network/trans_african/fastest_routes_graph_edges.csv", DataFrame)
edges.distance /= 1000    # to km
edges.border_dist /= 1000 # to km
edges.duration /= 60      # to hours
edges.total_cost = edges.ug_cost / 1e6  # to millions USD per edge

n = maximum([maximum(edges.from), maximum(edges.to)])

# Adjacency
adj_matrix = falses(n, n)
for i in 1:size(edges, 1)
    adj_matrix[edges.from[i], edges.to[i]] = adj_matrix[edges.to[i], edges.from[i]] = true
end

nodes = CSV.read("data/transport_network/trans_african/fastest_routes_graph_nodes.csv", DataFrame)
nodes.population /= 1000
nodes.outflows /= 1000

# Baseline infrastructure (speed) matrix
infra_matrix = zeros(n, n)
for i in 1:size(edges, 1)
    speed = edges.distance[i] / edges.duration[i]
    infra_matrix[edges.from[i], edges.to[i]] = infra_matrix[edges.to[i], edges.from[i]] = speed
end

# Iceberg trade cost matrix
iceberg_matrix = zeros(n, n)
for i in 1:size(edges, 1)
    iceberg_matrix[edges.from[i], edges.to[i]] = iceberg_matrix[edges.to[i], edges.from[i]] = 0.1158826 * log(edges.distance[i] / 1.609)
end
iceberg_matrix[iceberg_matrix .< 0] .= 0

# Infrastructure building costs (Collier et al. 2016) — total-cost matrix,
# then rescaled to "cost per unit speed increase" as in the reference file.
infra_building_matrix = zeros(n, n)
for i in 1:size(edges, 1)
    infra_building_matrix[edges.from[i], edges.to[i]] = infra_building_matrix[edges.to[i], edges.from[i]] = edges.total_cost[i]
end

# Population, productivity
population = nodes.population
population += (population .== 0) * 1e-6

sum(population .> 2000 .|| nodes.outflows .> 1000) == 47 || error("Expected 47 largest port cities")
productivity = zeros(n, maximum(nodes.product))
with_ports = true
for i in 1:n
    productivity[i, nodes.product[i]] = nodes.IWI[i]
    if with_ports && nodes.outflows[i] > 0
        productivity[i, nodes.product[i]] += (37 * nodes.outflows[i]) / population[i]
    end
end
all(sum(productivity .> 0, dims = 2) .== 1) || error("Each node must have exactly one product")

J = size(productivity, 1)
N = size(productivity, 2)

min_mask = infra_matrix
max_mask = max.(infra_matrix, adj_matrix .* 100)

# Rescale infra_building_matrix so `delta_i` gives cost per unit speed
# increase (matches the dual-loop file's setup).
infra_building_matrix ./= (max_mask - infra_matrix)
infra_building_matrix[isinf.(infra_building_matrix)] .= 0
infra_building_matrix[isnan.(infra_building_matrix)] .= 0
K_base = sum(infra_building_matrix .* infra_matrix) / 2

# -------------------------------------------------------------------
# 2. PIDA edge indices (row indices in `edges` / the CSV)
# -------------------------------------------------------------------
# Verbatim from PIDA_PE_analysis_trans_african.R L41-46: coverage-based
# match at buffer = 25km, threshold = 50%, computed in
# code/12_PIDA_evaluation/match_PIDA_to_trans_african.R.
PIDA_ind = [   7,  16,  17,  18,  23,  24,  26,  52,  53,  59,
              60,  69,  70,  71,  72,  73,  74,  86,  88,  89,
              98, 100, 106, 116, 117, 134, 138, 170, 177, 182,
             183, 185, 187, 189, 205, 212, 213, 216, 217, 218,
             220, 221, 232, 234, 235, 248, 249, 255, 259, 267,
             268, 279, 285, 293, 294, 309, 322, 330]
println("PIDA links: ", length(PIDA_ind), " of ", size(edges, 1),
        " trans-African edges (", round(length(PIDA_ind) / size(edges, 1) * 100, digits = 1), "% of network)")

# Post-PIDA infrastructure matrix: upgrade PIDA links to >= 100 km/h
# (mirrors `w_PIDA = ifelse(edges$is_PIDA, duration_imp, duration)` in
# the R sibling file, since duration_imp corresponds to speeds >= 100).
infra_matrix_pida = copy(infra_matrix)
for i in PIDA_ind
    f = edges.from[i]; t = edges.to[i]
    new_speed = max(infra_matrix_pida[f, t], 100.0)
    infra_matrix_pida[f, t] = infra_matrix_pida[t, f] = new_speed
end

# PIDA package cost (billions USD) — for reporting only; the static
# solver does not enforce a budget. Sanity check: must round to
# $6.759B (= Ki in the dual-loop file at L102), which is how the
# PIDA-equivalent budget was set in the optimally-allocated GE run.
pida_cost_M = sum(edges.total_cost[PIDA_ind])
println("PIDA package cost: \$", round(pida_cost_M / 1e3, digits = 3), "B")
println("PIDA package distance: ",
        round(Int, sum(edges.distance[PIDA_ind])), " km")
@assert round(pida_cost_M / 1e3, digits = 3) == 6.759 "PIDA cost check failed: expected \$6.759B (matches Ki = 6.759e3 in optimal_trans_african_networks_largest_pcities_dual_loop.jl L102), got \$$(round(pida_cost_M / 1e3, digits = 3))B — verify PIDA_ind matches the trans-African edge order."

# -------------------------------------------------------------------
# 3. Model parameters (identical to reference)
# -------------------------------------------------------------------
alpha = 0.7
a = 1
gamma_DRS = 0.946
gamma_IRS = 1.2 # (1.2446 * gamma_DRS)^2 / gamma_DRS  # = beta_TG^2 / gamma_DRS
beta = 1

# Output directory
mkpath("results/transport_network/GE_dual/trans_african/pida")

# -------------------------------------------------------------------
# 4. Static allocations across specs
# -------------------------------------------------------------------
# Same (gamma, rho, sigma) grid as the dual-loop file. For each spec we
# solve TWO static allocations:
#   res_stat_base -> uj_orig, Cj_orig, ... (I0 = infra_matrix)
#   res_stat_pida -> uj,       Cj,       ... (I0 = infra_matrix_pida)
# `solve_allocation = true` short-circuits the dual-loop optimisation
# of `Ijk` and just solves the FS equilibrium at the given I0.
# -------------------------------------------------------------------
for gamma in [gamma_IRS, gamma_DRS], rho in [0, 2], sigma in [3.8, 2]
    print("\n\n=== sigma = ", sigma,
          "  gamma = ", gamma == gamma_IRS ? "IRS" : "DRS",
          "  rho = ", rho, " ===\n")

    param = init_parameters(annealing = false, labor_mobility = false,
                            cross_good_congestion = true, duality = true,
                            a = a, sigma = sigma, N = N,
                            alpha = alpha, beta = beta,
                            gamma = gamma, rho = rho,
                            K = K_base * 2,   # placeholder; unused for solve_allocation
                            tol = 1e-5, min_iter = 20, max_iter = 60)

    graph = create_graph(param, type = "custom",
                         x = nodes.lon, y = nodes.lat,
                         adjacency = adj_matrix,
                         Lj = population, Zjn = productivity,
                         Hj = population .* (1 - alpha))
    graph[:delta_i]   = infra_building_matrix
    graph[:delta_tau] = iceberg_matrix

    param[:optimizer_attr] = Dict(:hsllib => HSL_jll.libhsl_path,
                                  :linear_solver => "ma57",
                                  :tol => 1e-5, :max_iter => 500)

    # Baseline (observed infrastructure)
    print("\n-- Baseline static allocation --\n")
    @time res_stat_base = optimal_network(param, graph, I0 = infra_matrix,
                                          verbose = true, solve_allocation = true)

    # PIDA-upgraded infrastructure
    print("\n-- PIDA static allocation --\n")
    @time res_stat_pida = optimal_network(param, graph, I0 = infra_matrix_pida,
                                          verbose = true, solve_allocation = true)

    # -------- Report gains --------------
    Cj_base = sum(res_stat_base[:Cj])
    Cj_pida = sum(res_stat_pida[:Cj])
    if rho == 2   # rho=2 utility is negative-reciprocal (see analyze_trans_african_results.R)
        u_base = (res_stat_base[:uj] .* -1) .^ -1
        u_pida = (res_stat_pida[:uj] .* -1) .^ -1
        wg = sum(res_stat_pida[:Lj] .* u_pida) /
             sum(res_stat_base[:Lj] .* u_base) - 1
    else
        wg = sum(res_stat_pida[:Lj] .* res_stat_pida[:uj]) /
             sum(res_stat_base[:Lj] .* res_stat_base[:uj]) - 1
    end
    cg = Cj_pida / Cj_base - 1
    println(" Welfare gain (util-weighted, pop-weighted): ",
            round(wg * 100, digits = 3), "%")
    println(" Consumption gain:                          ",
            round(cg * 100, digits = 3), "%")

    # -------- Save nodes (mirrors dual-loop schema) --------------
    res_nodes = deepcopy(nodes)
    res_nodes.uj_orig  = vec(res_stat_base[:uj])
    res_nodes.Lj_orig  = vec(res_stat_base[:Lj])
    res_nodes.Cj_orig  = vec(res_stat_base[:Cj])
    res_nodes.Dj_orig  = vec(res_stat_base[:Dj])
    res_nodes.PCj_orig = vec(res_stat_base[:PCj])
    res_nodes.uj   = vec(res_stat_pida[:uj])
    res_nodes.Lj   = vec(res_stat_pida[:Lj])
    res_nodes.Cj   = vec(res_stat_pida[:Cj])
    res_nodes.Dj   = vec(res_stat_pida[:Dj])
    res_nodes.PCj  = vec(res_stat_pida[:PCj])
    for k in 1:N
        res_nodes[!, Symbol("Lj_$(k)")] = res_stat_pida[:Ljn][:, k]
        res_nodes[!, Symbol("Dj_$(k)")] = res_stat_pida[:Djn][:, k]
        res_nodes[!, Symbol("Yj_$(k)")] = res_stat_pida[:Yjn][:, k]
        res_nodes[!, Symbol("Pj_$(k)")] = res_stat_pida[:Pjn][:, k]
    end

    tag = "pida_static_22g$(with_ports ? "" : "_noport")_fixed_cgc$(gamma == gamma_IRS ? "_irs_na" : "")_sigma$(sigma)_rho$(rho)_duality_julia"
    res_nodes |> CSV.write("results/transport_network/GE_dual/trans_african/pida/nodes_results_$(tag).csv")

    # -------- Save edges -----------------------
    function res_to_vec(Ijk, edges)
        m = size(edges, 1)
        rv = zeros(m)
        for i in 1:m
            rv[i] = (Ijk[edges.from[i], edges.to[i]] + Ijk[edges.to[i], edges.from[i]]) / 2
        end
        return rv
    end

    res_edges = deepcopy(edges)
    res_edges.Ijk_orig = res_to_vec(infra_matrix,      edges)
    res_edges.Ijk      = res_to_vec(infra_matrix_pida, edges)
    res_edges.is_PIDA  = [i in Set(PIDA_ind) for i in 1:size(edges, 1)]
    for k in 1:N
        res_edges[!, Symbol("Qjk_$(k)_orig")] = res_to_vec(res_stat_base[:Qjkn][:, :, k], edges)
        res_edges[!, Symbol("Qjk_$(k)")]      = res_to_vec(res_stat_pida[:Qjkn][:, :, k], edges)
    end
    res_edges |> CSV.write("results/transport_network/GE_dual/trans_african/pida/edges_results_$(tag).csv")
end

println("\nDone. Wrote nodes/edges CSVs for 8 specs to results/transport_network/GE_dual/trans_african/pida/")
