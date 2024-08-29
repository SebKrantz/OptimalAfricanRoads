################################
# Real Parameterization
################################

using DataFrames, CSV, LinearAlgebra, Statistics, StatsBase, Plots, NaNStatistics

# **************************
# **** Parameterization ****
# **************************

calib_data = CSV.read("data/QSE/QSE_model_calibration_data_ctry_min_imp.csv", DataFrame)
# TODO: Redding (2016) inly proves existence and uniqueness of equilibrium for symmetric distances. OSRM gives asymmetries. Check robustness to this.
durations = CSV.read("data/QSE/QSE_model_durations_matrix.csv", DataFrame)
# Check
string.(calib_data.cell) == names(durations)

# # Travel Time Across Borders based on World Bank Data
# border_time = CSV.read("data/QSE/QSE_model_border_time_mat_transit.csv", DataFrame)
# # Check
# iso3c = border_time.iso3c
# unique(calib_data.ISO3) == iso3c
# btime = Matrix{Float64}(border_time[:, 2:end])
# # Durations Border
# durations_borders = Matrix{Float64}(durations)
# for i in 1:length(iso3c), j in 1:length(iso3c) 
#     durations_borders[calib_data.ISO3 .== iso3c[i], calib_data.ISO3 .== iso3c[j]] .+= btime[i, j] 
# end
durations = Matrix{Float64}(durations)
# nanmean(durations_borders ./ durations)
# See also: https://truckingresearch.org/wp-content/uploads/2022/08/ATRI-Operational-Cost-of-Trucking-2022.pdf
# Xu and Yang 2021 for China: We find that an hour reduction in travel time between capital cities reduces domestic trade costs by 1–1.2%.
# From Ragjelbaum and Schall: 
# Hummels and Schaur (2013) quantifiedthat one additional day in transit is equivalent to 0.6 to 2.1 percent tariff, and 
# Djankov, Freund, and Pham (2010) argued that each additional day of delay is equivalent to a country distancing 70 km from its trade partner. 


# Descriptives & Plot
describe(calib_data)
# calib_data_continental = dropmissing(calib_data, :wage_gmean)
scatter(calib_data.lon, calib_data.lat; 
        zcolor= calib_data.wage_IWI_global, # calib_data.log_avg_rad, 
        markershape = :hexagon, ms = 2.5, msw = 0, mc = :turbo, 
        aspect_ratio=1, size=(700, 700), #  msc = :auto
        title = "Wage", legend = false, colorbar = true,
        xlabel = "Longitude", ylabel = "Latitude")


# **** Trade Costs ****

# Distance Matrix
dist = durations # Using travel time: considers road quality
# Trade costs are a power function of effective distance
dist = 0.13 * sqrt.(dist)    # Naive Solution: Works well! (0.13 was suggested by Claude 2 and gives good results here) 
# -> I'm more interested in spatial patterns than in the exact numbers, and I want to avoid autarky-like results.
# dist = 0.1 .* dist.^0.33   # Redding (2016) suggests 0.33
# dist = 0.0466 * log.(dist) # Tillman Graff estimate for Africa
extrema(dist)
# Own iceberg transport costs are one
dist .+= 1.0 


# **** Parameters ****

# Share of goods in consumption expenditure (1-housing share)
theta = 0.75
# Share of labour in production
alpha = 1.1 # should be greater than 1, labout is the only factor of production
# Elasticity of substitution
sigma = 5
# Congestion in amenity access
beta = -0.1

param = (alpha = alpha, theta = theta, sigma = sigma, beta = beta)


# **** Productivity = GDP per Capita ****

A = calib_data.wage_IWI_global

println("Summary statistics productivities")
println("mean(a) std(a) max(a) min(a)")
[mean(A) std(A) maximum(A) minimum(A)]

# **** Other Observables ****  

# Land area
H = calib_data.H
H /= sum(H)
# Population
L = calib_data.pop_gpw4

# **** Country indices for restricted mobility case ****

fund = (A=A, L=L, H=H, B = ones(length(A)))


iso3c = unique(calib_data.ISO3)
nctry = length(iso3c)

cindl = Vector{Vector{Int64}}(undef, nctry)
for i in 1:nctry
    cindl[i] = findall(calib_data.ISO3 .== iso3c[i])
end
# Test 
for i in 1:nctry
    print(all(calib_data.ISO3[cindl[i]] .== iso3c[i]))
end


# **************************
# **** Simulation ****
# **************************

include("code/5_redding_2016/CRS_model.jl")

w = calib_data.wage_IWI_global
w /= geomean(w) 


results_immob = solve_model_immob(param, fund, dist, 0.001, true)

# # Solve for Productivity and amenities for mobile labour..
# # -> Only immobile case seems to converge here. 
# results_ab = solve_model_ab(param, fund, w, dist, 0.01, true)
# fund.A = results_ab.A
# fund.B = results_ab.B
#
# results_mob = solve_model_free_mob(param, fund, dist, 0.001, true)
#
# results_mob_ctry = solve_model_free_mob_in_ctry(param, fund, dist, 0.001, true)

# Choose which one to analyze: 
res = results_immob

# **************************
# **** Checks ****
# **************************

ratio = res.income ./ res.expend
[mean(ratio) std(ratio) minimum(ratio) maximum(ratio)]

cor(res.w, fund.A)  # Pearsons correlation between productivity and wages

ots = res.dtradesh
[mean(ots) std(ots) minimum(ots) maximum(ots)]

# ****************************************
# **** Compute Quantities of Interest ****
# ****************************************

include("code/5_redding_2016/helpers.jl")

# Land price
r = landprice(param, fund, res)
# Price index
P = pindex(param, fund, res)
# Consumption
C = consumption(param, fund, res)
# Output
Y = output(param, fund, res)
# Real Wage
wR = realw(param, fund, res)
# Utility per Worker
U = utility_per_worker(param, fund, res)
# Welfare  
welf = welfare(param, fund, res)
welf == welfare2(param, fund, res.L, res.dtradesh) # Check welfare functions are equivalent
# Expected utility
EU = expectut(param, fund, res)
# Immobile Gains from Trade
IGFT = immobile_gains_from_trade(param, fund, res.dtradesh)
# Immobile Gains from Trade per Worker
IGFTPW = immobile_gains_from_trade_per_worker(param, fund, res.dtradesh)


# **************************
# **** Plots ****
# **************************


# create a scatter plot using x = calib_data.lon, y = calib_data.lat, and colour is ots

scatter(calib_data.lon, calib_data.lat; zcolor=res.dtradesh, 
        markershape = :hexagon, ms = 2.5, msw = 0, mc = :turbo, 
        aspect_ratio=1, size=(700, 700), #  msc = :auto
        title = "Own Trade Share", legend = false, colorbar = true,
        xlabel = "Longitude", ylabel = "Latitude")

lon = calib_data.lon
lat = calib_data.lat

# using ColorSchemes
function scatter_plot_map(z, title; colorbar = false)
    # zcolor = (z .- minimum(z)) ./ (maximum(z) .- minimum(z)) 
    scatter(lon, lat; zcolor= z, 
            markershape = :circle, ms = 2, msw = 0, 
            mc = :turbo, # turbo[zcolor],
            # markerstrokealpha = 0, # msc =:white, # msc=:auto, # https://github.com/JuliaPlots/Plots.jl/issues/2494
            aspect_ratio=1, # size=(700, 700), #  msc = :auto
            title = title, legend = false, colorbar = colorbar)
end


# MULTI-PANEL FIGURE
p1 = plot(
    scatter_plot_map(fund.A, "Panel A: Productivity"),
    # scatter_plot_map(fund.B, "Panel B: Amenities"),
    scatter_plot_map(log10.(res.L), "Panel B: Log10 Population"),
    scatter_plot_map(res.dtradesh, "Panel C: Own Trade Share"),
    scatter_plot_map(P, "Panel D: Price Index"),
    scatter_plot_map(log10.(wR), "Panel E: Log10 Real Wages"),   
    scatter_plot_map(r, "Panel F: Land Prices"),
    layout=(3, 2), size=(1000, 1500), titlefontsize = 13 # tickfontsize = 10
)
# savefig(p1, "graphs/initial_equil.pdf")

# MULTI-PANEL FIGURE
p2 = plot(
    scatter_plot_map(log10.(C), "Panel A: Log10 Consumption"),
    # scatter_plot_map(log10.(Y), "Panel B: Log10 Output"),
    scatter_plot_map(r, "Panel C: Land Prices"),
    scatter_plot_map(res.w, "Panel D: Wages"),
    scatter_plot_map(wR, "Panel D: Wages"),
    scatter_plot_map(U, "Panel E: Utility per Worker"),
    scatter_plot_map(IGFTPW, "Panel F: Immobile Gains from Trade"),
    layout=(3, 2), size=(1000, 1500), titlefontsize = 13 # tickfontsize = 10
)

Plots.gr_cbar_width[] = 0.01 # https://github.com/JuliaPlots/Plots.jl/issues/2345
p3 = plot(
    scatter_plot_map(res.dtradesh*100, "Panel A: Own Trade Share (%)"; colorbar = true),
    scatter_plot_map((IGFTPW.-1).*100, "Panel B: Immobile Gains from Trade (%)"; colorbar = true),
    layout=(1, 2), size=(1300, 600), titlefontsize = 13, # tickfontsize = 10
    plotattr = Dict(:margin => 0)
)

# rcParams["figure.dpi"] = 300 
savefig(p3, "figures/full_network/gains_from_trade.pdf")