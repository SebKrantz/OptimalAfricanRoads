#####################################################################
# African Transport Network: Partial Equilibrium Analysis
#####################################################################

library(fastverse)
set_collapse(mask = c("manip", "helper", "special"), nthreads = 4)
fastverse_extend(qs, sf, units, sfnetworks, tmap)
source("code/helpers/helpers.R")
fastverse_conflicts()

# All essential objects from previous sections
# Need only load("data/transport_network/trans_africa_network.RData"), the result of 7_get_transport_network.R
# The .RData file with the _param suffix is created at the end of this file and adds a 'parameterized' version 
# of the network using results computed in this file. It thus just has additional objects and also works as input. 
load("data/transport_network/trans_africa_network_param.RData")

# Average route efficiency per edge
aere <- unattrib(mean(edges$sp_distance / edges$distance)) 
print(aere)

# Adding distance to new edges based on average edge route efficiency
settfm(add_links, sp_distance = unattrib(st_length(geometry)))
settfm(add_links, distance = sp_distance / aere)


# Shortest Paths ----------------------------------------------------------

distances <- st_network_cost(net, weights = edges$distance) # distance in m
n_links <- igraph::distances(net)
range(n_links)

# Checks 
identical(st_geometry(net, "nodes"), nodes$geometry)
dist(unattrib(st_coordinates(nodes$geometry)), unattrib(qM(dist_ttime_mats$sources)))
sp_distances <- st_distance(nodes)

# Network Route Efficiency
nre <- mean(sp_distances / distances, na.rm = TRUE)
rnre <- mean(sp_distances / dist_ttime_mats$distances, na.rm = TRUE) # Real NRE

# Now adding edges
identical(st_geometry(net, "edges"), edges$geometry)
net_ext_data <- rbind(select(edges, distance, geometry), select(add_links, distance, geometry))
net_ext <- as_sfnetwork(net_ext_data, directed = FALSE)
plot(net_ext)
identical(st_geometry(net_ext, "nodes"), nodes$geometry) # Not the case, thus need to recalculate spherical distance as well
ind_ext <- ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net_ext, "nodes"))))
sp_distances_ext <- st_distance(st_geometry(net_ext, "nodes"))[ind_ext, ind_ext]
identical(sp_distances_ext, sp_distances)

distances_ext <- st_network_cost(net_ext, weights = "distance")[ind_ext, ind_ext]
nre_ext <- mean(sp_distances_ext / distances_ext, na.rm = TRUE)

nre_ext / nre
rnre * (nre_ext / nre) # Reported increase
mean(distances / distances_ext, na.rm = TRUE) 

# Per link gain in NRE: takes a few mins
add_links$nre_per_link <- sapply(seq_row(add_links), function(i) {
  net_extd = as_sfnetwork(rbind(select(edges, distance), 
                                subset(add_links, i, distance)), directed = FALSE)
  distances_extd = st_network_cost(net_extd, weights = "distance")
  ind = ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net_extd, "nodes"))))
  mean(sp_distances / distances_extd[ind, ind], na.rm = TRUE)
})
# Percent increase
add_links$nre_gain_perc <- (unattrib(add_links$nre_per_link / nre) - 1) * 100
descr(add_links$nre_gain_perc)

# Plot percent increase
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(add_links) + 
  tm_lines(col = "nre_gain_perc", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.002, 0.005, 0.025, 0.1, Inf)),   
           col.legend = tm_legend(expression(Delta~"%"~"NRE"), position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 

dev.copy(pdf, "figures/transport_network/trans_africa_network_NRE_gain_perc.pdf", width = 10, height = 10)
dev.off()

# Gravity weighted versions
gravity <- replace_inf(tcrossprod(nodes$population) / sp_distances, set = TRUE) |> unclass()
nre_wtd <- fmean(unattrib(sp_distances / distances), w = gravity)
rnre_wtd <- fmean(unattrib(sp_distances / dist_ttime_mats$distances), w = gravity)
nre_ext_wtd <- fmean(unattrib(sp_distances_ext / distances_ext), w = gravity)

rnre_wtd
nre_ext_wtd / nre_wtd
rnre_wtd * (nre_ext_wtd / nre_wtd) # Reported increase

# Per link gain in NRE: weighted: takes a few mins
add_links$nre_wtd_per_link <- sapply(seq_row(add_links), function(i) {
  net_extd = as_sfnetwork(rbind(select(edges, distance), 
                                subset(add_links, i, distance)), directed = FALSE)
  distances_extd = st_network_cost(net_extd, weights = "distance")
  ind = ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net_extd, "nodes"))))
  fmean(unattrib(sp_distances / distances_extd[ind, ind]), w = gravity)
})
# Percent increase
add_links$nre_wtd_gain_perc <- (unattrib(add_links$nre_wtd_per_link / nre_wtd) - 1) * 100
descr(add_links$nre_wtd_gain_perc)

# Plot percent increase
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(add_links) + 
  tm_lines(col = "nre_wtd_gain_perc", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.002, 0.005, 0.025, 0.1, Inf)),   
           col.legend = tm_legend(expression(Delta~"%"~"NRE WTD"), position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 

dev.copy(pdf, "figures/transport_network/trans_africa_network_NRE_wtd_gain_perc.pdf", width = 10, height = 10)
dev.off()



# Market Access -------------------------------------------------------------------------------

fastverse_extend(dggridR)

# Matching to grid, and choosing distance-weiighted nearest cell value within 30km
OUTCOMES <- qread("data/other_inputs/imputed_wealth_GDP_10km_hex.qs")
outcomes_coords <- OUTCOMES |> select(lon, lat) |> qM()
nodes_coord_mat <- st_coordinates(nodes) # graph_nodes |> select(lon, lat) |> qM() 
nodes$gdp_cap <- nodes$IWI <- NA_real_

for (i in seq_row(nodes)) {
  d = geodist::geodist(nodes_coord_mat[i, , drop = FALSE], outcomes_coords, measure = "haversine") 
  ind = which(d < 30e3)
  w = 1 / d[ind]
  nodes$gdp_cap[i] = fmean.default(OUTCOMES$GDP_per_capita_PPP[ind], w = w)
  nodes$IWI[i] = fmean.default(OUTCOMES$IWI[ind], w = w)
}
settfm(nodes, gdp = gdp_cap * population, wealth = IWI * population)
rm(d, ind, w, nodes_coord_mat, outcomes_coords, OUTCOMES); gc()
fndistinct(atomic_elem(nodes))

# Computing total market access
(MA_real <- total_MA(dist_ttime_mats$distances, nodes$gdp))
(MA <- total_MA(distances, nodes$gdp)) # ^3.8

# Total gain
(MA_ext <- total_MA(distances_ext, nodes$gdp)) # ^3.8 

MA_ext / MA

# Needed for later
ma_gain_per_km <- (MA_ext - MA) * 1000

# Compute change in MA from each link
add_links$MA_per_link <- sapply(seq_row(add_links), function(i) {
  nete = as_sfnetwork(rbind(select(edges, distance), 
                            subset(add_links, i, distance)), directed = FALSE)
  inv_dist = 1 / unclass(st_network_cost(nete, weights = "distance"))
  diag(inv_dist) = 0
  ind = ckmatch(mctl(st_coordinates(st_geometry(nete, "nodes"))), nodes_coord)
  sum(inv_dist %*% nodes$gdp[ind])
})
# Percent increase
add_links$MA_gain_perc <- (add_links$MA_per_link / MA - 1) * 100
descr(add_links$MA_gain_perc)

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(add_links) + 
  tm_lines(col = "MA_gain_perc", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.002, 0.005, 0.025, 0.1, Inf)),
           col.legend = tm_legend(expression(Delta~"%"~"MA"), position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 

dev.copy(pdf, "figures/transport_network/trans_africa_network_MA_gain_perc.pdf", width = 10, height = 10)
dev.off()



# Adding trade costs ----------------------------------------------------------------------

border_time <- fread("data/QSE/QSE_model_border_time_mat.csv") |> qM(1)
border_time_sym <- (border_time + t(border_time)) / 2
border_dist <- fread("data/QSE/QSE_model_border_dist_mat.csv") |> qM(1)
border_dist_sym <- (border_dist + t(border_dist)) / 2
border_dist_transit <- fread("data/QSE/QSE_model_border_dist_mat_transit.csv") |> qM(1)
border_time_transit <- fread("data/QSE/QSE_model_border_time_mat_transit.csv") |> qM(1)

# Adding Country Classification
GADM0_africa <- qread("data/other_inputs/GADM0_africa_simplified.qs")
# GADM0_africa <- st_read("/Users/sebastiankrantz/Documents/Data/GADM/gadm_410-levels.gpkg", layer = "ADM_0") %>%
#   subset(GID_0 %in% africamonitor::am_countries$ISO3) %>% st_make_valid()
# GADM0_africa %<>% rmapshaper::ms_simplify(keep = 0.2) %>% st_make_valid()
# mapview::mapView(GADM0_africa)
ctry <- st_within(nodes, GADM0_africa)
table(vlengths(ctry))
ctry[vlengths(ctry) == 0] <- NA
ctry <- as.integer(ctry)
# mapview::mapview(nodes[is.na(ctry), "geometry"])
nodes$iso3c <- GADM0_africa$GID_0[ctry]
nodes$iso3c[is.na(ctry)] <- c("GHA", "TUN")
anyNA(nodes$iso3c)
rm(GADM0_africa, ctry); gc()

# Adding trade costs (symmetric since graph is undirected)
edges$from_ctry <- nodes$iso3c[edges$from]
edges$to_ctry <- nodes$iso3c[edges$to]
edges$border_dist <- sapply(seq_row(edges), function(i) border_dist_sym[edges$from_ctry[i], edges$to_ctry[i]])
edges$border_time <- sapply(seq_row(edges), function(i) border_time_sym[edges$from_ctry[i], edges$to_ctry[i]])

# Checks: pretty large costs, but less for distance...
edges |> qDT() |> extract(border_dist > 0, distance / border_dist) |> descr()
edges |> qDT() |> extract(border_time > 0, duration / border_time) |> descr()

# Same for Additional Routes
add_links$from_ctry <- nodes$iso3c[add_links$from]
add_links$to_ctry <- nodes$iso3c[add_links$to]
add_links$border_dist <- sapply(seq_row(add_links), function(i) border_dist_sym[add_links$from_ctry[i], add_links$to_ctry[i]])
add_links$border_time <- sapply(seq_row(add_links), function(i) border_time_sym[add_links$from_ctry[i], add_links$to_ctry[i]])
# mapview::mapview(select(add_links, from_ctry, to_ctry))

# Now: Repeat Market Access Simulations

# Computing total real market access
bdt_nodes <- border_dist_transit[nodes$iso3c, nodes$iso3c]
MA_bc <- total_MA(distances + bdt_nodes, nodes$gdp) # dist_ttime_mats$distances

# Total gain
MA_ext_bc <- total_MA(distances_ext + bdt_nodes, nodes$gdp) 

MA_ext_bc / MA_bc

# Needed for later
ma_gain_per_km_bc <- (MA_ext_bc - MA_bc) * 1000

# Compute change in MA from each link, with border costs
add_links$MA_per_link_bc <- sapply(seq_row(add_links), function(i) {
  nete = as_sfnetwork(rbind(select(edges, distance), 
                            subset(add_links, i, distance)), directed = FALSE)
  ind = ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(nete, "nodes"))))
  inv_dist = 1 / (unclass(st_network_cost(nete, weights = "distance"))[ind, ind] + bdt_nodes)
  diag(inv_dist) = 0
  sum(inv_dist %*% nodes$gdp)
})
# Percent increase
add_links$MA_gain_perc_bc <- (add_links$MA_per_link_bc / MA_bc - 1) * 100
descr(add_links$MA_gain_perc_bc)

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(add_links) + 
  tm_lines(col = "MA_gain_perc_bc", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.002, 0.005, 0.025, 0.1, Inf)),
           col.legend = tm_legend(expression(Delta~"%"~"MA"), position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 

dev.copy(pdf, "figures/transport_network/trans_africa_network_MA_gain_perc_bc.pdf", width = 10, height = 10)
dev.off()

# Compute Ratios
perch_to_diff <- function(lev, perch) lev - (lev / (perch / 100 + 1))

settfm(add_links, 
       MA_gain_bc_ratio = perch_to_diff(MA_per_link_bc, MA_gain_perc_bc) / perch_to_diff(MA_per_link, MA_gain_perc), 
       MA_gain_perc_bc_ratio = MA_gain_perc_bc / MA_gain_perc)

add_links |> gvr("ratio") |> descr()
add_links$MA_gain_bc_ratio |> replace_outliers(c(0, 1), "clip", set = TRUE)

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(add_links) + 
  tm_lines(col = "MA_gain_bc_ratio", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = seq(0, 1, 0.2)),
           col.legend = tm_legend(expression(Delta~"MA Ratio"), position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 

dev.copy(pdf, "figures/transport_network/trans_africa_network_MA_gain_bc_ratio.pdf", width = 10, height = 10)
dev.off()


# Now doing the same with re-optimization of routing by Agents --------------------------------

# Adding border costs
settfm(add_links, total_dist = distance + border_dist)
settfm(edges, total_dist = distance + border_dist, total_time = duration + border_time)

# Recalculating distances
distances <- st_network_cost(net, weights = edges$distance)
distances_bc <- st_network_cost(net, weights = edges$total_dist) # distance in m
sum(distances_bc) / sum(distances)
sum(distances+bdt_nodes) / sum(distances)
mean(distances_bc / distances, na.rm = TRUE)
# Relative cost
mean(distances) / mean(edges$distance)
mean(distances_bc) / mean(edges$total_dist)

net_ext <- as_sfnetwork(rbind(select(edges, distance, total_dist), 
                              select(add_links, distance, total_dist)), directed = FALSE)
plot(net_ext)
ind_ext <- ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net_ext, "nodes"))))
sp_distances_ext <- st_distance(st_geometry(net_ext, "nodes"))[ind_ext, ind_ext]
identical(sp_distances_ext, sp_distances)

distances_ext <- st_network_cost(net_ext, weights = "distance")[ind_ext, ind_ext]
distances_ext_bc <- st_network_cost(net_ext, weights = "total_dist")[ind_ext, ind_ext]
sum(distances_ext_bc) / sum(distances_ext)
mean(distances_ext_bc / distances_ext, na.rm = TRUE)

# Computing total real market access
MA_bc_opt <- total_MA(distances_bc, nodes$gdp)

# Total gain
MA_ext_bc_opt <- total_MA(distances_ext_bc, nodes$gdp)

MA_ext_bc_opt / MA_bc_opt

# Needed for later
ma_gain_per_km_bc_opt <- (MA_ext_bc_opt - MA_bc_opt) * 1000

# Compute change in MA from each link, with border costs
add_links$MA_per_link_bc_opt <- sapply(seq_row(add_links), function(i) {
  nete = as_sfnetwork(rbind(select(edges, total_dist), 
                            subset(add_links, i, total_dist)), directed = FALSE)
  ind = ckmatch(mctl(st_coordinates(st_geometry(nete, "nodes"))), nodes_coord)
  inv_dist = 1 / unclass(st_network_cost(nete, weights = "total_dist"))
  diag(inv_dist) = 0
  sum(inv_dist %*% nodes$gdp[ind])
})
# Percent increase
add_links$MA_gain_perc_bc_opt <- (add_links$MA_per_link_bc_opt / MA_bc_opt - 1) * 100
descr(add_links$MA_gain_perc_bc_opt)

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(add_links) + 
  tm_lines(col = "MA_gain_perc_bc_opt", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.002, 0.005, 0.025, 0.1, Inf)),
           col.legend = tm_legend(expression(Delta~"%"~"MA"), position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 

dev.copy(pdf, "figures/transport_network/trans_africa_network_MA_gain_perc_bc_opt.pdf", width = 10, height = 10)
dev.off()

# Compute Ratios
settfm(add_links, 
       MA_gain_bc_opt_ratio = replace_outliers(perch_to_diff(MA_per_link_bc_opt, MA_gain_perc_bc_opt) / 
                                               perch_to_diff(MA_per_link, MA_gain_perc), c(0, 3), "clip"), 
       MA_gain_perc_bc_opt_ratio = replace_outliers(MA_gain_perc_bc_opt / MA_gain_perc, c(0, 4), "clip"))

add_links |> gvr("ratio") |> descr()

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(add_links) + 
  tm_lines(col = "MA_gain_bc_opt_ratio", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(seq(0, 1, 0.2), 2, 3)),
           col.legend = tm_legend(expression(Delta~"MA Ratio"), position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 

dev.copy(pdf, "figures/transport_network/trans_africa_network_MA_gain_bc_opt_ratio.pdf", width = 10, height = 10)
dev.off()


# Estimating Network Building Costs -------------------------------------------

# add_links_buff_3km <- st_buffer(add_links, as_units(3000, "m"))
# edges_buff_3km <- st_buffer(edges, as_units(3000, "m"))
#
# # Adding Ruggedness: https://diegopuga.org/data/rugged/
# rugg <- terra::rast("/Users/sebastiankrantz/Documents/Data/Ruggedness/tri.txt")
# # max(rugg)
# add_links$rugg <- exactextractr::exact_extract(rugg, add_links_buff_3km, fun = "mean")
# edges$rugg <- exactextractr::exact_extract(rugg, edges_buff_3km, fun = "mean")
# # Adding Population (WorldPop 2020 1km2 global)
# pop_wpop <- terra::rast("/Users/sebastiankrantz/Documents/Data/WorldPop/africa_pop_2020_1km.tif")
# # max(pop_wpop)
# add_links$pop_wpop <- exactextractr::exact_extract(pop_wpop, add_links_buff_3km, fun = "sum") 
# add_links$pop_wpop_km2 <- unattrib(add_links$pop_wpop / (st_area(add_links_buff_3km) / 1e6))
# edges$pop_wpop <- exactextractr::exact_extract(pop_wpop, edges_buff_3km, fun = "sum")
# edges$pop_wpop_km2 <- unattrib(edges$pop_wpop / (st_area(edges_buff_3km) / 1e6))

# Loading precomputed version
rugg_pop <- qread("data/transport_network/edges_rugg_pop.qs")
add_links %<>% join(rugg_pop$add_links)
edges %<>% join(rugg_pop$edges)
rm(rugg_pop)

# Plot Ruggedness
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(mutate(rbind(select(edges, rugg), select(add_links, rugg)), rugg = rugg / 1000)) +
  tm_lines(col = "rugg",
           col.scale = tm_scale_continuous_log1p(10, values = "turbo"),
           col.legend = tm_legend("Ruggedness", position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)

dev.copy(pdf, "figures/transport_network/trans_africa_network_rugg.pdf", width = 10, height = 10)
dev.off()

add_links |> gvr("_km2") |> descr()
edges |> gvr("_km2") |> descr()

# Plot Population Density
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(rbind(select(edges, pop_wpop_km2), select(add_links, pop_wpop_km2))) +
  tm_lines(col = "pop_wpop_km2",
           col.scale = tm_scale_continuous_log1p(10, values = "turbo"),
           col.legend = tm_legend(expression("Population/km"^2), position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)

dev.copy(pdf, "figures/transport_network/trans_africa_network_pop_wpop_km2.pdf", width = 10, height = 10)
dev.off()

# Loading Updated ROCKS Database: ----------------------------
# fastverse_extend(africamonitor)
# US_DEFL <- am_data("USA", "NY_GDP_DEFL_ZS", expand.date = TRUE, gen = "Year", keep.date = FALSE) |> rename(tolower)
# # G(c(US_DEFL[year == 2006, ny_gdp_defl_zs], US_DEFL[year == 2015, ny_gdp_defl_zs]))
# ROCKS <- readxl::read_xlsx("/Users/sebastiankrantz/Documents/Data/ROCKS/ROCKS-Update-June-2018.xlsx", skip = 2) |> 
#   janitor::clean_names() |> 
#   subset(cost_type == "Actual", -cost_type) |> 
#   join(US_DEFL, on = "year") |> 
#   transformv(c(cost_m_usd, unit_costs_m_usd_per_km), `/`, ny_gdp_defl_zs / 100)
# 
# ROCKS %<>% get_vars(varying(.))
# qsave(ROCKS, "data/ROCKS_2018.qs")

ROCKS <- qread("data/other_inputs/ROCKS_2018.qs")
continental_africa <- fread("data/other_inputs/continental_africa.csv")

options(scipen = 1000)
# Aggreaating for Continental Africa
ROCKS_AFR_AGG <- ROCKS |> 
  mutate(iso3c = countrycode::countryname(country, "iso3c")) |> 
  subset(iso3c %in% continental_africa$ISO3) |> 
  group_by(work_type) |> 
  summarise(N = n(), across(c(length_km, cost_m_usd, unit_costs_m_usd_per_km), fmedian)) 

ROCKS_AFR_table <- ROCKS_AFR_AGG |> 
  subset(N > 10 | startsWith(work_type, "New")) |> 
  roworder(-N) |> 
  mutate(weight = N * length_km) |> 
  colorder(length_km, weight, pos = "after") 

ROCKS_AFR_table |> xtable::xtable(digits = c(0, 0, 0, 1, 0, 2, 3)) |> print(booktabs = TRUE, include.r = FALSE)
ROCKS_AFR_table |> group_by(new_road = startsWith(work_type, "New")) |> 
  num_vars() |> fmean(weight, keep.w = FALSE) |> 
  xtable::xtable() |> print(booktabs = TRUE, include.r = FALSE)

# Aggreaating for all LMI Countries
LMI_ctry <- wlddev |> with(iso3c[income %in% c("Low income", "Lower middle income")]) |> unique() |> as.character()
ROCKS_LMI_AGG <- ROCKS |> 
  mutate(iso3c = countrycode::countryname(country, "iso3c")) |> 
  subset(iso3c %in% LMI_ctry) |> 
  group_by(work_type) |> 
  summarise(N = n(), across(c(length_km, cost_m_usd, unit_costs_m_usd_per_km), fmedian)) 

ROCKS_LMI_table <- ROCKS_LMI_AGG |> 
  subset(work_type %in% ROCKS_AFR_table$work_type | startsWith(work_type, "New")) |> 
  roworder(-N) |> 
  mutate(weight = N * length_km) |> 
  colorder(length_km, weight, pos = "after") 

ROCKS_LMI_table |> xtable::xtable(digits = c(0, 0, 0, 1, 0, 2, 3)) |> print(booktabs = TRUE, include.r = FALSE)
ROCKS_LMI_table |> group_by(new_road = startsWith(work_type, "New")) |> 
  num_vars() |> fmean(weight, keep.w = FALSE) |> 
  xtable::xtable() |> print(booktabs = TRUE, include.r = FALSE)

# All Countries
ROCKS_AGG <- ROCKS |> 
  group_by(work_type) |> 
  summarise(N = n(), across(c(length_km, cost_m_usd, unit_costs_m_usd_per_km), fmedian)) 

ROCKS_table <- ROCKS_AGG |> 
  subset(work_type %in% ROCKS_AFR_table$work_type | startsWith(work_type, "New")) |> 
  roworder(-N) |> 
  mutate(weight = N * length_km) |> 
  colorder(length_km, weight, pos = "after") 

ROCKS_table |> xtable::xtable(digits = c(0, 0, 0, 1, 0, 2, 3)) |> print(booktabs = TRUE, include.r = FALSE)
ROCKS_table |> group_by(new_road = startsWith(work_type, "New")) |> 
  num_vars() |> fmean(weight, keep.w = FALSE) |> 
  xtable::xtable() |> print(booktabs = TRUE, include.r = FALSE)


# Estimates from Collier, Kirchberger & Söderbom (2016)
# Pop Coef
mean(c(0.11, 0.088, 0.082,   # Table 4
       0.077, 0.083, 0.074)) # Taable 5

# Calibrating cost eauation to match median 2L Highway construction cost in Africa (611 million/km)
mean(with(edges, exp(log(120e3) - 0.11 * (distance > 50e3) + 0.12 * log(rugg) + 0.085 * log(pop_wpop_km2+1)))) / 1000 

# Applying Equation
settfm(edges, cost_km = exp(log(120e3) - 0.11 * (distance > 50e3) + 0.12 * log(rugg) + 0.085 * log(pop_wpop_km2+1))) # 0.00085 * pop_wpop_km2
descr(edges$cost_km)

settfm(add_links, cost_km = exp(log(120e3) - 0.11 * (distance > 50e3) + 0.12 * log(rugg) + 0.085 * log(pop_wpop_km2+1))) # 0.00085 * pop_wpop_km2
descr(add_links$cost_km)

descr(rbind(select(edges, cost_km), select(add_links, cost_km)))

# Total Network Length and Cost
sum(edges$distance / 1000) / 1e3
sum(edges$cost_km * edges$distance / 1000) / 1e9
sum(add_links$distance / 1000) / 1e3
sum(add_links$cost_km * add_links$distance / 1000) / 1e9

# Plots
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(rbind(select(edges, cost_km), select(add_links, cost_km))) +
  tm_lines(col = "cost_km",
           col.scale = tm_scale_continuous(10, values = "turbo"), 
           col.legend = tm_legend("USD'15/km", position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)

dev.copy(pdf, "figures/transport_network/trans_africa_network_cost_km.pdf", width = 10, height = 10)
dev.off()

# Algeria-Morocco roads partly exists -> I let them be 3x cheaper
alg_mor <- which((nodes$iso3c[add_links$from] == "DZA" & nodes$iso3c[add_links$to] == "MAR") | (nodes$iso3c[add_links$from] == "MAR" & nodes$iso3c[add_links$to] == "DZA"))
add_links$cost_km_adj <- add_links$cost_km
add_links$cost_km_adj[alg_mor] <- add_links$cost_km_adj[alg_mor] / 3
# rm(alg_mor)

tmap_mode("plot")

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(select(add_links, cost_km_adj)) +
  tm_lines(col = "cost_km_adj",
           col.scale = tm_scale_continuous(10, values = "turbo"), 
           col.legend = tm_legend("USD'15/km", position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)

dev.copy(pdf, "figures/transport_network/trans_africa_network_cost_km_adj.pdf", width = 10, height = 10)
dev.off()



# Cost-Benefit Analysis: Adding New Links -------------------------------------------

# Total Cost-Benefit Ratios
descr(add_links$cost_km)
# AFR_GDP <- africamonitor::am_data(series = "NY_GDP_MKTP_KD", expand.date = TRUE, gen = "Year", keep.date = FALSE) |> 
#            group_by(year = Year) |> summarise(gdp = fsum(NY_GDP_MKTP_KD))
# 
# # GDP growth
# G(AFR_GDP, t = ~year) |> subset(year >= 2000) |> fmedian()

# No Frictions
AFRGDP22 <- 2811259831806 # Africa GDP 2022 in constant 2015 USD
ma_gain_per_km / 1e9 # MA gain in billions
ma_gain_per_km / sum(with(add_links, cost_km_adj * distance / 1000)) # MA gain per investment
# With Frictions
ma_gain_per_km_bc / 1e9 # MA gain in billions
ma_gain_per_km_bc / sum(with(add_links, cost_km_adj * distance / 1000)) # MA gain per investment
# With Frictions and Optimizing Agents
ma_gain_per_km_bc_opt / 1e9 # MA gain in billions
ma_gain_per_km_bc_opt / sum(with(add_links, cost_km_adj * distance / 1000)) # MA gain per investment

# MA Gain per Dollar
settfm(add_links, MA_gain_pusd = perch_to_diff(MA_per_link, MA_gain_perc) * 1000 / (cost_km * distance / 1000)) # * 1216
descr(add_links$MA_gain_pusd)
proportions(table(add_links$MA_gain_pusd < 1))

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(add_links) + 
  tm_lines(col = "MA_gain_pusd", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 1, 2, 5, 10, 20, 50, 100)),
           col.legend = tm_legend(expression(Delta~"MA/USD"), position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 

dev.copy(pdf, "figures/transport_network/trans_africa_network_MA_gain_pusd.pdf", width = 10, height = 10)
dev.off()

# Under Frictions: Static
settfm(add_links, MA_gain_pusd_bc = perch_to_diff(MA_per_link_bc, MA_gain_perc_bc) * 1000 / (cost_km * distance / 1000)) # * 1216
descr(add_links$MA_gain_pusd_bc)

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(add_links) + 
  tm_lines(col = "MA_gain_pusd_bc", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 1, 2, 5, 10, 20, 50, 100)),
           col.legend = tm_legend(expression(Delta~"MA/USD"), position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 

dev.copy(pdf, "figures/transport_network/trans_africa_network_MA_gain_pusd_bc.pdf", width = 10, height = 10)
dev.off()

# Under Frictions: Optimizing Agents
settfm(add_links, MA_gain_pusd_bc_opt = perch_to_diff(MA_per_link_bc_opt, MA_gain_perc_bc_opt) * 1000 / (cost_km * distance / 1000)) # * 1216
descr(add_links$MA_gain_pusd_bc_opt)

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(add_links) + 
  tm_lines(col = "MA_gain_pusd_bc_opt", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 1, 2, 5, 10, 20, 50, 100)),
           col.legend = tm_legend(expression(Delta~"MA/USD"), position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 

dev.copy(pdf, "figures/transport_network/trans_africa_network_MA_gain_pusd_bc_opt.pdf", width = 10, height = 10)
dev.off()

# Consensus Package
settfm(add_links, 
       consensus = MA_gain_pusd > 1 & (MA_gain_pusd_bc > 1 | MA_gain_pusd_bc_opt > 1),
       MA_gain_pusd_cons = pmean(MA_gain_pusd, MA_gain_pusd_bc, MA_gain_pusd_bc_opt))

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(subset(add_links, consensus, MA_gain_pusd_cons)) + 
  tm_lines(col = "MA_gain_pusd_cons", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 1, 2, 5, 10, 20, 50, 100)),
           col.legend = tm_legend(expression(Delta~"MA/USD"), position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 

dev.copy(pdf, "figures/transport_network/trans_africa_network_MA_gain_pusd_cons.pdf", width = 10, height = 10)
dev.off()

# Consensus Gains
nrow(subset(add_links, consensus)) / nrow(add_links)
subset(add_links, consensus) |> with(sum(cost_km_adj * distance / 1000)) |> divide_by(1e9)

net_ext_cons <- as_sfnetwork(rbind(select(edges, distance, total_dist), 
                                   subset(add_links, consensus & seq_along(consensus) %!in% alg_mor, # Without Algeria-Morocco Links
                                          distance, total_dist)), directed = FALSE)
plot(net_ext_cons)
ind_ext_cons <- ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net_ext_cons, "nodes"))))
identical(st_distance(st_geometry(net_ext_cons, "nodes"))[ind_ext_cons, ind_ext_cons], sp_distances)

distances_ext_cons <- st_network_cost(net_ext_cons, weights = "distance")[ind_ext_cons, ind_ext_cons]
distances_ext_bc_cons <- st_network_cost(net_ext_cons, weights = "total_dist")[ind_ext_cons, ind_ext_cons]

# Here need to choose which matrices to use: with/without internal/added border costs: 
# distances_ext_bc_cons <- distances_ext_cons + bdt_nodes
sum(distances_ext_bc_cons) / sum(distances_ext_cons)
mean(distances_ext_bc_cons / distances_ext_cons, na.rm = TRUE)

# Total gain
(MA_ext_cons_bc_opt <- total_MA(distances_ext_bc_cons, nodes$gdp)) # + bdt_nodes

MA_ext_cons_bc_opt / MA_bc_opt

ma_gain_per_km_cons <- (MA_ext_cons_bc_opt - MA_bc_opt) * 1000

ma_gain_per_km_cons / 1e9 # MA gain in billions
ma_gain_per_km_cons / sum(with(subset(add_links, consensus), cost_km_adj * distance / 1000)) # MA gain per investment




# Now: Improving Existing Links ------------------------------------------------------------------------

settfm(edges, speed_kmh = (distance / 1000) / (duration / 60))
descr(edges$speed_kmh)
hist(edges$speed_kmh, breaks = 80, xlab = "Average Link Speed in km/h", main = NULL)
dev.copy(pdf, "figures/transport_network/trans_africa_network_average_link_speed_hist.pdf", width = 8, height = 8)
dev.off()

# Inspect
# edges |> select(speed_kmh) |> mapview::mapview() 

# Plot
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(col = "speed_kmh", 
           col.scale = tm_scale_continuous(values = "turbo"),
           col.legend = tm_legend("Avg. Speed in km/h", position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 

dev.copy(pdf, "figures/transport_network/trans_africa_network_average_link_speed.pdf", width = 10, height = 10)
dev.off()

# Computing Times 
times <- st_network_cost(net, weights = edges$duration)
edges %<>% mutate(speed_kmh_imp = iif(speed_kmh < 100, 100, speed_kmh),
                  duration_imp = duration * speed_kmh / speed_kmh_imp)
descr(edges, cols = .c(duration, duration_imp))
times_imp <- st_network_cost(net, weights = edges$duration_imp)

# Computing total real market access
(MA_real <- total_MA(dist_ttime_mats$durations, nodes$gdp)) # Original: 1748.128 billion USD/min
(MA <- total_MA(times, nodes$gdp)) #  

# Total gain
(MA_imp <- total_MA(times_imp, nodes$gdp))

MA_imp / MA
# Gain from original: 
(MA_imp / MA) * MA_real

# Needed for later
ma_gain_per_min <- MA_imp - MA

edges$MA_100_min_speed <- sapply(seq_row(edges), function(i) {
  w = copyv(edges$duration, i, edges$duration_imp, vind1 = TRUE)
  inv_dur = 1 / unclass(st_network_cost(net, weights = w))
  diag(inv_dur) = 0 
  sum(inv_dur %*% nodes$gdp) 
})
# Percent increase
edges$MA_100_min_speed_perc <- (edges$MA_100_min_speed / MA - 1) * 100
descr(edges$MA_100_min_speed_perc)

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(col = "MA_100_min_speed_perc", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.005, 0.01, 0.025, 0.1, Inf)),
           col.legend = tm_legend(expression(Delta~"%"~"MA [GDP/min]"), 
                                  position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)

dev.copy(pdf, "figures/transport_network/trans_africa_network_MA_100_min_speed_perc.pdf", width = 10, height = 10)
dev.off()

# Considering the addition of proposed links under 100km/h or 65km/h assumption

# Function to convert km/h to meters/min
kmh_to_mmin <- function(speed_kmh) speed_kmh * 1000 / 60 
# Function to convert meters/min to km/h
mmin_to_kmh <- function(speed_mmin) speed_mmin * 60 / 100
settfm(add_links, 
       duration_100kmh = distance / kmh_to_mmin(100), 
       duration_65kmh = distance / kmh_to_mmin(65))

# Temporary networks as needed
net_ext_tmp <- as_sfnetwork(rbind(select(edges, duration = duration_imp), 
                                  select(add_links, duration = duration_100kmh)), directed = FALSE)
ind_ext_tmp <- ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net_ext_tmp, "nodes"))))
times_ext_tmp <- st_network_cost(net_ext_tmp, weights = "duration")[ind_ext_tmp, ind_ext_tmp]

# Total gain
(MA_tmp <- total_MA(times_ext_tmp, nodes$gdp)) / 1e9

MA_tmp / MA # * MA_real
(MA_tmp - MA) / 1e9
rm(list = ls()[endsWith(ls(), "_tmp")]); gc()


# Simulating gains from new 100km/h links under existing and improved network ---------

# Existing Network
add_links$MA_per_link_100kmh <- sapply(seq_row(add_links), function(i) {
  nete = as_sfnetwork(rbind(select(edges, duration), 
                            subset(add_links, i, duration = duration_100kmh)), directed = FALSE)
  ind = ckmatch(mctl(st_coordinates(st_geometry(nete, "nodes"))), nodes_coord)
  inv_dur = 1 / unclass(st_network_cost(nete, weights = "duration"))
  diag(inv_dur) = 0
  sum(inv_dur %*% nodes$gdp[ind])
})
# Percent increase
add_links$MA_per_link_100kmh_perc <- (add_links$MA_per_link_100kmh / MA - 1) * 100
descr(add_links$MA_per_link_100kmh_perc)

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(add_links) + 
  tm_lines(col = "MA_per_link_100kmh_perc", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.005, 0.01, 0.025, 0.1, Inf)),
           col.legend = tm_legend(expression(Delta~"%"~"MA [GDP/min]"), 
                                  position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 

dev.copy(pdf, "figures/transport_network/trans_africa_network_MA_per_link_100kmh_perc.pdf", width = 10, height = 10)
dev.off()

# Improved network
add_links$MA_per_link_100kmh_imp <- sapply(seq_row(add_links), function(i) {
  nete = as_sfnetwork(rbind(select(edges, duration = duration_imp), 
                            subset(add_links, i, duration = duration_100kmh)), directed = FALSE)
  ind = ckmatch(mctl(st_coordinates(st_geometry(nete, "nodes"))), nodes_coord)
  inv_dur = 1 / unclass(st_network_cost(nete, weights = "duration"))
  diag(inv_dur) = 0
  sum(inv_dur %*% nodes$gdp[ind])
})
# Percent increase
add_links$MA_per_link_100kmh_imp_perc <- (add_links$MA_per_link_100kmh_imp / MA_imp - 1) * 100
descr(add_links$MA_per_link_100kmh_imp_perc)

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(add_links) + 
  tm_lines(col = "MA_per_link_100kmh_imp_perc", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.005, 0.01, 0.025, 0.1, Inf)),
           col.legend = tm_legend(expression(Delta~"%"~"MA [GDP/min]"), 
                                  position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 

dev.copy(pdf, "figures/transport_network/trans_africa_network_MA_per_link_100kmh_imp_perc.pdf", width = 10, height = 10)
dev.off()

# Compute Ratios
settfm(add_links, 
       MA_per_link_100kmh_ratio = replace_outliers(perch_to_diff(MA_per_link_100kmh_imp, MA_per_link_100kmh_imp_perc) / 
                                                   perch_to_diff(MA_per_link_100kmh, MA_per_link_100kmh_perc), c(0, 3), "clip"), 
       MA_per_link_100kmh_perc_ratio = replace_outliers(MA_per_link_100kmh_imp_perc / MA_per_link_100kmh_perc, c(0, 4), "clip"))

descr(add_links$MA_per_link_100kmh_perc_ratio)

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(add_links) + 
  tm_lines(col = "MA_per_link_100kmh_ratio", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(seq(0, 1, 0.2), 2, 3)),
           col.legend = tm_legend(expression(Delta~"MA Ratio"), 
                                  position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)

dev.copy(pdf, "figures/transport_network/trans_africa_network_MA_per_link_100kmh_ratio.pdf", width = 10, height = 10)
dev.off()

# Excursus: Check overlap with consensus extension
intersect(add_links |> subset(MA_per_link_100kmh_ratio > 1, id) |> extract2("id"),
          add_links |> subset(consensus, id) |> extract2("id"))



# Adding Border Frictions and Repeating ------------------------------------------

btt_nodes <- border_time_transit[nodes$iso3c, nodes$iso3c]

# Computing total real market access
(MA_bt <- total_MA(times + btt_nodes, nodes$gdp)) / 1e9 # dist_ttime_mats$durations 

MA_bt / MA * MA_real # Reported increase

# Total gain
(MA_imp_bt <- total_MA(times_imp + btt_nodes, nodes$gdp)) / 1e9
MA_imp_bt / MA * MA_real # Reported 
MA_imp_bt / MA_bt # 27% gains, vs. 42% without frictions

# Needed for later
ma_gain_per_min_bt <- MA_imp_bt - MA_bt 

# Per link
edges$MA_100_min_speed_bt <- sapply(seq_row(edges), function(i) {
  w = copyv(edges$duration, i, edges$duration_imp, vind1 = TRUE)
  inv_dur = 1 / (unclass(st_network_cost(net, weights = w)) + btt_nodes)
  diag(inv_dur) = 0 
  sum(inv_dur %*% nodes$gdp) 
})
# Percent increase
edges$MA_100_min_speed_bt_perc <- (edges$MA_100_min_speed_bt / MA_bt - 1) * 100
descr(edges$MA_100_min_speed_bt_perc)

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(col = "MA_100_min_speed_bt_perc", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.005, 0.01, 0.025, 0.1, Inf)),
           col.legend = tm_legend(expression(Delta~"%"~"MA [GDP/min]"), 
                                  position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)

dev.copy(pdf, "figures/transport_network/trans_africa_network_MA_100_min_speed_bt_perc.pdf", width = 10, height = 10)
dev.off()


# Compute Ratio
settfm(edges, 
       MA_100_min_speed_bt_ratio = replace_na(perch_to_diff(MA_100_min_speed_bt, MA_100_min_speed_bt_perc) / perch_to_diff(MA_100_min_speed, MA_100_min_speed_perc), 0), 
       MA_100_min_speed_bt_perc_ratio = MA_100_min_speed_bt_perc / MA_100_min_speed_perc)
descr(edges$MA_100_min_speed_bt_ratio)
edges$MA_100_min_speed_bt_ratio |> replace_outliers(c(0, 1), "clip", set = TRUE)

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(col = "MA_100_min_speed_bt_ratio", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = seq(0, 1, 0.2)), # breaks = c(0, 0.25, 0.5, 1, 1.5, 2),  # For Ratio
           col.legend = tm_legend(expression(Delta~"MA Ratio"), 
                                  position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)

dev.copy(pdf, "figures/transport_network/trans_africa_network_MA_100_min_speed_bt_ratio.pdf", width = 10, height = 10)
dev.off()


# Now doing the same with re-optimization of routing by Agents --------------------------------

# Adding border costs
settfm(edges, total_time = duration + border_time, total_time_imp = duration_imp + border_time)

# Recalculating distances
times_bt <- st_network_cost(net, weights = edges$total_time) # time in min
sum(times_bt) / sum(times)
mean(times_bt / times, na.rm = TRUE)
times_imp_bt <- st_network_cost(net, weights = edges$total_time_imp)
sum(times_imp_bt) / sum(times_imp)
mean(times_imp_bt / times_imp, na.rm = TRUE)

# Computing total real market access
(MA_bt_opt <- total_MA(times_bt, nodes$gdp)) / 1e9

MA_bt_opt / MA * MA_real

# Total gain
(MA_imp_bt_opt <- total_MA(times_imp_bt, nodes$gdp)) / 1e9
MA_imp_bt_opt / MA * MA_real
MA_imp_bt_opt / MA_bt_opt

# Needed for later
ma_gain_per_min_bt_opt <- MA_imp_bt_opt - MA_bt_opt

# Per link
edges$MA_100_min_speed_bt_opt <- sapply(seq_row(edges), function(i) {
  w = copyv(edges$total_time, i, edges$total_time_imp, vind1 = TRUE)
  inv_dur = 1 / unclass(st_network_cost(net, weights = w))
  diag(inv_dur) = 0 
  sum(inv_dur %*% nodes$gdp) 
})
# Percent increase
edges$MA_100_min_speed_bt_opt_perc <- (edges$MA_100_min_speed_bt_opt / MA_bt_opt - 1) * 100
descr(edges$MA_100_min_speed_bt_opt_perc)

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(col = "MA_100_min_speed_bt_opt_perc", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.005, 0.01, 0.025, 0.1, Inf)),
           col.legend = tm_legend(expression(Delta~"%"~"MA [GDP/min]"), 
                                  position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)

dev.copy(pdf, "figures/transport_network/trans_africa_network_MA_100_min_speed_bt_opt_perc.pdf", width = 10, height = 10)
dev.off()


# Compute Ratio
settfm(edges, 
       MA_100_min_speed_bt_opt_ratio = replace_na(perch_to_diff(MA_100_min_speed_bt_opt, MA_100_min_speed_bt_opt_perc) / perch_to_diff(MA_100_min_speed, MA_100_min_speed_perc), 0), 
       MA_100_min_speed_bt_opt_perc_ratio = replace_outliers(MA_100_min_speed_bt_opt_perc / MA_100_min_speed_perc, c(0,4), "clip"))

descr(edges$MA_100_min_speed_bt_opt_ratio)

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(col = "MA_100_min_speed_bt_opt_ratio", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(seq(0, 1, 0.2), 2, Inf)), # breaks = c(0, 0.25, 0.5, 1, 2, 4),  # For Ratio
           col.legend = tm_legend(expression(Delta~"MA Ratio"), 
                                  position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)

dev.copy(pdf, "figures/transport_network/trans_africa_network_MA_100_min_speed_bt_opt_ratio.pdf", width = 10, height = 10)
dev.off()




# Link Upgrading Costs ---------------------------------------------------

settfm(edges, upgrade_cat = nif(speed_kmh < 60, "Upgrade", speed_kmh >= 60 & speed_kmh < 80, "Mixed Works", 
                                speed_kmh >= 80 & speed_kmh < 100, "Asphalt Mix Resurfacing", speed_kmh >= 100, "Nothing") |> 
         factor(levels = c("Nothing", "Asphalt Mix Resurfacing", "Mixed Works", "Upgrade")))
table(edges$upgrade_cat)

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(col = "upgrade_cat", 
           col.scale = tm_scale_categorical(values = "turbo"),
           col.legend = tm_legend("Type of Work", 
                                  position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)

dev.copy(pdf, "figures/transport_network/trans_africa_network_type_of_work.pdf", width = 10, height = 10)
dev.off()

# Now costing the categories
descr(with(edges, # subset(edges, upgrade_cat == "Asphalt Mix Resurfacing"), 
           exp(log(28.4e3) - 0.11 * (distance > 50e3) + 0.12 * log(rugg) + 0.085 * log(pop_wpop_km2+1)))) 
descr(with(edges, # subset(edges, upgrade_cat == "Mixed Works"), 
           exp(log(64.6e3) - 0.11 * (distance > 50e3) + 0.12 * log(rugg) + 0.085 * log(pop_wpop_km2+1)))) 
descr(with(edges, # subset(edges, upgrade_cat == "Upgrade"), 
           exp(log(101.6e3) - 0.11 * (distance > 50e3) + 0.12 * log(rugg) + 0.085 * log(pop_wpop_km2+1)))) 

clip5perc <- function(x) replace_outliers(x, .quantile(x, c(0.025, 0.975), type = 8L), "clip")

edges %<>%
  mutate(ug_cost_km = -0.11 * (distance > 50e3) + 0.12 * log(rugg) + 0.085 * log(pop_wpop_km2+1), 
         ug_cost_km = nif(upgrade_cat == "Asphalt Mix Resurfacing", clip5perc(exp(ug_cost_km + log(28.4e3))),
                          upgrade_cat == "Mixed Works", clip5perc(exp(ug_cost_km + log(64.6e3))),
                          upgrade_cat == "Upgrade", clip5perc(exp(ug_cost_km + log(101.6e3))), 
                          upgrade_cat == "Nothing", 0)) 

descr(edges$ug_cost_km)
descr(edges, ug_cost_km ~ upgrade_cat)
hist(edges$ug_cost_km, breaks = 80)

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(col = "ug_cost_km", 
           col.scale = tm_scale_continuous(7, values = "turbo"), 
           col.legend = tm_legend("USD'15/km", position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)

dev.copy(pdf, "figures/transport_network/trans_africa_network_upgrading_costs.pdf", width = 10, height = 10)
dev.off()

# Total Costs and Breakdown
sum(edges$ug_cost_km * edges$distance / 1000) / 1e9
fsum(edges$ug_cost_km * edges$distance / 1000, edges$upgrade_cat) / 1e9



# Cost-Benefit Analysis ---------------------------------------------------

# No Frictions
ma_gain_per_min / 1e9 # MA gain in billions
ma_gain_per_min / sum(with(edges, ug_cost_km * distance / 1000)) # MA gain per investment
# With Frictions
ma_gain_per_min_bt / 1e9 # MA gain in billions
ma_gain_per_min_bt / sum(with(edges, ug_cost_km * distance / 1000)) # MA gain per investment
# With Frictions and Optimizing Agents
ma_gain_per_min_bt_opt / 1e9 # MA gain in billions
ma_gain_per_min_bt_opt / sum(with(edges, ug_cost_km * distance / 1000)) # MA gain per investment

# MA Gain per Dollar
settfm(edges, MA_gain_pusd = perch_to_diff(MA_100_min_speed, MA_100_min_speed_perc) / (ug_cost_km * distance / 1000)) # * 1216
descr(edges$MA_gain_pusd)
proportions(table(edges$MA_gain_pusd < 1))

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) + 
  tm_lines(col = "MA_gain_pusd", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 1, 2, 5, 10, 20, 50, 100, Inf)),
           col.legend = tm_legend(expression(Delta~"MA/USD"),
                                  position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)

dev.copy(pdf, "figures/transport_network/trans_africa_network_MA_gain_100_min_speed_pusd.pdf", width = 10, height = 10)
dev.off()

# Under Frictions: Static
settfm(edges, MA_gain_pusd_bt = perch_to_diff(MA_100_min_speed_bt, MA_100_min_speed_bt_perc) / (ug_cost_km * distance / 1000)) # * 1216
descr(edges$MA_gain_pusd_bt)
proportions(table(edges$MA_gain_pusd_bt < 1))

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(col = "MA_gain_pusd_bt", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 1, 2, 5, 10, 20, 50, 100, Inf)),
           col.legend = tm_legend(expression(Delta~"MA/USD"),
                                  position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)

dev.copy(pdf, "figures/transport_network/trans_africa_network_MA_gain_100_min_speed_pusd_bt.pdf", width = 10, height = 10)
dev.off()

# Under Frictions: Optimizing Agents
settfm(edges, MA_gain_pusd_bt_opt = perch_to_diff(MA_100_min_speed_bt_opt, MA_100_min_speed_bt_opt_perc) / (ug_cost_km * distance / 1000)) # * 1216
descr(edges$MA_gain_pusd_bt_opt)
proportions(table(edges$MA_gain_pusd_bt_opt < 1))

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(col = "MA_gain_pusd_bt_opt", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 1, 2, 5, 10, 20, 50, 100, Inf)),
           col.legend = tm_legend(expression(Delta~"MA/USD"),
                                  position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)

dev.copy(pdf, "figures/transport_network/trans_africa_network_MA_gain_100_min_speed_pusd_bt_opt.pdf", width = 10, height = 10)
dev.off()

# Consensus Package
settfm(edges, 
       consensus = MA_gain_pusd > 1 & (MA_gain_pusd_bt > 1 | MA_gain_pusd_bt_opt > 1), 
       MA_gain_pusd_cons = pmean(MA_gain_pusd, MA_gain_pusd_bt, MA_gain_pusd_bt_opt))

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(subset(edges, !consensus)) + tm_lines(lwd = 2, col = "grey70") +
  tm_shape(subset(edges, consensus, MA_gain_pusd_cons)) +
  tm_lines(col = "MA_gain_pusd_cons", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 1, 2, 5, 10, 20, 50, 100, Inf)),
           col.legend = tm_legend(expression(Delta~"MA/USD"),
                                  position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)

dev.copy(pdf, "figures/transport_network/trans_africa_network_MA_gain_100_min_speed_pusd_cons.pdf", width = 10, height = 10)
dev.off()

# Consensus Gains
nrow(subset(edges, consensus)) / nrow(edges)
subset(edges, consensus) |> with(sum(ug_cost_km * distance / 1000)) |> divide_by(1e9)

net_imp_cons <- as_sfnetwork(rbind(subset(edges, !consensus, duration, total_time), 
                                   subset(edges, consensus, duration = duration_imp, total_time = total_time_imp)), directed = FALSE)
plot(net_imp_cons)
ind_imp_cons <- ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net_imp_cons, "nodes"))))
identical(st_distance(st_geometry(net_imp_cons, "nodes"))[ind_imp_cons, ind_imp_cons], sp_distances)

times_imp_cons <- st_network_cost(net_imp_cons, weights = "duration")[ind_imp_cons, ind_imp_cons]
times_imp_bt_cons <- st_network_cost(net_imp_cons, weights = "total_time")[ind_imp_cons, ind_imp_cons]
# times_imp_bt_cons <- times_imp_cons + bdt_nodes
sum(times_imp_bt_cons) / sum(times_imp_cons)
mean(times_imp_bt_cons / times_imp_cons, na.rm = TRUE)

# Total gain
MA_imp_cons <- total_MA(times_imp_cons, nodes$gdp)

MA_imp_cons / MA

ma_gain_per_min_cons <- MA_imp_cons - MA

ma_gain_per_min_cons / 1e9 # MA gain in billions
ma_gain_per_min_cons / sum(with(subset(edges, consensus), ug_cost_km * distance / 1000)) # MA gain per investment



# Cost-Benefit Analysis: Joint Scenarios ---------------------------------------------------

settfm(add_links, total_time_100kmh = duration_100kmh + border_time, total_time_65kmh = duration_65kmh + border_time)

# Plot All Costs
all_costs <- rowbind(existing = select(edges, cost_km = ug_cost_km, distance, duration, duration_imp, border_time, total_time, total_time_imp), 
                     new = select(add_links, cost_km = cost_km_adj, distance, duration_imp = duration_100kmh, 
                                  border_time, total_time_imp = total_time_100kmh) |> 
                           transform(duration = duration_imp, total_time = total_time_imp), 
                     idcol = "type")

descr(all_costs$cost_km)
hist(all_costs$cost_km / 1000, breaks = 80, xlab = "Cost per Km in Thousands of 2015 USD", main = NULL)
dev.copy(pdf, "figures/transport_network/trans_africa_network_all_costs_hist.pdf", width = 8, height = 8)
dev.off()

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(all_costs) +
  tm_lines(col = "cost_km", 
           col.scale = tm_scale_continuous(7, values = "turbo"), 
           col.legend = tm_legend("USD'15/km", position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)

dev.copy(pdf, "figures/transport_network/trans_africa_network_all_costs.pdf", width = 10, height = 10)
dev.off()

# Total Costs
sum(all_costs$cost_km * all_costs$distance / 1000) / 1e9

# Need to repeat simulations for new links with MA denominated in time 

# Added Cost Scenario
add_links$MA_per_link_100kmh_bt <- sapply(seq_row(add_links), function(i) {
  nete = as_sfnetwork(rbind(select(edges, duration), 
                            subset(add_links, i, duration = duration_100kmh)), directed = FALSE)
  ind = ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(nete, "nodes"))))
  inv_dur = 1 / (unclass(st_network_cost(nete, weights = "duration"))[ind, ind] + btt_nodes)
  diag(inv_dur) = 0
  sum(inv_dur %*% nodes$gdp)
})
add_links$MA_per_link_100kmh_bt_perc <- (add_links$MA_per_link_100kmh_bt / MA_bt - 1) * 100
descr(add_links$MA_per_link_100kmh_bt_perc)

# Optimizing Agents Scenario
add_links$MA_per_link_100kmh_bt_opt <- sapply(seq_row(add_links), function(i) {
  nete = as_sfnetwork(rbind(select(edges, duration = total_time), 
                            subset(add_links, i, duration = total_time_100kmh)), directed = FALSE)
  ind = ckmatch(mctl(st_coordinates(st_geometry(nete, "nodes"))), nodes_coord)
  inv_dur = 1 / unclass(st_network_cost(nete, weights = "duration"))
  diag(inv_dur) = 0
  sum(inv_dur %*% nodes$gdp[ind])
})
add_links$MA_per_link_100kmh_bt_opt_perc <- (add_links$MA_per_link_100kmh_bt_opt / MA_bt_opt - 1) * 100
descr(add_links$MA_per_link_100kmh_bt_opt_perc)

# Computing Cost-Benefit Ratios
settfm(add_links, 
   MA_gain_100kmh_pusd = perch_to_diff(MA_per_link_100kmh, MA_per_link_100kmh_perc) / (cost_km_adj * distance / 1000),
   MA_gain_100kmh_pusd_bt = perch_to_diff(MA_per_link_100kmh_bt, MA_per_link_100kmh_bt_perc) / (cost_km_adj * distance / 1000),
   MA_gain_100kmh_pusd_bt_opt = perch_to_diff(MA_per_link_100kmh_bt_opt, MA_per_link_100kmh_bt_opt_perc) / (cost_km_adj * distance / 1000)
)

# Combining Datasets
all_cb_ratios <- rbind(edges |> select(MA_gain_pusd, MA_gain_pusd_bt, MA_gain_pusd_bt_opt),
                       add_links |> select(MA_gain_100kmh_pusd, MA_gain_100kmh_pusd_bt, MA_gain_100kmh_pusd_bt_opt) |> rm_stub("100kmh_", regex = TRUE))
tfm(all_cb_ratios) <- all_costs |> atomic_elem()
descr(all_cb_ratios)

for (v in .c(pusd, pusd_bt, pusd_bt_opt)) {
  print(v)
  pl <- tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
    tm_shape(all_cb_ratios) + 
    tm_lines(col = paste0("MA_gain_", v), 
            col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 1, 2, 5, 10, 20, 50, 100, Inf)),
            col.legend = tm_legend(expression(Delta~"MA/USD"),
                                    position = c("left", "bottom"), frame = FALSE, 
                                    text.size = 1.5, title.size = 2), lwd = 2) + 
    tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
    tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
    tm_layout(frame = FALSE)
  print(pl)
  dev.copy(pdf, sprintf("figures/transport_network/trans_africa_network_MA_gain_all_100kmh_%s.pdf", v), 
          width = 10, height = 10)
  dev.off()
}; rm(v)

for (i in c(1:2, 4L)) {
  cat("MA Gain greatern than ", i, fill = TRUE)

# Consensus Package
settfm(all_cb_ratios, 
       consensus = MA_gain_pusd > i & (MA_gain_pusd_bt > i | MA_gain_pusd_bt_opt > i), # 1, 2, or 4
       MA_gain_pusd_cons = pmean(MA_gain_pusd, MA_gain_pusd_bt, MA_gain_pusd_bt_opt))

pl <- tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(subset(all_cb_ratios, !consensus)) + tm_lines(lwd = 2, col = "grey70") +
  tm_shape(subset(all_cb_ratios, consensus, MA_gain_pusd_cons)) +
  tm_lines(col = "MA_gain_pusd_cons", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 1, 2, 5, 10, 20, 50, 100, Inf)),
           col.legend = tm_legend(expression(Delta~"MA/USD"), position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 
print(pl)

dev.copy(pdf, sprintf("figures/transport_network/trans_africa_network_MA_gain_all_100kmh_pusd_cons_MAg%d.pdf", i), width = 10, height = 10)
dev.off()

# Consensus Gains
nrow(subset(all_cb_ratios, consensus)) / nrow(all_cb_ratios)
all_cb_ratios %$% table(type, consensus) |> proportions(1)
# Cost
subset(all_cb_ratios, consensus) |> with(sum(cost_km * distance / 1000)) |> divide_by(1e9) |> print()

net_imp_cons <- as_sfnetwork(rbind(
  subset(all_cb_ratios, consensus, duration = duration_imp, total_time = total_time_imp),
  subset(all_cb_ratios, !consensus & type == "existing", duration, total_time)), directed = FALSE)

# plot(net_imp_cons)
ind_imp_cons <- ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net_imp_cons, "nodes"))))
identical(st_distance(st_geometry(net_imp_cons, "nodes"))[ind_imp_cons, ind_imp_cons], sp_distances)

times_imp_cons <- st_network_cost(net_imp_cons, weights = "duration")[ind_imp_cons, ind_imp_cons]
times_imp_bt_cons <- st_network_cost(net_imp_cons, weights = "total_time")[ind_imp_cons, ind_imp_cons]
# Again adjust frictions scenario. Default: cumulative frictions.
# times_imp_bt_cons <- times_imp_cons + btt_nodes
sum(times_imp_bt_cons) / sum(times_imp_cons)
mean(times_imp_bt_cons / times_imp_cons, na.rm = TRUE)

# Total gain
MA_imp_cons <- total_MA(times_imp_cons, nodes$gdp) # _bt
print(MA_imp_cons / MA)

ma_gain_per_min_cons <- MA_imp_cons - MA

print(ma_gain_per_min_cons / 1e9) # MA gain in billions
# MA gain per investment:
print(ma_gain_per_min_cons / sum(with(subset(all_cb_ratios, consensus), cost_km * distance / 1000))) 
}

# Macroeconomic Cost-Benefit Analysis (Minimal Required Growth Returns) --------------------------------------------

AFRGDP22 <- 2811259831806 # Affrica GDP 2022 in constant 2015 USD

my_PDV <- function(gdp = AFRGDP22, df = 1.1, gr_old_perc = 4.1, gr_new_perc = 4.11, max_years = 30) {
  gr_old = 1 + gr_old_perc / 100
  gr_new = 1 + gr_new_perc / 100
  years = 1:max_years
  FV = gdp * (gr_new^years - gr_old^years)
  PDV = sum(FV / df^years)
  return(PDV)
}

# What minimum new growth required 
inv_PDV <- function(PDV = 40e9, ...) {
  objective <- function(x) abs(PDV - my_PDV(gr_new_perc = x, ...))
  optimize(objective, c(0, 100), tol = .Machine$double.eps)$minimum
}

# Test
inv_PDV(my_PDV())

packages <- c(
  "Full Extension" = 56.1e9, # sum(add_links$cost_km * add_links$distance / 1000) 
  "Consensus Extension" = 12.1e9,
  "Full Upgrade" = 105.7e9,  # sum(edges$ug_cost_km * edges$distance / 1000) 
  "Consensus Upgrade" = 45e9,
  "All Links MA > 1" = 60.9e9, 
  "All Links MA > 2" = 36.6e9, 
  "All Links MA > 4" = 17.0e9
)

calc_rates <- function(x, bgr) {
  gr_new = inv_PDV(x, gr_old_perc = bgr)
  c("Rate" = gr_new, "Growth of Rate" = (gr_new / bgr - 1) * 100)
}

# Optimistic Growth Scenario
sapply(packages, calc_rates, 4.1) |> t() |> round(4)

# Pessimistic Growth Scenario
sapply(packages, calc_rates, 3) |> t() |> round(4)




#############################################
# Saving Graphs for OTN (GE Analysis)
#############################################

qsu(edges)

graph_orig <- edges |> qDT() |> 
  select(from, from_ctry, to, to_ctry, sp_distance, distance, duration, speed_kmh, 
         speed_kmh_imp, duration_imp, border_dist, border_time, total_time, total_time_imp, 
         rugg, pop_wpop, pop_wpop_km2, cost_km, cost_km_adj = cost_km, upgrade_cat, ug_cost_km)


settfm(add_links, total_time_100kmh = duration_100kmh + border_time, total_time_65kmh = duration_65kmh + border_time)

graph_add <- add_links |> qDT() |> 
  select(from, from_ctry, to, to_ctry, sp_distance, distance, duration_100kmh, 
         duration_65kmh, border_dist, border_time, total_dist, total_time_100kmh, total_time_65kmh,
         rugg, pop_wpop, pop_wpop_km2, cost_km, cost_km_adj)

anyNA(cities_ports$city_country)
any_duplicated(nodes$city_country)
sum(nodes$city_port)

graph_nodes <- nodes |> qDT() |> 
  transform(set_names(mctl(st_coordinates(geometry)), c("lon", "lat"))) |> 
  select(lon, lat, iso3c, city_country, city_port, population, gdp_cap, IWI, gdp, wealth) |> 
  join(compute(cities_ports, city_port = TRUE,
               keep = .c(city_country, port_locode, port_name, port_status, outflows)), 
       on = c("city_country", "city_port")) |> 
  mutate(outflows = replace_na(outflows))

# Consistency Checks
identical(graph_orig$from_ctry, graph_nodes$iso3c[graph_orig$from])
identical(graph_orig$to_ctry, graph_nodes$iso3c[graph_orig$to])
identical(graph_add$from_ctry, graph_nodes$iso3c[graph_add$from])
identical(graph_add$to_ctry, graph_nodes$iso3c[graph_add$to])

# Saving
for (name in .c(graph_orig, graph_add, graph_nodes)) {
  sprintf("data/transport_network/csv/%s.csv", name) |> 
    fwrite(x = get(name))
}

# Also Adding the information to the RData file
load("data/transport_network/trans_africa_network.RData")

# Load previous saved graphs
graphs <- sapply(.c(graph_orig, graph_add, graph_nodes), function(name)
  sprintf("data/transport_network/csv/%s.csv", name) |> fread())
graphs$graph_orig$add <- FALSE
graphs$graph_add$add <- TRUE

# Joining
nodes %<>% transform(qDF(round(st_coordinates(.), 5))) %>% 
  join(tfm(graphs$graph_nodes, X = round(lon, 5), Y = round(lat, 5)), 
       on = c("X", "Y", "population", "city_port"), drop = "x", overid = 2) %>% select(-X, -Y)
edges %<>% join(graphs$graph_orig, on = c("from", "to", "distance"), drop = "x", overid = 2)
add_links %<>% join(graphs$graph_add, on = c("from", "to"), drop = "x", overid = 2)

# Check that network aligns with nodes
allv(st_distance(st_geometry(net, "nodes"), st_geometry(nodes), by_element = TRUE), 0)
identical(st_geometry(net, "edges"), st_geometry(edges))

# Add to network
net %<>% activate("nodes") %>% dplyr::mutate(nodes |> atomic_elem() |> qDF())
net %<>% activate("edges") %>% dplyr::mutate(select(edges, -(from:gravity_dur)) |> atomic_elem() |> qDF())

# Save
TAN_env <- new.env()
load("data/transport_network/trans_africa_network.RData", envir = TAN_env)
TAN_env$nodes_param <- nodes
TAN_env$edges_param <- edges
TAN_env$add_links_param <- add_links
TAN_env$net_param <- net
save(list = ls(TAN_env), file = "data/transport_network/trans_africa_network_param.RData", envir = TAN_env)



# Plot population and productivity (for GE Calibration) ------------------------------------------------------

graph_nodes <- fread("data/transport_network/csv/graph_nodes.csv") 
graph_edges <- fread("data/transport_network/csv/graph_orig.csv") 

# Now: Plotting Productivity
graph_nodes %<>%
  mutate(citys = nif(outflows > 0, 4L, population >= 1e6, 3L, population >= 2e5, 2L, population >= 1, 1L, default = 0L), 
         prod_in = replace_na(37*outflows/population, 0),
         prod = IWI + prod_in, 
         total_prod = prod * population, 
         total_dom_prod = IWI * population) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Check: This is imports divided by African GDP: 26.5%
with(graph_nodes, sum(prod_in*population)/sum(IWI*population))
with(graph_nodes, sum(prod*population)/sum(IWI*population))

# Table
table(graph_nodes$citys)
graph_nodes |> qDT() |> 
  collap(prod_in ~ citys, list(fsum, fmean, fmedian), give.names = FALSE) |> 
  transpose(make.names = "citys", keep.names = "stat") |> 
  tfmv(is.numeric, scales::label_number(scale_cut = scales::cut_short_scale(), accuracy = 0.01)) |>
  xtable::xtable() |> print(include.r = FALSE, booktabs = TRUE)

# load("data/transport_network/trans_africa_network.RData")
edges_real <- qread("data/transport_network/edges_real_simplified.qs")

pdf("figures/transport_network/trans_africa_network_GE_parameterization_latest.pdf", width = 12, height = 12)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(mutate(edges_real, speed_kmh = (edges$distance/1000)/(edges$duration/60)) |> 
             rowbind(mutate(add_links, speed_kmh = 0) |> st_cast("MULTILINESTRING"), fill = TRUE)) +
  tm_lines(col = "speed_kmh", 
           col.legend = tm_legend("Speed (km/h)", position = tm_pos_in(0, 0.475), stack = "h", frame = FALSE, text.size = 1.3, title.size = 1.6),
           col.scale = tm_scale_continuous(values = "turbo"), # 7, 
           lwd = 2) + 
  tm_shape(subset(graph_nodes, citys == 4L) |>
             mutate(ofl = round(outflows / 1e6, 1))) + 
  tm_dots(size = "ofl", 
          size.scale = tm_scale_intervals(5, style = "jenks", values.scale = 2), # 
          size.legend = tm_legend("Port Outflows (M)", position = tm_pos_in(0, 0.475), frame = FALSE, text.size = 1.3, title.size = 1.6),
          fill = scales::alpha("black", 0.25)) +
  tm_shape(subset(graph_nodes, population > 0) |> mutate(pop = population / 1000)) + 
  tm_dots(size = "pop", 
          size.scale = tm_scale_intervals(breaks = c(0, 200, 1000, Inf),
                                          values = c(1, 3, 5)*0.12),
          size.legend = tm_legend("City Population (K)", position = tm_pos_in(0, 0.2), stack = "h", frame = FALSE, text.size = 1.3, title.size = 1.6),
          fill = "IWI",
          size.free = TRUE,
          fill.scale = tm_scale_intervals(values = "inferno"), #viridis::inferno(5, alpha = 0.5, direction = -1)),
          fill.legend = tm_legend("City Productivity (IWI)", position = tm_pos_in(0, 0.2), frame = FALSE, text.size = 1.3, title.size = 1.6)) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()

