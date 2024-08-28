#####################################################################
# Trans-African Network Between Cities > 100k and International Ports
#####################################################################

library(fastverse)
set_collapse(mask = c("manip", "helper", "special"), nthreads = 4)
fastverse_extend(qs, sf, s2, units, stplanr, sfnetworks, osrm, tmap, install = TRUE)
source("code/helpers.R")
fastverse_conflicts()

####################################
# Part 1: Important Cities and Ports
####################################

# Load Populated Cities from Simplemaps: https://simplemaps.com/data/world-cities
WC24_Africa <- fread("data/worldcities_continental_africa_2024.csv")
fsum(WC24_Africa$population)
WC24_Africa_g50k <- WC24_Africa |> subset(population > 5e4)

## Visualization
# WC24_Africa_g50k |> with(plot(lon, lat, cex = sqrt(population)/1e3, pch = 16, col = rgb(0, 0, 0, alpha=0.5)))
# WC24_Africa_g50k |> st_as_sf(coords = c("lon", "lat"), crs = st_crs(4326)) |> select(population) |> mapview::mapview()

# Weaving in administrative information through power weights
WC24_Africa_g50k[, pop_cap := 5^as.integer(factor(capital, levels = c("", "admin", "minor", "primary"))) * population]
# First Deduplication: extra weight for administrative function
WC24_Africa_g50k_red <- largest_within_radius(WC24_Africa_g50k, c("lon", "lat"), "pop_cap", radius_km = 30) 
WC24_Africa_g50k_red <- largest_within_radius(WC24_Africa_g50k_red, c("lon", "lat"), "pop_cap", radius_km = 50)  
## Test
# tmp = largest_within_radius(WC24_Africa_g50k_red, c("lon", "lat"), radius_km = 100, return_dist = TRUE)  
# diag(tmp) <- NA
# min(fmin(tmp))
# rm(tmp)

# Plotting
# WC24_Africa_g50k_red |> with(plot(lon, lat, cex = sqrt(population)/1e3, pch = 16, col = rgb(0, 0, 0, alpha=0.5)))
# WC24_Africa_g50k_red |> st_as_sf(coords = c("lon", "lat"), crs = st_crs(4326)) |> select(population, city) |> mapview::mapview()

# Now Adding populations within 30 km
WC24_Africa_g50k_red <- join_within_radius(WC24_Africa_g50k_red, WC24_Africa, c("lon", "lat"), radius_km = 30)

# Now a second deduplication using the 100k threshold
WC24_Africa_g100k_red <- largest_within_radius(WC24_Africa_g50k_red[population > 1e5], c("lon", "lat"), "population", radius_km = 70) 
# WC24_Africa_g100k_red |> with(plot(lon, lat, cex = sqrt(population)/1e3, pch = 16, col = rgb(0, 0, 0, alpha=0.5)))
# WC24_Africa_g100k_red |> st_as_sf(coords = c("lon", "lat"), crs = st_crs(4326)) |> select(population) |> mapview::mapview()

# Load International Ports from The World Bank: https://datacatalog.worldbank.org/search/dataset/0038118/Global---International-Ports -----

WBP_Africa <- fread("data/continental_african_ports.csv")
WBP_Africa_g100k <- WBP_Africa |> subset(outflows > 1e5)

## Visualization
# WBP_Africa_g100k |> with(plot(lon, lat, cex = sqrt(outflows)/1e3, pch = 16, col = rgb(0, 0, 0, alpha=0.5)))
# WBP_Africa_g100k |> st_as_sf(coords = c("lon", "lat"), crs = st_crs(4326)) |> select(outflows) |> mapview::mapview()

# Deduplication
WBP_Africa_g100k_red <- largest_within_radius(WBP_Africa_g100k, size = "outflows", radius_km = 100)
WBP_Africa_g100k_red <- join_within_radius(WBP_Africa_g100k_red, WBP_Africa, size = "outflows", radius_km = 99.9)

## Visualization
# WBP_Africa_g100k_red |> with(plot(lon, lat, cex = sqrt(outflows)/1e3, pch = 16, col = rgb(0, 0, 0, alpha=0.5)))
# WBP_Africa_g100k_red |> st_as_sf(coords = c("lon", "lat"), crs = st_crs(4326)) |> select(outflows) |> mapview::mapview()

# Now Joining Cities and Ports ---------------------------------

cities_ports_dmat <- st_distance(with(WC24_Africa_g100k_red, s2_lnglat(lon, lat)), with(WBP_Africa_g100k_red, s2_lnglat(lon, lat)))
fsum(cities_ports_dmat < as_units(30, "km"))

m_ind <- dapply(cities_ports_dmat < as_units(30, "km"), function(x) if (any(x)) which.max(x) else NA_integer_)

settransform(WBP_Africa_g100k_red, row = m_ind)
settransform(WC24_Africa_g100k_red, row = seq_along(city))

WC24_Africa_g100k_red <- WC24_Africa_g100k_red |> 
  join(subset(WBP_Africa_g100k_red, !is.na(row), row, port_locode = locode, port_name = name, port_status = status, outflows), on = "row") |> 
  select(-row)

# Matching remaining ports to cities
cities_ports_dmat <- st_distance(with(WC24_Africa_g50k_red, s2_lnglat(lon, lat)), with(WBP_Africa_g100k_red[is.na(row)], s2_lnglat(lon, lat)))
m_ind <- dapply(cities_ports_dmat < as_units(30, "km"), function(x) if (any(x)) which.max(x) else NA_integer_)
cities_ports_dmat <- st_distance(with(WC24_Africa_g50k, s2_lnglat(lon, lat)), with(WBP_Africa_g100k_red[is.na(row)], s2_lnglat(lon, lat)))
m_ind2 <- dapply(cities_ports_dmat < as_units(30, "km"), function(x) if (any(x)) which.max(x) else NA_integer_)

tmp <- Map(pfirst, ss(WC24_Africa_g50k_red, m_ind), ss(WC24_Africa_g50k, m_ind2)) |> qDT() |> 
  join_within_radius(WC24_Africa, c("lon", "lat"), radius_km = 39.9) |> 
  add_vars(subset(WBP_Africa_g100k_red, is.na(row), port_locode = locode, port_name = name, port_status = status, outflows))

tmp2 <- subset(WBP_Africa_g100k_red, is.na(row))

tmp[is.na(population), population := 1000]
setv(select(tmp, lon, lat, city), is.na(tmp$lat), select(tmp2, lon, lat, name))

# Final joining
cities_ports <- rowbind(WC24_Africa_g100k_red, tmp)

## Visualization
# WC24_Africa_g100k_red |> with(plot(lon, lat, cex = sqrt(population)/1e3, pch = 16, col = rgb(0, 0, 0, alpha=0.5)))
# WC24_Africa_g100k_red |> st_as_sf(coords = c("lon", "lat"), crs = st_crs(4326)) |> select(population) |> mapview::mapview()

# Cleanup
rm(list = setdiff(ls(), .c(cities_ports, split_large_dist_matrix)))


####################################
# Part 2: Transport Network
####################################

# Raw matrix of routes:
dist_ttime_mats <- split_large_dist_matrix(select(cities_ports, lon, lat))

# Checks
all.equal(dist_ttime_mats$sources, dist_ttime_mats$destinations)
allv(diag(dist_ttime_mats$durations), 0)
allv(diag(dist_ttime_mats$distances), 0)
diag(cor(dist_ttime_mats$sources, select(cities_ports, lon, lat)))
s2_distance(with(dist_ttime_mats$sources, s2_lnglat(lon, lat)),
            with(select(cities_ports, lon, lat), s2_lnglat(lon, lat))) |> descr()
pwcor(unattrib(dist_ttime_mats$distances), unattrib(dist_ttime_mats$durations))

# Now finding places that are on islands (e.g. Zanzibar City)
isl_ind <- which(fnobs(dist_ttime_mats$durations) < 100)

# -> Malabo and Zanzibar
cities_ports |> ss(isl_ind)

# Removing Island Cities
cities_ports <- ss(cities_ports, -isl_ind)
dist_ttime_mats$sources <- ss(dist_ttime_mats$sources, -isl_ind)
dist_ttime_mats$destinations <- ss(dist_ttime_mats$destinations, -isl_ind)
dist_ttime_mats$durations <- ss(dist_ttime_mats$durations, -isl_ind, -isl_ind)
dist_ttime_mats$distances <- ss(dist_ttime_mats$distances, -isl_ind, -isl_ind)

# Creating unique city identifier
fndistinct(cities_ports)
cities_ports[, city_country := stringi::stri_trans_general(paste(city, country, sep = " - "), "latin-ascii")]

# Saving
cities_ports |> fwrite("data/transport_network/cities_ports.csv")
dist_ttime_mats |> qsave("data/transport_network/cities_ports_dist_ttime_mats.qs")

# Reading again
cities_ports <- fread("data/transport_network/cities_ports.csv")
dist_ttime_mats <- qread("data/transport_network/cities_ports_dist_ttime_mats.qs")

# As spatial data frame
cities_ports_sf <- cities_ports |> st_as_sf(coords = c("lon", "lat"), crs = 4326)
# Distance between city centroids and route start/end points
descr(diag(st_distance(cities_ports_sf, st_as_sf(dist_ttime_mats$sources, coords = c("lon", "lat"), crs = 4326))))
# Create route-start-point version of the dataset
cities_ports_rsp_sf <- cities_ports_sf |> mutate(geometry = st_as_sf(dist_ttime_mats$sources, coords = c("lon", "lat"), crs = 4326)$geometry)

# Routes between all connections < 2000km apart ------------------------------------------

cities_ports_dmat <- dist_ttime_mats$sources |> with(s2_lnglat(lon, lat)) |> st_distance() |> set_units("km")
cities_ports_dmat[cities_ports_dmat > as_units(2000, "km")] <- NA
diag(cities_ports_dmat) <- NA
cities_ports_dmat[upper.tri(cities_ports_dmat)] <- NA

# Routes to be calculated
routes_ind <- which(!is.na(cities_ports_dmat), arr.ind = TRUE)
nrow(routes_ind)


# Determining Ideal Hypothetical (US-Grade) Network 

# See: https://christopherwolfram.com/projects/efficiency-of-road-networks/

# US Route efficeincy = 0.843, thus mr = 1/0.843 = 1.186 
keep_routes <- !intercepted_routes(routes_ind, dist_ttime_mats$sources, NULL, alpha = 20, mr = 1/0.843) 
sum(keep_routes)

# mapview::mapview(subset(routes, keep_routes), zcol = "distance") +
#   mapview::mapview(st_as_sf(cities_ports, coords = c("lon", "lat"), crs = st_crs(4326)))

# Plot ideal Network
with(cities_ports, {
  oldpar <- par(mar = c(0,0,0,0))
  plot(lon, lat, cex = sqrt(population)/1e3, pch = 16, col = rgb(0, 0, 0, alpha=0.5), 
       axes = FALSE, xlab = NA, ylab = NA, asp=1)
  for (r in mrtl(routes_ind[keep_routes, ])) { # & intersects == 2
    lines(lon[r], lat[r])
  }
  par(oldpar)
})

dev.copy(pdf, "figures/transport_network/trans_africa_network_US_20deg.pdf", width = 10, height = 10)
dev.off()

# EU Route Efficiency  
keep_routes <- !intercepted_routes(routes_ind, dist_ttime_mats$sources, NULL, alpha = 45, mr = 1/0.767) 
sum(keep_routes)

# Plot ideal Network
with(cities_ports, {
  oldpar <- par(mar = c(0,0,0,0))
  plot(lon, lat, cex = 1.3, pch = 16, axes = FALSE, xlab = NA, ylab = NA, asp=1)
  for (r in mrtl(routes_ind[keep_routes, ])) { # & intersects == 2
    lines(lon[r], lat[r])
  }
  par(oldpar)
})

dev.copy(pdf, "figures/transport_network/trans_africa_network_EU_45deg.pdf", width = 10, height = 10)
dev.off()


# Fetching all Routes and Simplifying -------------------------------------------
# Note: The routes file is about 2.6GB large! you can load it yourself by executing this code. 
# In this replication package I only provide final segments: the result of intersecting all routes.

routes <- data.table(from_city_country = cities_ports$city_country[routes_ind[, 1]], 
                     to_city_country = cities_ports$city_country[routes_ind[, 2]], 
                     duration = NA_real_, 
                     distance = NA_real_, 
                     geometry = list())
# Fetch Routes
i = 1L
for (r in mrtl(routes_ind)) {
  cat(i, " ")
  route <- osrmRoute(ss(cities_ports, r[1], c("lon", "lat")),
                     ss(cities_ports, r[2], c("lon", "lat")), overview = "full")
  set(routes, i, 3:5, select(route, duration, distance, geometry))
  i <- i + 1L
}
routes <- routes |> st_as_sf(crs = st_crs(route))

# Saving
routes |> qsave("data/transport_network/routes_raw.qs")

# Adding Gravity to Routes (https://en.wikipedia.org/wiki/Newton%27s_law_of_universal_gravitation)
dmat <- st_distance(cities_ports_sf)
diag(dmat) <- NA
frange(dmat) 

colnames(dmat) <- rownames(dmat) <- cities_ports$city_country
ddf <- pivot(qDF(dmat, "from"), "from", names = list("to", "sp_dist"), na.rm = TRUE) |> 
  join(select(cities_ports, from = city_country, from_pop = population)) |> 
  join(select(cities_ports, to = city_country, to_pop = population)) |> 
  mutate(sp_dist = sp_dist / 1000,
         gravity = from_pop * to_pop / sp_dist / 1e9) # Figure in TEU for port cities??

routes <- routes |> 
  join(ddf, on = c(from_city_country = "from", to_city_country = "to")) |>
  mutate(gravity_rd = from_pop * to_pop / distance / 1e9, 
         gravity_dur = from_pop * to_pop / duration / 1e9)

# Intersecting Routes
segments <- list()
n <- fnrow(routes)
for (i in seq(1, n, by = 1000)) {
  rows_i = i:min(i + 1000 - 1, n)
  segments[[as.character(i)]] <- overline2(routes[rows_i, ] |> mutate(passes = 1L), 
                                           attrib = c("passes", "gravity", "gravity_rd", "gravity_dur"))
}
rm(routes); gc()
segments <- rowbind(segments) |> st_make_valid()
segments <- overline2(segments, attrib = c("passes", "gravity", "gravity_rd", "gravity_dur"))
segments %<>% ss(!is_linepoint(.)) %>% st_make_valid()

# Saving
segments |> qsave("data/transport_network/segments.qs")


# Loading Segments --------------------------------------------------------------------

segments <- qread("data/transport_network/segments.qs")

# First Round of subdivision
segments <- rmapshaper::ms_simplify(segments, keep = 0.5, snap_interval = deg_m(500)) |> 
  subset(vlengths(geometry) > 0)
segments <- overline2(segments, attrib = c("passes", "gravity", "gravity_rd", "gravity_dur"))
if(any(is_linepoint(segments))) segments %<>% ss(!is_linepoint(.))
segments_simp <- segments |> rmapshaper::ms_simplify(keep = 0.25, snap_interval = deg_m(1500)) # For later plot
# plot(segments)
# mapview(segments) + mapview(rnet_get_nodes(segments))

# Note: The code bwloe can be executed to regenerate the network. For all practical purposes however I have 
# already provided the final network, so this code can also be skipped. 

# Creating Network
net <- as_sfnetwork(segments, directed = FALSE)
# plot(net)
summ_fun <- function(fun) list(passes = fun, gravity = fun, gravity_rd = fun, gravity_dur = fun, "ignore")
filter_smooth <- function(net) {
  net |> 
    tidygraph::activate("edges") |> 
    dplyr::filter(!tidygraph::edge_is_multiple()) |> 
    dplyr::filter(!tidygraph::edge_is_loop()) |> 
    tidygraph::convert(to_spatial_smooth, 
                       protect = cities_ports_rsp_sf$geometry,
                       summarise_attributes = summ_fun("mean"))
}


# Filtering and Smoothing
net <- filter_smooth(net)

# Second round of subdivision
segments <- net |> activate("edges") |> tidygraph::as_tibble() |> mutate(.tidygraph_edge_index = NULL)
segments <- rmapshaper::ms_simplify(segments, keep = 0.5, snap_interval = deg_m(1000)) |> 
  subset(vlengths(geometry) > 0)
segments <- overline2(segments, attrib = c("passes", "gravity", "gravity_rd", "gravity_dur"))
if(any(is_linepoint(segments))) segments %<>% ss(!is_linepoint(.))

# Create SF network
net <- as_sfnetwork(segments, directed = FALSE)

## Filtering, smoothing and spatial subdivision 
net <- filter_smooth(net)
net <- tidygraph::convert(net, to_spatial_subdivision)
net <- filter_smooth(net)
# plot(net)
# mapview(st_geometry(net, "edges")) + mapview(st_geometry(net, "nodes"))

# Saving Smoothed Version (pre contraction)
net |> qsave("data/transport_network/net_smoothed.qs")

## Contracting network: Manually 
segments <- net |> activate("edges") |> tidygraph::as_tibble() |> mutate(.tidygraph_edge_index = NULL)
nodes <- nodes_max_passes(segments, attrib = c("passes", "gravity", "gravity_rd", "gravity_dur"))
nodes_clustered <- cluster_nodes_by_cities(nodes, cities_ports_rsp_sf, 
                                           city_radius_km = 25, 
                                           cluster_radius_km = 15, 
                                           algo.dbscan = FALSE,
                                           weight = "gravity_rd")
segments_contracted <- contract_segments(segments, nodes_clustered, attrib = c("passes", "gravity", "gravity_rd", "gravity_dur"))
# plot(segments_contracted)
# mapview(segments_contracted) + mapview(rnet_get_nodes(segments_contracted))

net <- as_sfnetwork(segments_contracted, directed = FALSE)

# Filtering, smoothing and spatial subdivision 
net <- filter_smooth(net)
net <- tidygraph::convert(net, to_spatial_subdivision)
net <- filter_smooth(net)
# plot(net)
# mapview(st_geometry(net, "edges")) + mapview(st_geometry(net, "nodes"))

## Second contraction round
segments <- net |> activate("edges") |> tidygraph::as_tibble() |> mutate(.tidygraph_edge_index = NULL)
nodes <- nodes_max_passes(segments, attrib = c("passes", "gravity", "gravity_rd", "gravity_dur"))
nodes_clustered <- cluster_nodes_by_cities(nodes, cities_ports_rsp_sf, 
                                           city_radius_km = 30, 
                                           cluster_radius_km = 30, 
                                           algo.dbscan = FALSE,
                                           weight = "gravity_rd")
segments_contracted <- contract_segments(segments, nodes_clustered, attrib = c("passes", "gravity", "gravity_rd", "gravity_dur"))
# plot(segments_contracted)
# mapview(segments_contracted) + mapview(rnet_get_nodes(segments_contracted))

net <- as_sfnetwork(segments_contracted, directed = FALSE)
# plot(net)
# mapview(st_geometry(net, "edges")) + mapview(st_geometry(net, "nodes"))

# Filtering, smoothing and spatial subdivision 
net <- filter_smooth(net)
net <- tidygraph::convert(net, to_spatial_subdivision)
net <- filter_smooth(net)
# plot(net)
# mapview(st_geometry(net, "edges")) + mapview(st_geometry(net, "nodes")) # mapview(st_geometry(net_smoothed, "edges"))

## Saving 
net |> qsave("data/transport_network/net_discrete_final.qs")

## Loading Final Network -----------------------------------------

net <- qread("data/transport_network/net_discrete_final.qs")

## Plotting
plot(net)
nodes <- net |> activate("nodes") |> tidygraph::as_tibble() |> 
  mutate(.tidygraph_node_index = NULL)
edges <- net |> activate("edges") |> tidygraph::as_tibble() |> 
  mutate(.tidygraph_edge_index = NULL, log10_gravity = log10(gravity))

nodes$city_port <- rowSums(st_distance(nodes, cities_ports_rsp_sf) < as_units(20, "km")) > 0
sum(nodes$city_port)
descr(edges$gravity_rd)

# Needed throughout 
nodes_coord <- st_coordinates(nodes) |> qDF() |> setNames(c("lon", "lat"))

tmap_mode("plot")

# First a plot of just the routes
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(segments_simp) + tm_lines(col = "black") +
  tm_shape(subset(nodes, city_port)) + tm_dots(size = 0.2, fill = "orange2") +
  tm_layout(frame = FALSE)

dev.copy(pdf, "figures/transport_network/trans_africa_network_actual.pdf", width = 10, height = 10)
dev.off()

# Now the plot with the discretized representation
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(segments_simp) + tm_lines(col = "black") +
  tm_shape(edges) +
  tm_lines(col = "gravity_rd",
           col.scale = tm_scale_intervals(values = "inferno", breaks = c(0, 5, 25, 125, 625, Inf)),
           col.legend = tm_legend("Sum of Gravity", position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_shape(subset(nodes, city_port)) + tm_dots(size = 0.2) +
  tm_shape(subset(nodes, !city_port)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)

dev.copy(pdf, "figures/transport_network/trans_africa_network_actual_discretized_gravity_plus_orig.pdf", width = 10, height = 10)
dev.off()


# Recalculating Routes Matrix for Network ----------------------------------------------------------

nrow(nodes_coord)
all.equal(unattrib(nodes_coord), mctl(st_coordinates(nodes)))
dist_ttime_mats <- split_large_dist_matrix(nodes_coord, verbose = TRUE)

# Checks
all.equal(dist_ttime_mats$sources, dist_ttime_mats$destinations)
allv(diag(dist_ttime_mats$durations), 0)
allv(diag(dist_ttime_mats$distances), 0)
diag(cor(dist_ttime_mats$sources, nodes_coord))
s2_distance(with(dist_ttime_mats$sources, s2_lnglat(lon, lat)),
            with(nodes_coord, s2_lnglat(lon, lat))) |> descr()
pwcor(unattrib(dist_ttime_mats$distances), unattrib(dist_ttime_mats$durations))
# Now finding places that are on islands (e.g. Zanzibar City): should not exist here
if(any(which(fnobs(dist_ttime_mats$durations) < 200))) stop("Found Islands")

dist_ttime_mats |> qsave("data/transport_network/dist_ttime_mats_network_30km.qs")

# Loading again -----------------------------------
dist_ttime_mats <- qread("data/transport_network/dist_ttime_mats_network_30km.qs")

# Check
all.equal(st_as_sf(net, "nodes")$geometry, nodes$geometry)

# Making symmetric
sym_dist_mat <- (dist_ttime_mats$distances + t(dist_ttime_mats$distances)) / 2
sym_time_mat <- (dist_ttime_mats$durations + t(dist_ttime_mats$durations)) / 2

# Add average distance and travel time to edges
edges_ind <- edges |> qDF() |> select(from, to) |> qM()
edges$sp_distance <- st_length(edges)
edges$distance <- sym_dist_mat[edges_ind]
edges$duration <- sym_time_mat[edges_ind]


# Adding High-Value Links to Existing Network ----------------------------------------------------------
# Note: This is a self-contained section due to the use of specific indices for manual adjustments!! -> Skip it!
# It uses a previous version of the network, which is basically identical to the current version, but the nodes and edges are in a different order. 
# Because only the output of this section is relevant, I have not added the inputs to the replication package. 
# But users can still execute the section with the current network and distance matrices to verify that the network finding algorithm produces 
# sensible new links. Just the manual adjustment (removing links crossing water bodies and protected areas and adding a few links) won't work
# using the indices provided. The output is a spatial data frame 'add_routes' with the final proposed links. These can be added to the current network.

net <- qread("data/transport_network/old/net_discrete_final.qs")
dist_ttime_mats <- qread("~/Documents/IFW Kiel/Africa-Infrastructure/data/transport_network/old/dist_ttime_mats_network_30km.qs")
sym_dist_mat <- (dist_ttime_mats$distances + t(dist_ttime_mats$distances)) / 2

## Plotting
plot(net)
nodes <- net |> activate("nodes") |> tidygraph::as_tibble() |> 
  mutate(.tidygraph_node_index = NULL)

nodes_dmat <- st_distance(nodes) |> set_units("km")
nodes_dmat[nodes_dmat > as_units(2000, "km")] <- NA
diag(nodes_dmat) <- NA
nodes_dmat[upper.tri(nodes_dmat)] <- NA

# Routes to be calculated
routes_ind <- which(!is.na(nodes_dmat), arr.ind = TRUE)
nrow(routes_ind)

nodes_coord <- st_coordinates(nodes) |> qDF() |> setNames(c("lon", "lat"))

keep_routes <- !intercepted_routes(routes_ind, nodes_coord, sym_dist_mat, feasible = TRUE, alpha = 45, mr = 1/0.767, fmr = 1.5) # US: 0.843 EU: 0.767
sum(keep_routes)

add_routes <- with(nodes_coord, lapply(mrtl(routes_ind[keep_routes, ]), function(r) st_linestring(cbind(lon[r], lat[r])))) |> st_sfc(crs = 4326) |> st_as_sf()
add_routes_df <- line2df(add_routes)[-1] %*=% 1000 %>% dapply(as.integer)
edges_df <- line2df(edges)[-1] %*=% 1000 %>% dapply(as.integer)
m <- add_routes_df %in% edges_df | add_routes_df %in% gv(edges_df, c(3:4, 1:2))
nrow(add_routes) - sum(m)
add_routes <- add_routes[!m, ]
rm(add_routes_df, edges_df)

add_routes |> qsave("~/Documents/IFW Kiel/Africa-Infrastructure/data/transport_network/add_routes_network_30km_alpha45_mrEU_fmr15.qs")
add_routes <- qread("~/Documents/IFW Kiel/Africa-Infrastructure/data/transport_network/add_routes_network_30km_alpha45_mrEU_fmr15.qs")

# Plot
# mapview(edges, map.types = c(mapviewGetOption("basemaps"), "Esri.WorldStreetMap", "Esri.WorldTerrain")) + mapview(nodes) + mapview(add_routes, color = "green") # [remove, ]
remove <- c(305, 396, 399, 400, 402, 512, 513, 514, 516, 403, 409, 410, 390, 386, 457, 420, 429, 428, 426, 427, 434,
            358, 357, 363, 364, 344, 345, 165, 175, 106, 27, 10, 69, 68, 433, 490, 499, 444, 326, 324, 419, 279, 317, 178, 331, 270)
# link 59 and 57 should have cost of 0 -> political
# Links 32, 33, 39, 62, 67, 118, 128, 129, 141, 412, 487, 489 could have better roads
# add_routes[remove, "x"] <- st_linestring()
add <- list(c(63, 105), c(862, 717), c(499, 519), c(499, 455), c(312, 379), c(1184, 1271), c(46, 37), c(1071, 842),
            c(901, 751), c(751, 710), c(593, 514))
add_routes <- rbind(add_routes[-remove, ],
                    with(nodes_coord, lapply(add, function(r) st_linestring(cbind(lon[r], lat[r])))) |> st_sfc(crs = 4326) |> st_as_sf())

add_routes |> qsave("~/Documents/IFW Kiel/Africa-Infrastructure/data/transport_network/add_routes_network_30km_alpha45_mrEU_fmr15_adjusted.qs")


# Loading, Integrating and Plotting ---------------------------------------------------

add_routes <- qread("data/transport_network/add_routes_network_30km_alpha45_mrEU_fmr15_adjusted.qs")
add_routes_df <- line2points(add_routes)
dmat <- st_distance(nodes$geometry, add_routes_df$geometry)
add_routes_df$node <- dapply(dmat, which.min)
add_routes_df$geometry <- nodes$geometry[add_routes_df$node]
add_routes <- add_routes_df |> group_by(id) |> 
  summarise(from = ffirst(node),
            to = flast(node),
            geometry = st_combine(geometry)) |> st_cast("LINESTRING")
# Checks
all(line2df(add_routes) %>% select(fx, fy) %in% qDF(st_coordinates(nodes)))
all(line2df(add_routes) %>% select(tx, ty) %in% qDF(st_coordinates(nodes)))
all(select(line2df(add_routes), fx:ty) %!in% select(line2df(edges), fx:ty))
rm(dmat, add_routes_df)

tmap_mode("plot")

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(col = "gravity_rd", 
           col.scale = tm_scale_intervals(values = "inferno", breaks = c(0, 5, 25, 125, 625, Inf)),
           col.legend = tm_legend("Sum of Gravity", position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_shape(add_routes) + tm_lines(col = "green4", lwd = 1) + 
  tm_shape(subset(nodes, city_port)) + tm_dots(size = 0.2) +
  tm_shape(subset(nodes, !city_port)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 

dev.copy(pdf, "figures/transport_network/trans_africa_network_actual_discretized_gravity_new_roads.pdf", 
         width = 10, height = 10)
dev.off()

# Now Adding Populations ----------------------------------------------------------

# min(nodes_dmat, na.rm = TRUE) # Need to assign to nearest node

# Distance Matrix
dmat <- nodes |> st_distance(st_as_sf(WC24_Africa, coords = c("lon", "lat"), crs = 4326)) 
# Remove cities part of port cities
used <- dapply(dmat[nodes$city_port, ] < as_units(30, "km"), any)
table(used)
dmat <- dmat[!nodes$city_port, !used]
dim(dmat)
dmat[dmat >= as_units(30, "km")] <- NA
col <- numeric(nrow(dmat))
m <- dapply(dmat, function(x) if (allNA(x)) col else replace(col, which.min(x), 1)) 
nodes$population <- NA
nodes$population[!nodes$city_port] <- rowSums(m %r*% WC24_Africa$population[!used])
WC24_Africa[, pop_cap := 5^as.integer(factor(capital, levels = c("", "admin", "minor", "primary"))) * population]
nodes$city_country <- NA_character_
nodes$city_country[!nodes$city_port] <- WC24_Africa |> 
  with(paste(city, country, sep = " - ")) |> 
  extract(dapply(m %r*% WC24_Africa$pop_cap[!used], 
                 function(x) if(any(x > 0.1)) which.max(x) else NA_integer_, MARGIN = 1))
# Now adding port cities (closest matches)
ind <- dapply(st_distance(cities_ports_rsp_sf, nodes), function(x) if (any(x < 20e3)) which.min(x) else NA_integer_)
nodes$population[!is.na(ind)] <- cities_ports_rsp_sf$population[na_rm(ind)]
nodes$city_country[!is.na(ind)] <- cities_ports_rsp_sf$city_country[na_rm(ind)]
# Cleanup
rm(col, m, ind, used, dmat)

# Ratios
sum(nodes$population) / sum(cities_ports_rsp_sf$population)
sum(nodes$population) / sum(WC24_Africa$population)
(sum(nodes$population > 0) - nrow(cities_ports_rsp_sf)) / (nrow(nodes) - nrow(cities_ports_rsp_sf))

# Saving all necessary objects in an RData file ---------------------------------------------------

# First adding info to edges
tmp <- net |> st_as_sf("edges") |> atomic_elem() |> qDT() |> 
  join(atomic_elem(edges), overid = 2L) |> 
  select(-from, -to, -.tidygraph_edge_index)

net %<>% activate("edges") %>% dplyr::mutate(tmp)
rm(tmp)

save(nodes, edges, edges_ind, nodes_coord, net, add_routes, WC24_Africa, 
     cities_ports, cities_ports_sf, cities_ports_rsp_sf, 
     dist_ttime_mats, sym_dist_mat, sym_time_mat, 
     file = "data/transport_network/trans_africa_network.RData")




########################################################
# Part 2.1: Recalculating Routes for Visual Analysis
########################################################

load("data/transport_network/trans_africa_network.RData")

# This is the previous Plot
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) +
  tm_lines(col = "gravity_rd",
           col.scale = tm_scale_intervals(values = "inferno", breaks = c(0, 5, 25, 125, 625, Inf)),
           col.legend = tm_legend("Sum of Gravity", position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_shape(add_routes) + tm_lines(col = "green4", lwd = 1) + 
  tm_shape(subset(nodes, city_port)) + tm_dots(size = 0.2) +
  tm_shape(subset(nodes, !city_port)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 

# dev.copy(pdf, "figures/transport_network/trans_africa_network_actual_discretized_gravity_new_roads.pdf", 
#          width = 10, height = 10)
# dev.off()

# Now: Load Simplified Routes
edges_real <- edges |> qDT() |> 
  transform(geometry = list(NULL), distance = NA_real_, duration = NA_real_)
nodes_coord <- st_coordinates(nodes) |> qDF() |> setNames(c("lon", "lat"))

# Fetch Routes
i <- 1L
for (r in mrtl(edges_ind)) {
  cat(i, " ")
  route <- osrmRoute(ss(nodes_coord, r[1]),
                     ss(nodes_coord, r[2]), overview = "simplified")
  set(edges_real, i, c("duration", "distance", "geometry"), 
      select(route, duration, distance, geometry))
  i <- i + 1L
}
edges_real <- edges_real |> st_as_sf(crs = st_crs(route)) |> st_make_valid()
edges_real |> qsave("data/transport_network/edges_real_simplified.qs")
rm(route, i, r)

# Now the updated plot
edges_real <- qread("data/transport_network/edges_real_simplified.qs")

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges_real) +
  tm_lines(col = "gravity_rd",
           col.scale = tm_scale_intervals(values = "inferno", breaks = c(0, 5, 25, 125, 625, Inf)),
           col.legend = tm_legend("Sum of Gravity", position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_shape(add_routes) + tm_lines(col = "green4", lwd = 1) + # limegreen # , lty = "twodash"
  tm_shape(subset(nodes, city_port)) + tm_dots(size = 0.15) +
  tm_shape(subset(nodes, !city_port & population > 0)) + tm_dots(size = 0.1, fill = "grey20") +
  tm_shape(subset(nodes, !city_port & population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) #, inner.margins = c(0.1, 0.1, 0.1, 0.1))

dev.copy(pdf, "figures/transport_network/trans_africa_network_actual_discretized_gravity_new_roads_real_edges.pdf", 
         width = 10, height = 10)
dev.off()


########################################################
# Part 2.2: Identifying Major Transport Routes
########################################################

# TODO: Keep here? Or separate file??

load("data/transport_network/trans_africa_network.RData")
edges_real <- qread("data/transport_network/edges_real_simplified.qs")
nodes %<>% join(atomic_elem(cities_ports_rsp_sf), on = c("population", "city_country"), drop = "y")


# Plot high gravity roads
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(subset(edges_real, gravity_rd >= 25)) +
  tm_lines(col = "gravity_rd",
           col.scale = tm_scale_intervals(values = "inferno", breaks = c(0, 5, 25, 125, 625, Inf)),
           col.legend = tm_legend("Sum of Gravity", position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_shape(subset(edges_real, gravity_rd < 25)) +
  tm_lines(col = "grey50", lwd = 2) +
  tm_shape(subset(nodes, city_port)) + tm_dots(size = 0.15) +
  tm_shape(subset(nodes, !city_port & population > 0)) + tm_dots(size = 0.1, fill = "grey20") +
  tm_shape(subset(nodes, !city_port & population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 

# 47 largest port-cities
large_cities <- nodes %$% which(population > 2e6 | outflows > 1e6)
length(large_cities)

# Fastest Routes between them
# igraph::all_shortest_paths(net, large_cities[1], large_cities)
large_city_paths <- lapply(large_cities, function(i) 
  st_network_paths(net, from = i, to = large_cities, weights = "duration", mode = "all")) |>
  rowbind()

large_city_paths <- list(nodes = unique(unlist(large_city_paths$node_paths)), 
                         edges = unique(unlist(large_city_paths$edge_paths)))

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(subset(edges, large_city_paths$edges)) +
  tm_lines(col = "gravity_rd",
           col.scale = tm_scale_intervals(values = "inferno", breaks = c(0, 5, 25, 125, 625, Inf)),
           col.legend = tm_legend("Sum of Gravity", position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_shape(subset(edges, -large_city_paths$edges)) +
  tm_lines(col = "grey70", lwd = 2) +
  tm_shape(subset(nodes, large_city_paths$nodes)) + tm_dots(size = 0.2) +
  tm_shape(subset(nodes, -large_city_paths$nodes)) + tm_dots(size = 0.1, fill = "grey50") +
  tm_shape(subset(nodes, large_cities)) + tm_dots(size = 0.2, fill = "red") +
  tm_layout(frame = FALSE) 

net_main_routes <- as_sfnetwork(subset(edges, large_city_paths$edges), 
                                directed = FALSE, edges_as_lines = TRUE)

filter_smooth2 <- function(net, funs = list(passes = "mean", gravity = "mean", gravity_rd = "mean", gravity_dur = "mean", 
                                            distance = "sum", duration = "sum", "ignore")) {
  net |> 
    tidygraph::activate("edges") |> 
    dplyr::filter(!tidygraph::edge_is_multiple()) |> 
    dplyr::filter(!tidygraph::edge_is_loop()) |> 
    tidygraph::convert(to_spatial_smooth, 
                       protect = subset(nodes, large_cities, geometry)$geometry,
                       summarise_attributes = funs)
}

net_main_routes <- filter_smooth2(net_main_routes)
net_main_routes <- filter_smooth2(net_main_routes)

plot(net_main_routes)
st_as_sf(net_main_routes, "edges") |> names()

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(subset(edges, -large_city_paths$edges)) + tm_lines(col = "grey70", lwd = 2) +
  tm_shape(nodes) + tm_dots(size = 0.1, fill = "grey50") +
  tm_shape(st_as_sf(activate(net_main_routes, "edges"))) +
  tm_lines(col = "gravity_rd",
           col.scale = tm_scale_intervals(values = "inferno", breaks = c(0, 5, 25, 125, 625, Inf)),
           col.legend = tm_legend("Sum of Gravity", position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_shape(st_as_sf(activate(net_main_routes, "nodes"))) + tm_dots(size = 0.2) +
  tm_shape(subset(nodes, large_cities)) + tm_dots(size = 0.3, fill = "red") +
  tm_layout(frame = FALSE) #, inner.margins = c(0.1, 0.1, 0.1, 0.1))

dev.copy(pdf, "figures/transport_network/trans_africa_network_reduced_47_duration.pdf", 
         width = 10, height = 10)
dev.off()

# Join nodes
nodes_sf <- st_as_sf(net_main_routes, "nodes")
nodes_sf %<>% transform(qDF(st_coordinates(.))) %>% 
  join(transform(nodes, qDF(st_coordinates(geometry))), on  = c("X", "Y"), drop = TRUE) %>%
  select(-X, -Y, -.tidygraph_node_index)
net_main_routes %<>% activate("nodes") %>% dplyr::mutate(qDF(atomic_elem(nodes_sf)))

results <- list(
  fastest_routes = c(large_city_paths, list(network = net_main_routes))
)

# Now With All Routes, incl. new ones ---------------------------------------

identical(st_geometry(net, "edges"), edges$geometry)

net_ext_data <- rowbind(select(edges_param, from, to, add, gravity_rd, distance, duration, duration_imp, 
                               border_dist, border_time, total_time, cost_per_km = ug_cost_km, geometry) |>
                          mutate(total_cost = distance / 1000 * cost_per_km) , 
                        select(add_routes_param, from, to, add, distance, duration = duration_100kmh,  
                               border_dist, border_time, total_time = total_time_100kmh, cost_per_km = cost_km_adj, geometry) |>
                          transform(total_cost = distance / 1000 * cost_per_km, 
                                    gravity_rd = NA_real_, duration_imp = NA_real_))
net_ext <- sfnetwork(nodes_param, net_ext_data, directed = FALSE)
# Check
identical(st_geometry(net_ext, "nodes"), st_geometry(net_param, "nodes"))

large_city_paths <- lapply(large_cities, function(i) 
  rowbind(st_network_paths(net_param, from = i, to = large_cities, weights = "distance", mode = "all"),
          st_network_paths(net_ext, from = i, to = large_cities, weights = "distance", mode = "all"),
          st_network_paths(net_param, from = i, to = large_cities, weights = "duration", mode = "all"))
) |>
  rowbind()

large_city_paths <- list(nodes = unique(unlist(large_city_paths$node_paths)), 
                         edges = unique(unlist(large_city_paths$edge_paths)))

tmp = subset(st_as_sf(net_ext, "edges"), large_city_paths$edges)
settfm(tmp, 
       duration_imp = iif(add, duration, duration_imp),
       duration = iif(add, (distance/1000)*60, duration)) # Setting to 1 km/h
net_main_routes <- as_sfnetwork(tmp, directed = FALSE, edges_as_lines = TRUE)

net_main_routes <- filter_smooth2(net_main_routes, 
                                  funs = list(add = "mean", gravity_rd = "mean", 
                                              distance = "sum", duration = "sum", duration_imp = "sum", 
                                              border_dist = "sum", border_time = "sum", total_time = "sum",
                                              total_cost = "sum", "ignore"))
# net_main_routes <- tidygraph::convert(net_main_routes, to_spatial_subdivision)
# net_main_routes <- filter_smooth2(net_main_routes)

plot(net_main_routes)
st_as_sf(net_main_routes, "edges") |> fnobs()

tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) + tm_lines(col = "grey70", lwd = 2) +
  tm_shape(nodes) + tm_dots(size = 0.1, fill = "grey50") +
  tm_shape(st_as_sf(net_main_routes, "edges")) +
  tm_lines(col = "gravity_rd",
           col.scale = tm_scale_intervals(values = "inferno", breaks = c(0, 5, 25, 125, 625, Inf), 
                                          value.na = "green4", label.na = "Proposed or Mixed"),
           col.legend = tm_legend("Sum of Gravity"), lwd = 2) +
  tm_shape(st_as_sf(activate(net_main_routes, "nodes"))) + tm_dots(size = 0.2) +
  tm_shape(subset(nodes, large_cities)) + tm_dots(size = 0.3, fill = "red") +
  tm_legend(position = c("left", "bottom"), frame = FALSE, 
            text.size = 1.5, title.size = 2) +
  tm_layout(frame = FALSE) #, inner.margins = c(0.1, 0.1, 0.1, 0.1))

dev.copy(pdf, "figures/transport_network/trans_africa_network_reduced_47_all.pdf", 
         width = 10, height = 10)
dev.off()

# Join nodes
nodes_sf <- st_as_sf(net_main_routes, "nodes")
nodes_sf %<>% transform(qDF(st_coordinates(.))) %>% 
  join(transform(nodes_param, qDF(st_coordinates(geometry))), on  = c("X", "Y"), drop = TRUE) %>%
  select(-X, -Y, -.tidygraph_node_index)
net_main_routes %<>% activate("nodes") %>% dplyr::mutate(qDF(atomic_elem(nodes_sf)))

# Combining results
results$all_routes <- c(large_city_paths, list(network = net_main_routes))

qsave(results, "data/transport_network/largest_pcities/trans_africa_network_47_largest.qs")

# Now saving as CSV ---------------------

results <- qread("data/transport_network/largest_pcities/trans_africa_network_47_largest.qs")

for (nname in .c(fastest_routes, all_routes)) {
  
  results[[nname]]$network %>% st_as_sf("nodes") %>% atomic_elem() %>% qDT() %>% 
    select(-.tidygraph_node_index) %>% 
    fwrite(sprintf("data/transport_network/largest_pcities/%s_graph_nodes.csv", nname))
  
  results[[nname]]$network %>% st_as_sf("edges") %>% atomic_elem() %>% qDT() %>% 
    fwrite(sprintf("data/transport_network/largest_pcities/%s_graph_edges.csv", nname))
}

# Plotting --------------------------------

nname <- "fastest_routes"
# edges_real <- qread("~/Documents/IFW Kiel/Africa-Infrastructure/data/transport_network/edges_real_simplified.qs") 
nodes <- results[[nname]]$network |> st_as_sf("nodes") 
edges <- results[[nname]]$network |> st_as_sf("edges") 

# Classification
settfm(edges, speed_kmh = (distance / 1000) / (duration / 60))
# settfm(edges_real, speed_kmh = (edges_param$distance/1000)/(edges_param$duration/60))
# edges_real %<>% subset(results$fastest_routes$edges)

# nodes <- fread(sprintf("data/transport_network/largest_pcities/%s_graph_nodes.csv", nname))
settfm(nodes, major_city_port = population > 2e6 | outflows > 1e6)
sum(nodes$major_city_port)
largest <- c("Dakar - Senegal", "Casablanca - Morocco", "Abidjan - Côte d’Ivoire", 
             "Kumasi - Ghana", "Algiers - Algeria", "Lagos - Nigeria", "Kano - Nigeria", "Yaoundé - Cameroon", 
             "Luanda - Angola", "Kinshasa - Congo (Kinshasa)", "Johannesburg - South Africa", "Cape Town - South Africa",
             "Cairo - Egypt", "Khartoum - Sudan", "Nairobi - Kenya", "Addis Ababa - Ethiopia", 
             "Dar es Salaam - Tanzania")

settfm(nodes, product = nif(major_city_port & base::match(city_country, largest, 0L) > 0L, NA_integer_, # Heterogeneous products
                            population > 1e6 & outflows > 1e6, 5L, # Large Port-City
                            population > 2e6, 4L,   # Large City
                            outflows > 0, 3L,       # Port
                            population > 2e5, 2L,   # Medium-Sized City
                            default = 1L))          # Town/Node
table(nodes$product, na.exclude = FALSE)
setv(nodes$product, whichNA(nodes$product), seq_along(largest) + 5L)
# fwrite(nodes, sprintf("data/transport_network/largest_pcities/%s_graph_nodes.csv", nname))
attr(nodes$product, "levels") <- c("Small City/Node", "City > 200K", "Port", "City > 2M", "Large Port-City", paste("Megacity", seq_along(largest)))
class(nodes$product) <- "factor"

# Plotting
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges) + # _real
  tm_lines(col = "speed_kmh", 
           col.legend = tm_legend("Speed", position = c("left", "bottom"), stack = "h", 
                                  frame = FALSE, text.size = 1.5, title.size = 2),
           col.scale = tm_scale_continuous(6, values = "turbo"), # 7, 
           lwd = 2) + 
  # tm_shape(subset(nodes, !major_city_port)) + tm_dots(size = 0.1, fill = "black") +
  # tm_shape(subset(nodes, major_city_port)) + tm_dots(size = 0.3, fill = "red") +
  tm_shape(droplevels(mutate(nodes, product = fifelse(unclass(product) > 5L, NA, product)))) + 
  tm_dots(fill = "product", size = 0.25, 
          fill.scale = tm_scale_categorical(values = "turbo", value.na = "purple3", label.na = "Megacity (Own Product)"),
          fill.legend = tm_legend("Product", position = c("left", "bottom"), # stack = "h", 
                                  frame = FALSE, text.size = 1.5, title.size = 2)) +
  tm_shape(subset(nodes, unclass(product) > 5L)) + tm_dots(size = 0.5, fill = "purple3") +
  tm_layout(frame = FALSE)

dev.copy(pdf, "figures/transport_network/trans_africa_network_reduced_22_products.pdf", 
         width = 10, height = 10)
dev.off()         

