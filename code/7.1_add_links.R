#############################################
# Adding High-Value Links to Existing Network
#############################################

# Note: This is a self-contained section due to the use of specific indices for manual adjustments. 
# It uses a previous version of the network identical to the current/generated version but where the nodes and edges are in a different (random) order. 
# The output is a spatial data frame 'add_links' with the final proposed links. These can be added to the current/generated version.

library(fastverse)
set_collapse(mask = c("manip", "helper", "special"), nthreads = 4)
fastverse_extend(qs, sf, s2, units, sfnetworks, stplanr, install = TRUE)
source("code/helpers/helpers.R")
fastverse_conflicts()

net <- qread("data/transport_network/old/net_discrete_final.qs")
dist_ttime_mats <- qread("data/transport_network/old/net_dist_ttime_mats.qs")
sym_dist_mat <- (dist_ttime_mats$distances + t(dist_ttime_mats$distances)) / 2

## Plotting
plot(net)
nodes <- net |> activate("nodes") |> tidygraph::as_tibble() |> 
  mutate(.tidygraph_node_index = NULL)
edges <- net |> activate("edges") |> tidygraph::as_tibble() |> 
  mutate(.tidygraph_edge_index = NULL)

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

add_links <- with(nodes_coord, lapply(mrtl(routes_ind[keep_routes, ]), function(r) st_linestring(cbind(lon[r], lat[r])))) |> 
              st_sfc(crs = 4326) |> st_as_sf()
add_links_df <- line2df(add_links)[-1] %*=% 1000 %>% dapply(as.integer)
edges_df <- line2df(edges)[-1] %*=% 1000 %>% dapply(as.integer)
m <- add_links_df %in% edges_df | add_links_df %in% gv(edges_df, c(3:4, 1:2))
nrow(add_links) - sum(m)
add_links <- add_links[!m, ]
rm(add_links_df, edges_df)

add_links |> qsave("data/transport_network/add_links_network_30km_alpha45_mrEU_fmr15.qs")
add_links <- qread("data/transport_network/add_links_network_30km_alpha45_mrEU_fmr15.qs")

# Manual Adjustments: Removing links crossing waterbodies or protected areas
# mapview(edges, map.types = c(mapviewGetOption("basemaps"), "Esri.WorldStreetMap", "Esri.WorldTerrain")) + mapview(nodes) + mapview(add_links, color = "green") # [remove, ]
remove <- c(305, 396, 399, 400, 402, 512, 513, 514, 516, 403, 409, 410, 390, 386, 457, 420, 429, 428, 426, 427, 434,
            358, 357, 363, 364, 344, 345, 165, 175, 106, 27, 10, 69, 68, 433, 490, 499, 444, 326, 324, 419, 279, 317, 178, 331, 270)
add <- list(c(63, 105), c(862, 717), c(499, 519), c(499, 455), c(312, 379), c(1184, 1271), c(46, 37), c(1071, 842),
            c(901, 751), c(751, 710), c(593, 514))
add_links <- rbind(add_links[-remove, ],
                    with(nodes_coord, lapply(add, function(r) st_linestring(cbind(lon[r], lat[r])))) |> st_sfc(crs = 4326) |> st_as_sf())

nrow(add_links)
if(nrow(add_links) != 481) stop("File should have 481 rows")
add_links |> qsave("data/transport_network/add_links_network_30km_alpha45_mrEU_fmr15_adjusted.qs")
