###############################################################################
# Identifying Major Transport Routes (for GE analysis of Trans-African Network)
###############################################################################

library(fastverse)
set_collapse(mask = c("manip", "helper", "special"))
fastverse_extend(qs, sf, sfnetworks, tmap, install = TRUE)
fastverse_conflicts()

load("data/transport_network/trans_africa_network_param.RData")
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
  st_network_paths(net_param, from = i, to = large_cities, weights = "duration", mode = "all")) |>
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

net_main_routes <- as_sfnetwork(subset(edges_param, large_city_paths$edges) |> 
                                  mutate(ug_cost = distance / 1000 * ug_cost_km), 
                                directed = FALSE, edges_as_lines = TRUE)

filter_smooth2 <- function(net, funs = list(passes = "mean", gravity = "mean", gravity_rd = "mean", gravity_dur = "mean", 
                                            distance = "sum", duration = "sum", duration_imp = "sum", 
                                            border_dist = "sum", border_time = "sum", total_dist = "sum", total_time = "sum",
                                            ug_cost = "sum", "ignore")) {
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
  tm_layout(frame = FALSE) 

dev.copy(pdf, "figures/transport_network/trans_africa_network_reduced_47_duration.pdf", 
         width = 10, height = 10)
dev.off()

# Join nodes
nodes_sf <- st_as_sf(net_main_routes, "nodes")
nodes_sf %<>% transform(qDF(st_coordinates(.))) %>% 
  join(transform(nodes_param, qDF(st_coordinates(geometry))), on  = c("X", "Y"), drop = TRUE) %>%
  select(-X, -Y, -.tidygraph_node_index)
net_main_routes %<>% activate("nodes") %>% dplyr::mutate(qDF(atomic_elem(nodes_sf)))

results <- list(
  fastest_routes = c(large_city_paths, list(network = net_main_routes))
)

# Now With All Routes, incl. new ones ---------------------------------------

identical(st_geometry(net_param, "edges"), edges_param$geometry)

net_ext_data <- rowbind(select(edges_param, from, to, add, gravity_rd, distance, duration, duration_imp, 
                               border_dist, border_time, total_time, cost_per_km = ug_cost_km, geometry) |>
                          mutate(total_cost = distance / 1000 * cost_per_km) , 
                        select(add_links_param, from, to, add, distance, duration = duration_100kmh,  
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
           col.legend = tm_legend("Sum of Gravity", position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_shape(st_as_sf(activate(net_main_routes, "nodes"))) + tm_dots(size = 0.2) +
  tm_shape(subset(nodes, large_cities)) + tm_dots(size = 0.3, fill = "red") +
  tm_layout(frame = FALSE) 

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
# edges_real <- qread("data/transport_network/edges_real_simplified.qs") 
nodes <- results[[nname]]$network |> st_as_sf("nodes") 
edges <- results[[nname]]$network |> st_as_sf("edges") 

# Classification
settfm(edges, speed_kmh = (distance / 1000) / (duration / 60))
# settfm(edges_real, speed_kmh = (edges_param$distance/1000)/(edges_param$duration/60))
# edges_real %<>% subset(results$fastest_routes$edges)

# nodes <- fread(sprintf("data/transport_network/largest_pcities/%s_graph_nodes.csv", nname))
settfm(nodes, major_city_port = population > 2e6 | outflows > 1e6)
sum(nodes$major_city_port)
largest <- c("Dakar - Senegal", "Casablanca - Morocco", "Abidjan - Cote d'Ivoire", 
             "Kumasi - Ghana", "Algiers - Algeria", "Lagos - Nigeria", "Kano - Nigeria", 
             "Yaounde - Cameroon", "Luanda - Angola", "Kinshasa - Congo (Kinshasa)", 
             "Johannesburg - South Africa", "Cape Town - South Africa", "Cairo - Egypt", 
             "Khartoum - Sudan", "Nairobi - Kenya", "Addis Ababa - Ethiopia", 
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

