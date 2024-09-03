#######################################################
# Optimal Graph Representation of Africa's Road Network
#######################################################

library(fastverse)
fastverse_extend(dggridR, s2, cppRouting, qs, install = TRUE)
set_collapse(mask = c("manip", "helper", "special"), nthreads = 4, sort = FALSE)
fastverse_conflicts()

# Finding nearest cells and removing some cells --------------------------------------------------------

spherical_dist <- qread("data/full_network/africa_full_spherical_distance_matrix_r9.qs")
sp_distances <- spherical_dist$distances
centroids <- spherical_dist$centroids
diag(sp_distances) <- NA
range(sp_distances)
rm(spherical_dist); gc()

# Cells entirely above water: eliminate
# mapview::mapview(dgcellstogrid(dgconstruct(res = 9), centroids$cell), col.regions = NA)
cells_above_water <- as.integer(c(55878, 55796, 55877, 55958, 55876, 56120, 56202, 56283, 56203, 56284, 56122, 56040, 55959, 56121, 56039, # Lake Victoria
                                  55709, 56195, # Lake Taganyika
                                  # Coastline cells
                                  39789, 39627, 39546, 	39464, 	39384, 35398, 35074, 31830, 	31506, 31102, 	39577, 	39578, 	32722,
                                  40472, 	42591, 42835, 43563, 	46241, 46244, 46650, 59094, 59410, 	59488, 59566, 	59401, 	59482,
                                  59644, 	59563, 	59725, 	59643, 59642, 	59641, 	59722, 	59803, 	59640, 	59721, 59883, 59802,
                                  59964, 	59882, 	60204, 	60285, 60365, 60283, 	60930, 	61661, 	61823, 	61986, 	62392, 	62717, 	62716,
                                  62226, 61573, 59295, 59049, 	58885, 	58803, 58721, 	58959, 	58878, 	58957, 	138255, 138336, 138415,
                                  138328, 138002, 	137835, 59016, 	59011, 58929, 	58518, 	58273, 56726, 56402, 56078, 54298, 53813,
                                  53733, 	53652, 53571, 53491, 	53248, 	52924, 	52681, 	52438, 51712, 51473, 50988, 50421, 49935, 49530,
                                  48476, 	48478, 43196, 	40760, 40518, 39951,
                                  # Sinai Peninsula
                                  59336, 59498, 59660, 	59741, 	59579, 	59416, 	59335, 59497, 59659, 59822,
                                  59740, 59578, 59415, 59253, 	59334, 	59496, 	59658, 	59821, 	59739, 	59577,
                                  59414, 	59252, 59333, 59495, 	59657, 	59576, 59413, 	59251, 	59332, 59494,	59412)) 
cow_ind <- centroids$cell %in% cells_above_water

# Madagascar: eliminate
mdg_ind <- centroids$ISO3 == "MDG"

# Subsetting
centroids %<>% subset(!(cow_ind | mdg_ind))
sp_distances %<>% extract(!(cow_ind | mdg_ind), !(cow_ind | mdg_ind))

# mapview::mapview(dgcellstogrid(dgconstruct(res = 9), centroids$cell), col.regions = NA)

# Testing: 70 or 80 thousand gives them same
tmp = fsum(sp_distances < 70000)
sum(tmp)
qsu(tmp)
rm(tmp)

# Nearest neighbours: rows are sources
nn_ind_list <- lapply(mrtl(sp_distances < 70000), which)
nn_cell_list <- nn_ind_list
names(nn_cell_list) <- centroids$cell
nn_cell_list %<>% lapply(function(i) centroids$cell[i])

# Adjusting

# Simple graph based on distances
africa_dist <- qread("data/full_network/africa_full_distance_matrix_r9_adjusted.qs")
africa_dist[c("sources", "centroids")] %<>% lapply(ss, !(cow_ind | mdg_ind))
africa_dist[c("distances", "durations", "distances_nosphere")] %<>% lapply(ss, !(cow_ind | mdg_ind), !(cow_ind | mdg_ind))

if(!identical(centroids, africa_dist$centroids)) stop("Mismatch")

# Creating Basic Full Graphs --------------------------------------------------------

graphs <- list(
  sp_distance = Map(function(z, ind, cell) setNames(z[ind], cell),
                    mrtl(sp_distances), nn_ind_list, nn_cell_list),
  distance = Map(function(z, ind, cell) setNames(z[ind], cell),
                 mrtl(africa_dist$distances), nn_ind_list, nn_cell_list),
  duration = Map(function(z, ind, cell) setNames(z[ind], cell),
                 mrtl(africa_dist$durations), nn_ind_list, nn_cell_list)
) |> lapply(setNames, names(nn_cell_list))

# Graph Data Frames
graph_dfs <- lapply(graphs, function(x) {
  lapply(x, qDF, "to") |> 
    rowbind(idcol = "from", id.factor = FALSE) |> 
    rename(i = cost) |> 
    transformv(c("from", "to"), as.integer)
})

# Combine
graphs_df <- Reduce(join, Map(setNames, graph_dfs, lapply(names(graph_dfs), \(x) c("from", "to", x)))) |> 
             roworder(from, to) 

# Saving
graphs_df |> 
  transform(from = group(from), to = ckmatch(to, unique(from)), 
            from_cell = from, to_cell = to) |> 
  colorder(from, from_cell, to, to_cell) |> 
  fwrite("data/full_network/full_graph_df.csv")

# Also saving cell-data
calib_data <- fread("data/QSE/QSE_model_calibration_data_ctry_min_imp.csv") |> 
  subset(cell %in% unique(graphs_df$from)) |> 
  roworder(cell) |> 
  mutate(id = ckmatch(cell, unique(graphs_df$from))) |> 
  roworder(id) |> 
  colorder(id) 

calib_data |> fwrite("data/QSE/QSE_model_calibration_data_ctry_min_imp_full_graph.csv")


# Creating Optimal Full Graphs --------------------------------------------------------

# Read graphs data!
graphs_df <- fread("data/full_network/full_graph_df.csv")
# Check
descr(africa_dist$durations[with(graphs_df, cbind(from, to))] - graphs_df$duration)
# Add original distances
graphs_df$distance_nosphere <- africa_dist$distances_nosphere[with(graphs_df, cbind(from, to))]
descr(with(graphs_df, distance / distance_nosphere))
graphs_df$distance_nosphere |> replace_na(Inf, set = TRUE)

# Calculating Shortest path 
dist_mat_from_df <- function(df, ...) {
  graph <- makegraph(df) # directed = FALSE # cpp_simplify()
  nodes <- graph$dict$ref
  get_distance_matrix(graph, from = nodes, to = nodes, ...)
}

# RcppParallel::setThreadOptions(numThreads = 8)
system.time({
  distance_mat <- dist_mat_from_df(select(graphs_df, from_cell, to_cell, duration), algorithm = "mch")
})

all_identical(dimnames(distance_mat))

# Getting true matrix
cell_ind <- ckmatch(rownames(distance_mat), africa_dist$centroids$cell)
true_distance_mat <- africa_dist$durations[cell_ind, cell_ind]
diag(true_distance_mat) <- NA
# Check: How much worse is the matrix of all shortest paths on the graph
descr(unattrib(distance_mat / true_distance_mat))
rm(cell_ind, true_distance_mat); gc()

# Now the actual implementation  --------------------

identical(centroids, africa_dist$centroids)

# Can use distance_nosphere or duration here
road_distances <- africa_dist$durations
diag(road_distances) <- 0
dimnames(road_distances) <- NULL
diag(sp_distances) <- 0
dimnames(sp_distances) <- NULL
gc()

# Removing cells
world_res9 <- dgconstruct(res = 9)
all.equal(centroids$cell, dgGEO_to_SEQNUM(world_res9, centroids$lon, centroids$lat)$seqnum)
starts_in_cell <- as.integer(centroids$cell) == as.integer(with(africa_dist$sources, dgGEO_to_SEQNUM(world_res9, lon, lat)$seqnum))
table(starts_in_cell)

all.equal(as.integer(names(nn_cell_list)), centroids$cell)
nn_cell_df <- data.table(from = rep(as.integer(names(nn_cell_list)), vlengths(nn_cell_list)), 
                         to = unlist(nn_cell_list, use.names = FALSE)) 

all.equal(unattrib(nn_cell_df), unattrib(select(graphs_df, from_cell, to_cell)))

compute_weighted_graph <- function(alpha = 30, gamma = 2) {
  oldopts <- options(warn = -1)
  on.exit(options(oldopts))
  raw_graph <- Map(function(x, i, sd, rd) { # x are the nearest neighbours, i the cell index, dst the sperical distances..
    b <- sd[-i]  # Spherical distance from i to all others # -c(i, x_i) # x_i can be part of it
    br <- rd[-i] # Road distance (or travel time)
    vapply(x, function(x_i) { # This loops over nearest neighbours indices 
      a = sd[x_i] # Distance from i to its nearest neighbour x_i
      c = sp_distances[x_i, -i] # Distance from x_i to all other cells except i # -c(i, x_i)
      theta = acos((a^2 + b^2 - c^2)/(2*a*b)) * 180 / pi # Angle between a and b
      ind = which(theta < alpha) # filtering cells sufficiently in the direction of x_i
      b_ind = b[ind]
      br_ind = br[ind] # -c(i, x_i)
      kappa = fmean.default(br_ind/b_ind, w = 1e11/b_ind^gamma) # Average network route/time efficiency in the direction of i
      kappa * a # Final distance/duration measure is average network efficiency in the direction times spherical distance
    }, numeric(1), USE.NAMES = FALSE)
  }, nn_ind_list, seq_along(nn_ind_list), mrtl(sp_distances), mrtl(road_distances))
  add_vars(nn_cell_df, cost = unlist(raw_graph, use.names = FALSE))
}

# Objective function only considers cells where routes start/end within the cell
objective <- function(param, weights = NULL) {
  graph_weighted_distance <- compute_weighted_graph(param[1], param[2])
  distance_mat <- dist_mat_from_df(graph_weighted_distance, algorithm = "mch")
  cell_ind <- ckmatch(africa_dist$centroids$cell, rownames(distance_mat))
  if(!identical(cell_ind, seq_row(distance_mat))) {
    distance_mat <- distance_mat[cell_ind, cell_ind]
  }
  diag(distance_mat) <- NA_real_
  fmean.default((distance_mat %/=% road_distances)[starts_in_cell, starts_in_cell], w = weights)
}

# objective(c(30, 2))
print(objective(c(30, 1)))

# Now finding best parameters for optimal graph (may take a while)
best_param <- optim(c(30, 2), objective, NULL, method = "L-BFGS-B",
                    lower = c(1, 1), upper = c(40, 5))
# Result for distance:
best_param <- list(par = c(27.7176642781289, 1), value = 1.16517890938921, 
                   counts = c(`function` = 23L, gradient = 23L), convergence = 0L, 
                   message = "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH")
# Result for duration:
best_param <- list(par = c(29.8834997902231, 1), value = 1.37535288169957, 
                   counts = c(`function` = 18L, gradient = 18L), convergence = 0L, 
                   message = "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH")

# # Result Weighted: favouring short routes in the objective -> no big difference in optimal parameters
# weights <- (1e9 / sp_distances^2)[starts_in_cell, starts_in_cell]
# diag(weights) <- NA
# objective(c(30, 1), weights)
# 
# best_param_w <- optim(c(30, 2), objective, NULL, weights, method = "L-BFGS-B",
#                       lower = c(1, 1), upper = c(40, 5))
# rm(weights)
# # Result for duration:
# best_param_w <- list(par = c(30.0673814764975, 1), value = 0.930531476795317, 
#                      counts = c(`function` = 28L, gradient = 28L), convergence = 0L, 
#                      message = "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH")

# Market Access Optimization ----------------------------------------------------------

calib_data <- fread("data/QSE/QSE_model_calibration_data_ctry_min_imp_full_graph.csv")

centroids %<>% join(calib_data, on = "cell", drop = "y")

# best_param <- list()
# best_param$par <- c(30, 1)

optimized_graph <- compute_weighted_graph(best_param$par[1], best_param$par[2]) |>
  join(graphs_df |> select(from = from_cell, to = to_cell, sp_distance:distance_nosphere))

# Check
identical(unique(optimized_graph$from), centroids$cell)

# See goodness of fit
objective_graph <- function(graph, weights = NULL) {
  distance_mat <- dist_mat_from_df(graph, algorithm = "mch")
  cell_ind <- ckmatch(africa_dist$centroids$cell, rownames(distance_mat))
  if(!identical(cell_ind, seq_row(distance_mat))) {
    distance_mat <- distance_mat[cell_ind, cell_ind]
  }
  diag(distance_mat) <- NA_real_
  fmean.default((distance_mat %/=% road_distances)[starts_in_cell, starts_in_cell], w = weights)
}

# How many times longer are the shortest paths versus the real routes?
objective_graph(select(optimized_graph, from, to, cost))
# objective_graph(select(optimized_graph, from, to, cost), weights)
# Distance
objective_graph(compute(optimized_graph, cost = distance_nosphere, keep = .c(from, to)))
# objective_graph(compute(optimized_graph, cost = distance_nosphere, keep = .c(from, to)), weights)
# Times
objective_graph(select(optimized_graph, from, to, duration))
# objective_graph(select(optimized_graph, from, to, duration), weights)

# Plot Optimal Graph
settfm(optimized_graph, from_ind = ckmatch(from, calib_data$cell), to_ind = ckmatch(to, calib_data$cell))

map_to_palette <- function(vector, option = "turbo") {
  # Normalize the vector to the range [0, 1]
  normalized_vector <- (vector - min(vector)) / (max(vector) - min(vector))
  # Get the Inferno colors corresponding to the normalized vector
  colors <- viridis::viridis(length(vector), option = "turbo")[as.integer(cut(normalized_vector, breaks = length(vector), labels = FALSE))]
  return(colors)
}

settfm(optimized_graph, cost_colour = map_to_palette(cost))

pdf("figures/full_network/full_network_optimal_duration_graph.pdf", width = 10, height = 10)
with(calib_data, {
  oldpar <- par(mar = c(0,0,0,0))
  plot(lon, lat, cex = NA, pch = 16, axes = FALSE, xlab = NA, ylab = NA, asp=1)
  for (r in seq_row(optimized_graph)) { # & intersects == 2
    ind = c(optimized_graph$from_ind[r], optimized_graph$to_ind[r])
    lines(lon[ind], lat[ind], cex = 0.1, col = optimized_graph$cost_colour[r])
  }
  points(lon, lat, cex = 0.2, pch = 16)
  par(oldpar)
})
dev.off()

# Now compute Market Access
compute_MA <- function(graph) {
  inv_duration <- 1 / dist_mat_from_df(graph, algorithm = "mch")
  diag(inv_duration) <- 0
  setv(inv_duration, NA_real_, 0)
  if(!all_identical(dimnames(inv_duration))) stop("Not Symmetric")
  if(!identical(as.integer(rownames(inv_duration)), centroids$cell)) stop("Mismatch")
  sum(inv_duration %*% centroids$GDP_PPP)
}
system.time(MA <- compute_MA(select(optimized_graph, from, to, cost)))
print(MA)

# Route Efficiency: use only if cost is distance (road_distances variable on line 141 is a distance measure)
descr(optimized_graph$sp_distance)
local_nre <- optimized_graph |> with(sp_distance/cost)
descr(local_nre) # Max US route efficiency is 0.843. I use 0.85
quantile(local_nre, seq(0.99, 1, 0.0005))
optimized_graph %$% all.equal(sp_distance/local_nre, cost)
optimized_graph |> with(pmin(sp_distance/0.85, cost)) |> descr()
rm(local_nre)

# This can take ~10-24h
MA_new <- sapply(centroids$cell, function(i) {
  graph <- compute(optimized_graph, 
                   cost = iif(from == i | to == i, pmin(sp_distance/0.85, cost), cost), 
                   keep = .c(from, to))
  c(N_Edges = sum(graph$from == i | graph$to == i),
    MA = compute_MA(graph), 
    Reduction = sum(optimized_graph$cost) - sum(graph$cost))
})

MA_new %<>% t() %>% qDT() %>% add_vars(select(centroids, cell:pop_wpop), pos = "front")
attr(MA_new, "MA_orig") <- MA
MA_new$MA_ratio = MA_new$MA / MA
MA_new$MA_growth = (MA_new$MA_ratio - 1)*100
descr(MA_new$Reduction)
descr(MA_new$MA_growth * 100)

qsave(MA_new, "results/full_network/africa_full_MA_0.85_NRE.qs")

# Time efficiency: use only if cost is duration (road_distances variable on line 141 is a travel time measure)
descr(optimized_graph$sp_distance)
speed_kmh <- optimized_graph |> with(sp_distance/cost * 60/1000)
descr(speed_kmh) # Max US time efficiency is 50kmh
quantile(speed_kmh, seq(0.99, 1, 0.0005))
optimized_graph %$% all.equal(sp_distance/speed_kmh * 60/1000, cost)
optimized_graph |> with(pmin(sp_distance/50 * 60/1000, cost)) |> descr()
rm(speed_kmh)

# This can take ~10-24h
MA_new <- sapply(centroids$cell, function(i) {
  graph <- compute(optimized_graph, 
                   cost = iif(from == i | to == i, pmin(sp_distance/50 * 60/1000, cost), cost), 
                   keep = .c(from, to))
  c(N_Edges = sum(graph$from == i | graph$to == i),
    MA = compute_MA(graph), 
    Reduction = sum(optimized_graph$cost) - sum(graph$cost))
})

MA_new %<>% t() %>% qDT() %>% add_vars(select(centroids, cell:pop_wpop), pos = "front")
attr(MA_new, "MA_orig") <- MA
MA_new$MA_ratio = MA_new$MA / MA
MA_new$MA_growth = (MA_new$MA_ratio - 1)*100

descr(MA_new$Reduction)
descr(MA_new$MA_growth * 100)

qsave(MA_new, "results/full_network/africa_full_MA_50kmh_NTE.qs")



# Visual Exploration: Travel Time ---------------------------------

MA_new <- qread("results/full_network/africa_full_MA_50kmh_NTE.qs")
graph_count <- fread("data/full_network/full_graph_df.csv") %>%
  subset(from_cell %in% MA_new$cell & to_cell %in% MA_new$cell) %>% {
    join(count(., cell = from_cell), count(., cell = to_cell), on = "cell") %>% 
      transform(N = N + N_y, N_y = NULL)
  }
MA_new %<>% join(graph_count, on = "cell", drop = TRUE)

qsu(MA_new)

fastverse_extend(ggplot2, viridis, install = TRUE)

# Average Travel Time Reduction (Investment)
MA_new %>% 
  ggplot(aes(x = lon, y = lat, fill = Reduction/60/N)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  scale_fill_viridis(option = "H", n.breaks = 10, # trans = "sqrt", 
                     na.value = "#30123BFF") +
  labs(title = "Average Travel Time Reduction (Hours, All Directions)", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/average_travel_time_reduction_50kmh.pdf", width = 6, height = 6)

# Total MA Gain (%)
MA_new %>% 
  ggplot(aes(x = lon, y = lat, fill = MA_ratio-1)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  scale_fill_viridis(option = "H", breaks = c(0.1, 1, 2, 5, 10, 15, 20) * 1e-4,  # c(0.1, 1, 2, 4, 6, 8, 10, 12, 14, 16) * 1e-4, #  n.breaks = 10, 
                     trans = "sqrt", 
                     # limits = c(1e-6, NA),
                     na.value = "#30123BFF",
                     labels = scales::label_percent()) +
  labs(title = "Continental Market Access Gain (%)", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/MA_50kmh_NTE_growth.pdf", width = 6, height = 6)


# MA Gain per hour reduction
MA_new %>% 
  mutate(MA_gain_per_minute = (MA - MA/MA_ratio)/Reduction) %>% # descr()
  ggplot(aes(x = lon, y = lat, fill = MA_gain_per_minute/1e6)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  scale_fill_viridis(option = "H", 
                     # limits = c(NA, 50),
                     breaks = c(0.1, 1, 2, 5, seq(10, 80, 10)), # c(0.1, 1, 2, seq(5, 50, 5)), #  n.breaks = 10, 
                     labels = function(x) paste0("$", signif(x, 2), "M"),
                     trans = "sqrt", 
                     na.value = "#30123BFF") +
  labs(title = "MA Gain Per Minute Reduction (2015 GDP PPP per Minute)", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/MA_50kmh_NTE_gain_per_minute.pdf", width = 6, height = 6)


# Visual Exploration: Road Distance ---------------------------------

MA_new <- qread("results/full_network/africa_full_MA_0.85_NRE.qs")
qsu(MA_new)

fastverse_extend(ggplot2, viridis, install = TRUE)

# Average Road Distance Reduction (Investment)
MA_new %>% 
  # st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% select(Reduction) %>% mapview::mapview()
  # subset(Reduction > 0) %>%
  mutate(Reduction = iif(Reduction <= 0, fnth(Reduction, 0.98), Reduction)) %>%
  ggplot(aes(x = lon, y = lat, fill = Reduction/1000/N)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  scale_fill_viridis(option = "H", n.breaks = 10, trans = "log10", 
                     na.value = "#30123BFF") +
  labs(title = "Average Road Distance Reduction (Km, All Directions)", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/average_road_distance_reduction_0.85.pdf", width = 6, height = 6)

# Total MA Gain (%)
MA_new %>% 
  ggplot(aes(x = lon, y = lat, fill = MA_ratio-1)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  scale_fill_viridis(option = "H", breaks = c(0.1, 1, 2, 5, 10, 15) * 1e-4,  # c(0.1, 1, 2, 4, 6, 8, 10, 12, 14, 16) * 1e-4, #  n.breaks = 10, 
                     trans = "sqrt", 
                     # limits = c(1e-6, NA),
                     na.value = "#30123BFF",
                     labels = scales::label_percent()) +
  labs(title = "Continental Market Access Gain (%)", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/MA_0.85_NRE_growth.pdf", width = 6, height = 6)


# MA Gain per Km reduction
MA_new %>% 
  mutate(MA_gain_per_km = (MA - MA/MA_ratio)/Reduction*1e6) %>% # descr()
  ggplot(aes(x = lon, y = lat, fill = MA_gain_per_km/1e6)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  scale_fill_viridis(option = "H", 
                     # limits = c(NA, 50),
                     breaks = c(0.1, 1, 2, 5, seq(10, 70, 10)), # c(0.1, 1, 2, seq(5, 50, 5)), # n.breaks = 10, 
                     labels = function(x) paste0("$", signif(x, 2), "M"),
                     trans = "sqrt", 
                     na.value = "#30123BFF") +
  labs(title = "MA Gain Per Km Reduction (2015 GDP PPP per Km)", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/MA_0.85_NRE_gain_per_km.pdf", width = 6, height = 6)




