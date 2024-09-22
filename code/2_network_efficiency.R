#######################################################
# Local and Global Estimates of Road Network Efficiency
#######################################################

library(fastverse)
fastverse_extend(qs, ggplot2, install = TRUE)

africa_dist <- qread("data/full_network/africa_full_distance_matrix_r9_adjusted.qs")

spherical_dist <- qread("data/full_network/africa_full_spherical_distance_matrix_r9.qs")
identical(spherical_dist$centroids, africa_dist$centroids)

# TRA(spherical_dist$distances, africa_dist$centroids$pop_wpop, "fill") %/=%
weights <- 1e6 / (pi * spherical_dist$distances^2)
diag(weights) <- 0
weights <- weights %/=% rowSums(weights, na.rm = TRUE)
frange(weights)
fquantile(weights, seq(0.95, 1, 0.001))
rowSums(weights)
fsum(weights)

# Plot WorldPop Population
# <Figure A1: RHS>
africa_dist$centroids %>% 
  ggplot(aes(x = lon, y = lat, fill = pop_wpop)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  viridis::scale_fill_viridis(option = "H", n.breaks = 8, transform = "log10", 
                              labels = scales::label_log()) +
  labs(title = "WorldPop 2020 Total Population Estimate", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/WorldPop2020_hex9_grid.pdf", width = 6, height = 6)


################################
# Average Network Speed per Cell
################################

# Calculating average network speed in the cell
distances <- africa_dist$distances / 1000
durations <- africa_dist$durations / 60
diag(distances) <- 1
diag(durations) <- 1

# weights <- 1e6 / (pi * distances^2)
# diag(weights) <- 0
# weights <- weights %/=% rowSums(weights, na.rm = TRUE)
# frange(weights)

dd_ratio <- distances / durations
rdd <- frange(dd_ratio) %T>% print()
ans <- rowSums(dd_ratio * weights, na.rm = TRUE)
qsu(ans)
qsu(ans, w = africa_dist$centroids$pop_wpop)
# # ans[ans < rdd[1L] | ans > rdd[2L]] <- NA
# # frange(ans)
rm(dd_ratio); gc() # , weights

# <Figure A1: LHS>
africa_dist$centroids %>% 
  add_vars(ans) %>% 
  ggplot(aes(x = lon, y = lat, fill = ans)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  viridis::scale_fill_viridis(option = "H", limits = c(10, 70)) +
  labs(title = "Average Network Speed (km/h)", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/average_network_speed.pdf", width = 6, height = 6)

# # Interactive Plot
# africa_dist$centroids %>% 
#   add_vars(ans) %>% 
#   fselect(ans, cell, lon, lat) %>% 
#   ftransform(longitude = lon, latitude = lat) %>% 
#   sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
#   mapview::mapview(zcol = "ans", lwd = 0, alpha = 2)


################################
# Route Efficiency of Network
################################

identical(spherical_dist$centroids, africa_dist$centroids)

spherical_ratio <- spherical_dist$distances 
diag(spherical_ratio) <- NA
spherical_ratio %/=% africa_dist$distances_nosphere

nre <- rowSums(spherical_ratio * weights, na.rm = TRUE)
qsu(nre)
qsu(nre, w = africa_dist$centroids$pop_wpop)
rm(spherical_ratio)

# <Figure 4: RHS>
spherical_dist$centroids %>% 
  add_vars(nre) %>% 
  ggplot(aes(x = lon, y = lat, fill = nre)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  viridis::scale_fill_viridis(option = "H", n.breaks = 7) +
  labs(title = "Network Route Efficiency (Geodesic/Road Distance)", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/network_route_efficiency.pdf", width = 6, height = 6)


################################
# Time Efficiency of Network
################################

# https://christopherwolfram.com/projects/efficiency-of-road-networks/
# Idea: geodesic distance divided by travel time (or length): road network efficiency...

identical(spherical_dist$centroids, africa_dist$centroids)

spherical_ratio <- spherical_dist$distances 
diag(spherical_ratio) <- NA
spherical_ratio %/=% 1000
spherical_ratio %/=% (africa_dist$durations / 60)

nte <- rowSums(spherical_ratio * weights, na.rm = TRUE)
qsu(nte)
qsu(nte, w = africa_dist$centroids$pop_wpop)
rm(spherical_ratio)

# <Figure 4: LHS>
spherical_dist$centroids %>% 
  add_vars(nte) %>% 
  ggplot(aes(x = lon, y = lat, fill = nte)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  viridis::scale_fill_viridis(option = "H", n.breaks = 8) +
  labs(title = "Network Time Efficiency (Geodesic km/Route h)", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/network_time_efficiency.pdf", width = 6, height = 6)

# Saving Estimates --------------------

qsave(add_vars(spherical_dist$centroids, ANS = ans, NRE = nre, NTE = nte), 
      file = "results/full_network/network_efficiency.qs")


##################################################
# Global Efficiency Estimates using Law of Gravity
##################################################

# Excluding Madagascar
MDG_ind <- africa_dist$centroids$ISO3 %==% "MDG"

weights <- tcrossprod(africa_dist$centroids$pop_gpw4) %/=% (spherical_dist$distances^2)
replace_inf(weights, set = TRUE)

# Average Network Speed
fmean.default(distances / durations, w = weights)
fmean.default((distances / durations)[-MDG_ind, -MDG_ind], w = weights[-MDG_ind, -MDG_ind])

# Network Route Efficiency
fmean.default(spherical_dist$distances / africa_dist$distances_nosphere, w = weights)
fmean.default((spherical_dist$distances / africa_dist$distances_nosphere)[-MDG_ind, -MDG_ind], w = weights[-MDG_ind, -MDG_ind])

# Network Time Efficiency
fmean.default(spherical_dist$distances / africa_dist$durations, w = weights) * 60/1000
fmean.default((spherical_dist$distances / africa_dist$durations)[-MDG_ind, -MDG_ind], w = weights[-MDG_ind, -MDG_ind]) * 60/1000

