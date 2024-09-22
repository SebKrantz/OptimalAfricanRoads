############################################
# Processing Distance and Travel Time Matrix
############################################

library(fastverse)
fastverse_extend(qs, s2, install = TRUE)

# Note: need to generate using '0_get_distance_matrix.R' or download from https://drive.google.com/file/d/1oE_9i3SdqvYKcdl880uS9q774dD9KXmS/view?usp=sharing
africa_dist <- qread("data/full_network/africa_full_distance_matrix_r9.qs")
fnobs.default(africa_dist$durations)
fnobs.default(africa_dist$distances)
N <- fnobs(africa_dist$durations)
identical(N, fnobs(t(africa_dist$durations)))
sort(qtab(N), decreasing = TRUE)
names(N) <- NULL
identical(africa_dist$sources, africa_dist$destinations)
identical(unattrib(africa_dist$sources), slt(africa_dist$centroids, lon, lat))
# The coordinates are snapped to the grid...
print(pwcor(africa_dist$sources, slt(africa_dist$centroids, lon, lat)), digits = 5)

# <Figure 3>
par(mfrow = c(1, 2))
# This shows the connected points
africa_dist$centroids[N >= 200, ] %$% plot(lon, lat, pch = 16, cex = 0.25, main = "Grid Centroids")
#points(-6.5, 23, col = "red", pch = "x") # This is in the middle of nowhere, no connections...
# Non-reachable continental areas: 
africa_dist$centroids %>% 
  sbt(N < 11500 & between(lon, -13, -2) & between(lat, 19, 28)) %$%
  points(lon, lat, col = "grey80", pch = 16, cex = 0.25) 
africa_dist$centroids %>% 
  sbt(N < 11500 & between(lon, 9.2, 36) & between(lat, -10, 15)) %$%
  points(lon, lat, col = "grey80", pch = 16, cex = 0.25) 
# Non-reachable areas in Madagascar
africa_dist$centroids %>% 
  sbt(N < 200 & between(lon, 45.5, 50) & between(lat, -20, -11)) %$%
  points(lon, lat, col = "grey80", pch = 16, cex = 0.25) 
# This shows the actual rout start and and points
africa_dist$sources[N > 200, ] %$% plot(lon, lat, pch = 16, cex = 0.25, main = "Route Starting / Destination Points")
# points(-6.5, 23, col = "red", pch = "x") # This is in the middle of nowhere, no connections...
par(mfrow = c(1, 1))

dev.copy(pdf, "figures/full_network/route_starting_points.pdf", width = 14, height = 8.27)
dev.off()

# # Plot real stating positions 
# africa_dist$sources[N >= 20, ] %>% 
#   ftransform(longitude = lon, latitude = lat) %>% 
#   sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
#   mapview::mapview()

missing_centroids <- africa_dist$centroids %>% 
  fsubset((N < 11500 & between(lon, -13, -2) & between(lat, 19, 28)) |
          (N < 11500 & between(lon, 9.2, 36) & between(lat, -10, 15)) |
          (N < 200 & between(lon, 45.5, 50) & between(lat, -20, -11)))

# Subsetting data
ind <- which(N > 200 | africa_dist$centroids$cell %in% missing_centroids$cell)
africa_dist %<>% lapply(function(x) if(is.matrix(x)) x[ind, ind] else ss(x, ind))
gc()
# All continental cells + madagascar
africa_dist$centroids %$% plot(lon, lat, pch = ".", main = "Grid Centroids")

# Now some processing: first basic imputation with max for zero distances
distances_adj <- africa_dist$distances
durations_adj <- africa_dist$durations
diag(durations_adj) <- diag(distances_adj) <- NA
# fquantile(distances_adj, probs = seq(0,0.03,0.001))
# fquantile(durations_adj, probs = seq(0,0.03,0.001))
# Cells are 50 km apart (dggridR::dgconstruct(area = 5000))
# distances_adj[distances_adj < 20000] <- NA # less than 30 should not be, even with snapping
# durations_adj[durations_adj < 20] <- NA # going at 120km'h means 25 min 
# Summary Stats
qsu.default(distances_adj, stable.algo = FALSE)
qsu.default(durations_adj, stable.algo = FALSE)
# Examining zero distances
kit::count(distances_adj, 0)
w0 <- which(distances_adj == 0, arr.ind = TRUE)
# These are areas where the starting and end point of the route (snapped to the grid) is the same
africa_dist$sources[c(w0[1,]), ]

# Now following Graff (2019) adjusting matrices
africa_dist %$% identical(sources, destinations)
starting_points <- africa_dist$sources %$% s2_lnglat(lon, lat)
centroid_points <- africa_dist$centroids %$% s2_lnglat(lon, lat)

# First imputing missing cells with nearest neighbours value
miss_ind <- which(africa_dist$centroids$cell %in% missing_centroids$cell)
qtab(fnobs(t(distances_adj[miss_ind, ])))
kit::count(distances_adj[miss_ind, ], 0) # fortunately few zero routes
# Replacing starting points
good_starting_points <- starting_points[-miss_ind]
repl_ind <- sapply(miss_ind, function(i) which.min(s2_distance(centroid_points[i], good_starting_points)))
sources_adj <- africa_dist$sources
sources_adj[miss_ind, ] <- sources_adj[-miss_ind, ][repl_ind, ]  
distances_adj[miss_ind, ] <- distances_adj[-miss_ind, ][repl_ind, ] 
durations_adj[miss_ind, ] <- durations_adj[-miss_ind, ][repl_ind, ] 
distances_adj[, miss_ind] <- distances_adj[, -miss_ind][, repl_ind] 
durations_adj[, miss_ind] <- durations_adj[, -miss_ind][, repl_ind] 
qtab(fnobs(durations_adj))
starting_points <- sources_adj %$% s2_lnglat(lon, lat)
identical(good_starting_points, starting_points[-miss_ind])
start_cen_dist <- s2_distance(starting_points, centroid_points)
descr(start_cen_dist)
rm(good_starting_points); gc()

# Adding walking distance
descr(start_cen_dist / fmean(distances_adj))
setop(distances_adj, "+", start_cen_dist)
setop(distances_adj, "+", start_cen_dist, rowwise = TRUE)
# On average 0.53% increase
fmean.default(replace_inf(distances_adj / africa_dist$distances, set = TRUE))

# Travel time walking following Graff (2019) is 4km/h, thus 1/15 km/min, so 15min per km (seems too much though)
# Better: 10km/h, so 1/6 km/min so 6 min per km
descr((start_cen_dist / 1000 * 6) / fmean(durations_adj))
setop(durations_adj, "+", (start_cen_dist / 1000 * 6))
setop(durations_adj, "+", (start_cen_dist / 1000 * 6), rowwise = TRUE)
# On average 8.6% increase walking at 4km/h, and 3.4% moving at 10km/h
# Now creating spherical distance matrix, takes a minute
system.time(
  spherical_distances <- geodist::geodist(cbind(x = s2_x(centroid_points), y = s2_y(centroid_points)), measure = "geodesic")
  # lapply(seq_along(centroid_points), function(i) s2_distance(centroid_points[i], centroid_points)) %>% qM()
)
# Some checks
isSymmetric(spherical_distances)
allv(diag(spherical_distances), 0)
# Check road network efficiency
fquantile(spherical_distances / distances_adj, seq(0, 1, 0.1))

# Creating durations matrix at 10km/h
spherical_durations <- (spherical_distances / 1000) %*=% 6
# Check: for some routes the spherical duration is shorter than the road distance
fquantile(spherical_durations / durations_adj, seq(0, 1, 0.1))

# Imputation with spherical distances and durations
diag(distances_adj) <- diag(durations_adj) <- 0
distances_adj_nosphere <- copy(distances_adj)
# Checking remaining NA's (apart from Madagascar)
MDG_ind <- africa_dist$centroids$ISO3 %==% "MDG"
DD_na <- is.na(distances_adj) | is.na(durations_adj)
DD_na[MDG_ind, ] <- DD_na[, MDG_ind] <- FALSE
rem_na <- funique(c(which(DD_na, arr.ind = TRUE)))
length(rem_na)
africa_dist$centroids[rem_na, ] # Seems mainly distances between cells in the sahara
qtab(rem_na %in% miss_ind)
africa_dist$centroids %$% plot(lon, lat, pch = ".")
africa_dist$centroids[rem_na, ] %$% points(lon, lat, pch = "x", col = "red")

# Replacing missing values
setv(distances_adj, DD_na, spherical_distances)
setv(durations_adj, DD_na, spherical_durations)

# Now replacing accordingly 
ind <- which(durations_adj > spherical_durations)
length(ind) / length(durations_adj) * 100 # 0.67% of the data
setv(distances_adj, ind, spherical_distances)
setv(durations_adj, ind, spherical_durations)
rm(ind)

# Finally: dealing with Madagascar: use spherical distances and durations as well
allNA(distances_adj[MDG_ind, -MDG_ind])
distances_adj[MDG_ind, -MDG_ind] <- spherical_distances[MDG_ind, -MDG_ind]
allNA(distances_adj[-MDG_ind, MDG_ind])
distances_adj[-MDG_ind, MDG_ind] <- spherical_distances[-MDG_ind, MDG_ind]
allNA(durations_adj[MDG_ind, -MDG_ind])
durations_adj[MDG_ind, -MDG_ind] <- spherical_durations[MDG_ind, -MDG_ind]
allNA(durations_adj[-MDG_ind, MDG_ind])
durations_adj[-MDG_ind, MDG_ind] <- spherical_durations[-MDG_ind, MDG_ind]

# Finally: set NA's inside Madagascar to spherical values
MDG_nas <- unique(c(whichNA(durations_adj), whichNA(distances_adj)))
setv(distances_adj, MDG_nas, spherical_distances)
setv(durations_adj, MDG_nas, spherical_durations)
rm(MDG_nas)

# Final checks
allv(diag(durations_adj), 0)
anyNA(durations_adj)
qsu.default(durations_adj)

# Same for distances
allv(diag(distances_adj), 0)
anyNA(distances_adj)
qsu.default(distances_adj)

# Travel speed in km/h
kmh = (distances_adj / 1000) / (durations_adj / 60)
.range(kmh)

# Saving:
list(durations = durations_adj,
     distances = distances_adj, 
     distances_nosphere = distances_adj_nosphere, 
     sources = sources_adj, 
     centroids = africa_dist$centroids) %>% 
  qsave("data/full_network/africa_full_distance_matrix_r9_adjusted.qs")

list(distances = spherical_distances, 
     centroids = africa_dist$centroids) %>% 
  qsave("data/full_network/africa_full_spherical_distance_matrix_r9.qs")

# Save as csv 
durations_adj %>% 
  setColnames(africa_dist$centroids$cell) %>% 
  qDT() %>% 
  fwrite("data/QSE/model_durations_matrix.csv")


