
split_large_dist_matrix <- function(data, chunk_size = 100, verbose = FALSE) {
  n = nrow(data)
  res_list = list()
  
  # Loop over each chunk to compute the pairwise distances and travel times
  count = 0
  for (i in seq(1, n, by = chunk_size)) {
    for (j in seq(1, n, by = chunk_size)) {
      # Define the row indices for the current chunks
      rows_i = i:min(i + chunk_size - 1, n)
      rows_j = j:min(j + chunk_size - 1, n)
      
      # Extract the data for the current chunks
      ds_i = data[rows_i, ]
      ds_j = data[rows_j, ]
      
      if(verbose) {
        count = count + 1L
        cat(count," ")
      }
      
      # Perform the API call for the current chunks
      r_ij = osrmTable(src = ds_i, dst = ds_j, measure = c('duration', 'distance'))
      
      # Store the result in a list for later combination
      res_list[[paste(i, j, sep = "_")]] = r_ij
    }
  }
  
  # Combine the results from the list into one large matrix for durations and distances
  res_sources = matrix(NA, n, 2)
  res_destinations = matrix(NA, n, 2)
  res_durations = matrix(NA, n, n)
  res_distances = matrix(NA, n, n)
  for (i in seq(1, n, by = chunk_size)) {
    for (j in seq(1, n, by = chunk_size)) {
      rows_i = i:min(i + chunk_size - 1, n)
      rows_j = j:min(j + chunk_size - 1, n)
      
      # Retrieve the result from the list
      r_ij = res_list[[paste(i, j, sep = "_")]]
      
      # Place the result into the corresponding location in the matrix
      res_sources[rows_i, ] = qM(r_ij$sources)
      res_destinations[rows_j, ] = qM(r_ij$destinations)
      res_durations[rows_i, rows_j] = r_ij$durations
      res_distances[rows_i, rows_j] = r_ij$distances
    }
  }
  
  # Create a result list to return
  res = list(
    sources = qDF(copyAttrib(mctl(res_sources), data)),
    destinations = qDF(copyAttrib(mctl(res_destinations), data)),
    durations = res_durations,
    distances = res_distances
  )
  
  rn = rownames(res$sources)
  if(length(rn) && suppressWarnings(!identical(as.integer(rn), seq_along(rn)))) {
    dimnames(res$durations) <- dimnames(res$distances) <- list(rn, rn)
  }
  
  return(res)
}

# Function for deduplication
largest_within_radius <- function(data, coords = c("lon", "lat"), size = "population", radius_km = 100, 
                                  iter = 5, return_dist = FALSE) {
  for (i in 1:iter) {
    dmat <- s2_lnglat(data[[coords[1]]], data[[coords[2]]]) |> st_distance() |> set_units("km")
    if (return_dist) return(dmat)
    # Getting the unique largest items in each column less than radius_km apart
    keep_ind <- ((dmat < as_units(radius_km, "km")) * data[[size]]) |> dapply(which.max) |> unique(sort = TRUE)
    data <- data |> ss(keep_ind)  
  }
  return(data)
}

join_within_radius <- function(d1, d2, coords = c("lon", "lat"), size = "population", radius_km = 100) {
  d1_p <- s2_lnglat(d1[[coords[1]]], d1[[coords[2]]])
  d2_p <- s2_lnglat(d2[[coords[1]]], d2[[coords[2]]])
  dmat <- st_distance(d2_p, d1_p) |> set_units("km")
  d1[[size]] <- fsum((dmat < as_units(radius_km, "km")) * d2[[size]])
  d1
}

max_ratio <- function(alpha = 20) { # Angle between Hypothenuse and Adjacent (a)
  # Pythagoras:
  h = 100    # Hypothenuse    
  g = sin(alpha * pi/180) * h # Length of "gegenüber"
  a = cos(alpha * pi/180) * h # Length of adjacent
  # sqrt(g^2 + a^2)             # Pythagoras Theorem
  # (g+a-h)/h*100 # Percent delay
  (g+a)/h # Max ratio
}

# # https://www.calculator.net/triangle-calculator.html
# a <- s2_distance(p1, p2) 
# b <- s2_distance(p1, p3) 
# c <- s2_distance(p2, p3) 
# # Angle between a and b: In general, c is the side opposite the angle!!
# acos((a^2 + b^2 - c^2)/(2*a*b)) * 180 / pi
intercepted_routes <- function(routes_ind, points_df, travel_dist_mat = NULL, 
                               alpha = 20, mr = NULL, two_sided = TRUE, strict = FALSE, 
                               feasible = FALSE, fmr = mr) {
  if(!two_sided) strict = FALSE
  # Maximum delay based on angle
  if(is.null(mr)) mr = max_ratio(alpha) 
  # Spherical distance matrix
  sdmat = points_df |> with(s2_lnglat(lon, lat)) |> st_distance() |> unclass()
  diag(sdmat) = Inf
  # Function to calculate interception
  is_intercepted <- function(from, to) {
    sd = sdmat[from, to] # Spherical Distance
    orsd = sdmat[from, ] # Other routes from origin (row)
    if (any(orsd < sd)) { # If any smaller than current route
      ind = which(orsd < sd) # < not <= to avoid including the route itself
      orsd_ind = orsd[ind]
      c = sdmat[ind, to] # Routes from these points to the destination 
      angle = unattrib(acos((sd^2 + orsd_ind^2 - c^2)/(2*sd*orsd_ind)) * 180 / pi) # Triangle Equation: Angle between sd and orsd
      check = is.finite(angle) & angle < alpha & (orsd_ind + c) / sd < mr
      if (any(check)) { # If any intercepting points, also check travel distance
        if (feasible || is.null(travel_dist_mat)) return(if(strict) ind[check] else TRUE)
        ind = ind[check]
        rd = travel_dist_mat[from, to] # Note: this is in m!!‚
        orrd_ind = travel_dist_mat[from, ind]
        c = travel_dist_mat[ind, to]
        if (any((orrd_ind + c) / rd < mr)) return(if(strict) ind else TRUE)
      }
    }
    return(if(strict) integer(0) else FALSE)
  } 
  non_feasible <- function(from, to) travel_dist_mat[from, to] / sdmat[from, to] < fmr
  # Intercepted routes
  interc_routes = logical(nrow(routes_ind))
  # Loop to check if routes are intercepted by other points
  for (i in seq_row(routes_ind)) {
    from = routes_ind[i, 1L]
    to = routes_ind[i, 2L]
    interc_routes[i] = if(two_sided && strict) any(is_intercepted(from, to) %in% is_intercepted(to, from)) else if(two_sided) 
      is_intercepted(from, to) || is_intercepted(to, from) else is_intercepted(from, to)
    if(feasible) interc_routes[i] = interc_routes[i] | non_feasible(from, to)
  }
  return(interc_routes)
}


proj_crs_from_bbox <- function(bbox) {
  # Calculate the bounding box and its centroid
  centroid <- st_centroid(st_as_sfc(bbox))
  # Extract longitude from the centroid to calculate UTM zone
  lon <- st_coordinates(centroid)[1] # assuming the first column is longitude
  utm_zone <- floor((lon + 180) / 6) + 1
  # Construct the UTM CRS string. Northern hemisphere zones are positive; southern are negative.
  hemisphere <- ifelse(st_coordinates(centroid)[2] >= 0, 'north', 'south')
  epsg_code <- ifelse(hemisphere == 'north',
                      paste0('326', stringr::str_pad(utm_zone, width = 2, pad = "0")),
                      paste0('327', stringr::str_pad(utm_zone, width = 2, pad = "0")))
  paste0('EPSG:', epsg_code)
}

deg_m <- function(m) m / (40075017 / 360)  

m_to_rad <- function(m) m / 6371000 # Average Earth radius in meters

nodes_max_passes <- function(segments, attrib = "passes") {
  points = line2points(segments)
  passes = segments |> fmutate(id = seq_along(geometry)) |> qDF() |> get_vars(c("id", attrib))
  points = points |> join(passes, on = "id", verbose = 0)
  geom = points |> st_geometry()
  .subset(points, attrib) |> fmax(group(mctl(st_coordinates(geom))), "fill", set = TRUE)
  points[!base::duplicated(geom), attrib]
}

cluster_nodes_by_cities <- function(nodes, cities, city_radius_km = 7, cluster_radius_km = 20, 
                                    weight = "passes", algo.dbscan = TRUE) {
  nodes <- copy(nodes)
  # Cluster nodes close to cities
  dmat <- st_distance(nodes, cities) < units::as_units(city_radius_km, "km")
  # fsum(dmat)
  city_nodes <- dapply(dmat, function(x) if(any(x)) which.max(x) else NA_integer_, MARGIN = 1)
  # Cluster non-city nodes
  non_city_nodes <- st_coordinates(nodes)[is.na(city_nodes), ]
  clusters <- if(algo.dbscan) dbscan::dbscan(non_city_nodes, eps = deg_m(cluster_radius_km * 1000), minPts = 1)$cluster else
    leaderCluster::leaderCluster(non_city_nodes, cluster_radius_km, distance = "haversine")$cluster_id 
  # Add geometry to city nodes
  tfm(nodes) <- cities |> qDF() |> tfm(qDF(st_coordinates(geometry))) |> tfm(geometry = NULL) |> 
    ss(city_nodes) |> rnm(X = lon, Y = lat)
  ## Combining
  set(nodes, whichNA(city_nodes), c("lon", "lat"), 
      mctl(non_city_nodes[fmode(seq_row(non_city_nodes), clusters, nodes[[weight]][is.na(city_nodes)], "fill"), ]))
  # Adding cluster variable
  city_nodes[is.na(city_nodes)] <- clusters + fmax(city_nodes)
  nodes$cluster <- city_nodes
  return(nodes)
}

contract_segments <- function(segments, nodes_clustered, attrib = "passes") {
  
  # Now Manually contracting graph
  graph_straight <- line2df(segments) |> fselect(-L1) |> 
    add_vars(.subset(segments, attrib)) |>  
    # Removing singleton edges
    fsubset(!(as.integer(fx*8) == as.integer(tx*8) & as.integer(fy*8) == as.integer(ty*8))) |>  
    # Aggregating (if applicable)
    fgroup_by(1:4) |> fsum()
  
  # Creating undirected graph (removing identical directional edges)
  ind <- fmatch(fselect(graph_straight, tx, ty, fx, fy), gv(graph_straight, 1:4), overid = 2)
  if(!allNA(ind)) graph_straight <- graph_straight[is.na(ind), ] 
  
  # Contracting graph 
  graph_straight |> 
    join(nodes_clustered |> qDF() |> 
           ftransform(qDF(st_coordinates(geometry))) |> 
           fselect(fx = X, fy = Y, lon, lat)) |> 
    ftransform(fx = lon, fy = lat, lon = NULL, lat = NULL) |> 
    join(nodes_clustered |> qDF() |> 
           ftransform(qDF(st_coordinates(geometry))) |> 
           fselect(tx = X, ty = Y, lon, lat)) |>
    ftransform(tx = lon, ty = lat, lon = NULL, lat = NULL) |> 
    fgroup_by(1:4) |> fsum() |> na_omit() |> 
    st_as_sf(coords = c("fx", "fy"), crs = 4326, sf_column_name = "from") |> qDF() |> 
    st_as_sf(coords = c("tx", "ty"), crs = 4326, sf_column_name = "to") |> qDF() |> 
    fmutate(row = seq_along(from)) |> 
    pivot(c("row", attrib), c("from", "to")) |> 
    collap(~ row, custom = list(ffirst = attrib, st_combine = c("geometry" = "value")), 
           keep.col.order = FALSE) |> 
    st_as_sf(crs = 4326) |> 
    st_cast("LINESTRING") |> 
    fmutate(row = NULL) 
}

total_MA <- function(distances, weights) {
  inv_distance <- 1 / unclass(distances) 
  diag(inv_distance) <- 0 
  market_access <- inv_distance %*% weights
  sum(market_access)
}
