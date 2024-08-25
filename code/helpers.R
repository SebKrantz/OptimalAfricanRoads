
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
  if(length(rn) && (is.character(rn) || !identical(rn, seq_along(rn)))) {
    dimnames(res$durations) <- dimnames(res$distances) <- list(rn, rn)
  }
  
  return(res)
}
