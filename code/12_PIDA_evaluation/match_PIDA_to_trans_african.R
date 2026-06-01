#####################################################################
# Improved PIDA -> trans-African network matching
# -------------------------------------------------------------------
# The original matching (in `code/9_major_transport_routes.R`, section
# "Match to PIDA") used:
#   ind <- st_nearest_feature(st_centroid(PIDA_edges),
#                             st_centroid(real_paths_simplified)) |> unique()
# which considers only centroid proximity. As a result, some matched
# trans-African edges are PERPENDICULAR to the actual PIDA road, and
# many PIDA corridors that are clearly represented by a parallel
# trans-African link are missed.
#
# This script replaces that with a coverage-based match:
#   For each trans-African edge T_i, compute the fraction of T_i's
#   length that falls within a buffer around the union of all PIDA
#   edges (real road geometry). A trans-African edge is flagged
#   PIDA-matched when:
#     coverage_fraction >= COVER_THRESHOLD
#
# Output: an updated `ind` vector, saved to
#   data/PIDA/PIDA_trans_african_indices.csv
# plus diagnostic comparison maps for visual iteration.
#####################################################################

library(fastverse)
set_collapse(mask = c("manip", "helper", "special"), nthreads = 4)
fastverse_extend(qs2, sf, sfnetworks, tmap, install = TRUE)
source("code/helpers/helpers.R")
fastverse_conflicts()

dir.create("figures/transport_network/PIDA/match",
           recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# 1. Inputs
# -------------------------------------------------------------------
load("data/transport_network/trans_africa_network_param.RData")
edges_real_full <- qs2::qs_read("data/transport_network/edges_real_simplified.qs2")
PIDA_bridge     <- fread("data/PIDA/PIDA_edges_bridge.csv")

# Trans-African 47-largest network (330 edges) + simplified OSRM real paths
real_paths_simplified <- qs2::qs_read(
  "data/transport_network/trans_african/trans_africa_network_47_largest_fastest_real_edges.qs2"
) |> rmapshaper::ms_simplify(keep = 0.1) |> st_make_valid()
network_obj <- qs2::qs_read(
  "data/transport_network/trans_african/trans_africa_network_47_largest.qs2"
) |> extract2("fastest_routes") |> extract2("network")

# -------------------------------------------------------------------
# 2. Build PIDA edge geometries on the full grid
#    (use the actual real-road geometry, not the straight-line
#    edges_param geometry, so the buffer matches real corridor shape)
# -------------------------------------------------------------------
add_renamed <- rename(add_links_param,
                      duration_65kmh   = duration,
                      duration_100kmh  = duration_imp,
                      total_time_65kmh = total_time,
                      total_time_100kmh= total_time_imp)
drop_cols <- intersect(c("id", "total_dist"), names(add_renamed))
if (length(drop_cols))
  add_renamed <- add_renamed[, setdiff(names(add_renamed), drop_cols)]
edges_full_param <- rowbind(edges_param, add_renamed, fill = TRUE)

# `edges_real_full` aligns to `edges_param`; for add_links use straight-
# line geometry (acceptable since add_links are mostly short candidate
# new links). Use union for the buffer.
PIDA_full_idx   <- which(join(edges_full_param, PIDA_bridge,
                              on = c("from", "to"),
                              verbose = 0, overid = 2)$PIDA == "Yes")
cat("PIDA edges on full grid:", length(PIDA_full_idx), "\n")

# For existing edges use the simplified real geometry; for new links
# use the straight-line geometry
n_existing      <- nrow(edges_param)
existing_PIDA   <- PIDA_full_idx[PIDA_full_idx <= n_existing]
new_PIDA        <- PIDA_full_idx[PIDA_full_idx >  n_existing] - n_existing

PIDA_geom_existing <- st_geometry(edges_real_full)[existing_PIDA]
PIDA_geom_new      <- if (length(new_PIDA))
  st_geometry(add_links_param)[new_PIDA] else NULL
PIDA_geom <- c(PIDA_geom_existing, PIDA_geom_new) |>
  st_sfc(crs = st_crs(real_paths_simplified))

cat("PIDA real-geometry edges: ", length(PIDA_geom), "\n",
    "(", length(PIDA_geom_existing), " existing + ",
    length(PIDA_geom_new), " new)\n", sep = "")

# -------------------------------------------------------------------
# 3. Coverage-based matching
# -------------------------------------------------------------------
# Switch to a meter-based projected CRS for accurate buffer/length
# operations.  EPSG:3395 (World Mercator) is fine at continental scale
# for this purpose.
proj_crs <- 3395
PIDA_geom_m <- st_transform(PIDA_geom, proj_crs)
# IMPORTANT: extract geometry as sfc; seq_along(sf) returns col-count not row-count
ta_geom_m   <- st_transform(st_geometry(real_paths_simplified), proj_crs)
stopifnot(length(ta_geom_m) == nrow(real_paths_simplified))

# Make a single buffered polygon around all PIDA edges
make_match <- function(buffer_km, cover_thresh) {
  cat(sprintf("\n--- buffer=%dkm, threshold=%.2f ---\n", buffer_km, cover_thresh))
  PIDA_buf <- st_buffer(PIDA_geom_m, dist = buffer_km * 1000) |>
              st_union() |> st_make_valid()

  # For each trans-African edge, compute fraction of length inside buffer
  ta_len <- as.numeric(st_length(ta_geom_m))
  coverage <- sapply(seq_along(ta_geom_m), function(i) {
    inside <- suppressWarnings(st_intersection(ta_geom_m[i], PIDA_buf))
    if (length(inside) == 0) return(0)
    as.numeric(sum(st_length(inside))) / ta_len[i]
  })
  matched <- which(coverage >= cover_thresh)
  cat(sprintf("Matched %d / %d trans-African edges (%.1f%%)\n",
              length(matched), length(ta_geom_m),
              length(matched) / length(ta_geom_m) * 100))
  list(matched = matched, coverage = coverage,
       buffer_km = buffer_km, cover_thresh = cover_thresh)
}

# -------------------------------------------------------------------
# 4. Try several parameter combos, save diagnostic maps
# -------------------------------------------------------------------
ind_original <- c(158L, 232L, 250L, 171L,  88L, 248L, 169L, 170L, 182L, 249L,
                  177L, 234L, 100L,   7L, 235L, 236L, 183L, 259L, 205L,  16L,
                   17L,  18L, 187L, 188L,  24L,  26L, 285L, 211L, 284L, 238L,
                  279L, 116L,  98L, 189L, 216L,  27L,  37L, 212L, 255L, 192L,
                  286L, 217L, 218L, 213L, 138L, 301L, 200L, 133L,  52L,  54L,
                  308L, 294L, 322L, 293L,  53L, 221L, 310L,  71L,  70L, 268L,
                  143L, 309L, 267L, 270L, 222L, 109L,  61L,  74L, 223L, 220L,
                  244L, 135L,  86L, 330L,  72L, 303L,  73L, 296L,  76L, 136L,
                   77L, 247L,  60L)

plot_match <- function(matched, PIDA_geom, real_paths_simplified, fname, label) {
  pl <- tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
    # All trans-African edges in grey
    tm_shape(real_paths_simplified) +
    tm_lines(col = "grey80", lwd = 0.8) +
    # PIDA real geometry in orange (semi-transparent)
    tm_shape(PIDA_geom) +
    tm_lines(col = "darkorange", lwd = 1.8) +
    # Matched trans-African edges in navy (on top)
    tm_shape(real_paths_simplified[matched, ]) +
    tm_lines(col = "navy", lwd = 1.6) +
    tm_add_legend(type = "lines",
                  labels = c("Trans-African (all)",
                             "PIDA real geometry",
                             paste0("Matched trans-African (", label, ")")),
                  col = c("grey80", "darkorange", "navy"),
                  lwd = c(0.8, 1.8, 1.6),
                  title = "Layers",
                  position = c("left", "bottom"),
                  frame = FALSE, text.size = 1.0, title.size = 1.3) +
    tm_layout(frame = FALSE)
  tmap_save(pl, fname, width = 10, height = 10)
}

# Run a grid
combos <- list(
  list(buf = 15, thr = 0.40),
  list(buf = 15, thr = 0.50),
  list(buf = 20, thr = 0.40),
  list(buf = 20, thr = 0.50),
  list(buf = 25, thr = 0.40),
  list(buf = 25, thr = 0.50),
  list(buf = 30, thr = 0.50)
)
results <- list()
for (c_ in combos) {
  k <- sprintf("buf%d_thr%.2f", c_$buf, c_$thr)
  r <- make_match(c_$buf, c_$thr)
  results[[k]] <- r
  plot_match(r$matched, PIDA_geom, real_paths_simplified,
             sprintf("figures/transport_network/PIDA/match/match_%s.pdf", k),
             label = sprintf("buf=%dkm, cov>=%.0f%%", c_$buf, c_$thr * 100))
}

# Also plot the ORIGINAL (centroid-nearest) match for comparison
plot_match(ind_original, PIDA_geom, real_paths_simplified,
           "figures/transport_network/PIDA/match/match_original_centroid.pdf",
           label = "centroid st_nearest_feature (original)")

# -------------------------------------------------------------------
# 5. Summary table
# -------------------------------------------------------------------
cat("\n=== Summary ===\n")
cat(sprintf("Original centroid match: %d edges\n", length(ind_original)))
for (k in names(results)) {
  r <- results[[k]]
  overlap_with_orig <- length(intersect(r$matched, ind_original))
  cat(sprintf("%-15s  matched=%3d  overlap-with-original=%3d (%.0f%%)\n",
              k, length(r$matched), overlap_with_orig,
              overlap_with_orig / length(ind_original) * 100))
}

# Save all combos for later inspection
qs2::qs_save(list(combos_results = results, ind_original = ind_original,
                  PIDA_geom = PIDA_geom),
             "data/PIDA/PIDA_trans_african_match_candidates.qs2")

cat("\nDiagnostic maps written to figures/transport_network/PIDA/match/\n")
cat("Inspect them visually and pick the best (buffer, threshold) combo,\n")
cat("then save the chosen ind to data/PIDA/PIDA_trans_african_indices.csv.\n")
