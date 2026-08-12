#####################################################################
# PIDA Road Projects: Trans-African GE Evaluation
# -------------------------------------------------------------------
# Visualises the GE-optimal allocations from
#   code/11_GE_simulation_trans_african/optimal_trans_african_networks_largest_pcities_dual_loop.jl
# at the PIDA-equivalent budget K_i = 10.089B 2015 USD, and overlays
# the subset of trans-African links that correspond to PIDA projects.
#
# Network used: 47-largest-port-cities fastest-routes graph
# (330 condensed edges, 212 nodes); see
# `code/9_major_transport_routes.R` for construction. The PIDA-overlap
# indices were computed in that same script via
#   ind <- st_nearest_feature(st_centroid(PIDA_edges),
#                             st_centroid(real_paths_simplified)) |> unique()
# and pasted here verbatim.
#
# Frictions intentionally omitted: §4 of the paper shows that border
# frictions have only minor effects on the *shape* of the GE-optimal
# allocation. OSBP effects are precisely identified in the PE section.
#
# Inputs:
#   data/transport_network/trans_african/trans_africa_network_47_largest.qs2
#   data/transport_network/trans_african/trans_africa_network_47_largest_fastest_real_edges.qs2
#   results/transport_network/GE_dual/trans_african/{nodes,edges}_results_22g_10089m_fixed_cgc[_irs_na]_sigma3.8_rho{0,2}_duality_julia.csv
#
# Outputs:
#   figures/transport_network/PIDA/GE_trans_african/
#     trans_africa_network_GE_<spec>_perc_ug_with_PIDA.pdf  (4)
#     trans_africa_network_GE_<spec>_upw_gain.pdf            (4)
#   results/transport_network/PE/PIDA_GE_trans_african_results.qs2
#####################################################################

library(fastverse)
set_collapse(mask = c("manip", "helper", "special"), nthreads = 4)
fastverse_extend(qs2, sf, units, sfnetworks, tmap, install = TRUE)
source("code/helpers/helpers.R")
fastverse_conflicts()

dir.create("figures/transport_network/PIDA/GE_trans_african",
           recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# 1. PIDA edge indices in the 47-largest fastest-routes graph
# -------------------------------------------------------------------
# Computed in code/12_PIDA_evaluation/match_PIDA_to_trans_african.R via
# coverage-based matching: for each trans-African edge T, the fraction
# of T's length falling within a 25km buffer of the union of all PIDA
# (real-geometry) edges must be >= 50% for T to be flagged PIDA.
# (This replaces the earlier centroid-only st_nearest_feature match,
# which produced many perpendicular / spurious matches.)
PIDA_ind <- c(  7L,  16L,  17L,  18L,  23L,  24L,  26L,  52L,  53L,  59L,
               60L,  69L,  70L,  71L,  72L,  73L,  74L,  86L,  88L,  89L,
               98L, 100L, 106L, 116L, 117L, 134L, 138L, 170L, 177L, 182L,
              183L, 185L, 187L, 189L, 205L, 212L, 213L, 216L, 217L, 218L,
              220L, 221L, 232L, 234L, 235L, 248L, 249L, 255L, 259L, 267L,
              268L, 279L, 285L, 293L, 294L, 309L, 322L, 330L)
# 58 distinct trans-African edges flagged as PIDA-overlapping.

# -------------------------------------------------------------------
# 2. Load network + real-paths geometry
# -------------------------------------------------------------------
network_obj <- qs2::qs_read(
  "data/transport_network/trans_african/trans_africa_network_47_largest.qs2"
) |> extract2("fastest_routes") |> extract2("network")

nodes_sf  <- network_obj |> st_as_sf("nodes")
edges_sf  <- network_obj |> st_as_sf("edges")

# Pretty geometry for maps: real routed OSRM paths simplified to 10%
edges_real <- qs2::qs_read(
  "data/transport_network/trans_african/trans_africa_network_47_largest_fastest_real_edges.qs2"
) |> rmapshaper::ms_simplify(keep = 0.1) |> st_make_valid()
stopifnot(nrow(edges_real) == nrow(edges_sf))            # row alignment
tfm(edges_real) <- atomic_elem(
  select(edges_sf, from, to, distance, duration)
)

# Bounds check on PIDA indices
stopifnot(all(PIDA_ind >= 1L & PIDA_ind <= nrow(edges_real)))
cat("PIDA indices: ", length(PIDA_ind), " of ", nrow(edges_real),
    " trans-African edges (", round(length(PIDA_ind) /
    nrow(edges_real) * 100, 1), "% of network)\n", sep = "")

# -------------------------------------------------------------------
# 3. Helpers
# -------------------------------------------------------------------
load_spec <- function(spec) {
  list(
    nodes = fread(sprintf(
      "results/transport_network/GE_dual/trans_african/nodes_results_%s.csv", spec)),
    edges = fread(sprintf(
      "results/transport_network/GE_dual/trans_african/edges_results_%s.csv", spec))
  )
}

compute_stats <- function(res, spec) {
  e <- res$edges
  n <- res$nodes
  if (grepl("rho2", spec, fixed = TRUE)) {
    # rho=2 utility is negative-reciprocal scaled
    n <- mutate(n, uj      = (uj      * (-1))^(-1),
                   uj_orig = (uj_orig * (-1))^(-1))
  }
  perc_ug <- with(e, replace_inf(pmax(
    (Ijk - Ijk_orig) / (pmax(Ijk_orig, 100) - Ijk_orig), 0), 0))

  cost_vec   <- perc_ug * e$distance * (e$total_cost / e$distance) # = perc_ug * total_cost
  dist_vec   <- perc_ug * e$distance
  in_PIDA    <- seq_len(nrow(e)) %in% PIDA_ind

  budget_total   <- sum(cost_vec, na.rm = TRUE)
  budget_pida    <- sum(cost_vec[ in_PIDA], na.rm = TRUE)
  budget_nonpida <- sum(cost_vec[!in_PIDA], na.rm = TRUE)
  km_total       <- sum(dist_vec, na.rm = TRUE)
  km_pida        <- sum(dist_vec[ in_PIDA], na.rm = TRUE)
  km_nonpida     <- sum(dist_vec[!in_PIDA], na.rm = TRUE)
  share_pida_bud <- budget_pida / budget_total
  share_pida_km  <- km_pida     / km_total

  # km / budget split by add type (FALSE/NA = upgrade, TRUE = new construction).
  # In the trans-African fastest-routes network `add` is logical with only
  # FALSE (154) and NA (176) values; treat both as upgrade.
  if ("add" %in% names(e)) {
    is_new        <- !is.na(e$add) & as.logical(e$add)
    km_upgrade    <- sum(dist_vec[!is_new], na.rm = TRUE)
    km_new        <- sum(dist_vec[ is_new], na.rm = TRUE)
    bud_upgrade   <- sum(cost_vec[!is_new], na.rm = TRUE)
    bud_new       <- sum(cost_vec[ is_new], na.rm = TRUE)
  } else {
    km_upgrade  <- km_total ; km_new  <- 0
    bud_upgrade <- budget_total ; bud_new <- 0
  }

  # Welfare gain (utility-weighted, population-weighted)
  wg <- sum(n$uj * n$Lj, na.rm = TRUE) /
        sum(n$uj_orig * n$Lj_orig, na.rm = TRUE) - 1
  # Consumption gain
  cg <- sum(n$Cj, na.rm = TRUE) /
        sum(n$Cj_orig, na.rm = TRUE) - 1

  list(
    spec           = spec,
    budget_total   = budget_total,
    budget_pida    = budget_pida,
    budget_nonpida = budget_nonpida,
    km_total       = km_total,
    km_pida        = km_pida,
    km_nonpida     = km_nonpida,
    km_upgrade     = km_upgrade,
    km_new         = km_new,
    bud_upgrade    = bud_upgrade,
    bud_new        = bud_new,
    share_pida_bud = share_pida_bud,
    share_pida_km  = share_pida_km,
    wg_perc        = wg * 100,
    cg_perc        = cg * 100
  )
}

compute_MA_gain <- function(res) {
  e <- res$edges
  e <- mutate(e,
              duration     = distance / Ijk_orig,
              duration_new = distance / Ijk)
  # Build sfnetwork via sf-on-left join
  geom_union <- select(edges_sf, from, to)
  e_sf <- join(geom_union, e, on = c("from", "to"), verbose = 0, overid = 2)
  net_e <- as_sfnetwork(e_sf, directed = FALSE)
  ind   <- ckmatch(round(select(qDF(res$nodes), lon, lat), 5),
                   mctl(round(st_coordinates(st_geometry(net_e, "nodes")), 5)))
  times     <- st_network_cost(net_e, weights = e_sf$duration     * 60)[ind, ind]
  times_new <- st_network_cost(net_e, weights = e_sf$duration_new * 60)[ind, ind]
  MA     <- total_MA(times,     res$nodes$gdp)
  MA_new <- total_MA(times_new, res$nodes$gdp)
  (MA_new / MA - 1) * 100
}

# -------------------------------------------------------------------
# 4. Map helper: GE-optimal investments with PIDA overlay
# -------------------------------------------------------------------
nice_label <- function(spec) {
  irs <- grepl("irs_na", spec, fixed = TRUE)
  rho2 <- grepl("rho2",  spec, fixed = TRUE)
  paste0(if (irs) "IRS planner" else "Standard planner",
         " (", if (rho2) "rho=2, IA" else "rho=0", ")")
}

make_perc_ug_map_PIDA <- function(res, spec, st) {
  # perc_ug per edge (intensive margin)
  e <- res$edges |> mutate(perc_ug = pmin(pmax(
    (Ijk - Ijk_orig) / (100 - Ijk_orig) * 100, 0), 100))
  # Attach perc_ug to the pretty real-edges geometry by row position
  e_real <- edges_real
  e_real$perc_ug <- e$perc_ug
  # Mark PIDA edges
  e_real$is_PIDA <- seq_len(nrow(e_real)) %in% PIDA_ind

  # Nodes: rebuild `product` factor (levels come from analyze_trans_african_results.R L52-54)
  n <- res$nodes
  largest <- with(subset(qDT(n), unclass(product) > 5L),
                  set_names(product, city_country)) |> sort() |> names()
  attr(n$product, "levels") <- c("Small City/Node", "City > 200K", "Port",
                                 "City > 2M", "Large Port-City", largest)
  class(n$product) <- "factor"
  # Compact display factor: strip "/Node" and "Large ", drop megacities to NA
  # so the legend labels them "Megacity (Own)" via label.na.
  n <- n |> mutate(
    prod2 = set_attr(product, "levels",
                     gsub("/Node|Large ", "", levels(product))),
    prod2 = droplevels(fifelse(unclass(prod2) > 5L, NA, prod2))
  )
  n_sf <- st_as_sf(n, coords = c("lon", "lat"), crs = 4326)

  # Statistics box (mirrors stats_legend in analyze_trans_african_results.R L157-160)
  # Budget & km are millions & km respectively; gains are already in %.
  stats_legend <- c(
    sprintf("Budget: $%.2fB",     st$budget_total / 1e3),
    sprintf("Upgraded: %s km (%.2fB)",
            format(round(st$km_upgrade), big.mark = ","),
            st$bud_upgrade / 1e3)
  )
  if (st$km_new > 0.5) {
    stats_legend <- c(stats_legend,
      sprintf("New: %s km (%.2fB)",
              format(round(st$km_new), big.mark = ","),
              st$bud_new / 1e3))
  }
  stats_legend <- c(stats_legend,
    sprintf("PIDA share: %.1f%%b, %.1f%%km",
            st$share_pida_bud * 100, st$share_pida_km * 100),
    sprintf("Gains: %.1f%%C, %.1f%%W, %.1f%%MA",
            st$cg_perc, st$wg_perc, st$ma_perc))

  pl <- tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
    tm_shape(e_real) +
    tm_lines(col = "perc_ug",
             col.scale  = tm_scale_continuous(ticks = seq(0, 100, 20),
                                              values = "brewer.yl_or_rd"),
             col.legend = tm_legend(expression(Delta~"%"~"UG"),
                                    position = c("left", "bottom"),
                                    stack = "h", frame = FALSE, bg.alpha = 0,
                                    height = 16, item.width = 0.5,
                                    text.size = 1.2, title.size = 1.5,
                                    title.padding = c(-0.5, 0, 0, 0)),
             lwd = 2) +
    # PIDA overlay: dotted navy
    tm_shape(subset(e_real, is_PIDA)) +
    tm_lines(col = "navy", lwd = 1.1, lty = "dotted") +
    # Populated nodes (sized by population, colored by product type)
    tm_shape(subset(n_sf, population > 0)) +
    tm_dots(size = "population",
            size.scale = tm_scale_intervals(
              breaks = c(0, 0.5e3, 2e3, Inf),
              values = c(1, 2, 3) * 0.1,
              labels = c("0 to 500", "500 to 2,000", "2,000 or more")),
            size.legend = tm_legend("Population (K)",
                                    position = tm_pos_in(0.79, 0.15),
                                    frame = FALSE, bg.alpha = 0,
                                    text.size = 1.1,
                                    item.width = 0.4, title.size = 1.1),
            fill = "prod2",
            fill.scale = tm_scale_categorical(values = "turbo",
                                              value.na = "purple3",
                                              label.na = "Megacity (Own)"),
            fill.legend = tm_legend("Product",
                                    position = c("left", "bottom"),
                                    size = 0.8, frame = FALSE, bg.alpha = 0,
                                    text.size = 1.2, title.size = 1.5,
                                    title.padding = c(0, 0, 0, 0),
                                    item.width = 1)) +
    # Small unpopulated nodes (from the graph, not res$nodes)
    tm_shape(subset(nodes_sf, population <= 0)) +
    tm_dots(size = 0.07, fill = "grey70") +
    # PIDA overlay legend - top-right corner. Use tm_pos_in() with
    # just.h/just.v to anchor the top-right corner of the legend at
    # (0.99, 0.99) of the map area, so text sits tight to the edge
    # instead of getting clipped or dropped low.
    tm_add_legend(type = "lines", labels = "PIDA-overlap link",
                  col = "navy", lty = "dotted", lwd = 2,
                  title = "Overlay",
                  position = tm_pos_in(0.99, 0.99,
                                       just.h = "right", just.v = "top"),
                  text.size = 1, title.size = 1.3,
                  frame = FALSE, bg.alpha = 0) +
    # Statistics legend
    tm_add_legend(title = "Statistics", type = "lines",
                  labels = stats_legend,
                  position = tm_pos_in(0.133, 0.23),
                  text.size = 1, title.size = 1.5,
                  item.width = 0.2, item.space = 0.2,
                  frame = FALSE, bg.alpha = 0) +
    tm_layout(frame = FALSE)

  out <- sprintf(
    "figures/transport_network/PIDA/GE_trans_african/trans_africa_network_GE_%s_perc_ug_with_PIDA.pdf",
    spec)
  tmap_save(pl, out, width = 8, height = 8.3)
  invisible(out)
}

make_welfare_map <- function(res, spec) {
  n <- res$nodes
  if (grepl("rho2", spec, fixed = TRUE)) {
    n <- mutate(n, uj      = (uj      * (-1))^(-1),
                   uj_orig = (uj_orig * (-1))^(-1))
  }
  n_sf <- n |> mutate(ugain = (uj / uj_orig - 1) * 100) |>
    st_as_sf(coords = c("lon", "lat"), crs = 4326)

  pl <- tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
    tm_shape(edges_real) +
    tm_lines(col = "grey80", lwd = 0.6) +
    tm_shape(n_sf) +
    tm_dots(col = NULL, fill = "ugain",
            fill.scale  = tm_scale_continuous(
              limits = c(-30, 30), outliers.trunc = rep(TRUE, 2),
              values = "-tableau.classic_orange_blue"),
            fill.legend = tm_legend("Welfare Gain (%)",
                                    position = c("left", "bottom"),
                                    frame = FALSE, height = 16,
                                    item.width = 0.5,
                                    text.size = 1.1, title.size = 1.4),
            size = 0.3) +
    tm_layout(frame = FALSE)
  out <- sprintf(
    "figures/transport_network/PIDA/GE_trans_african/trans_africa_network_GE_%s_upw_gain.pdf",
    spec)
  tmap_save(pl, out, width = 10, height = 10)
  invisible(out)
}

# -------------------------------------------------------------------
# 5. Main loop: 4 specs at sigma=3.8, frictionless
# -------------------------------------------------------------------
specs <- c(
  # sigma = 3.8 (Armington baseline; used in main paper figure)
  "22g_10089m_fixed_cgc_sigma3.8_rho0_duality_julia",        # standard, rho=0
  "22g_10089m_fixed_cgc_sigma3.8_rho2_duality_julia",        # standard, rho=2
  "22g_10089m_fixed_cgc_irs_na_sigma3.8_rho0_duality_julia", # IRS,      rho=0
  "22g_10089m_fixed_cgc_irs_na_sigma3.8_rho2_duality_julia", # IRS,      rho=2
  # sigma = 2.0 (lower elasticity = more differentiated goods, less substitutable)
  "22g_10089m_fixed_cgc_sigma2.0_rho0_duality_julia",        # standard, rho=0
  "22g_10089m_fixed_cgc_sigma2.0_rho2_duality_julia",        # standard, rho=2
  "22g_10089m_fixed_cgc_irs_na_sigma2.0_rho0_duality_julia", # IRS,      rho=0
  "22g_10089m_fixed_cgc_irs_na_sigma2.0_rho2_duality_julia"  # IRS,      rho=2
)

all_stats <- list()
for (spec in specs) {
  cat("\n=== ", spec, " (", nice_label(spec), ") ===\n", sep = "")
  res <- load_spec(spec)
  stopifnot(nrow(res$edges) == nrow(edges_sf))
  st  <- compute_stats(res, spec)
  st$ma_perc <- compute_MA_gain(res)
  all_stats[[spec]] <- st

  cat(sprintf("Budget spent: $%.3fB on %.0f km (target $10.089B)\n",
              st$budget_total / 1e3, st$km_total))
  cat(sprintf("PIDA-overlap share of budget: %.1f%% ($%.2fB on %.0f km)\n",
              st$share_pida_bud * 100, st$budget_pida / 1e3, st$km_pida))
  cat(sprintf("PIDA-overlap share of km: %.1f%%\n", st$share_pida_km * 100))
  cat(sprintf("Welfare gain: %.3f%% | Consumption gain: %.3f%% | MA gain: %.2f%%\n",
              st$wg_perc, st$cg_perc, st$ma_perc))

  make_perc_ug_map_PIDA(res, spec, st)
  make_welfare_map(res, spec)
}

# -------------------------------------------------------------------
# 6. Summary
# -------------------------------------------------------------------
stats_tbl <- rowbind(lapply(all_stats, as.data.frame), idcol = "spec")
print(stats_tbl)

qs2::qs_save(list(
  stats_tbl = stats_tbl,
  all_stats = all_stats,
  PIDA_ind  = PIDA_ind,
  network_n_edges = nrow(edges_sf),
  network_n_nodes = nrow(nodes_sf)
), "results/transport_network/PE/PIDA_GE_trans_african_results.qs2")

cat("\nDone. Wrote PIDA_GE_trans_african_results.qs2 and 8 figures",
    "to figures/transport_network/PIDA/GE_trans_african/\n")
