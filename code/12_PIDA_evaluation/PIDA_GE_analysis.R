#####################################################################
# PIDA Road Projects: General Equilibrium Evaluation
# -------------------------------------------------------------------
# Adapts `code/10_GE_simulation_regional/analyze_total_results.R` to
# compare the welfare-maximising GE planner's investments at the
# PIDA-equivalent budget ($10B, sigma=3.8) against the PIDA package.
# PIDA links are overlaid on each optimal-investment map in a
# contrasting colour so that the two allocations can be inspected
# side-by-side.
#
# Frictions are intentionally omitted: as established earlier in the
# paper (and verified in the GE simulations there), border frictions
# have only minimal effects on the *shape* of the GE-optimal
# allocation, so all four runs below are frictionless and the OSBP
# component is set aside (it is treated in the PE section, where its
# effect is well identified).
#
# Inputs:
#   data/transport_network/trans_africa_network_param.RData
#   data/transport_network/edges_real_simplified.qs2
#   data/PIDA/PIDA_edges_bridge.csv
#   results/transport_network/GE/total/{nodes,edges}_results_<spec>.csv
#     for spec in:
#       add_22g_10b_fixed_sigma3.8_rho0_julia        (standard, rho=0)
#       add_22g_10b_fixed_sigma3.8_rho2_julia        (standard, rho=2)
#       add_22g_10b_fixed_irs_na_sigma3.8_rho0_julia (IRS,      rho=0)
#       add_22g_10b_fixed_irs_na_sigma3.8_rho2_julia (IRS,      rho=2)
#
# Outputs:
#   figures/transport_network/PIDA/GE/trans_africa_network_GE_<spec>_perc_ug_with_PIDA.pdf  (4)
#   figures/transport_network/PIDA/GE/trans_africa_network_GE_<spec>_upw_gain.pdf            (4)
#   results/transport_network/PE/PIDA_GE_results.qs2  (per-spec stats table)
#####################################################################

library(fastverse)
set_collapse(mask = c("manip", "helper", "special"), nthreads = 4)
fastverse_extend(qs2, sf, units, sfnetworks, tmap, install = TRUE)
source("code/helpers/helpers.R")
fastverse_conflicts()

dir.create("figures/transport_network/PIDA/GE", recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# Load network + PIDA bridge
# -------------------------------------------------------------------
load("data/transport_network/trans_africa_network_param.RData")
edges_real <- qs2::qs_read("data/transport_network/edges_real_simplified.qs2") |>
  select(from, to) |> rmapshaper::ms_simplify(keep = 0.06) |> st_make_valid()

PIDA_bridge <- fread("data/PIDA/PIDA_edges_bridge.csv")  # from, to, PIDA, OSBP

# -------------------------------------------------------------------
# Helper: load a single GE result spec
# -------------------------------------------------------------------
load_spec <- function(spec) {
  list(
    nodes = fread(sprintf("results/transport_network/GE/total/nodes_results_%s.csv", spec)),
    edges = fread(sprintf("results/transport_network/GE/total/edges_results_%s.csv", spec))
  )
}

attach_geom <- function(res) {
  # Attach PIDA flags. Geometry is *not* attached here (the join would
  # strip sf class via collapse's join). For map plotting we rebuild
  # geometry from edges_real + add_links by row order; for the MA
  # computation in `compute_MA_gain` we use the sf union directly.
  res$edges <- join(res$edges, PIDA_bridge, on = c("from", "to"),
                    verbose = 0, overid = 2)
  res
}

compute_stats <- function(res, spec) {
  e <- res$edges
  n <- res$nodes
  if (grepl("rho2", spec, fixed = TRUE)) {
    n <- mutate(n, uj = (uj * (-1))^(-1), uj_orig = (uj_orig * (-1))^(-1))
  }
  perc_ug <- with(e, replace_inf(pmax((Ijk - Ijk_orig) / (pmax(Ijk_orig, 100) - Ijk_orig), 0), 0))
  cost_vec <- perc_ug * e$distance * e$cost_per_km   # millions USD spent per edge
  dist_vec <- perc_ug * e$distance                   # km worked per edge
  is_new   <- as.logical(e$add)
  is_pida  <- e$PIDA == "Yes"
  # Total budget actually spent (in millions; cost_per_km is USD/km, distance is km)
  budget_total   <- sum(cost_vec, na.rm = TRUE)
  budget_upgrade <- sum(cost_vec[!is_new], na.rm = TRUE)
  budget_new     <- sum(cost_vec[ is_new], na.rm = TRUE)
  km_upgrade     <- sum(dist_vec[!is_new], na.rm = TRUE)
  km_new         <- sum(dist_vec[ is_new], na.rm = TRUE)
  budget_pida    <- sum(cost_vec[ is_pida], na.rm = TRUE)
  budget_nonpida <- sum(cost_vec[!is_pida], na.rm = TRUE)
  km_pida        <- sum(dist_vec[ is_pida], na.rm = TRUE)
  km_nonpida     <- sum(dist_vec[!is_pida], na.rm = TRUE)
  # Share of optimal budget falling on PIDA links
  share_pida    <- budget_pida / budget_total
  # Welfare gain (utility-weighted)
  wg <- sum(n$uj * n$Lj, na.rm = TRUE) / sum(n$uj_orig * n$Lj_orig, na.rm = TRUE) - 1
  # Consumption gain
  cg <- sum(n$Cj, na.rm = TRUE) / sum(n$Cj_orig, na.rm = TRUE) - 1
  list(
    spec          = spec,
    budget_total  = budget_total,
    budget_upgrade= budget_upgrade,
    budget_new    = budget_new,
    km_upgrade    = km_upgrade,
    km_new        = km_new,
    budget_pida   = budget_pida,
    budget_nonpida= budget_nonpida,
    km_pida       = km_pida,
    km_nonpida    = km_nonpida,
    share_pida    = share_pida,
    wg_perc       = wg * 100,
    cg_perc       = cg * 100
  )
}

# -------------------------------------------------------------------
# Helper: MA gain on the optimal network (matches analyze_total_results)
# -------------------------------------------------------------------
compute_MA_gain <- function(res) {
  e <- res$edges
  e <- mutate(e, duration = iif(add == 1L, Inf, duration),
                 duration_new = distance / Ijk)
  # collapse::join preserves the sf class only when the sf object is on
  # the LEFT.  Put the geometry table first; the GE-result columns from
  # `e` are appended via a left-join (matched on from/to).
  geom_union <- rbind(select(edges, from, to),
                      select(add_links, from, to))
  e_sf <- join(geom_union, e, on = c("from", "to"), verbose = 0, overid = 2)
  net_e <- as_sfnetwork(e_sf, directed = FALSE)
  ind <- ckmatch(round(select(res$nodes, lon, lat), 5),
                 mctl(round(st_coordinates(st_geometry(net_e, "nodes")), 5)))
  times     <- st_network_cost(net_e, weights = e_sf$duration     * 60)[ind, ind]
  times_new <- st_network_cost(net_e, weights = e_sf$duration_new * 60)[ind, ind]
  MA     <- total_MA(times,     res$nodes$gdp)
  MA_new <- total_MA(times_new, res$nodes$gdp)
  (MA_new / MA - 1) * 100
}

# -------------------------------------------------------------------
# Helper: produce the optimal-investments + PIDA-overlay map
# -------------------------------------------------------------------
make_perc_ug_map <- function(res, spec, st) {
  e <- res$edges |>
    mutate(perc_ug = pmin(pmax((Ijk - Ijk_orig)/(100 - Ijk_orig)*100, 0), 100))
  # Need sf for tmap; rebuild geometry from canonical edges_real + add_links
  e_sf <- st_sf(qDT(e)[, c("from", "to", "perc_ug", "PIDA", "add"), with = FALSE],
                geometry = c(st_geometry(edges_real), st_geometry(add_links)))

  pl <- tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
    # All optimal investments, coloured by upgrade intensity
    tm_shape(e_sf) +
    tm_lines(col = "perc_ug",
             col.scale  = tm_scale_continuous(ticks = seq(0, 100, 20),
                                              values = "brewer.yl_or_rd"),
             col.legend = tm_legend(expression(Delta~"%"~"UG"),
                                    position = c("left", "bottom"),
                                    frame = FALSE, height = 16, item.width = 0.5,
                                    text.size = 1.1, title.size = 1.4),
             lwd = 1.5) +
    # PIDA links overlay - dashed navy
    tm_shape(subset(e_sf, PIDA == "Yes")) +
    tm_lines(col = "navy", lwd = 0.7, lty = "dashed") +
    tm_add_legend(type = "lines", labels = "PIDA link",
                  col = "navy", lty = "dashed", lwd = 1.5,
                  title = "Overlay",
                  position = c("right", "bottom"),
                  text.size = 1.0, title.size = 1.3, frame = FALSE) +
    tm_layout(frame = FALSE)

  out <- sprintf("figures/transport_network/PIDA/GE/trans_africa_network_GE_%s_perc_ug_with_PIDA.pdf", spec)
  tmap_save(pl, out, width = 10, height = 10)
  invisible(out)
}

# -------------------------------------------------------------------
# Helper: local welfare gain (per worker) map
# -------------------------------------------------------------------
make_welfare_map <- function(res, spec) {
  n <- res$nodes
  if (grepl("rho2", spec, fixed = TRUE)) {
    n <- mutate(n, uj = (uj * (-1))^(-1), uj_orig = (uj_orig * (-1))^(-1))
  }
  n_sf <- n |> mutate(ugain = (uj / uj_orig - 1) * 100) |>
    st_as_sf(coords = c("lon", "lat"), crs = 4326)

  pl <- tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
    tm_shape(n_sf) +
    tm_dots(col = NULL,
            fill = "ugain",
            fill.scale  = tm_scale_continuous(limits = c(-30, 30),
                                              outliers.trunc = rep(TRUE, 2),
                                              values = "-tableau.classic_orange_blue"),
            fill.legend = tm_legend("Welfare Gain (%)",
                                    position = c("left", "bottom"), frame = FALSE,
                                    height = 16, item.width = 0.5,
                                    text.size = 1.1, title.size = 1.4),
            size = 0.2) +
    tm_layout(frame = FALSE)
  out <- sprintf("figures/transport_network/PIDA/GE/trans_africa_network_GE_%s_upw_gain.pdf", spec)
  tmap_save(pl, out, width = 10, height = 10)
  invisible(out)
}


# -------------------------------------------------------------------
# Main loop over the 4 GE configurations
# -------------------------------------------------------------------
specs <- c(
  "add_22g_10b_fixed_sigma3.8_rho0_julia",          # standard, rho=0
  "add_22g_10b_fixed_sigma3.8_rho2_julia",          # standard, rho=2
  "add_22g_10b_fixed_irs_na_sigma3.8_rho0_julia",   # IRS,      rho=0
  "add_22g_10b_fixed_irs_na_sigma3.8_rho2_julia"    # IRS,      rho=2
)

all_stats <- list()
for (spec in specs) {
  cat("\n=== ", spec, " ===\n", sep = "")
  res <- attach_geom(load_spec(spec))
  st  <- compute_stats(res, spec)
  st$ma_perc <- compute_MA_gain(res)
  all_stats[[spec]] <- st

  cat(sprintf("Budget spent: $%.2fB  (upgrade $%.2fB / new $%.2fB)\n",
              st$budget_total/1e3, st$budget_upgrade/1e3, st$budget_new/1e3))
  cat(sprintf("Km worked:    %.0f km upgrade + %.0f km new\n", st$km_upgrade, st$km_new))
  cat(sprintf("Share of budget falling on PIDA links: %.2f%%  ($%.2fB on %.0fkm of PIDA)\n",
              st$share_pida * 100, st$budget_pida/1e3, st$km_pida))
  cat(sprintf("Welfare gain: %.3f%% | Consumption gain: %.3f%% | MA gain: %.2f%%\n",
              st$wg_perc, st$cg_perc, st$ma_perc))

  make_perc_ug_map(res, spec, st)
  make_welfare_map(res, spec)
}

# -------------------------------------------------------------------
# Compile summary table for the paper
# -------------------------------------------------------------------
stats_tbl <- rowbind(lapply(all_stats, as.data.frame), idcol = "spec")
print(stats_tbl)

qs2::qs_save(list(
  stats_tbl = stats_tbl,
  all_stats = all_stats
), "results/transport_network/PE/PIDA_GE_results.qs2")

cat("\nDone. Wrote PIDA_GE_results.qs2 and 8 figures to figures/transport_network/PIDA/GE/\n")
