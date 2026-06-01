#####################################################################
# PIDA Road Projects: Partial Equilibrium Evaluation
# -------------------------------------------------------------------
# This script jointly evaluates the PIDA package of road upgrades and
# new links + One-Stop Border Posts (OSBPs) in PE terms, mirroring the
# CEMAC analogue in `OptimalCEMACRoads/code/2_PE_analysis.R` L1944-2262
# ("Evaluation Regional Road Projects").
#
# Inputs (read-only):
#   - data/transport_network/trans_africa_network_param.RData
#         (provides edges_param, add_links_param, nodes_param, net_param,
#          dist_ttime_mats, cities_ports — parameterised network state
#          produced at the end of 8_PE_analysis.R)
#   - data/PIDA/PIDA_edges_bridge.csv      (from PIDA_edges_osbp_export.R)
#   - data/PIDA/PIDA_OSBP_border_merged.qs2 (OSBP point geometries)
#   - results/transport_network/PE/PE_results.qs2 (cached MA columns)
#   - data/QSE/model_border_time_mat_transit.csv (frictions matrix)
#   - data/transport_network/edges_real_simplified.qs2 (pretty geometry)
#
# Outputs:
#   - figures/transport_network/PIDA/*.pdf  (8 maps)
#   - results/transport_network/PE/PIDA_PE_results.qs2
#   - LaTeX-ready PDV table printed to console
#
# Note: uses `qs2::qs_read()` / `qs2::qs_save()` throughout. The `qs`
# package is deprecated.
#####################################################################

library(fastverse)
set_collapse(mask = c("manip", "helper", "special"), nthreads = 4)
fastverse_extend(qs2, sf, units, sfnetworks, tmap, install = TRUE)
source("code/helpers/helpers.R")
fastverse_conflicts()

dir.create("figures/transport_network/PIDA", recursive = TRUE, showWarnings = FALSE)
dir.create("results/transport_network/PE", recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# Step 2 - Load parameterised network + cached PE results + PIDA bridge
# -------------------------------------------------------------------

# Parameterised network state from 8_PE_analysis.R: carries network
# topology (`net_param`) and `nodes_param`, but NOT the per-link MA
# columns (which live only in the PE_results cache below).
load("data/transport_network/trans_africa_network_param.RData")
nodes <- nodes_param
net   <- net_param

# Cached per-link PE results (qs2 format).  `pe_cache$edges` and
# `pe_cache$add_links` are sf objects carrying every MA_* column
# produced by `8_PE_analysis.R` (in the same row order as `net_param`).
pe_cache  <- qs2::qs_read("results/transport_network/PE/PE_results.qs2")
edges     <- pe_cache$edges       # sf, 2344 rows
add_links <- pe_cache$add_links   # sf, 481 rows
stopifnot(nrow(edges) + nrow(add_links) == nrow(fread("data/PIDA/PIDA_edges_bridge.csv")))

# Pretty edge geometry for maps
edges_real <- qs2::qs_read("data/transport_network/edges_real_simplified.qs2") |>
  rmapshaper::ms_simplify(keep = 0.06) |> st_make_valid()

# PIDA flags
PIDA_bridge   <- fread("data/PIDA/PIDA_edges_bridge.csv")  # from, to, PIDA, OSBP
PIDA_OSBP_pts <- qs2::qs_read("data/PIDA/PIDA_OSBP_border_merged.qs2")

# Frictions matrix (transit cumulative across countries)
border_time_transit <- fread("data/QSE/model_border_time_mat_transit.csv") |> qM(1)
btt_nodes <- border_time_transit[nodes$iso3c, nodes$iso3c]

# Baseline MAs (recompute from the parameterised network)
# Equivalent to 8_PE_analysis.R L763 / L993 (optimising-agents).
times    <- st_network_cost(net, weights = edges$duration)
times_bt <- st_network_cost(net, weights = edges$total_time)
nodes_coord <- mctl(st_coordinates(nodes))
ind <- ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net, "nodes"))))
times    <- times[ind, ind]
times_bt <- times_bt[ind, ind]

MA        <- total_MA(times,    nodes$gdp)   # frictionless baseline
MA_bt_opt <- total_MA(times_bt, nodes$gdp)   # frictions baseline (opt-agents)
cat(sprintf("Baseline MA = %.3e ; MA_bt_opt = %.3e ; FR/NoFR = %.3f\n",
            MA, MA_bt_opt, MA_bt_opt / MA))


# -------------------------------------------------------------------
# Step 3 - Build edges_all (existing + new links, PIDA flag joined)
# -------------------------------------------------------------------
# Mirrors `edges_all_param` in PIDA_edges_osbp_export.R but uses the
# parameterised in-memory objects so every MA_* column comes along.

# Defensive name alignment: ensure add_links has the schema used below
# (duration / duration_imp / total_time / total_time_imp).
# In some saves add_links keeps `duration_100kmh` etc. instead.
.rn_if <- function(df, from, to) {
  if (from %in% names(df) && !(to %in% names(df))) {
    setnames(df, from, to)
  }
  df
}

# Capture geometry BEFORE flattening to data.table (qDT drops sf class).
geom_edges     <- st_geometry(edges_real)         # simplified existing-edge geometry
geom_add_links <- st_geometry(add_links)          # new-link geometry
stopifnot(length(geom_edges) == nrow(edges))
stopifnot(length(geom_add_links) == nrow(add_links))

add_links_dt <- qDT(add_links)
.rn_if(add_links_dt, "duration_100kmh",  "duration_imp")
.rn_if(add_links_dt, "duration_65kmh",   "duration")
.rn_if(add_links_dt, "total_time_100kmh","total_time_imp")
.rn_if(add_links_dt, "total_time_65kmh", "total_time")
# add_links uses cost_km_adj (Algeria-Morocco 1/3 cost adjust); use as ug_cost_km proxy
if (!("ug_cost_km" %in% names(add_links_dt))) add_links_dt$ug_cost_km <- add_links_dt$cost_km_adj

edges_all <- rowbind(
  mutate(qDT(edges),    type = "existing"),
  mutate(add_links_dt,  type = "new"),
  fill = TRUE
)
edges_all <- join(edges_all, PIDA_bridge, on = c("from", "to"), verbose = 0, overid = 2)
stopifnot(!anyNA(edges_all$PIDA), !anyNA(edges_all$OSBP))

cat(sprintf("PIDA flag set on %d of %d edges (%d existing + %d new).\n",
            sum(edges_all$PIDA == "Yes"), nrow(edges_all),
            sum(edges_all$PIDA == "Yes" & edges_all$type == "existing"),
            sum(edges_all$PIDA == "Yes" & edges_all$type == "new")))


# -------------------------------------------------------------------
# Step 4 - Cost sanity-check (~$10.1B)
# -------------------------------------------------------------------
cost_per_edge <- with(edges_all, pfirst(ug_cost_km, cost_km) * distance / 1000)
cost_PIDA <- sum(cost_per_edge[edges_all$PIDA == "Yes"], na.rm = TRUE)
cat(sprintf("\ncost_PIDA = %.3e USD'15  (expected ~1.009e10)\n", cost_PIDA))
stopifnot(abs(cost_PIDA - 1.009e10) / 1.009e10 < 0.05)

cat("\nPIDA cost breakdown by type:\n")
fsum(cost_per_edge, edges_all[, .(PIDA, type)]) |> print()
cat("\nPIDA km breakdown by type:\n")
fsum(edges_all$distance/1e3, edges_all[, .(PIDA, type)]) |> print()


# -------------------------------------------------------------------
# Step 5 - PIDA per-link gains: maps (use cached MA_* columns)
# -------------------------------------------------------------------
# Mirrors CEMAC `2_PE_analysis.R` L2108-2173. We reuse the per-link
# columns already computed in 8_PE_analysis.R: for existing edges these
# are `MA_100_min_speed*` and `MA_gain_pusd*`; for new links they are
# `MA_per_link_100kmh*` and `MA_gain_100kmh_pusd*`. To facilitate a
# unified map, harmonise the names into a single `all_cb_ratios` table.

# Existing: MA_*_perc + MA_gain_pusd*
existing_cols <- intersect(c("MA_100_min_speed_perc", "MA_100_min_speed_bt_perc",
                             "MA_100_min_speed_bt_opt_perc",
                             "MA_gain_pusd", "MA_gain_pusd_bt", "MA_gain_pusd_bt_opt"),
                           names(edges))
new_cols      <- intersect(c("MA_per_link_100kmh_perc", "MA_per_link_100kmh_bt_perc",
                             "MA_per_link_100kmh_bt_opt_perc",
                             "MA_gain_100kmh_pusd", "MA_gain_100kmh_pusd_bt",
                             "MA_gain_100kmh_pusd_bt_opt"),
                           names(add_links_dt))
# Harmonised names (existing schema)
.map <- c("MA_per_link_100kmh_perc"        = "MA_100_min_speed_perc",
          "MA_per_link_100kmh_bt_perc"     = "MA_100_min_speed_bt_perc",
          "MA_per_link_100kmh_bt_opt_perc" = "MA_100_min_speed_bt_opt_perc",
          "MA_gain_100kmh_pusd"            = "MA_gain_pusd",
          "MA_gain_100kmh_pusd_bt"         = "MA_gain_pusd_bt",
          "MA_gain_100kmh_pusd_bt_opt"     = "MA_gain_pusd_bt_opt")

all_cb_ratios <- rowbind(
  existing = qDT(edges)[, c("from", "to", existing_cols), with = FALSE],
  new      = setnames(add_links_dt[, c("from", "to", new_cols), with = FALSE],
                       new_cols, unname(.map[new_cols])),
  idcol = "type",
  fill = TRUE
)
all_cb_ratios <- join(all_cb_ratios, PIDA_bridge, on = c("from", "to"), verbose = 0, overid = 2)
all_cb_ratios$cost_total <- cost_per_edge
all_cb_ratios$distance   <- edges_all$distance
# Geometry for plotting (existing uses edges_real, new uses add_links geom)
all_cb_ratios_sf <- st_sf(all_cb_ratios,
  geometry = c(geom_edges, geom_add_links))

# Consensus pusd metric for ranking (used in Step 7 too)
all_cb_ratios$MA_gain_pusd_cons <- with(all_cb_ratios,
  pmean(MA_gain_pusd, MA_gain_pusd_bt, MA_gain_pusd_bt_opt))

# Six maps (3 percent + 3 $/min/$, NoFR / FR-opt / Ratio)
make_pida_map <- function(col, breaks, legend_label, file_suffix) {
  pl <- tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
    tm_shape(all_cb_ratios_sf) + tm_lines(lwd = 1, col = "grey80") +
    tm_shape(subset(all_cb_ratios_sf, PIDA == "Yes")) +
    tm_lines(col = col,
             col.scale  = tm_scale_intervals(values = "turbo", breaks = breaks),
             col.legend = tm_legend(legend_label, position = c("left", "bottom"),
                                    frame = FALSE, text.size = 1.3, title.size = 1.6),
             lwd = 2) +
    tm_shape(subset(nodes, population > 0))  + tm_dots(size = 0.1) +
    tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
    tm_layout(frame = FALSE)
  tmap_save(pl, sprintf("figures/transport_network/PIDA/trans_africa_network_PIDA_%s.pdf", file_suffix),
            width = 10, height = 10)
  invisible(pl)
}

# Δ%-MA breaks (continent-scale gains per link are small)
ma_pct_breaks  <- c(0, 0.005, 0.01, 0.025, 0.1, 0.25, Inf)
# $/min/$ breaks copied from CEMAC L2150
ma_usd_breaks  <- c(0, 0.1, 0.2, 0.5, 1, 2, 5, Inf)

make_pida_map("MA_100_min_speed_perc",        ma_pct_breaks,
              expression(Delta~"%"~"MA [GDP/min]"), "MA_100_min_speed_perc")
make_pida_map("MA_100_min_speed_bt_opt_perc", ma_pct_breaks,
              expression(Delta~"%"~"MA [GDP/min]"), "MA_100_min_speed_bt_opt_perc")
# Ratio of friction-to-frictionless gain
all_cb_ratios_sf <- mutate(all_cb_ratios_sf,
  MA_100_min_speed_bt_opt_ratio = replace_na(replace_outliers(
    perch_to_diff(MA_100_min_speed_bt_opt_perc + MA, MA_100_min_speed_bt_opt_perc) /
    perch_to_diff(MA_100_min_speed_perc        + MA, MA_100_min_speed_perc),
    c(0, 1.5), "clip"), 0))
make_pida_map("MA_100_min_speed_bt_opt_ratio", seq(0, 1.5, 0.25),
              expression(Delta~"MA Ratio (FR/NoFR)"), "MA_100_min_speed_bt_opt_ratio")
make_pida_map("MA_gain_pusd",        ma_usd_breaks,
              expression(Delta~"MA/USD"), "MA_gain_pusd")
make_pida_map("MA_gain_pusd_bt_opt", ma_usd_breaks,
              expression(Delta~"MA/USD"), "MA_gain_pusd_bt_opt")
all_cb_ratios_sf <- mutate(all_cb_ratios_sf,
  MA_gain_pusd_bt_opt_ratio = replace_na(replace_outliers(
    MA_gain_pusd_bt_opt / MA_gain_pusd, c(0, 1.5), "clip"), 0))
make_pida_map("MA_gain_pusd_bt_opt_ratio", seq(0, 1.5, 0.25),
              expression(Delta~"MA Ratio (FR/NoFR)"), "MA_gain_pusd_bt_opt_ratio")

cat("\nWrote 6 PIDA per-link maps to figures/transport_network/PIDA/\n")


# -------------------------------------------------------------------
# Step 6 - Joint PIDA package MA gain (3 scenarios)
# -------------------------------------------------------------------
# Build sfnetwork with: all existing edges + PIDA-new add_links.
# Replace duration on PIDA-Yes edges with the 100km/h equivalent.
# Scenarios:
#   A. Frictionless        weight = duration
#   B. Frictions (opt)     weight = total_time (= duration + border_time)
#   C. Frictions + OSBP-50 weight = duration + 0.5*border_time on OSBP=Yes, else total_time

is_PIDA <- edges_all$PIDA == "Yes"
is_OSBP <- edges_all$OSBP == "Yes"

# Edge set for the PIDA-evaluated network:
#   - all existing edges
#   - only PIDA=Yes new add_links (other add_links are not in the network)
keep_for_net <- edges_all$type == "existing" | (edges_all$type == "new" & is_PIDA)

# Duration weights per scenario for kept edges
dur_A <- with(edges_all, fifelse(is_PIDA, duration_imp, duration))
dur_B <- with(edges_all, fifelse(is_PIDA, total_time_imp, total_time))
# Scenario C: OSBP edges get duration + 0.5 * border_time (whether PIDA-Yes or not,
# because halving the friction at an OSBP affects every traveller through it -
# but the OSBP flag was set only on candidate cross-border edges; non-OSBP border
# edges keep full friction).
osbp_dur <- with(edges_all, fifelse(is_PIDA, duration_imp, duration) + 0.5 * border_time)
dur_C    <- with(edges_all, fifelse(is_OSBP, osbp_dur, fifelse(is_PIDA, total_time_imp, total_time)))

# Use the canonical sfnetwork-building pattern from 8_PE_analysis.R L909:
# attach a `weight` column to the original sf objects, then `as_sfnetwork`
# the row-bound union. This guarantees CRS/coordinate alignment with the
# pre-built `nodes` (and `net_param`).
n_edges <- nrow(edges)
n_add   <- nrow(add_links)
stopifnot(length(is_PIDA) == n_edges + n_add)

joint_MA <- function(weights, label, base) {
  # Split weights back into existing-edge and new-link halves
  w_edges <- weights[seq_len(n_edges)]
  w_new   <- weights[n_edges + seq_len(n_add)]
  # Build network: full existing edge set + PIDA-Yes new links only
  edges_tmp <- edges
  edges_tmp$weight <- w_edges
  add_subset <- add_links[is_PIDA[n_edges + seq_len(n_add)], ]
  add_subset$weight <- w_new[is_PIDA[n_edges + seq_len(n_add)]]
  net_j <- as_sfnetwork(rbind(select(edges_tmp, weight),
                              select(add_subset, weight)),
                        directed = FALSE)
  ind   <- ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net_j, "nodes"))))
  times_j <- st_network_cost(net_j, weights = "weight")[ind, ind]
  MA_j <- total_MA(times_j, nodes$gdp)
  delta <- MA_j - base
  cat(sprintf("Scenario %s: MA = %.3e ; %%-gain vs base = %.3f%% ; $/min/$ = %.3f\n",
              label, MA_j, (MA_j / base - 1) * 100, delta / cost_PIDA))
  list(MA = MA_j, delta = delta, perc_gain = (MA_j / base - 1) * 100,
       usd_min_usd = delta / cost_PIDA)
}

cat("\n=== Joint PIDA package MA gains ===\n")
res_PIDA <- list(
  NoFR    = joint_MA(dur_A, "NoFR",    MA),
  FR      = joint_MA(dur_B, "FR",      MA_bt_opt),
  OSBP_50 = joint_MA(dur_C, "OSBP_50", MA_bt_opt)
)

# Sum-of-marginal (for comparison with joint, CEMAC framing)
PIDA_marg_NoFR <- sum(all_cb_ratios$MA_100_min_speed_perc[all_cb_ratios$PIDA == "Yes"],
                      na.rm = TRUE)
PIDA_marg_FR   <- sum(all_cb_ratios$MA_100_min_speed_bt_opt_perc[all_cb_ratios$PIDA == "Yes"],
                      na.rm = TRUE)
cat(sprintf("\nSum of marginal %%-MA gains: NoFR = %.3f%% ; FR = %.3f%%\n",
            PIDA_marg_NoFR, PIDA_marg_FR))
cat(sprintf("Joint > Sum?  NoFR: %s ; FR: %s\n",
            res_PIDA$NoFR$perc_gain > PIDA_marg_NoFR,
            res_PIDA$FR$perc_gain   > PIDA_marg_FR))


# -------------------------------------------------------------------
# Step 7 - Top-N consensus benchmark at ~$10 B
# -------------------------------------------------------------------
# Sort all candidate links by `MA_gain_pusd_cons` desc, accumulate cost
# until we hit cost_PIDA. Then compute joint MA on this optimal package.

ord <- order(-all_cb_ratios$MA_gain_pusd_cons,
             na.last = TRUE)
cum_cost <- cumsum(all_cb_ratios$cost_total[ord])
N_top    <- which.min(abs(cum_cost - cost_PIDA))
top_idx  <- ord[seq_len(N_top)]
all_cb_ratios$top10B <- FALSE
all_cb_ratios$top10B[top_idx] <- TRUE
cost_Top10B <- sum(all_cb_ratios$cost_total[top_idx], na.rm = TRUE)
cat(sprintf("\nTop-N consensus benchmark: N = %d links, cost = %.3e (PIDA = %.3e)\n",
            N_top, cost_Top10B, cost_PIDA))

# Joint MA for the Top-N package (same canonical pattern as PIDA above).
is_Top10B <- all_cb_ratios$top10B
dur_A_top <- with(edges_all, fifelse(is_Top10B, duration_imp, duration))
dur_B_top <- with(edges_all, fifelse(is_Top10B, total_time_imp, total_time))

joint_MA_top <- function(weights, base, label) {
  w_edges <- weights[seq_len(n_edges)]
  w_new   <- weights[n_edges + seq_len(n_add)]
  edges_tmp <- edges
  edges_tmp$weight <- w_edges
  add_subset <- add_links[is_Top10B[n_edges + seq_len(n_add)], ]
  add_subset$weight <- w_new[is_Top10B[n_edges + seq_len(n_add)]]
  net_j <- as_sfnetwork(rbind(select(edges_tmp, weight),
                              select(add_subset, weight)),
                        directed = FALSE)
  ind   <- ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net_j, "nodes"))))
  times_j <- st_network_cost(net_j, weights = "weight")[ind, ind]
  MA_j <- total_MA(times_j, nodes$gdp)
  delta <- MA_j - base
  cat(sprintf("Top-N %s: MA = %.3e ; %%-gain = %.3f%% ; $/min/$ = %.3f\n",
              label, MA_j, (MA_j / base - 1) * 100, delta / cost_Top10B))
  list(MA = MA_j, delta = delta, perc_gain = (MA_j / base - 1) * 100,
       usd_min_usd = delta / cost_Top10B)
}

cat("\n=== Joint Top-N (~$10B optimal) MA gains ===\n")
res_Top10B <- list(
  NoFR = joint_MA_top(dur_A_top, MA,        "NoFR"),
  FR   = joint_MA_top(dur_B_top, MA_bt_opt, "FR")
)


# -------------------------------------------------------------------
# Step 8 - Side-by-side comparison map (PIDA vs Top-N)
# -------------------------------------------------------------------
# Propagate the top10B flag from the data.table to the sf companion
all_cb_ratios_sf$top10B   <- all_cb_ratios$top10B
all_cb_ratios_sf$which_pkg <- with(all_cb_ratios_sf,
  fifelse(PIDA == "Yes" & top10B, "Both",
  fifelse(PIDA == "Yes",           "PIDA only",
  fifelse(top10B,                  "Top-N only", "Neither"))))

pl_cmp <- tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(subset(all_cb_ratios_sf, which_pkg == "Neither")) +
  tm_lines(lwd = 1, col = "grey85") +
  tm_shape(subset(all_cb_ratios_sf, which_pkg != "Neither")) +
  tm_lines(col = "which_pkg",
           col.scale  = tm_scale_categorical(
             values = c("Both" = "purple3", "PIDA only" = "orange", "Top-N only" = "steelblue")),
           col.legend = tm_legend("Package", position = c("left", "bottom"),
                                  frame = FALSE, text.size = 1.3, title.size = 1.6),
           lwd = 2) +
  tm_shape(subset(nodes, population > 0))  + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
tmap_save(pl_cmp, "figures/transport_network/PIDA/trans_africa_network_PIDA_vs_top10B_consensus.pdf",
          width = 10, height = 10)


# -------------------------------------------------------------------
# Step 9 - Reverse-PDV macro CBA (extend Table tab:MACAB)
# -------------------------------------------------------------------
# Reuses the my_PDV / inv_PDV / calc_rates pattern from
# 8_PE_analysis.R L1497-1529.

AFRGDP22 <- 2811259831806  # Africa GDP 2022 in constant 2015 USD

my_PDV <- function(gdp = AFRGDP22, df = 1.1, gr_old_perc = 4.1,
                   gr_new_perc = 4.11, max_years = 30) {
  gr_old <- 1 + gr_old_perc / 100
  gr_new <- 1 + gr_new_perc / 100
  years  <- 1:max_years
  FV     <- gdp * (gr_new^years - gr_old^years)
  sum(FV / df^years)
}

inv_PDV <- function(PDV = 40e9, ...) {
  objective <- function(x) abs(PDV - my_PDV(gr_new_perc = x, ...))
  optimize(objective, c(0, 100), tol = .Machine$double.eps)$minimum
}

calc_rates <- function(x, bgr) {
  gr_new <- inv_PDV(x, gr_old_perc = bgr)
  c(Rate = gr_new, `Growth of Rate` = (gr_new / bgr - 1) * 100)
}

packages_PIDA <- c(
  "PIDA Projects"        = cost_PIDA,
  "Top-N Consensus ~10B" = cost_Top10B,
  "All Links MA > 4"     = 17.0e9   # existing reference row from tab:MACAB
)

PDV_table_41 <- sapply(packages_PIDA, calc_rates, 4.1) |> t() |> round(4)
PDV_table_30 <- sapply(packages_PIDA, calc_rates, 3.0) |> t() |> round(4)

cat("\n=== Reverse-PDV: required growth bump at 4.1% baseline ===\n")
print(PDV_table_41)
cat("\n=== Reverse-PDV: required growth bump at 3.0% baseline ===\n")
print(PDV_table_30)


# -------------------------------------------------------------------
# Step 10 - Save results
# -------------------------------------------------------------------
PIDA_summary <- list(
  cost_PIDA           = cost_PIDA,
  cost_Top10B         = cost_Top10B,
  N_Top10B            = N_top,
  baseline            = list(MA = MA, MA_bt_opt = MA_bt_opt),
  joint_results_PIDA  = res_PIDA,
  joint_results_Top10B = res_Top10B,
  sum_marginal        = list(NoFR = PIDA_marg_NoFR, FR = PIDA_marg_FR),
  PDV_table_41        = PDV_table_41,
  PDV_table_30        = PDV_table_30,
  km_breakdown        = fsum(edges_all$distance/1e3, edges_all[, .(PIDA, type)]),
  cost_breakdown      = fsum(cost_per_edge, edges_all[, .(PIDA, type)])
)

qs2::qs_save(list(
  summary       = PIDA_summary,
  edges_all     = atomic_elem(edges_all),
  all_cb_ratios = atomic_elem(qDT(all_cb_ratios))
), "results/transport_network/PE/PIDA_PE_results.qs2")

cat("\nDone. Wrote results/transport_network/PE/PIDA_PE_results.qs2\n")
