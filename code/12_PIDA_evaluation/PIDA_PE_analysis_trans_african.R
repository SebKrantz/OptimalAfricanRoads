#####################################################################
# PIDA Road Projects: Trans-African PE Evaluation
# -------------------------------------------------------------------
# Mirrors `code/12_PIDA_evaluation/PIDA_PE_analysis.R` on the SAME
# reduced network used in `PIDA_GE_analysis_trans_african.R`:
# the 47-largest-port-cities fastest-routes graph (330 condensed
# corridor segments, 212 nodes).
#
# This gives an apples-to-apples PE complement to the GE evaluation
# in §7.1 of the paper: both PE and GE now compare PIDA against an
# alternative allocation of the SAME budget on the SAME strategic
# trans-African network.
#
# Frictionless only (consistent with the GE counterpart; §4 of the
# paper shows that frictions move levels but not the spatial
# allocation shape on this strategic network).
#
# Outputs:
#   results/transport_network/PE/PIDA_PE_trans_african_results.qs2
#   figures/transport_network/PIDA/PE_trans_african/
#     trans_africa_network_PE_TA_PIDA_MA_100_min_speed_perc.pdf
#     trans_africa_network_PE_TA_PIDA_MA_gain_pusd.pdf
#     trans_africa_network_PE_TA_PIDA_vs_topN_consensus.pdf
#####################################################################

library(fastverse)
set_collapse(mask = c("manip", "helper", "special"), nthreads = 4)
fastverse_extend(qs2, sf, units, sfnetworks, tmap, install = TRUE)
source("code/helpers/helpers.R")
fastverse_conflicts()

dir.create("figures/transport_network/PIDA/PE_trans_african",
           recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# 1. PIDA edge indices (same as in PIDA_GE_analysis_trans_african.R)
# -------------------------------------------------------------------
# Coverage-based match (buffer = 25km, threshold = 50%), computed in
# code/12_PIDA_evaluation/match_PIDA_to_trans_african.R.
# Replaces the earlier centroid-only match.
PIDA_ind <- c(  7L,  16L,  17L,  18L,  23L,  24L,  26L,  52L,  53L,  59L,
               60L,  69L,  70L,  71L,  72L,  73L,  74L,  86L,  88L,  89L,
               98L, 100L, 106L, 116L, 117L, 134L, 138L, 170L, 177L, 182L,
              183L, 185L, 187L, 189L, 205L, 212L, 213L, 216L, 217L, 218L,
              220L, 221L, 232L, 234L, 235L, 248L, 249L, 255L, 259L, 267L,
              268L, 279L, 285L, 293L, 294L, 309L, 322L, 330L)

# -------------------------------------------------------------------
# 2. Load network + simplified real-routes geometry
# -------------------------------------------------------------------
network_obj <- qs2::qs_read(
  "data/transport_network/trans_african/trans_africa_network_47_largest.qs2"
) |> extract2("fastest_routes") |> extract2("network")

nodes <- network_obj |> st_as_sf("nodes")
edges <- network_obj |> st_as_sf("edges")

edges_real <- qs2::qs_read(
  "data/transport_network/trans_african/trans_africa_network_47_largest_fastest_real_edges.qs2"
) |> rmapshaper::ms_simplify(keep = 0.1) |> st_make_valid()
stopifnot(nrow(edges_real) == nrow(edges))

# Units recap (verified directly from the qs2 edges):
#   distance      [m]       -> divide by 1000 for km
#   duration      [min]     -> divide by 60   for hours
#   duration_imp  [min]     == duration at >=100 km/h (otherwise distance / 100kmh)
#   ug_cost_km    [USD/km]
#   ug_cost       [USD]     == (distance / 1000) * ug_cost_km

edges$cost      <- edges$ug_cost                        # USD per edge
edges$is_PIDA   <- seq_len(nrow(edges)) %in% PIDA_ind

cost_PIDA      <- sum(edges$cost[PIDA_ind])
cost_PIDA_B    <- cost_PIDA / 1e9
n_PIDA         <- length(PIDA_ind)

cat("Trans-African network: ", nrow(edges), " edges, ", nrow(nodes), " nodes\n", sep = "")
cat("PIDA-overlap edges: ", n_PIDA, " (", round(n_PIDA / nrow(edges) * 100, 1),
    "% of network)\n", sep = "")
cat("PIDA package cost on this network: $", round(cost_PIDA_B, 2), "B\n", sep = "")
cat("PIDA package distance: ",
    round(sum(edges$distance[PIDA_ind]) / 1e3, 0), " km\n", sep = "")

# -------------------------------------------------------------------
# 3. Baseline MA (frictionless)
# -------------------------------------------------------------------
# Map sfnetwork nodes back to the `nodes` table (lon/lat match)
nodes_coord <- round(st_coordinates(nodes), 5)
net_node_coord <- round(st_coordinates(st_geometry(network_obj, "nodes")), 5)
ind_nodes <- ckmatch(mctl(nodes_coord), mctl(net_node_coord))
stopifnot(!anyDuplicated(ind_nodes))

times_base <- st_network_cost(network_obj, weights = edges$duration)[ind_nodes, ind_nodes]
MA_base    <- total_MA(times_base, nodes$gdp)
cat(sprintf("Baseline MA on trans-African network: %.3e\n", MA_base))

# -------------------------------------------------------------------
# 4. Per-link MA gain (single-edge upgrade to >=100 km/h)
# -------------------------------------------------------------------
cat("Computing per-link MA gains (330 link upgrades)...\n")
edges$MA_per_link <- sapply(seq_len(nrow(edges)), function(i) {
  w <- copyv(edges$duration, i, edges$duration_imp, vind1 = TRUE)
  inv_dur <- 1 / unclass(st_network_cost(network_obj, weights = w)[ind_nodes, ind_nodes])
  diag(inv_dur) <- 0
  sum(inv_dur %*% nodes$gdp)
})
edges$MA_perc       <- (edges$MA_per_link / MA_base - 1) * 100
# $/min/$ : marginal MA gain (USD-of-GDP per minute) per USD invested
edges$MA_gain_pusd  <- perch_to_diff(edges$MA_per_link, edges$MA_perc) / edges$cost

q1 <- quantile(edges$MA_perc, c(0.5, 0.75, 0.9, 0.99))
q2 <- quantile(edges$MA_gain_pusd, c(0.5, 0.75, 0.9, 0.99))
cat(sprintf("Per-link MA gain percentiles (50/75/90/99): %.4f%% / %.4f%% / %.4f%% / %.4f%%\n",
            q1[1], q1[2], q1[3], q1[4]))
cat(sprintf("Per-link $/min/$ percentiles (50/75/90/99): %.2e / %.2e / %.2e / %.2e\n",
            q2[1], q2[2], q2[3], q2[4]))

# -------------------------------------------------------------------
# 5. Joint PIDA package MA gain
# -------------------------------------------------------------------
w_PIDA      <- ifelse(edges$is_PIDA, edges$duration_imp, edges$duration)
times_PIDA  <- st_network_cost(network_obj, weights = w_PIDA)[ind_nodes, ind_nodes]
MA_PIDA     <- total_MA(times_PIDA, nodes$gdp)
PIDA_perc   <- (MA_PIDA / MA_base - 1) * 100
PIDA_pusd   <- (MA_PIDA - MA_base) / cost_PIDA

cat(sprintf("\nJoint PIDA: MA gain = %.2f%%, $/min/$ = %.2e\n",
            PIDA_perc, PIDA_pusd))
cat(sprintf("Sum of marginal PIDA gains: %.2f%%\n", sum(edges$MA_perc[PIDA_ind])))

# -------------------------------------------------------------------
# 6. Top-N consensus benchmark at PIDA-equivalent cost
# -------------------------------------------------------------------
# Rank edges by $/min/$ (marginal return), accumulate until cost ≈ PIDA cost
ord       <- order(edges$MA_gain_pusd, decreasing = TRUE)
cum_cost  <- cumsum(edges$cost[ord])
n_keep    <- max(which(cum_cost <= cost_PIDA))
top_inds  <- ord[seq_len(n_keep)]
cost_top  <- sum(edges$cost[top_inds])
cat(sprintf("\nTop-N benchmark: N = %d links, cost = $%.2fB (PIDA = $%.2fB)\n",
            length(top_inds), cost_top / 1e9, cost_PIDA_B))

w_top      <- ifelse(seq_len(nrow(edges)) %in% top_inds,
                     edges$duration_imp, edges$duration)
times_top  <- st_network_cost(network_obj, weights = w_top)[ind_nodes, ind_nodes]
MA_top     <- total_MA(times_top, nodes$gdp)
top_perc   <- (MA_top / MA_base - 1) * 100
top_pusd   <- (MA_top - MA_base) / cost_top

cat(sprintf("Joint Top-N: MA gain = %.2f%%, $/min/$ = %.2e\n",
            top_perc, top_pusd))

# Overlap
both_inds <- intersect(top_inds, PIDA_ind)
pida_only <- setdiff(PIDA_ind, top_inds)
top_only  <- setdiff(top_inds, PIDA_ind)
cat(sprintf("\nOverlap: Both = %d links | PIDA only = %d | Top-N only = %d\n",
            length(both_inds), length(pida_only), length(top_only)))
cat(sprintf("Top-N km = %d | PIDA km = %d\n",
            round(sum(edges$distance[top_inds]) / 1e3),
            round(sum(edges$distance[PIDA_ind])  / 1e3)))

# -------------------------------------------------------------------
# 7. Maps
# -------------------------------------------------------------------
# Attach per-link results to the simplified real-edges geometry
e_real <- edges_real
e_real$MA_perc       <- edges$MA_perc
e_real$MA_gain_pusd  <- edges$MA_gain_pusd
e_real$is_PIDA       <- edges$is_PIDA
e_real$which_pkg     <- factor(
  ifelse(seq_len(nrow(e_real)) %in% both_inds, "Both (PIDA & Top-N)",
  ifelse(e_real$is_PIDA,                                    "PIDA only",
  ifelse(seq_len(nrow(e_real)) %in% top_inds,               "Top-N only",
                                                            "Neither"))),
  levels = c("Both (PIDA & Top-N)", "PIDA only", "Top-N only", "Neither"))

# --- Map 1: per-link MA gain (%) on PIDA edges only ----------------
pl <- tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(e_real) + tm_lines(col = "grey80", lwd = 0.6) +
  tm_shape(subset(e_real, is_PIDA)) +
  tm_lines(col = "MA_perc",
           col.scale  = tm_scale_intervals(values = "turbo",
                                           breaks = c(0, 0.05, 0.1, 0.25, 0.5, 1, 2, Inf)),
           col.legend = tm_legend(expression(Delta~"%"~"MA"),
                                  position = c("left", "bottom"), frame = FALSE,
                                  text.size = 1.0, title.size = 1.3), lwd = 2.5) +
  tm_layout(frame = FALSE)
tmap_save(pl,
  "figures/transport_network/PIDA/PE_trans_african/trans_africa_network_PE_TA_PIDA_MA_100_min_speed_perc.pdf",
  width = 10, height = 10)

# --- Map 2: per-link $/min/$ on PIDA edges -------------------------
pl <- tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(e_real) + tm_lines(col = "grey80", lwd = 0.6) +
  tm_shape(subset(e_real, is_PIDA)) +
  tm_lines(col = "MA_gain_pusd",
           col.scale  = tm_scale_intervals(values = "turbo",
                                           breaks = c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10, Inf)),
           col.legend = tm_legend(expression(Delta~"MA"/"USD"),
                                  position = c("left", "bottom"), frame = FALSE,
                                  text.size = 1.0, title.size = 1.3), lwd = 2.5) +
  tm_layout(frame = FALSE)
tmap_save(pl,
  "figures/transport_network/PIDA/PE_trans_african/trans_africa_network_PE_TA_PIDA_MA_gain_pusd.pdf",
  width = 10, height = 10)

# --- Map 3: PIDA vs Top-N comparison -------------------------------
# Overlap statistics (link counts + km) for the on-figure text block
km_both_TA   <- sum(edges$distance[both_inds]) / 1e3
km_pida_o_TA <- sum(edges$distance[pida_only]) / 1e3
km_top_o_TA  <- sum(edges$distance[top_only])  / 1e3

overlap_text_TA <- paste0(
  "Overlap Statistics\n",
  sprintf("Both:        %d links (%s km)\n",
          length(both_inds),
          format(round(km_both_TA),   big.mark = ",")),
  sprintf("PIDA only:   %d links (%s km)\n",
          length(pida_only),
          format(round(km_pida_o_TA), big.mark = ",")),
  sprintf("Top-N only:  %d links (%s km)\n",
          length(top_only),
          format(round(km_top_o_TA),  big.mark = ",")),
  sprintf("PIDA MA gain:   %.1f%%  |  $%.2fB\n",
          PIDA_perc, cost_PIDA_B),
  sprintf("Top-N MA gain:  %.1f%%  |  $%.2fB",
          top_perc, cost_top / 1e9)
)

pl <- tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(subset(e_real, which_pkg == "Neither")) +
    tm_lines(col = "grey80", lwd = 0.6) +
  tm_shape(subset(e_real, which_pkg != "Neither")) +
  tm_lines(col = "which_pkg",
           col.scale = tm_scale_categorical(
             values = c("Both (PIDA & Top-N)" = "purple4",
                        "PIDA only"           = "darkorange",
                        "Top-N only"          = "royalblue3",
                        "Neither"             = "grey80")),
           col.legend = tm_legend("Package",
                                  position = c("left", "bottom"),
                                  frame = FALSE, bg.alpha = 0,
                                  text.size = 1.0, title.size = 1.3),
           lwd = 2.5) +
  # Place overlap text ABOVE the Package legend (tm_credits is top-anchored).
  # Larger font: this panel has less clutter than the full-grid panel.
  tm_credits(overlap_text_TA,
             position = tm_pos_in(0.02, 0.45),
             size = 1.25, fontface = "plain",
             bg.color = "white", bg.alpha = 0) +
  tm_layout(frame = FALSE)
tmap_save(pl,
  "figures/transport_network/PIDA/PE_trans_african/trans_africa_network_PE_TA_PIDA_vs_topN_consensus.pdf",
  width = 10, height = 10)

# -------------------------------------------------------------------
# 8. Reverse-PDV macro CBA (mirrors 8_PE_analysis.R L1497–1535)
# -------------------------------------------------------------------
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
  c(Rate = gr_new, GrowthOfRate = (gr_new / bgr - 1) * 100)
}

packages_TA <- c("PIDA (TA)" = cost_PIDA,
                 "Top-N TA"  = cost_top)
cat("\n=== Reverse-PDV (4.1% baseline) ===\n")
print(round(t(sapply(packages_TA, calc_rates, 4.1)), 4))
cat("\n=== Reverse-PDV (3.0% baseline) ===\n")
print(round(t(sapply(packages_TA, calc_rates, 3.0)), 4))

# -------------------------------------------------------------------
# 9. Save summary
# -------------------------------------------------------------------
summary_list <- list(
  PIDA = list(
    n_links     = n_PIDA,
    cost_USD    = cost_PIDA,
    distance_km = sum(edges$distance[PIDA_ind]) / 1e3,
    MA_gain_perc= PIDA_perc,
    MA_gain_pusd= PIDA_pusd,
    sum_of_marginal_perc = sum(edges$MA_perc[PIDA_ind])
  ),
  TopN = list(
    n_links     = length(top_inds),
    cost_USD    = cost_top,
    distance_km = sum(edges$distance[top_inds]) / 1e3,
    MA_gain_perc= top_perc,
    MA_gain_pusd= top_pusd
  ),
  overlap = list(
    both      = length(both_inds),
    PIDA_only = length(pida_only),
    TopN_only = length(top_only)
  ),
  baseline_MA = MA_base,
  PIDA_ind    = PIDA_ind,
  top_inds    = top_inds,
  edges       = edges |> qDF()
)

qs2::qs_save(summary_list,
  "results/transport_network/PE/PIDA_PE_trans_african_results.qs2")

cat("\nDone. Wrote PIDA_PE_trans_african_results.qs2 and 3 figures",
    "to figures/transport_network/PIDA/PE_trans_african/\n")
