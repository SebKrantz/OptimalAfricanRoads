#####################################################################
# Transport Network: Analyze Optimal Regional GE Investments
#####################################################################

library(fastverse)
set_collapse(mask = c("manip", "helper", "special"), nthreads = 4)
fastverse_extend(qs, sf, units, sfnetworks, tmap)
source("code/helpers/helpers.R")
fastverse_conflicts()

# Load essential network data
load("data/transport_network/trans_africa_network_param.RData")

res_name <- "4g_50b_fixed_cgc_sigma15"

results <- list(
  nodes = fread(sprintf("results/transport_network/regional/nodes_results_%s.csv", res_name)),
  edges = fread(sprintf("results/transport_network/regional/edges_results_%s.csv", res_name))
)

results$edges %<>% join(x = rowbind(select(edges, from, to), select(add_links, from, to)), 
                        on = c("from", "to"))

if(res_name %ilike% "_bc") { # With frictions: add original distance variables
  results$edges %<>% join(rowbind(select(qDF(edges_param), from, to, distance), 
                                  select(qDF(add_links_param), from, to, distance)), 
                          on = c("from", "to"), drop = "x") %>% 
                     mutate(distance = distance / 1000)
}

# Check utility and consumption correspondence
alpha <- if(res_name %ilike% "_alpha01") 0.1 else 0.7
cor(with(results$nodes, utility(Cj*10/copyv(population, 0, 1e-6), alpha = alpha)), results$nodes$uj)
with(results$nodes, cor(inv_utility(uj, alpha = alpha),Cj*10/copyv(population, 0, 1e-6)))
with(results$nodes, all.equal(utility(inv_utility(uj, alpha = alpha), alpha = 0.7), uj))

# Statistics on the work extent -----------------------------------------------------------------

# Cost of all possible work (should be 162 billion in millions)
results$edges |> with(sum(distance*cost_per_km))
# Budget spent (should give 50 billion in millions)
results$edges |> with(sum(replace_inf(pmax((Ijk-Ijk_orig)/(pmax(Ijk_orig, 100)-Ijk_orig), 0), 0)*distance*cost_per_km))
# Amount spent on different types of work (1 = new construction, 0 = upgrade)
results$edges |> with(fsum(replace_inf(pmax((Ijk-Ijk_orig)/(pmax(Ijk_orig, 100)-Ijk_orig), 0), 0)*distance*cost_per_km, add)) # |> proportions()

# Road km built/upgraded
results$edges |> with(sum(replace_inf(pmax((Ijk-Ijk_orig)/(pmax(Ijk_orig, 100)-Ijk_orig), 0), 0)*distance))
# Km on different types of work
results$edges |> with(fsum(replace_inf(pmax((Ijk-Ijk_orig)/(pmax(Ijk_orig, 100)-Ijk_orig), 0), 0)*distance, add)) # |> proportions()

# Number of roads worked on (TRUE) (extensive margin)
results$edges |> with(Ijk-Ijk_orig > 1) |> table()  # |> proportions()
results$edges |> with(table(Ijk-Ijk_orig > 1, add)) # proportions() |> addmargins() # By type

# Work Intensity (intensive margin = km/h added)
results$edges |> with(descr((Ijk-Ijk_orig)[Ijk-Ijk_orig > 1]))
results$edges |> subset(Ijk-Ijk_orig > 1) |> with(descr(Ijk-Ijk_orig, add)) # By type


# Statistics on the economic gains ----------------------------------------------------------------------

results$nodes %<>% mutate(type = nif(outflows > 0, 4L, population >= 1000, 3L, population >= 200, 2L, population > 0, 1L, default = 0L))

# Welfare Gains on Original Scale if Needed (but better resolve model with alpha = 0.7)
if(res_name %ilike% "_alpha01" && !(res_name %ilike% "res_alpha07")) {
  results$nodes %<>% transform(uj = utility(inv_utility(uj, alpha = 0.1)), 
                               uj_orig = utility(inv_utility(uj_orig, alpha = 0.1)))
}

# Global Welfare Gains (Ratio)
results$nodes |> with(sum(uj * Lj) / sum(uj_orig * Lj_orig))
# Local Welfare Gains (Ratio)
results$nodes |> with(descr(uj / uj_orig))
# Correlates of Local Welfare Gains (Ratio)
results$nodes |> mutate(ugain = uj / uj_orig) |>
  select(ugain, population, gdp_cap, IWI, gdp, wealth) |> pwcor()

# Consumption Gains
results$nodes |> with(sum(Cj) / sum(Cj_orig)) # Global
results$nodes |> with(descr(Cj / Cj_orig))    # Local

# Consumption Percent by City Type
# Overall
results$nodes |> group_by(type) |> gvr("Dj_") |> select(-Dj_orig) |> fsum() |> 
  subset(type > 0) |> tfmv(-1, fsum, TRA = "%", apply = FALSE)
# At the city level + median aggregation 
results$nodes |> group_by(type) |> gvr("Dj_") |> select(-Dj_orig) %>%
  dapply(`/`, psum(.)/100) |> fmedian() |> subset(type > 0)
# Per capita
results$nodes |> group_by(type) |> gvr("Dj_|pop") |> select(-Dj_orig) |> 
  tfmv(-population, `/`, population+1) |>
  fmedian() |> subset(type > 0) |> tfmv(-(1:2), fsum, TRA = "%", apply = FALSE)


# Market Access Gains --------------------------------------

results$edges %<>% mutate(duration = iif(add == 1L, Inf, duration))

# Check
results$edges |> with(duration / (distance / Ijk_orig)) |> descr()
results$edges |> with(duration / (distance / Ijk)) |> replace_inf() |> descr()

# New duration
results$edges %<>% mutate(duration_new = distance / Ijk)

# Market Access
net_res <- suppressWarnings(as_sfnetwork(results$edges, directed = FALSE))
ind_res <- ckmatch(round(select(results$nodes, lon, lat), 5), mctl(round(st_coordinates(st_geometry(net_res, "nodes")), 5)))
any_duplicated(ind_res)

# Computing Times 
times <- st_network_cost(net_res, weights = results$edges$duration*60)[ind_res, ind_res]
times_new <- st_network_cost(net_res, weights = results$edges$duration_new*60)[ind_res, ind_res]

# Computing total real market access
(MA <- total_MA(times, results$nodes$gdp)) / 1e9

# Total gain
(MA_new <- total_MA(times_new, results$nodes$gdp)) / 1e9

(MA_new / MA - 1) * 100 # Percent gain
MA_new / MA * 1748.128 - 1748.128 # Gain in billion USD'15/minute


# Plots of Final Network and Optimal Investments ---------------------------------------------------

tmap_options(raster.max.cells = 1e6)

# Final Network
pdf(sprintf("figures/transport_network/GE/regional/trans_africa_network_GE_%s.pdf", res_name), width = 8, height = 8)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(results$edges) +
  tm_lines(col = "Ijk", 
           col.scale = tm_scale_continuous(7, values = "inferno", limits = c(0, 130)),
           col.legend = tm_legend("Speed (km/h)", position = c("left", "bottom"), frame = FALSE,
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 
dev.off()

# Difference (Intensive Margin)
pdf(sprintf("figures/transport_network/GE/regional/trans_africa_network_GE_%s_diff.pdf", res_name), width = 8, height = 8)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(mutate(results$edges, diff = pmin(pmax(Ijk - Ijk_orig, 0), 100))) +
  tm_lines(col = "diff", 
           col.scale = tm_scale_continuous(ticks = seq(0, 100, 20), values = "yl_or_rd", limits = c(0, 100)),
           col.legend = tm_legend("Difference (km/h)", position = c("left", "bottom"), frame = FALSE,
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 
dev.off()

# Upgrade Percent (Extensive Margin)
pdf(sprintf("figures/transport_network/GE/regional/trans_africa_network_GE_%s_perc_ug.pdf", res_name), width = 8, height = 8)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(mutate(results$edges, perc_ug = pmin(pmax((Ijk - Ijk_orig)/(100 - Ijk_orig)*100, 0), 100))) +
  tm_lines(col = "perc_ug", 
           col.scale = tm_scale_continuous(ticks = seq(0, 100, 20), values = "yl_or_rd"),
           col.legend = tm_legend("Upgrade Percent", position = c("left", "bottom"), frame = FALSE,
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 
dev.off()

# Local Welfare (Utility Per Worker) Gains
pdf(sprintf("figures/transport_network/GE/regional/trans_africa_network_GE_%s_upw_gain.pdf", res_name), width = 8, height = 8)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(mutate(results$nodes, ugain = (uj / uj_orig - 1) * 100) |>
             st_as_sf(coords = .c(lon, lat), crs = 4326)) +
  tm_dots(fill = "ugain", 
          fill.scale = tm_scale_intervals(7, breaks = c(-Inf, -25, 0, 25, 50, 100, Inf), values = "turbo"),
          fill.legend = tm_legend("Welfare Gain (%)", position = c("left", "bottom"), frame = FALSE,
                                  text.size = 1.5, title.size = 2), size = 0.2) +
  tm_layout(frame = FALSE) 
dev.off()


# Flows of Goods --------------------------------------------------------------

tmap_options(raster.max.cells = 1e7)

# Flow of Goods
pdf(sprintf("figures/transport_network/GE/regional/trans_africa_network_GE_%s_good_flows.pdf", res_name), width = 10, height = 10)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(results$edges |> gvr("Qjk_") |> 
             rename(Qjk_1 = "Good of Cities Sizes 1-200K",
                    Qjk_2 = "Good of Cities Sizes 200K-1M",
                    Qjk_3 = "Good of Cities Sizes 1M+",
                    Qjk_4 = "International Good of Port Cities"
             ) |> pivot("geometry")) +
  tm_facets_wrap("variable", nrows = 2) + 
  tm_lines(col = "value", 
           col.scale = tm_scale_continuous(8, trans = "log1p", values = "yl_or_rd"), # ticks = c(0, 10, 25, 50, 100, 250, 1000, 5000)
           col.legend = tm_legend("Flow of Good", position = c("left", "bottom"), frame = FALSE), lwd = 2) +
  # tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  # tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE, 
            panel.label.bg.color = "white",
            panel.label.frame = FALSE) 
dev.off()

# Flow of Goods and Local Consumption 
pdf(sprintf("figures/transport_network/GE/regional/trans_africa_network_GE_%s_good_flows_cons.pdf", res_name), width = 10, height = 10)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(results$edges |> gvr("Qjk_") |> 
             rename(Qjk_1 = "Good of Cities Sizes 1-200K",
                    Qjk_2 = "Good of Cities Sizes 200K-1M",
                    Qjk_3 = "Good of Cities Sizes 1M+",
                    Qjk_4 = "International Good of Port Cities") |> pivot("geometry")) +
  tm_lines(col = "value", 
           col.scale = tm_scale_continuous(8, trans = "log1p", values = "yl_or_rd"), # ticks = c(0, 10, 25, 50, 100, 250, 1000, 5000)
           col.legend = tm_legend("Flow", position = c("left", "bottom"), stack = "h", frame = FALSE), lwd = 2) +
  tm_facets_wrap("variable", nrows = 2) + 
  tm_shape(results$nodes |> 
             tfmv(c(Dj_1, Dj_2, Dj_3, Dj_4), `/`, population+1) |>
             st_as_sf(coords = c("lon", "lat"), crs = 4326) |> gvr("Dj_") |> slt(-Dj_orig) |>
             rename(Dj_1 = "Good of Cities Sizes 1-200K",
                    Dj_2 = "Good of Cities Sizes 200K-1M",
                    Dj_3 = "Good of Cities Sizes 1M+",
                    Dj_4 = "International Good of Port Cities") |> pivot("geometry") |>
             mutate(value = value / 1)) +
  tm_dots(fill = "value", 
          fill.scale = tm_scale_intervals(breaks = c(0, 10, 25, 50, 100, 250, 500, 1000, 2500, Inf), values = "inferno"), # tm_scale_intervals(7, style = "fisher", values = "inferno")
          fill.legend = tm_legend("Cons./1KP", position = c("left", "bottom"), frame = FALSE), size = 0.1) + 
  tm_facets_wrap("variable", nrows = 2) + 
  tm_layout(frame = FALSE, 
            panel.label.bg.color = "white",
            panel.label.frame = FALSE) 
dev.off()


# Impact of Frictions on Optimal Investments ----------------------------------------------------------------

tmap_options(raster.max.cells = 1e6)

# res_name <- "4g_50b_fixed_cgc_sigma15"
edges_res <- list(NoFR = fread(sprintf("results/transport_network/regional/edges_results_%s.csv", res_name)),
                  FR = fread(sprintf("results/transport_network/regional/edges_results_%s_bc.csv", res_name)))
edges_res %<>% lapply(select, from, to, Ijk, Ijk_orig) %>% 
  rowbind(idcol = "data") %>% 
  mutate(perc_ug = pmin(pmax((Ijk - Ijk_orig)/(100 - Ijk_orig)*100, 0), 100)) %>% 
  pivot(c("from", "to"), "perc_ug", "data", how = "w") %>% 
  mutate(perc_ug_diff = replace_outliers(FR - NoFR, c(-100, 100), "clip"))
edges_res %<>% join(x = rowbind(select(edges, from, to), select(add_links, from, to)), on = c("from", "to"))

descr(edges_res$perc_ug_diff)

# Different in % Upgraded (Extensive Margin)
pdf(sprintf("figures/transport_network/GE/regional/trans_africa_network_GE_%s_Ijk_bc_perc_ug_diff.pdf", res_name), width = 8, height = 8)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges_res) +
  tm_lines(col = "perc_ug_diff", 
           col.scale = tm_scale_continuous(ticks = c(-100, -50, -25, 0, 25, 50, 100), 
                                           limits = c(-100, 100), midpoint = 0, values = "-classic_red_blue"), # "-spectral"
           col.legend = tm_legend("Diff in % UG", position = c("left", "bottom"), frame = FALSE,
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_layout(frame = FALSE, 
            panel.label.bg.color = "white",
            panel.label.frame = FALSE) 
dev.off()


# Impact of Frictions on Trade Flows ----------------------------------------------------------------

tmap_options(raster.max.cells = 1e7)

# res_name <- "4g_50b_fixed_cgc_sigma15"
edges_res <- list(NoFR = fread(sprintf("results/transport_network/regional/edges_results_%s.csv", res_name)),
                  FR = fread(sprintf("results/transport_network/regional/edges_results_%s_bc.csv", res_name)))
edges_res %<>% lapply(gvr, "^from$|^to$|^Qjk_") # %>% rowbind(idcol = "data")
edges_res$Ratio <- edges_res$FR %>% tfm(slt(., Qjk_1:Qjk_4) %c/% slt(join(slt(edges_res$FR, from, to), edges_res$NoFR), Qjk_1:Qjk_4) %>% 
                                          replace_outliers(c(0, 100), "clip"))
edges_res %<>% lapply(join, x = rowbind(select(edges, from, to), select(add_links, from, to)), on = c("from", "to"))

descr(edges_res$Ratio)
edges_res %>% lapply(. %>% atomic_elem() %>% num_vars() %>% fsum())

pdf(sprintf("figures/transport_network/GE/regional/trans_africa_network_GE_%s_good_flows_bc_ratio.pdf", res_name), width = 10, height = 10)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges_res$Ratio |> gvr("Qjk_") |> 
             rename(Qjk_1 = "Good of Cities Sizes 1-200K",
                    Qjk_2 = "Good of Cities Sizes 200K-1M",
                    Qjk_3 = "Good of Cities Sizes 1M+",
                    Qjk_4 = "International Good of Port Cities") |> pivot("geometry") |> na_omit()) +
  tm_facets_wrap("variable", nrows = 2) + 
  tm_lines(col = "value", 
           col.scale = tm_scale_intervals(breaks = c(0, 0.25, 0.5, 0.75, 1, 1.33, 2, 4, Inf), 
                                          midpoint = 1, values = "-rd_bu"),
           col.legend = tm_legend("Flows Ratio", position = c("left", "bottom"), frame = FALSE), lwd = 2) +
  tm_layout(frame = FALSE, 
            panel.label.bg.color = "white",
            panel.label.frame = FALSE) 
dev.off()

