#####################################################################
# Transport Network: Analyze Optimal Regional GE Investments
#####################################################################

library(fastverse)
set_collapse(mask = c("manip", "helper", "special"))
fastverse_extend(qs, sf, units, sfnetworks, tmap, install = TRUE)
source("code/helpers/helpers.R")
fastverse_conflicts()

# Load essential network data
load("data/transport_network/trans_africa_network_param.RData")
edges_real <- qread("data/transport_network/edges_real_simplified.qs") |> 
  select(from, to) |> rmapshaper::ms_simplify(keep = 0.06) |> st_make_valid()

res_name <- "add_22g_35b_fixed_irs1.2_na_sigma3.8_rho2_julia"

# nodes_results_add_22g_20b_fixed_irs_na_sigma2.0_rho2_julia

results <- list(
  nodes = fread(sprintf("results/transport_network/GE_dual/total/nodes_results_%s.csv", res_name)),
  edges = fread(sprintf("results/transport_network/GE_dual/total/edges_results_%s.csv", res_name))
)

results$edges %<>% join(x = rbind(select(add_links, from, to), select(edges_real, from, to)), 
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

stats <- list()

# Cost of all possible work (should be 162 billion in millions)
results$edges |> with(sum(distance*cost_per_km))
# Budget spent (should give 50 billion in millions)
(stats$b <- results$edges |> with(sum(replace_inf(pmax((Ijk-Ijk_orig)/(pmax(Ijk_orig, 100)-Ijk_orig), 0), 0)*distance*cost_per_km)))
# Amount spent on different types of work (1 = new construction, 0 = upgrade)
(stats$bt <- results$edges |> with(fsum(replace_inf(pmax((Ijk-Ijk_orig)/(pmax(Ijk_orig, 100)-Ijk_orig), 0), 0)*distance*cost_per_km, add))) # |> proportions()

# Road km built/upgraded
results$edges |> with(sum(replace_inf(pmax((Ijk-Ijk_orig)/(pmax(Ijk_orig, 100)-Ijk_orig), 0), 0)*distance))
# Km on different types of work
(stats$wt <- results$edges |> with(fsum(replace_inf(pmax((Ijk-Ijk_orig)/(pmax(Ijk_orig, 100)-Ijk_orig), 0), 0)*distance, add))) # |> proportions()
(stats$wtp <-  stats$wt / (results$edges |> with(fsum(distance, add))) * 100) # |> proportions()

# Number of roads worked on (TRUE) (extensive margin)
results$edges |> with(Ijk-Ijk_orig > 1) |> table()  # |> proportions()
results$edges |> with(table(Ijk-Ijk_orig > 1, add)) # proportions() |> addmargins() # By type

# Work Intensity (intensive margin = km/h added)
results$edges |> with(descr((Ijk-Ijk_orig)[Ijk-Ijk_orig > 1]))
results$edges |> subset(Ijk-Ijk_orig > 1) |> with(descr(Ijk-Ijk_orig, add)) # By type


# Statistics on the economic gains ----------------------------------------------------------------------

largest <- results$nodes %>% subset(unclass(product) > 5L) %$% set_names(product, city_country) %>% sort() %>% names()
attr(results$nodes$product, "levels") <- c("Small City/Node", "City > 200K", "Port", "City > 2M", "Large Port-City", largest)
class(results$nodes$product) <- "factor"

if(res_name %ilike% "rho2") results$nodes %<>% mutate(uj = (uj*(-1))^(-1), uj_orig = (uj_orig*(-1))^(-1))

# Global Welfare Gains (Ratio)
(stats$wg <- results$nodes |> with(sum(uj * Lj) / sum(uj_orig * Lj_orig)))
# Local Welfare Gains (Ratio)
results$nodes |> with(descr(uj / uj_orig))
# Correlates of Local Welfare Gains (Ratio)
qDT(results$nodes) |> mutate(ugain = uj / uj_orig) |>
  select(ugain, population, gdp_cap, IWI, gdp, wealth) |> pwcor()

# Consumption Gains
(stats$cg <- results$nodes |> with(sum(Cj) / sum(Cj_orig))) # Global
results$nodes |> with(descr(Cj / Cj_orig))    # Local

# Consumption Percent by City Type: May need to change Dj -> Cj if cgc is off (all_routes/dual solutions)
# Overall
qDT(results$nodes) |> group_by(product) |> gvr("^Cj_") |> select(-Cj_orig) |> fsum() |> 
  tfmv(-1, fsum, TRA = "%", apply = FALSE) |> qM(1) %>% {set_names(diag(.), rownames(.))}
# At the city level + median aggregation 
qDT(results$nodes) |> group_by(product) |> gvr("^Cj_") |> select(-Cj_orig) %>%
  dapply(`/`, psum(.)) |> fmedian() |> qM(1) %>% {set_names(diag(.), rownames(.))}
# Per capita
qDT(results$nodes) |> group_by(product) |> gvr("^Cj_|pop") |> select(-Cj_orig) |> 
  tfmv(-population, `/`, population+1) |>
  fmedian() |> tfmv(-(1:2), fsum, TRA = "%", apply = FALSE) |> select(-population) |>
  qM(1) %>% {set_names(diag(.), rownames(.))}


# Market Access Gains --------------------------------------

results$edges %<>% mutate(duration = iif(add == 1L, Inf, duration))

# Check
results$edges |> with(duration / (distance / Ijk_orig)) |> descr()
results$edges |> with(duration / (distance / Ijk)) |> replace_inf() |> descr()

# New duration
results$edges %<>% mutate(duration_new = distance / Ijk)

# Market Access
net_res <- results$edges |> join(rbind(select(edges, from, to), select(add_links, from, to)), 
                                 on = c("from", "to"), drop = "x") |> as_sfnetwork(directed = FALSE)
ind_res <- ckmatch(round(select(results$nodes, lon, lat), 5), mctl(round(st_coordinates(st_geometry(net_res, "nodes")), 5)))
any_duplicated(ind_res)

# Computing Times 
times <- st_network_cost(net_res, weights = results$edges$duration*60)[ind_res, ind_res]
times_new <- st_network_cost(net_res, weights = results$edges$duration_new*60)[ind_res, ind_res]

# Computing total real market access
(MA <- total_MA(times, results$nodes$gdp)) / 1e9

# Total gain
(MA_new <- total_MA(times_new, results$nodes$gdp)) / 1e9

(stats$magp <- (MA_new / MA - 1) * 100) # Percent gain
MA_new / MA * 1748.128 - 1748.128 # Gain in billion USD'15/minute


# Plots of Final Network and Optimal Investments ---------------------------------------------------

results$nodes %<>% mutate(prod2 = set_attr(product, "levels", gsub("/Node|Large ", "", levels(product))), 
                          prod2 = droplevels(fifelse(unclass(prod2) > 5L, NA, prod2))) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

stats_legend <- c(paste0("Budget: ", round(stats$b/1e3, 1), "B"), 
                  paste0("Upgraded: ", round(stats$wt[1]), "km (", round(stats$bt[1]/1e3, 1), "B)"),
                  paste0("New Roads: ", round(stats$wt[2]), "km (", round(stats$bt[2]/1e3, 1), "B)"),
                  paste0("Gains: ", paste(paste0(round(c((stats$cg-1)*100, (stats$wg-1)*100, stats$magp), 1), "%", c("C", "W", "MA")), collapse = ", ")))

# Upgrade Percent (Extensive Margin)
# <Figures 34 and 38: RHS> 
# <Figures 36, 39, A14, and A18: Top Panels> 
# <Figures A15 and A19: LHS>
# pdf(sprintf("figures/transport_network/GE_dual/total/trans_africa_network_GE_%s_perc_ug.pdf", res_name), width = 8, height = 8)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(mutate(results$edges, perc_ug = pmin(pmax((Ijk - Ijk_orig)/(100 - Ijk_orig)*100, 0), 100))) +
  tm_lines(col = "perc_ug", 
           col.scale = tm_scale_continuous(ticks = seq(0, 100, 20), values = "brewer.yl_or_rd"),
           col.legend = tm_legend(expression(Delta~"%"~"UG"), position = c("left", "bottom"), 
                                  stack = "h", frame = FALSE, height = 16, item.width = 0.5, 
                                  text.size = 1.2, title.size = 1.5, title.padding = c(-0.5, 0, 0, 0)), lwd = 1.5) +
  # tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(results$nodes, population > 0)) + 
  tm_dots(size = "population",
          size.scale = tm_scale_intervals(breaks = c(0, 0.5e3, 2e3, Inf), values = c(1, 2, 3)*0.07),     
          size.legend = tm_legend("Population (K)", position = c("right", "bottom"), # bg.color = "white", bg.alpha = 0.5,
                                  frame = FALSE, text.size = 1.1, item.width = 0.4, title.size = 1.1),
          fill = "prod2", # size = 0.1, 
          fill.scale = tm_scale_categorical(values = "turbo", value.na = "purple3", label.na = "Megacity (Own)"),
          fill.legend = tm_legend("Product", position = c("left", "bottom"), size = 0.8,
                                  frame = FALSE, text.size = 1.2, title.size = 1.5, title.padding = c(0, 0, 0, 0), item.width = 1)) +
  # tm_shape(subset(results$nodes, is.na(prod2))) + tm_dots(size = 0.3, fill = "purple3") +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.07, fill = "grey70") +
  tm_add_legend(title = "Statistics", type = "lines", labels = stats_legend,
                position = tm_pos_in(0.133, 0.23), text.size = 1, title.size = 1.5, 
                item.width = 0.2, item.space = 0.2, frame = FALSE) +
  tm_layout(frame = FALSE) 
# dev.off()

# Local Welfare (Utility Per Worker) Gains
# <Figures 37 and 40: All Panels>
# <Figures A14 and A18 : Bottom Panels>
pdf(sprintf("figures/transport_network/GE_dual/total/trans_africa_network_GE_%s_upw_gain.pdf", res_name), width = 8, height = 8)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(mutate(results$nodes, ugain = (uj / uj_orig - 1) * 100) |>
             st_as_sf(coords = .c(lon, lat), crs = 4326)) +
  tm_dots(col = NULL,
          fill = "ugain", 
          fill.scale = tm_scale_continuous(limits = c(-30, 30), outliers.trunc = rep(TRUE, 2), 
                                           values = "-tableau.classic_orange_blue"), # "matplotlib.coolwarm",  -brewer.rd_yl_bu
          fill.legend = tm_legend("Welfare Gain (%)", position = c("left", "bottom"), frame = FALSE,
                                  height = 16, item.width = 0.5, 
                                  text.size = 1.2, title.size = 1.5, title.padding = c(-0.5, 0, 0, 0)), size = 0.2) +
  tm_layout(frame = FALSE) 
dev.off()


# Flows of Goods --------------------------------------------------------------

# Flow of Goods: 8 Cities
# <Figures A23, A24, A27, A28, and A30: Bottom Panels>
# <Figure A35: Top Panels>
pdf(sprintf("figures/transport_network/GE_dual/total/trans_africa_network_GE_%s_good_flows_8_city.pdf", res_name), width = 9, height = 5)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(results$edges |> gvr("Qjk_") |> pivot("geometry") |> 
             mutate(variable = set_attr(variable, "levels", sub(" (Kinshasa)", "", stringi::stri_trans_general(levels(results$nodes$product), "latin-ascii"), fixed = TRUE))) |>
             subset(variable %ilike% "Khartoum|Kano|Kinshasa|Nairobi|Dakar|Cairo|Johannesburg|Casablanca")) +
  tm_facets_wrap("variable", ncol = 4) + 
  tm_lines(col = "value", 
           col.scale = tm_scale_continuous(8, trans = "log1p", values = "brewer.yl_or_rd"),
           col.legend = tm_legend("Flow", position = tm_pos_in(0.02, 0.54), height = 11.5, 
                                  title.size = 0.75, text.size = 0.5, frame = FALSE), lwd = 1.5) +
  tm_layout(frame = FALSE, 
            panel.label.bg.color = "white",
            panel.label.frame = FALSE) 
dev.off()








# Final Network
# <Figures 34 and 38: LHS>
pdf(sprintf("figures/transport_network/GE_dual/total/trans_africa_network_GE_%s.pdf", res_name), width = 8, height = 8)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(results$edges) +
  tm_lines(col = "Ijk", 
           col.scale = tm_scale_continuous(7, values = "-inferno", limits = c(0, 130)),
           col.legend = tm_legend("Km/h", position = c("left", "bottom"), 
                                  stack = "h", frame = FALSE, height = 16, item.width = 0.5, 
                                  text.size = 1.2, title.size = 1.5, title.padding = c(-0.5, 0, 0, 0)), 
           lwd = 1.5) +
  # tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(results$nodes, population > 0)) + 
  tm_dots(fill = "prod2", size = 0.1, 
          fill.scale = tm_scale_categorical(values = "turbo", value.na = "purple3", label.na = "Megacity (Own)"),
          fill.legend = tm_legend("Product", position = c("left", "bottom"), size = 0.8,
                                  frame = FALSE, text.size = 1.2, title.size = 1.5, title.padding = c(0, 0, 0, 0), item.width = 1)) +
  tm_shape(subset(results$nodes, is.na(prod2))) + tm_dots(size = 0.3, fill = "purple3") +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_add_legend(title = "Statistics", type = "lines", labels = stats_legend,
                position = tm_pos_in(0.133, 0.23), text.size = 1, title.size = 1.5, 
                item.width = 0.2, item.space = 0.2, frame = FALSE) +
  tm_layout(frame = FALSE) 
dev.off()

# Difference (Intensive Margin)
# <Figures 34 and 38: Middle>
pdf(sprintf("figures/transport_network/GE_dual/total/trans_africa_network_GE_%s_diff.pdf", res_name), width = 8, height = 8)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(mutate(results$edges, diff = pmin(pmax(Ijk - Ijk_orig, 0), 100))) +
  tm_lines(col = "diff", 
           col.scale = tm_scale_continuous(ticks = seq(0, 100, 20), values = "brewer.yl_or_rd", limits = c(0, 100)),
           col.legend = tm_legend(expression(Delta~"km/h"), position = c("left", "bottom"), 
                                  stack = "h", frame = FALSE, height = 16, item.width = 0.5, 
                                  text.size = 1.2, title.size = 1.5, title.padding = c(-0.5, 0, 0, 0)), lwd = 1.5) +
  # tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(results$nodes, population > 0)) + 
  tm_dots(fill = "prod2", size = 0.1, 
          fill.scale = tm_scale_categorical(values = "turbo", value.na = "purple3", label.na = "Megacity (Own)"),
          fill.legend = tm_legend("Product", position = c("left", "bottom"), size = 0.8,
                                  frame = FALSE, text.size = 1.2, title.size = 1.5, title.padding = c(0, 0, 0, 0), item.width = 1)) +
  tm_shape(subset(results$nodes, is.na(prod2))) + tm_dots(size = 0.3, fill = "purple3") +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_add_legend(title = "Statistics", type = "lines", labels = stats_legend,
                position = tm_pos_in(0.133, 0.23), text.size = 1, title.size = 1.5, 
                item.width = 0.2, item.space = 0.2, frame = FALSE) +
  tm_layout(frame = FALSE) 
dev.off()










tmap_options(raster.max.cells = 1e7)

# Flow of Goods: All Cities
# <Figures A22, A25, and A29: All Panels>
pdf(sprintf("figures/transport_network/GE_dual/total/trans_africa_network_GE_%s_good_flows.pdf", res_name), width = 10, height = 12)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(results$edges |> gvr("Qjk_") |> pivot("geometry") |> 
             mutate(variable = set_attr(variable, "levels", stringi::stri_trans_general(levels(results$nodes$product), "latin-ascii")))) +
  tm_facets_wrap("variable", ncols = 5) + 
  tm_lines(col = "value", 
           col.scale = tm_scale_continuous(8, trans = "log1p", values = "brewer.yl_or_rd"),
           col.legend = tm_legend("Flow", position = tm_pos_in(0.02, 0.54), height = 10, 
                                  title.size = 0.75, text.size = 0.55, frame = FALSE), lwd = 2) +
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

res_name <- "add_22g_35b_fixed_irs1.2_na_sigma3.8_rho0_julia"
edges_res <- list(NoFR = fread(sprintf("results/transport_network/GE_dual/total/edges_results_%s.csv", res_name)),
                  FR = fread(sprintf("results/transport_network/GE_dual/total/edges_results_%s.csv", sub("julia", "bc_julia", res_name))))
edges_res %<>% lapply(select, from, to, Ijk, Ijk_orig) %>% 
  rowbind(idcol = "data") %>% 
  mutate(perc_ug = pmin(pmax((Ijk - Ijk_orig)/(100 - Ijk_orig)*100, 0), 100)) %>% 
  pivot(c("from", "to"), "perc_ug", "data", how = "w") %>% 
  mutate(perc_ug_diff = replace_outliers(FR - NoFR, c(-100, 100), "clip"))
edges_res %<>% join(x = rbind(select(add_links, from, to), select(edges_real, from, to)), on = c("from", "to"))

descr(edges_res$perc_ug_diff)

# Different in % Upgraded (Extensive Margin)
# <Figures 36 and 39: Bottom Panels> 
# <Figures A15 and A19: RHS> 
pdf(sprintf("figures/transport_network/GE_dual/total/trans_africa_network_GE_%s_Ijk_bc_perc_ug_diff.pdf", res_name), width = 8, height = 8)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges_res) +
  tm_lines(col = "perc_ug_diff", 
           col.scale = tm_scale_continuous(ticks = c(-100, -50, -25, 0, 25, 50, 100), 
                                           limits = c(-100, 100), midpoint = 0, values = "-classic_red_blue"), # "-spectral"
           col.legend = tm_legend(expression(Delta~"%"~"UG"), position = c("left", "bottom"), 
                                  stack = "h", frame = FALSE, height = 16, item.width = 0.5, 
                                  text.size = 1.2, title.size = 1.5, title.padding = c(-0.5, 0, 0, 0)), lwd = 1.5) +
  # tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.07, fill = "grey30") +
  # tm_shape(subset(results$nodes, population > 0)) + 
  # tm_dots(size = "population",
  #         size.scale = tm_scale_intervals(breaks = c(0, 0.5e3, 2e3, Inf), values = c(1, 2, 3)*0.07),     
  #         size.legend = tm_legend("Population (K)", position = c("right", "bottom"), # bg.color = "white", bg.alpha = 0.5,
  #                                 frame = FALSE, text.size = 1.1, item.width = 0.4, title.size = 1.1),
  #         fill = "prod2", # size = 0.1, 
  #         fill.scale = tm_scale_categorical(values = "turbo", value.na = "purple3", label.na = "Megacity (Own)"),
  #         fill.legend = tm_legend("Product", position = c("left", "bottom"), size = 0.8,
  #                                 frame = FALSE, text.size = 1.2, title.size = 1.5, title.padding = c(0, 0, 0, 0), item.width = 1)) +
  # tm_shape(subset(results$nodes, is.na(prod2))) + tm_dots(size = 0.3*0.7, fill = "purple3") +
  # tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.07, fill = "grey70") +
  tm_layout(frame = FALSE) 
dev.off()

# Impact on Welfare
nodes_res <- list(NoFR = fread(sprintf("results/transport_network/GE_dual/total/nodes_results_%s.csv", res_name)),
                  FR = fread(sprintf("results/transport_network/GE_dual/total/nodes_results_%s.csv", sub("julia", "bc_julia", res_name))))
# Aggregate Welfare
nodes_res |> lapply(with, fsum(uj_orig * Lj_orig)) |> with((FR/NoFR-1)*100) # Global
nodes_res %<>% lapply(select, city_country, lon, lat, uj, uj_orig) %>% 
  rowbind(idcol = "data") %>% 
  mutate(wg = (uj/uj_orig-1)*100) %>% 
  pivot(c("city_country", "lon", "lat"), "uj_orig", "data", how = "w") %>% 
  mutate(loss = replace_outliers((FR / NoFR - 1)*100, c(-100, 100), "clip"),
         diff = replace_outliers(FR - NoFR, c(-100, 100), "clip"))

descr(nodes_res$diff)

# Different in Welfare
pdf(sprintf("figures/transport_network/GE_dual/total/trans_africa_network_GE_%s_upw_gain_bc_loss.pdf", res_name), width = 8, height = 8)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(st_as_sf(nodes_res, coords = c("lon", "lat"), crs = 4326)) +
  tm_dots(col = NULL,
          fill = "loss", 
          fill.scale = tm_scale_continuous(# limits = c(-30, 30), outliers.trunc = rep(TRUE, 2), 
                                           midpoint = 1,
                                           values = "-tableau.classic_orange_blue"), # "matplotlib.coolwarm",  -brewer.rd_yl_bu
          fill.legend = tm_legend("Welfare Gain (%)", position = c("left", "bottom"), frame = FALSE,
                                  height = 16, item.width = 0.5, 
                                  text.size = 1.2, title.size = 1.5, title.padding = c(0, 0, 0, 0)), size = 0.2) +
  tm_layout(frame = FALSE) 
dev.off()


# Impact of Frictions on Trade Flows ----------------------------------------------------------------

tmap_options(raster.max.cells = 1e7)

# res_name <- "4g_50b_fixed_cgc_sigma15"
edges_res <- list(NoFR = fread(sprintf("results/transport_network/GE/regional/edges_results_%s.csv", res_name)),
                  FR = fread(sprintf("results/transport_network/GE/regional/edges_results_%s_bc.csv", res_name)))
edges_res %<>% lapply(gvr, "^from$|^to$|^Qjk_") # %>% rowbind(idcol = "data")
edges_res$Ratio <- edges_res$FR %>% tfm(slt(., Qjk_1:Qjk_4) %c/% slt(join(slt(edges_res$FR, from, to), edges_res$NoFR), Qjk_1:Qjk_4) %>% 
                                          replace_outliers(c(0, 100), "clip"))
edges_res %<>% lapply(join, x = rowbind(select(edges, from, to), select(add_links, from, to)), on = c("from", "to"))

descr(edges_res$Ratio)
edges_res %>% lapply(. %>% atomic_elem() %>% num_vars() %>% fsum())

# <Figures A16 and A20: All Panels> 
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

