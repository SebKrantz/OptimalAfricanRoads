#####################################################################
# Transport Network: Analyze Optimal Trans-African GE Investments
#####################################################################

library(fastverse)
set_collapse(mask = c("manip", "helper", "special"))
fastverse_extend(qs2, sf, units, sfnetworks, tmap, install = TRUE)
source("code/helpers/helpers.R")
fastverse_conflicts()

net <- "fastest_routes" # Or 'all_routes' to also include shortest routes
dir <- "trans_african"  # Or 'trans_african_add' to get shortest routes results  
res_name <- "22g_15b_fixed_cgc_irs1.2_na_sigma2.0_rho2_duality_julia" # '22g_add_10b_fixed_duality_sigma3'

results <- list(
  nodes = fread(sprintf("results/transport_network/GE_dual/%s/nodes_results_%s.csv", dir, res_name)),
  edges = fread(sprintf("results/transport_network/GE_dual/%s/edges_results_%s.csv", dir, res_name))
)

network <- qs_read("data/transport_network/trans_african/trans_africa_network_47_largest.qs2") |> 
           extract2(net) |> extract2("network") 
nodes <- network |> st_as_sf("nodes")
edges <- network |> st_as_sf("edges")
if(net == "fastest_routes") {
  edges_real <- qs_read("data/transport_network/trans_african/trans_africa_network_47_largest_fastest_real_edges.qs2") |>
    rmapshaper::ms_simplify(keep = 0.1) |> st_make_valid()
  tfm(edges_real) <- atomic_elem(select(edges, from, to, distance, duration))
}

## ----------------------------------------------------------
# # This information is already part of the results csv files
# settfm(nodes, major_city_port = population > 2e6 | outflows > 1e6)
# sum(nodes$major_city_port)
# largest <- c("Dakar - Senegal", "Casablanca - Morocco", "Abidjan - Cote d'Ivoire", 
#              "Kumasi - Ghana", "Algiers - Algeria", "Lagos - Nigeria", "Kano - Nigeria", 
#              "Yaounde - Cameroon", "Luanda - Angola", "Kinshasa - Congo (Kinshasa)", 
#              "Johannesburg - South Africa", "Cape Town - South Africa", "Cairo - Egypt", 
#              "Khartoum - Sudan", "Nairobi - Kenya", "Addis Ababa - Ethiopia", 
#              "Dar es Salaam - Tanzania")
# settfm(nodes, product = nif(major_city_port & base::match(city_country, largest, 0L) > 0L, NA_integer_, # Heterogeneous products
#                             population > 1e6 & outflows > 1e6, 5L, # Large Port-City
#                             population > 2e6, 4L,   # Large City
#                             outflows > 0, 3L,       # Port
#                             population > 2e5, 2L,   # Medium-Sized City
#                             default = 1L))          # Town/Node
# table(nodes$product, na.exclude = FALSE)
# setv(nodes$product, whichNA(nodes$product), seq_along(largest) + 5L)
# attr(nodes$product, "levels") <- c("Small City/Node", "City > 200K", "Port", "City > 2M", "Large Port-City", paste("Megacity", seq_along(largest)))
# class(nodes$product) <- "factor"
# -----------------------------------------------------------

largest <- results$nodes %>% subset(unclass(product) > 5L) %$% set_names(product, city_country) %>% sort() %>% names()
attr(results$nodes$product, "levels") <- c("Small City/Node", "City > 200K", "Port", "City > 2M", "Large Port-City", largest)
class(results$nodes$product) <- "factor"

results$edges %<>% join(x = edges, on = c("from", "to"), how = "inner", drop = "x") #, if(res_name %ilike% "_bc") "y" else "x")
results$nodes %<>% join(x = nodes, on = c("lon", "lat"), how = "inner", drop = "x") #  if(res_name %ilike% "_bc") "y" else "x")
results$edges %<>% mutate(# distance = distance / 1000,
                          cost_per_km = total_cost / distance)

# Check utility and consumption correspondence
alpha <- if(res_name %ilike% "_alpha01") 0.1 else 0.7
cor(with(results$nodes, utility(Cj*10/copyv(population, 0, 1e-6), alpha = alpha)), results$nodes$uj)
with(results$nodes, cor(inv_utility(uj, alpha = alpha),Cj*10/copyv(population, 0, 1e-6)))
with(results$nodes, all.equal(utility(inv_utility(uj, alpha = alpha), alpha = 0.7), uj))

# Statistics on the upgrade extent -----------------------------------------------------------------

stats <- list()

# Cost of all possible work (should be 33 billion in millions)
results$edges |> with(sum(distance*cost_per_km))
# Budget spent (should give 10 or 20 billion in millions)
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


# Statistics on the economic gains -----------------------------------------------------------------

if(res_name %ilike% "rho2") results$nodes %<>% mutate(uj = (uj*(-1))^(-1), uj_orig = (uj_orig*(-1))^(-1))

# Global Welfare Gains (Ratio)
(stats$wg <- results$nodes |> with(sum(uj * Lj) / sum(uj_orig * Lj_orig)))
# Local Welfare Gains (Ratio)
results$nodes |> with(descr(uj / uj_orig))
# Correlates of Local Welfare Gains (Ratio)
qDF(results$nodes) |> mutate(ugain = uj / uj_orig) |>
  select(ugain, population, gdp_cap, IWI, gdp, wealth) |> pwcor()

# Consumption Gains
(stats$cg <- results$nodes |> with(sum(Cj) / sum(Cj_orig))) # Global
results$nodes |> with(descr(Cj / Cj_orig))    # Local

# Consumption Percent by City Type: May need to change Dj -> Cj if cgc is off (all_routes/dual solutions)
# Overall
qDT(results$nodes) |> group_by(product) |> gvr("^Dj_") |> select(-Dj_orig) |> fsum() |> 
  tfmv(-1, fsum, TRA = "%", apply = FALSE) |> qM(1) %>% {set_names(diag(.), rownames(.))}
# At the city level + median aggregation 
qDT(results$nodes) |> group_by(product) |> gvr("^Dj_") |> select(-Dj_orig) %>%
  dapply(`/`, psum(.)) |> fmedian() |> qM(1) %>% {set_names(diag(.), rownames(.))}
# Per capita
qDT(results$nodes) |> group_by(product) |> gvr("^Dj_|pop") |> select(-Dj_orig) |> 
  tfmv(-population, `/`, population+1) |>
  fmedian() |> tfmv(-(1:2), fsum, TRA = "%", apply = FALSE) |> select(-population) |>
  qM(1) %>% {set_names(diag(.), rownames(.))}

# Market Access Gains --------------------------------------

results$edges %<>% mutate(duration = iif(is.finite(add) & add == 1L, Inf, distance / Ijk_orig))

# Check
results$edges |> with(duration / (distance / Ijk_orig)) |> descr()
results$edges |> with(duration / (distance / Ijk)) |> replace_inf() |> descr()

# New duration
results$edges %<>% mutate(duration_new = distance / Ijk)

# Market Access
net_res <- results$edges |> join(select(edges, from, to), on = c("from", "to"), drop = "x") |> 
  as_sfnetwork(directed = FALSE)
ind_res <- ckmatch(round(select(qDF(results$nodes), lon, lat), 5), mctl(round(st_coordinates(st_geometry(net_res, "nodes")), 5)))
any_duplicated(ind_res)

# Computing Times 
times <- st_network_cost(net_res, weights = results$edges$duration*60)[ind_res, ind_res]
times_new <- st_network_cost(net_res, weights = results$edges$duration_new*60)[ind_res, ind_res]

# Computing total real market access
(MA <- total_MA(times, results$nodes$gdp)) / 1e9

# Total gain
(MA_new <- total_MA(times_new, results$nodes$gdp)) / 1e9

(stats$magp <- (MA_new / MA - 1) * 100) # Percent gain


# Plots of Final Network and Optimal Investments ---------------------------------------------------

results$nodes %<>% mutate(prod2 = set_attr(product, "levels", gsub("/Node|Large ", "", levels(product))), 
                          prod2 = droplevels(fifelse(unclass(prod2) > 5L, NA, prod2)))

stats_legend <- c(paste0("Budget: ", round(stats$b/1e3, 1), "B"), 
                  paste0("Upgraded: ", round(stats$wt[1]), "km (", round(stats$bt[1]/1e3, 1), "B)"),
                  paste0("Mixed: ", round(stats$wt[2]), "km (", round(stats$bt[2]/1e3, 1), "B)"),
                  paste0("Gains: ", paste(paste0(round(c((stats$cg-1)*100, (stats$wg-1)*100, stats$magp), 1), "%", c("C", "W", "MA")), collapse = ", ")))

edges_real %<>% join(results$edges, on = c("from", "to"), drop = "y")
# Upgrade Percent (Extensive Margin)
# <Figures 42, 45, and A26: RHS>
# <Figures 43, 44, A23, A24, A27, A28, A30, A32, A33: Top Panels>
# <Figures 46 and A34: All Panels>
pdf(sprintf("figures/transport_network/GE_dual/%s/trans_africa_network_GE_%s_perc_ug.pdf", dir, res_name), width = 8, height = 8.3)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(mutate(edges_real, perc_ug = pmin(pmax((Ijk - Ijk_orig)/(100 - Ijk_orig)*100, 0), 100))) +
  tm_lines(col = "perc_ug", 
           col.scale = tm_scale_continuous(ticks = seq(0, 100, 20), values = "brewer.yl_or_rd"),
           col.legend = tm_legend(expression(Delta~"%"~"UG"), position = c("left", "bottom"), 
                                  stack = "h", frame = FALSE, height = 16, item.width = 0.5, 
                                  text.size = 1.2, title.size = 1.5, title.padding = c(-0.5, 0, 0, 0)), lwd = 2) +
  # tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(results$nodes, population > 0)) + 
  tm_dots(size = "population",
          size.scale = tm_scale_intervals(breaks = c(0, 0.5e3, 2e3, Inf), values = c(1, 2, 3)*0.1),     
          size.legend = tm_legend("Population (K)", position = tm_pos_in(0.79, 0.15), # bg.color = "white", bg.alpha = 0.5,
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
dev.off()


# Flows of Goods -----------------------------------------------------------------------------------


# Flow of Goods: 8 Cities
# <Figures A23, A24, A27, A28, and A30: Bottom Panels>
# <Figure A35: Top Panels>
pdf(sprintf("figures/transport_network/GE_dual/%s/trans_africa_network_GE_%s_good_flows_8_city.pdf", dir, res_name), width = 9, height = 5)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges_real |> gvr("Qjk_") |> pivot("geometry") |> 
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




# Flow of Goods: All Cities
# <Figures A22, A25, and A29: All Panels>
pdf(sprintf("figures/transport_network/GE/%s/trans_africa_network_GE_%s_good_flows.pdf", dir, res_name), width = 10, height = 12)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
tm_shape(results$edges |> gvr("Qjk_") |> pivot("geometry") |> 
           mutate(variable = set_attr(variable, "levels", stringi::stri_trans_general(levels(results$nodes$product), "latin-ascii")))) +
  tm_facets_wrap("variable", ncols = 5) + 
  tm_lines(col = "value", 
           col.scale = tm_scale_continuous(8, trans = "log1p", values = "yl_or_rd"),
           col.legend = tm_legend("Flow", position = tm_pos_in(0.02, 0.54), height = 10, 
                                  title.size = 0.75, text.size = 0.55, frame = FALSE), lwd = 2) +
  tm_layout(frame = FALSE, 
            panel.label.bg.color = "white",
            panel.label.frame = FALSE) 
dev.off()



# Final Network
# <Figures 42, 45, and A26: LHS>
pdf(sprintf("figures/transport_network/GE/%s/trans_africa_network_GE_%s.pdf", dir, res_name), width = 8, height = 8.3)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(results$edges) +
  tm_lines(col = "Ijk", 
           col.scale = tm_scale_continuous(7, values = "-inferno", limits = c(0, 130)),
           col.legend = tm_legend("Km/h", position = c("left", "bottom"), 
                                  stack = "h", frame = FALSE, height = 16, item.width = 0.5, 
                                  text.size = 1.2, title.size = 1.5, title.padding = c(-0.5, 0, 0, 0)), 
           lwd = 2) +
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
# <Figures 42, 45, and A26: Middle>
pdf(sprintf("figures/transport_network/GE/%s/trans_africa_network_GE_%s_diff.pdf", dir, res_name), width = 8, height = 8.3)
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

# Local Welfare (Utility Per Worker) Gains
# <Figure A31: All Panels>
# <Figures A32, A33, and A45: Bottom Panels>
pdf(sprintf("figures/transport_network/GE/%s/trans_africa_network_GE_%s_upw_gain.pdf", dir, res_name), width = 8, height = 8)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(mutate(results$nodes, ugain = (uj / uj_orig - 1) * 100) |>
             st_as_sf(coords = .c(lon, lat), crs = 4326)) +
  tm_dots(fill = "ugain", 
          fill.scale = tm_scale_intervals(7, breaks = c(-Inf, -25, 0, 25, 50, 100, Inf), values = "turbo"),
          fill.legend = tm_legend("Welfare Gain (%)", position = c("left", "bottom"), frame = FALSE,
                                  text.size = 1.5, title.size = 2), size = 0.2) +
  tm_layout(frame = FALSE) 
dev.off()

# Impact of Frictions on Optimal Investments -------------------------------------------------------

tmap_options(raster.max.cells = 1e6)

# res_name <- "22g_10b_fixed_cgc_sigma3.8_rho0_julia"
res_name <- "22g_20b_fixed_cgc_irs1.2_na_sigma1.5_rho2_duality_julia" # '22g_add_10b_fixed_duality_sigma3'
edges_res <- list(NoFR = fread(sprintf("results/transport_network/GE_dual/%s/edges_results_%s.csv", dir, res_name)),
                  FR = fread(sprintf("results/transport_network/GE_dual/%s/edges_results_%s.csv", dir, sub("_duality", "_bc_duality", res_name))))
edges_res %<>% lapply(select, from, to, Ijk, Ijk_orig) %>% 
  rowbind(idcol = "data") %>% 
  mutate(perc_ug = pmin(pmax((Ijk - Ijk_orig)/(100 - Ijk_orig)*100, 0), 100)) %>% 
  pivot(c("from", "to"), "perc_ug", "data", how = "w") %>% 
  mutate(perc_ug_diff = replace_outliers(FR - NoFR, c(-100, 100), "clip"))
edges_res %<>% join(x = edges, on = c("from", "to"))

descr(edges_res$perc_ug_diff)

# Different in % Upgraded (Extensive Margin)
# <Figures 43 and 44: Bottom Panels>
# <Figures A23, A24, A27, A28, and A30: Middle Panels>
pdf(sprintf("figures/transport_network/GE_dual/%s/trans_africa_network_GE_%s_Ijk_bc_perc_ug_diff.pdf", dir, res_name), width = 8, height = 8)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges_res) +
  tm_lines(col = "perc_ug_diff", 
           col.scale = tm_scale_continuous(ticks = c(-100, -50, -25, 0, 25, 50, 100), 
                                           limits = c(-100, 100), midpoint = 0, values = "-classic_red_blue"), # "-spectral"
           col.legend = tm_legend(expression(Delta~"%"~"UG"), position = c("left", "bottom"), 
                                  stack = "h", frame = FALSE, height = 16, item.width = 0.5, 
                                  text.size = 1.2, title.size = 1.5, title.padding = c(-0.5, 0, 0, 0)), lwd = 1.5) +
  tm_layout(frame = FALSE) 
dev.off()


# Impact of Frictions on Trade Flows ---------------------------------------------------------------

tmap_options(raster.max.cells = 1e7)

# res_name <- "22g_10b_fixed_cgc_sigma3.8_rho0_julia"
edges_res <- list(NoFR = fread(sprintf("results/transport_network/GE/%s/edges_results_%s.csv", dir, res_name)),
                  FR = fread(sprintf("results/transport_network/GE/%s/edges_results_%s.csv", dir, sub("_julia", "_bc_julia", res_name))))
edges_res %<>% lapply(gvr, "^from$|^to$|^Qjk_") # %>% rowbind(idcol = "data")
edges_res$Ratio <- edges_res$FR %>% tfm(slt(., Qjk_1:Qjk_22) %c/% slt(join(slt(edges_res$FR, from, to), edges_res$NoFR), Qjk_1:Qjk_22) %>% 
                                          replace_outliers(c(0, 100), "clip"))
edges_res %<>% lapply(join, x = edges, on = c("from", "to"))

descr(atomic_elem(edges_res$Ratio))
edges_res %>% lapply(. %>% atomic_elem() %>% num_vars() %>% fsum())

pdf(sprintf("figures/transport_network/GE/%s/trans_africa_network_GE_%s_good_flows_bc_ratio.pdf", dir, res_name), width = 10, height = 10)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges_res$Ratio |> gvr("Qjk_") |> pivot("geometry") |> 
           mutate(variable = set_attr(variable, "levels", stringi::stri_trans_general(levels(results$nodes$product), "latin-ascii"))) |> na_omit()) +
  tm_facets_wrap("variable", nrows = 5) + 
  tm_lines(col = "value", 
           col.scale = tm_scale_intervals(breaks = c(0, 0.25, 0.5, 0.75, 1, 1.33, 2, 4, Inf), 
                                          midpoint = 1, values = "-rd_bu"),
           col.legend = tm_legend("Flows Ratio", position = tm_pos_in(0.02, 0.54), height = 10, 
                                  title.size = 0.75, text.size = 0.55, frame = FALSE), lwd = 2) +
  tm_layout(frame = FALSE, 
            panel.label.bg.color = "white",
            panel.label.frame = FALSE) 
dev.off()

