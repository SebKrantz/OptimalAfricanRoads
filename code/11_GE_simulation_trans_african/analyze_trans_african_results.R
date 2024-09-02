#####################################################################
# Transport Network: Analyze Optimal Trans-African GE Investments
#####################################################################

library(fastverse)
set_collapse(mask = c("manip", "helper", "special"))
fastverse_extend(qs, sf, units, sfnetworks, tmap)
source("code/helpers/helpers.R")
fastverse_conflicts()

net <- "fastest_routes" # Or 'all_routes' to also include shortest routes
dir <- "trans_african"  # Or 'trans_african_add' to get shortest routes results  
res_name <- "22g_10b_fixed_cgc_sigma3.8_rho0_julia" 

results <- list(
  nodes = fread(sprintf("results/transport_network/GE/%s/nodes_results_%s.csv", dir, res_name)),
  edges = fread(sprintf("results/transport_network/GE/%s/edges_results_%s.csv", dir, res_name))
)

network <- qread("data/transport_network/trans_african/trans_africa_network_47_largest.qs") |> 
           extract2(net) |> extract2("network") 
nodes <- network |> st_as_sf("nodes")
edges <- network |> st_as_sf("edges")

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

results$edges %<>% join(x = edges, on = c("from", "to"), how = "inner", drop = if(res_name %ilike% "_bc") "y" else "x")
results$nodes %<>% join(x = nodes, on = c("lon", "lat"), how = "inner", drop = if(res_name %ilike% "_bc") "y" else "x")
results$edges %<>% mutate(distance = distance / 1000,
                          cost_per_km = total_cost / distance)

# Check utility and consumption correspondence
alpha <- if(res_name %ilike% "_alpha01") 0.1 else 0.7
cor(with(results$nodes, utility(Cj*10/copyv(population, 0, 1e-6), alpha = alpha)), results$nodes$uj)
with(results$nodes, cor(inv_utility(uj, alpha = alpha),Cj*10/copyv(population, 0, 1e-6)))
with(results$nodes, all.equal(utility(inv_utility(uj, alpha = alpha), alpha = 0.7), uj))

# Statistics on the upgrade extent -----------------------------------------------------------------

# Cost of all possible work (should be 33 billion in millions)
results$edges |> with(sum(distance*cost_per_km))
# Budget spent (should give 10 or 20 billion in millions)
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


# Statistics on the economic gains -----------------------------------------------------------------

if(res_name %ilike% "rho2") results$nodes %<>% mutate(uj = (uj*(-1))^(-1), uj_orig = (uj_orig*(-1))^(-1))

# Global Welfare Gains (Ratio)
results$nodes |> with(sum(uj * Lj) / sum(uj_orig * Lj_orig))
# Local Welfare Gains (Ratio)
results$nodes |> with(descr(uj / uj_orig))
# Correlates of Local Welfare Gains (Ratio)
qDF(results$nodes) |> mutate(ugain = uj / uj_orig) |>
  select(ugain, population, gdp_cap, IWI, gdp, wealth) |> pwcor()

# Consumption Gains
results$nodes |> with(sum(Cj) / sum(Cj_orig)) # Global
results$nodes |> with(descr(Cj / Cj_orig))    # Local

# Consumption Percent by City Type
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


# Plots of Final Network and Optimal Investments ---------------------------------------------------

results$nodes %<>% mutate(prod2 = set_attr(product, "levels", gsub("/Node|Large ", "", levels(product))), 
                          prod2 = droplevels(fifelse(unclass(prod2) > 5L, NA, prod2)))

tmap_options(raster.max.cells = 1e6)

# Final Network
pdf(sprintf("figures/transport_network/GE/%s/trans_africa_network_GE_%s.pdf", dir, res_name), width = 8, height = 8.3)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(results$edges) +
  tm_lines(col = "Ijk", 
           col.scale = tm_scale_continuous(7, values = "inferno", limits = c(0, 130)),
           col.legend = tm_legend("km/h", position = c("left", "bottom"), stack = "h", frame = FALSE,
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_shape(results$nodes) + 
  tm_dots(fill = "prod2", size = if(res_name %ilike% "_add") 0.15 else 0.24, 
          fill.scale = tm_scale_categorical(values = "turbo", value.na = "purple3", label.na = "Megacity (Own)"),
          fill.legend = tm_legend("Product", position = c("left", "bottom"), 
                                  frame = FALSE, text.size = 1.5, title.size = 2)) +
  tm_shape(subset(results$nodes, is.na(prod2))) + tm_dots(size = if(res_name %ilike% "_add") 0.3 else 0.5, fill = "purple3") +
  tm_layout(frame = FALSE)
dev.off()

# Difference (Intensive Margin)
pdf(sprintf("figures/transport_network/GE/%s/trans_africa_network_GE_%s_diff.pdf", dir, res_name), width = 8, height = 8.3)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(mutate(results$edges, diff = pmin(pmax(Ijk - Ijk_orig, 0), 100))) +
  tm_lines(col = "diff", 
           col.scale = tm_scale_continuous(ticks = seq(0, 100, 20), values = "yl_or_rd", limits = c(0, 100)),
           col.legend = tm_legend(expression(Delta~"km/h"), position = c("left", "bottom"), stack = "h", frame = FALSE,
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_shape(results$nodes) + 
  tm_dots(fill = "prod2", size = if(res_name %ilike% "_add") 0.15 else 0.24, 
          fill.scale = tm_scale_categorical(values = "turbo", value.na = "purple3", label.na = "Megacity (Own)"),
          fill.legend = tm_legend("Product", position = c("left", "bottom"), 
                                  frame = FALSE, text.size = 1.5, title.size = 2)) +
  tm_shape(subset(results$nodes, is.na(prod2))) + tm_dots(size = if(res_name %ilike% "_add") 0.3 else 0.5, fill = "purple3") +
  tm_layout(frame = FALSE)
dev.off()

# Upgrade Percent (Extensive Margin)
pdf(sprintf("figures/transport_network/GE/%s/trans_africa_network_GE_%s_perc_ug.pdf", dir, res_name), width = 8, height = 8.3)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(mutate(results$edges, perc_ug = pmin(pmax((Ijk - Ijk_orig)/(100 - Ijk_orig)*100, 0), 100))) +
  tm_lines(col = "perc_ug", 
           col.scale = tm_scale_continuous(ticks = seq(0, 100, 20), values = "yl_or_rd"),
           col.legend = tm_legend("% UG", position = c("left", "bottom"), stack = "h", frame = FALSE,
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_shape(results$nodes) + 
  tm_dots(fill = "prod2", size = if(res_name %ilike% "_add") 0.15 else 0.24, 
          fill.scale = tm_scale_categorical(values = "turbo", value.na = "purple3", label.na = "Megacity (Own)"),
          fill.legend = tm_legend("Product", position = c("left", "bottom"), 
                                  frame = FALSE, text.size = 1.5, title.size = 2)) +
  tm_shape(subset(results$nodes, is.na(prod2))) + tm_dots(size = if(res_name %ilike% "_add") 0.3 else 0.5, fill = "purple3") +
  tm_layout(frame = FALSE)
dev.off()

# Local Welfare (Utility Per Worker) Gains
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


# Flows of Goods -----------------------------------------------------------------------------------

tmap_options(raster.max.cells = 1e7)

# Flow of Goods: All Cities
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

# Flow of Goods: 4 Cities
pdf(sprintf("figures/transport_network/GE/%s/trans_africa_network_GE_%s_good_flows_4_city.pdf", dir, res_name), width = 5, height = 5)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(results$edges |> gvr("Qjk_") |> pivot("geometry") |> 
             mutate(variable = set_attr(variable, "levels", sub(" (Kinshasa)", "", stringi::stri_trans_general(levels(results$nodes$product), "latin-ascii"), fixed = TRUE))) |>
             subset(variable %ilike% "Khartoum|Kano|Kinshasa|Nairobi")) +
  tm_facets_wrap("variable", ncols = 2) + 
  tm_lines(col = "value", 
           col.scale = tm_scale_continuous(8, trans = "log1p", values = "yl_or_rd"),
           col.legend = tm_legend("Flow", position = tm_pos_in(0.02, 0.54), height = 10.5, 
                                  title.size = 0.75, text.size = 0.55, frame = FALSE), lwd = 2) +
  tm_layout(frame = FALSE, 
            panel.label.bg.color = "white",
            panel.label.frame = FALSE) 
dev.off()


# Impact of Frictions on Optimal Investments -------------------------------------------------------

tmap_options(raster.max.cells = 1e6)

# res_name <- "22g_10b_fixed_cgc_sigma3.8_rho0_julia"
edges_res <- list(NoFR = fread(sprintf("results/transport_network/GE/%s/edges_results_%s.csv", dir, res_name)),
                  FR = fread(sprintf("results/transport_network/GE/%s/edges_results_%s.csv", dir, sub("_julia", "_bc_julia", res_name))))
edges_res %<>% lapply(select, from, to, Ijk, Ijk_orig) %>% 
  rowbind(idcol = "data") %>% 
  mutate(perc_ug = pmin(pmax((Ijk - Ijk_orig)/(100 - Ijk_orig)*100, 0), 100)) %>% 
  pivot(c("from", "to"), "perc_ug", "data", how = "w") %>% 
  mutate(perc_ug_diff = replace_outliers(FR - NoFR, c(-100, 100), "clip"))
edges_res %<>% join(x = edges, on = c("from", "to"))

descr(edges_res$perc_ug_diff)

# Different in % Upgraded (Extensive Margin)
pdf(sprintf("figures/transport_network/GE/%s/trans_africa_network_GE_%s_Ijk_bc_perc_ug_diff.pdf", dir, res_name), width = 8, height = 8)
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

