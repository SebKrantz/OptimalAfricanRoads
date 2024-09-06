#########################
# Real Market Access Maps
#########################

library(fastverse)
fastverse_extend(qs, dggridR, ggplot2, viridis, install = TRUE)

# Loading Distance Matrix
africa_dist <- qread("data/africa_full_distance_matrix_r9_adjusted.qs")
names(africa_dist)

# Loading Network Efficiency Estimates 
NE <- qread("results/full_network/network_efficiency.qs")

OUTCOMES <- fread("data/QSE/model_calibration_data_ctry_min_imp.csv") %>%
            fmutate(total_wealth = IWI * pop_wpop)

centroids <- africa_dist$centroids
centroids %<>% 
  ftransform(ss(OUTCOMES, match(centroids$cell, OUTCOMES$cell), c("GDP_PPP", "total_wealth")))
replace_na(centroids[, .(GDP_PPP, total_wealth)], set = TRUE)

# Save all MA estimates
market_access_all <- list()

###################
# Distance Matrices
###################

durations <- africa_dist$durations / 60
diag(durations) <- 50.276890494384 / 2 / NE$NTE # Time efficiency
frange(durations)

distances <- replace_na(africa_dist$distances_nosphere, Inf) / 1000
diag(distances) <- 50.276890494384 / 2 / NE$NRE
frange(distances)
fquantile(distances)

################
# 2015 GDP PPP
################

# Using Travel Time  -------------------------------

# Theta == 1 (Peng & Chen, 2021)
market_access <- (1/durations) %*% centroids$GDP_PPP
frange(market_access)
market_access_all$baseline$ttime$GDP_PPP$theta_1 <- drop(market_access)
cor(market_access, centroids$GDP_PPP)

centroids %>% 
  add_vars(MA = replace_outliers(drop(market_access) / 1e9, 1, NA, "min")) %>% 
  ggplot(aes(x = lon, y = lat, fill = MA)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  scale_fill_viridis(option = "H", n.breaks = 10, transform = "log10", na.value = "#30123BFF",
                     labels = function(x) paste0("$", signif(x, 2), "B")) +
  labs(title = "Market Access (Total 2015 GDP PPP per Road Hour)", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/gdp_per_hour.pdf", width = 6, height = 6)

# Theta == 3.8 (Jedwab & Storeygard, 2022)
market_access <- (1/durations^3.8) %*% centroids$GDP_PPP
frange(market_access)
market_access_all$baseline$ttime$GDP_PPP$theta_3.8 <- drop(market_access)
cor(market_access, centroids$GDP_PPP)

centroids %>% 
  add_vars(MA = drop(market_access) / 1e9) %>% 
  ggplot(aes(x = lon, y = lat, fill = MA)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  scale_fill_viridis(option = "H", n.breaks = 9, transform = "log10", limits = c(NA, 100),
                     labels = function(x) paste0("$", signif(x, 2), "B")) +
  labs(title = "Market Access (Total 2015 GDP PPP per Road Hour to the Power 3.8)", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = -0.1), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/gdp_per_hour_theta_3.8.pdf", width = 6, height = 6)


# Using Road Distance  -------------------------------

# Theta == 1 (Peng & Chen, 2021)
market_access <- (1/distances) %*% centroids$GDP_PPP
frange(market_access)
market_access_all$baseline$tdist$GDP_PPP$theta_1 <- drop(market_access)
cor(market_access, centroids$GDP_PPP)

centroids %>% 
  add_vars(MA = replace_outliers(drop(market_access) / 1e6, 10, NA, "min")) %>% 
  ggplot(aes(x = lon, y = lat, fill = MA)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  scale_fill_viridis(option = "H", n.breaks = 10, transform = "log10", na.value = "#30123BFF",
                     labels = function(x) paste0("$", signif(x, 2), "M")) +
  labs(title = "Market Access (Total 2015 GDP PPP per Road Km)", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/gdp_per_km.pdf", width = 6, height = 6)


# Theta == 3.8 (Jedwab & Storeygard, 2022)
market_access <- (1/distances^3.8) %*% centroids$GDP_PPP
frange(market_access)
market_access_all$baseline$tdist$GDP_PPP$theta_3.8 <- drop(market_access)
cor(market_access, centroids$GDP_PPP)

centroids %>% 
  add_vars(MA = drop(market_access)) %>% 
  ggplot(aes(x = lon, y = lat, fill = MA)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  scale_fill_viridis(option = "H", n.breaks = 8, transform = "log10", limits = c(NA, 1e4),
                     labels = function(x) paste0("$", signif(x, 2))) +
  labs(title = "Market Access (Total 2015 GDP PPP per Road Km to the Power 3.8)", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = -0.1), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/gdp_per_km_theta_3.8.pdf", width = 6, height = 6)

# # Interactive Plot
# data <- centroids %>% 
#   add_vars(log10MA = drop(log10(market_access))) %>% 
#   fselect(log10MA, cell, lon, lat) %>% 
#   ftransform(longitude = lon, latitude = lat) %>% 
#   replace_NA() %>% 
#   st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
#   st_make_valid()
# 
# mapview::mapview(data["log10MA"], lwd = 0)


########################
# Total IWI-Based Wealth
########################

# Using Travel Time  -------------------------------

# Theta == 1 (Peng & Chen, 2021)
market_access <- (1/durations) %*% centroids$total_wealth
frange(market_access)
market_access_all$baseline$ttime$total_wealth$theta_1 <- drop(market_access)
cor(market_access, centroids$total_wealth)

centroids %>% 
  add_vars(MA = replace_outliers(drop(market_access) / 1e6, 10, NA, "min")) %>% 
  ggplot(aes(x = lon, y = lat, fill = MA)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  scale_fill_viridis(option = "H", n.breaks = 10, limits = c(NA, 200), transform = "log10", na.value = "#30123BFF",
                     labels = function(x) paste0("$", signif(x, 2), "M")) +
  labs(title = "Market Access (Total IWI-Based Wealth per Road Hour)", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/wealth_per_hour.pdf", width = 6, height = 6)

# Theta == 3.8 (Jedwab & Storeygard, 2022)
market_access <- (1/durations^3.8) %*% centroids$total_wealth
frange(market_access)
market_access_all$baseline$ttime$total_wealth$theta_3.8 <- drop(market_access)
cor(market_access, centroids$total_wealth)

centroids %>% 
  add_vars(MA = drop(market_access) / 1e6) %>% 
  ggplot(aes(x = lon, y = lat, fill = MA)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  scale_fill_viridis(option = "H", n.breaks = 9, transform = "log10", limits = c(NA, 600),
                     labels = function(x) paste0("$", signif(x, 2), "M")) +
  labs(title = "Market Access (Total IWI-Based Wealth per Road Hour to the Power 3.8)", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = -0.1), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/wealth_per_hour_theta_3.8.pdf", width = 6, height = 6)


# Using Road Distance  -------------------------------

# Theta == 1 (Peng & Chen, 2021)
market_access <- (1/distances) %*% centroids$total_wealth
frange(market_access)
market_access_all$baseline$tdist$total_wealth$theta_1 <- drop(market_access)
cor(market_access, centroids$total_wealth)

centroids %>% 
  add_vars(MA = replace_outliers(drop(market_access) / 1e6, 0.1, NA, "min")) %>% 
  ggplot(aes(x = lon, y = lat, fill = MA)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  scale_fill_viridis(option = "H", n.breaks = 10, transform = "log10", na.value = "#30123BFF", limits = c(NA, 3),
                     labels = function(x) paste0("$", signif(x, 2), "M")) +
  labs(title = "Market Access (Total IWI-Based Wealth per Road Km)", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/wealth_per_km.pdf", width = 6, height = 6)


# Theta == 3.8 (Jedwab & Storeygard, 2022)
market_access <- (1/distances^3.8) %*% centroids$total_wealth
frange(market_access)
market_access_all$baseline$tdist$total_wealth$theta_3.8 <- drop(market_access)
cor(market_access, centroids$total_wealth)

centroids %>% 
  add_vars(MA = drop(market_access)) %>% 
  ggplot(aes(x = lon, y = lat, fill = MA)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  scale_fill_viridis(option = "H", n.breaks = 8, transform = "log10", 
                     labels = function(x) paste0("$", signif(x, 2))) +
  labs(title = "Market Access (Total IWI-Based Wealth per Road Km to the Power 3.8)", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = -0.1), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/wealth_per_km_theta_3.8.pdf", width = 6, height = 6)



#############################
# Now Adding Border Frictions
#############################

# Using Travel Time

border_time <- fread("data/QSE/model_border_time_mat_transit.csv") |> qM(1)

durations_borders <- durations

for (i in rownames(border_time)) {
  for (j in colnames(border_time)) {
    durations_borders[centroids$ISO3 %==% i, centroids$ISO3 %==% j] <- 
      durations_borders[centroids$ISO3 %==% i, centroids$ISO3 %==% j] + border_time[i, j] / 60
  }
}

# Using Travel Distance

border_dist <- fread("data/QSE/model_border_dist_mat_transit.csv") |> qM(1)

distances_borders <- distances

for (i in rownames(border_dist)) {
  for (j in colnames(border_dist)) {
    distances_borders[centroids$ISO3 %==% i, centroids$ISO3 %==% j] <- 
      distances_borders[centroids$ISO3 %==% i, centroids$ISO3 %==% j] + border_dist[i, j] / 1000
  }
}


##############################
# 2015 GDP PPP: With Frictions
##############################

# Using Travel Time  -------------------------------

# Theta == 1 (Peng & Chen, 2021)
market_access <- (1/durations_borders) %*% centroids$GDP_PPP
frange(market_access)
market_access_all$frictions$ttime$GDP_PPP$theta_1 <- drop(market_access)
cor(market_access, centroids$GDP_PPP)

centroids %>% 
  add_vars(MA = replace_outliers(drop(market_access) / 1e9, 1, NA, "min")) %>% 
  ggplot(aes(x = lon, y = lat, fill = MA)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  scale_fill_viridis(option = "H", n.breaks = 10, transform = "log10", na.value = "#30123BFF",
                     labels = function(x) paste0("$", signif(x, 2), "B")) +
  labs(title = "Market Access (Total 2015 GDP PPP per Road Hour)", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/gdp_per_hour_btime.pdf", width = 6, height = 6)

# Theta == 3.8 (Jedwab & Storeygard, 2022)
market_access <- (1/durations_borders^3.8) %*% centroids$GDP_PPP
frange(market_access)
market_access_all$frictions$ttime$GDP_PPP$theta_3.8 <- drop(market_access)
cor(market_access, centroids$GDP_PPP)

centroids %>% 
  add_vars(MA = drop(market_access) / 1e9) %>% 
  ggplot(aes(x = lon, y = lat, fill = MA)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  scale_fill_viridis(option = "H", n.breaks = 9, transform = "log10", limits = c(NA, 100),
                     labels = function(x) paste0("$", signif(x, 2), "B")) +
  labs(title = "Market Access (Total 2015 GDP PPP per Road Hour to the Power 3.8)", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = -0.1), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/gdp_per_hour_btime_theta_3.8.pdf", width = 6, height = 6)


# Using Road Distance  -------------------------------

# Theta == 1 (Peng & Chen, 2021)
market_access <- (1/distances_borders) %*% centroids$GDP_PPP
frange(market_access)
market_access_all$frictions$tdist$GDP_PPP$theta_1 <- drop(market_access)
cor(market_access, centroids$GDP_PPP)

centroids %>% 
  add_vars(MA = replace_outliers(drop(market_access) / 1e6, 10, NA, "min")) %>% 
  ggplot(aes(x = lon, y = lat, fill = MA)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  scale_fill_viridis(option = "H", n.breaks = 10, transform = "log10", na.value = "#30123BFF",
                     labels = function(x) paste0("$", signif(x, 2), "M")) +
  labs(title = "Market Access (Total 2015 GDP PPP per Road Km)", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/gdp_per_km_bdist.pdf", width = 6, height = 6)


# Theta == 3.8 (Jedwab & Storeygard, 2022)
market_access <- (1/distances_borders^3.8) %*% centroids$GDP_PPP
frange(market_access)
market_access_all$frictions$tdist$GDP_PPP$theta_3.8 <- drop(market_access)
cor(market_access, centroids$GDP_PPP)

centroids %>% 
  add_vars(MA = drop(market_access)) %>% 
  ggplot(aes(x = lon, y = lat, fill = MA)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  scale_fill_viridis(option = "H", n.breaks = 8, transform = "log10", limits = c(NA, 1e4),
                     labels = function(x) paste0("$", signif(x, 2))) +
  labs(title = "Market Access (Total 2015 GDP PPP per Road Km to the Power 3.8)", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = -0.1), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/gdp_per_km_bdist_theta_3.8.pdf", width = 6, height = 6)


########################################
# Total IWI-Based Wealth: With Frictions
########################################

# Using Travel Time  -------------------------------

# Theta == 1 (Peng & Chen, 2021)
market_access <- (1/durations_borders) %*% centroids$total_wealth
frange(market_access)
market_access_all$frictions$ttime$total_wealth$theta_1 <- drop(market_access)
cor(market_access, centroids$total_wealth)

centroids %>% 
  add_vars(MA = replace_outliers(drop(market_access) / 1e6, 10, NA, "min")) %>% 
  ggplot(aes(x = lon, y = lat, fill = MA)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  scale_fill_viridis(option = "H", n.breaks = 10, limits = c(NA, 200), 
                     transform = "log10", na.value = "#30123BFF",
                     labels = function(x) paste0("$", signif(x, 2), "M")) +
  labs(title = "Market Access (Total IWI-Based Wealth per Road Hour)", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/wealth_per_hour_btime.pdf", width = 6, height = 6)

# Theta == 3.8 (Jedwab & Storeygard, 2022)
market_access <- (1/durations_borders^3.8) %*% centroids$total_wealth
frange(market_access)
market_access_all$frictions$ttime$total_wealth$theta_3.8 <- drop(market_access)
cor(market_access, centroids$total_wealth)

centroids %>% 
  add_vars(MA = drop(market_access) / 1e6) %>% 
  ggplot(aes(x = lon, y = lat, fill = MA)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  scale_fill_viridis(option = "H", n.breaks = 9, transform = "log10", limits = c(NA, 600),
                     labels = function(x) paste0("$", signif(x, 2), "M")) +
  labs(title = "Market Access (Total IWI-Based Wealth per Road Hour to the Power 3.8)", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = -0.1), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/wealth_per_hour_btime_theta_3.8.pdf", width = 6, height = 6)


# Using Road Distance  -------------------------------


# Theta == 1 (Peng & Chen, 2021)
market_access <- (1/distances_borders) %*% centroids$total_wealth
frange(market_access)
market_access_all$frictions$tdist$total_wealth$theta_1 <- drop(market_access)
cor(market_access, centroids$total_wealth)

centroids %>% 
  add_vars(MA = replace_outliers(drop(market_access) / 1e6, 0.1, NA, "min")) %>% 
  ggplot(aes(x = lon, y = lat, fill = MA)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  scale_fill_viridis(option = "H", n.breaks = 10, transform = "log10", na.value = "#30123BFF", limits = c(NA, 3),
                     labels = function(x) paste0("$", signif(x, 2), "M")) +
  labs(title = "Market Access (Total IWI-Based Wealth per Road Km)", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/wealth_per_km_bdist.pdf", width = 6, height = 6)


# Theta == 3.8 (Jedwab & Storeygard, 2022)
market_access <- (1/distances_borders^3.8) %*% centroids$total_wealth
frange(market_access)
market_access_all$frictions$tdist$total_wealth$theta_3.8 <- drop(market_access)
cor(market_access, centroids$total_wealth)

centroids %>% 
  add_vars(MA = drop(market_access)) %>% 
  ggplot(aes(x = lon, y = lat, fill = MA)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  # geom_point(size = 0.5) + 
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  scale_fill_viridis(option = "H", n.breaks = 8, transform = "log10", 
                     labels = function(x) paste0("$", signif(x, 2))) +
  labs(title = "Market Access (Total IWI-Based Wealth per Road Km to the Power 3.8)", fill = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = -0.1), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/wealth_per_km_bdist_theta_3.8.pdf", width = 6, height = 6)



########################################
# Analysis of All Scenarios
########################################

str(market_access_all)

market_access_all_dt <- market_access_all |>
  rapply2d(function(x) data.table(value = x, cell = centroids$cell)) |> 
  unlist2d(idcols = .c(scenario, metric, measure, theta), 
           DT = TRUE, id.factor = TRUE) |> # get_vars(1:5) |> any_duplicated()
  pivot(names = "measure", how = "w") 

market_access_all_sum <- market_access_all_dt |>
  collap(GDP_PPP + total_wealth ~ scenario + metric + theta, fsum)

market_access_all_cor <- market_access_all |>
  rapply2d(cor, centroids |> with(cbind(GDP_PPP, total_wealth))) |>
  unlist2d(idcols = .c(scenario, metric, measure, theta), 
           DT = TRUE, id.factor = TRUE) |> pivot(1:4) |> 
  fsubset(measure == variable, -variable) |>
  pivot(names = "measure", how = "w") 

# Saving
list(RAW = market_access_all_dt, 
     SUM = market_access_all_sum, 
     COR = market_access_all_cor) |>
  qsave("results/full_network/real_market_access.qs")

# MA in billions of min or km
options(scipen = 1000)
market_access_all_sum |>
  tfmv(is.numeric, `/`, iif(metric == "ttime", 60, 1) * 1e9)

# MA Decline
market_access_all_sum |>
  roworder(metric, theta) |>
  tfmv(is.numeric, fgrowth, g = list(metric, theta)) |>
  na_omit()

# Cell-level MA decline
market_access_all_dt |>
  roworder(metric, theta, cell) |>
  tfmv(c(GDP_PPP, total_wealth), fgrowth, 
       g = list(metric, theta, cell), apply = FALSE) |>
  na_omit() |> # descr( ~ metric + theta)
  join(centroids, on = "cell", drop = TRUE) |>
  fsubset(theta == "theta_1") |>
  fmutate(metric = set_attr(metric, "levels", recode_char(levels(metric), ttime = "Travel Time", tdist = "Road Distance"))) |>
  
  ggplot(aes(x = lon, y = lat, fill = GDP_PPP)) +
  ggstar::geom_star(starshape = "hexagon", size = 1.1, color = NA, angle = 30) +
  coord_fixed() +
  facet_wrap(~ metric) + # facet_grid(metric ~ theta) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 1)) +
  scale_fill_viridis(option = "magma", n.breaks = 8) +
  labs(fill = NULL) +
  theme_void() +
  theme(strip.text = element_text(size = 14), 
        legend.key.height = unit(2, "cm"), 
        legend.key.width = unit(2, "mm"), 
        plot.margin = margin())

ggsave("figures/full_network/gdp_per_km_MA_loss_perc_theta_1.pdf", width = 11, height = 6)




