################################
# Cost of Trading Across Borders
################################

library(fastverse)
set_collapse(mask = "manip")
fastverse_extend(africamonitor, tmap, sf, qs, install = TRUE) # tmap >= 3.99

africa_ctry_shp <- qread("data/other_inputs/africa_countries.qs")
  
################################
# Doing Business 
################################

# Using Doing Business Estimates
SERIES <- am_series()
exc <- SERIES[startsWith(Topic, "25"), Series]
exc_data <- am_data(series = exc, from = 2019, to = 2020) |> group_by(ISO3) |> flast()
exc_data %>% num_vars() %>% pwcor()

# Add transit time whenever more than 2 countries involved
time_transit <- 48 * (60 / 4) # time_border_compliance %>% num_vars() %>% fmedian() %>% fmean() %>% divide_by(2)
dist_transit <- 48 * 4 * (1000 / 1.6) # Multiply by 1000 since meters

# Getting export and import time vectors
africa_dist <- qread("data/full_network/africa_full_distance_matrix_r9_adjusted.qs")
iso3c <- unique(africa_dist$centroids$ISO3)
export_time_cost <- exc_data %$% set_names(IC_EXP_TMBC * (60 / 4) + IC_EXP_TMDC * (60 / 10), ISO3) %>% extract(iso3c) %>% set_names(iso3c)
import_time_cost <- exc_data %$% set_names(IC_IMP_TMBC * (60 / 4) + IC_IMP_TMDC * (60 / 10), ISO3) %>% extract(iso3c) %>% set_names(iso3c)
export_dist_cost <- exc_data %$% set_names(IC_EXP_CSBC_CD * (1000 / 1.6) + IC_EXP_CSDC_CD * (1000 / 1.6), ISO3) %>% extract(iso3c) %>% set_names(iso3c)
import_dist_cost <- exc_data %$% set_names(IC_IMP_CSBC_CD * (1000 / 1.6) + IC_IMP_CSDC_CD * (1000 / 1.6), ISO3) %>% extract(iso3c) %>% set_names(iso3c)

# Imputation
export_time_cost["ESH"] <- export_time_cost["MAR"]
import_time_cost["ESH"] <- import_time_cost["MAR"]
export_time_cost["ERI"] <- mean(export_time_cost[c("ETH", "SDN")]) 
import_time_cost["ERI"] <- mean(import_time_cost[c("ETH", "SDN")]) 

export_dist_cost["ESH"] <- export_dist_cost["MAR"]
import_dist_cost["ESH"] <- import_dist_cost["MAR"]
export_dist_cost["ERI"] <- mean(export_dist_cost[c("ETH", "SDN")]) 
import_dist_cost["ERI"] <- mean(import_dist_cost[c("ETH", "SDN")]) 

# CEPII GeoDist: Do countries have contiguous border?
geodist <- fread("data/other_inputs/africa_contig.csv")
geodist_mat <- psmat(geodist, contig ~ iso3_o, ~ iso3_d)[iso3c, iso3c]
diag(geodist_mat) <- 1

# Now creating trade border time matrix
border_time_mat <- copyv(geodist_mat, 1, 0) %c+% export_time_cost %r+% import_time_cost
diag(border_time_mat) <- 0
qsu(border_time_mat[border_time_mat > 0])
border_time_mat_transit <- border_time_mat
border_time_mat_transit[geodist_mat == 0] <- border_time_mat_transit[geodist_mat == 0] + time_transit
qsu(border_time_mat_transit[border_time_mat_transit > 0])

# Save as csv 
border_time_mat %>% qDF("iso3c") %>% 
  fwrite("data/QSE/model_border_time_mat.csv")
border_time_mat_transit %>% qDF("iso3c") %>% 
  fwrite("data/QSE/model_border_time_mat_transit.csv")

# Now creating trade border distance matrix
border_dist_mat <- copyv(geodist_mat, 1, 0) %c+% export_dist_cost %r+% import_dist_cost
diag(border_dist_mat) <- 0
qsu(border_dist_mat[border_dist_mat > 0])
border_dist_mat_transit <- border_dist_mat
border_dist_mat_transit[geodist_mat == 0] <- border_dist_mat_transit[geodist_mat == 0] + dist_transit
qsu(border_dist_mat_transit[border_dist_mat_transit > 0])

# Save as csv 
border_dist_mat %>% qDF("iso3c") %>% 
  fwrite("data/QSE/model_border_dist_mat.csv")
border_dist_mat_transit %>% qDF("iso3c") %>% 
  fwrite("data/QSE/model_border_dist_mat_transit.csv")


#####################################################
# Trading Across Borders Summary Table (Tables 3 & 4)
#####################################################

SERIES = am_series()
exc = SERIES[startsWith(Topic, "25"), Series]
exc_data = am_data(ctry = NULL, series = exc, from = 2019, to = 2020)

exc_table = exc_data |> select(-Date) |> 
  join(wlddev, on = c("ISO3" = "iso3c")) |> 
  group_by(region) |> gvr("IC_") |> fmean() 

# By Region (selected)
exc_table |> subset(region %ilike% "Europe|South|Latin|Middle|Sub-") |> 
  xtable::xtable(digits = 1) |> print(include.rownames = FALSE)

# By OECD Membership
exc_data |> select(-Date) |> 
  join(wlddev, on = c("ISO3" = "iso3c")) |> 
  group_by(OECD) |> gvr("IC_") |> fmean() |> 
  xtable::xtable(digits = 1) |> print(include.rownames = FALSE)

# By Income Status
exc_data |> select(-Date) |> 
  join(wlddev, on = c("ISO3" = "iso3c")) |> 
  group_by(income) |> gvr("IC_") |> fmean() |> 
  xtable::xtable(digits = 1) |> print(include.rownames = FALSE)

# GDP weighted all Africa Average
gdp_ppp = am_data(ctry = NULL, series = "NY_GDP_MKTP_PP_CD", from = 2017, to = 2019) |> 
          group_by(ISO3) |> flast()

exc_data |> 
  join(gdp_ppp, on = "ISO3", drop = TRUE) |> 
  subset(ISO3 %in% am_countries$ISO3) |>
  mutate(Africa = 1L) |> 
  group_by(Africa) |> gvr("IC_|_PP_") |> 
  fmean(NY_GDP_MKTP_PP_CD, keep.w = FALSE) |> 
  xtable::xtable(digits = 1) |> print(include.rownames = FALSE)

# Manually getting domestic transport information for 2020 from: 
# https://archive.doingbusiness.org/en/data/exploreeconomies

dom_trans <- fread("data/other_inputs/DB20_domestic_transport.csv")

# Average travel speed
dom_trans |> 
  transform(speed = distance_km / time_hours) |> 
  collap(speed ~ variable, list(mean = fmean, median = fmedian))

subset(dom_trans, variable == "Import") %$% plot(x = time_hours, y = cost_usd)

lm(cost_usd ~ time_hours + variable, dom_trans)
lm(cost_usd ~ distance_km + variable, dom_trans)

# Slightly better
lm(cost_usd ~ time_hours + log(divide_by(distance_km, time_hours)) + variable, dom_trans) |> summary()
lm(cost_usd ~ distance_km + log(divide_by(distance_km, time_hours)) + variable, dom_trans) |> summary()

# Adding Interaction term (heterogeneous slopes)
lm(cost_usd ~ time_hours * variable + log(divide_by(distance_km, time_hours)), dom_trans) |> summary()
lm(cost_usd ~ distance_km * variable + log(divide_by(distance_km, time_hours)), dom_trans) |> summary()
# -> insignificant

# Best fit but interpretation is worse
lm(log(cost_usd) ~ log(time_hours) + log(divide_by(distance_km, time_hours)) + variable, dom_trans) |> summary()

# Both
model <- lm(cost_usd ~ time_hours + distance_km + variable, dom_trans)
summary(model)
# plot(model)
# This is bad: holding distance fixed is not useful...
model <- robustbase::lmrob(cost_usd ~ time_hours + distance_km + variable, dom_trans)
coef(model)
model <- robustbase::lmrob(cost_usd ~ time_hours + variable, dom_trans)
coef(model)
model <- robustbase::lmrob(cost_usd ~ distance_km + variable, dom_trans)
coef(model)
model <- robustbase::lmrob(cost_usd ~ time_hours * variable + log(divide_by(distance_km, time_hours)), dom_trans)
coef(model)
model <- robustbase::lmrob(cost_usd ~ distance_km * variable + log(divide_by(distance_km, time_hours)), dom_trans)
coef(model)

# Summary: Price is 15 USD per hour and 1.5 USD per km. 
library(fixest)
m1 <- feols(cost_usd ~ time_hours + log(divide_by(distance_km, time_hours)) + variable, dom_trans, vcov = "hetero")
m2 <- feols(cost_usd ~ distance_km + log(divide_by(distance_km, time_hours)) + variable, dom_trans, vcov = "hetero")
esttable(m1, m2)
esttex(m1, m2)

# Now running similar regressions for border/documentary compliance
exc_data |> namlab()

exc_data_long <- exc_data |> pivot(1:2, labels = TRUE) |> 
  transform(type = iif(label %like% "Export", "Export", "Import"),
      compliance = iif(label %like% "Border", "Border", "Documentary"),
         measure = iif(label %like% "Time", "Time", "Cost")) |> 
  pivot(c("ISO3", "Date", "type"), "value", names = c("compliance", "measure"), how = "w") |> 
  rename(type = variable)

exc_data_long |> 
  subset(ISO3 %in% am_countries$ISO3) |> 
  transform(Border_Ratio = Border_Cost / Border_Time, 
            Documentary_Ratio = Documentary_Cost / Documentary_Time) |> 
  group_by(variable) |> num_vars() |> fmean()

m3 <- feols(Border_Cost ~ Border_Time + variable, subset(exc_data_long, ISO3 %in% am_countries$ISO3), vcov = "hetero")
m4 <- feols(Documentary_Cost ~ Documentary_Time + variable, subset(exc_data_long, ISO3 %in% am_countries$ISO3), vcov = "hetero")

# This is Table 4
esttable(m1, m2, m3, m4)
esttex(m1, m2, m3, m4) 


######################################
# Trade tariffs (Table 2, Figure 7)
######################################

tt = SERIES[startsWith(Topic, "24"), Series]
tt_data = am_data(ctry = NULL, series = tt, from = 2018, to = 2021) |> group_by(ISO3) |> flast()

tt_table = tt_data |> select(-Date) |> 
  join(wlddev, on = c("ISO3" = "iso3c")) |> 
  group_by(region) |> gvr("_MRCH_|_MANF_") |> fmean() 

tt_table |> namlab()

# By Region (selected)
tt_table  |> subset(region %ilike% "Europe|South|Latin|Middle|Sub-") |> 
  xtable::xtable(digits = 1) |> print(include.rownames = FALSE)

# By OECD Membership
tt_data |> select(-Date) |> 
  join(wlddev, on = c("ISO3" = "iso3c")) |> 
  group_by(OECD) |> gvr("_MRCH_|_MANF_") |> fmean() |> 
  xtable::xtable(digits = 1) |> print(include.rownames = FALSE)

# By Income Status
tt_data |> select(-Date) |> 
  join(wlddev, on = c("ISO3" = "iso3c")) |> 
  group_by(income) |> gvr("_MRCH_|_MANF_") |> fmean() |> 
  xtable::xtable(digits = 1) |> print(include.rownames = FALSE)

# GDP weighted all Africa Average
gdp_ppp = am_data(ctry = NULL, series = "NY_GDP_MKTP_PP_CD", from = 2018, to = 2021) |> group_by(ISO3) |> flast()

tt_data |> 
  join(gdp_ppp, on = "ISO3", drop = TRUE) |> 
  subset(ISO3 %in% am_countries$ISO3) |>
  mutate(Africa = 1L) |> 
  group_by(Africa) |> gvr("_MRCH_|_MANF_|_PP_") |> 
  fmean(NY_GDP_MKTP_PP_CD, keep.w = FALSE) |> 
  xtable::xtable(digits = 1) |> print(include.rownames = FALSE)

# Plot
africa_ctry_shp %<>% join(
  tt_data |> 
    subset(order(ISO3, Date), ISO3, Date, TM_TAX_MRCH_WM_AR_ZS) |>
    subset(Date == flast(Date, ISO3, "fill") & ISO3 %in% am_countries$ISO3), 
  on = "ISO3"
)

pdf("figures/WDI_weighted_average_tariff_map.pdf", width = 8, height = 8)
# tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
tm_shape(africa_ctry_shp) + 
  tm_polygons(fill = "TM_TAX_MRCH_WM_AR_ZS", 
              fill.scale = tm_scale_intervals(5, values = "yl_or_rd", label.na = FALSE), # , value.na = "grey40"
              fill.legend = tm_legend("Tariff (%)", position = c("left", "bottom"), frame = FALSE, bg.color = NA,
                                      text.size = 1.5, title.size = 2),
              lwd = 2) +
  tm_layout(frame = FALSE) 
dev.off()


######################################
# Non-Tariff Measures
######################################

# WIIW NTM Database
NTM_AGG <- qread("data/other_inputs/WIIW_NTM_AGG.qs")
frange(NTM_AGG$NTM_AGG$year)

NTM_AGG$NTM_AGG |> subset(year == flast(year, imp_iso3, "fill") & aff_iso3 == "WTO") |>
  group_by(Africa = imp_iso3 %in% am_countries$ISO3) |> select(N) |> fmean()

NTM_AGG$NTM_AGG |> subset(year == flast(year, imp_iso3, "fill") & aff_iso3 == "WTO") |>
  group_by(Africa = imp_iso3 %in% am_countries$ISO3) |> 
  summarise(.data[topn(N)])

# Plot
africa_ctry_shp %<>% join(
  NTM_AGG$NTM_AGG |> 
    subset(year == flast(year, imp_iso3, "fill")) |>
    fcount(imp_iso3, w = N), 
  on = c("ISO3" = "imp_iso3")
)

pdf("figures/WIIW_NTMs_total_map.pdf", width = 8, height = 8)
# tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
tm_shape(africa_ctry_shp) + 
  tm_polygons(fill = "N", 
              fill.scale = tm_scale_intervals(5, breaks = c(0, 2, 10, 50, 250, 1000, Inf), values = "yl_or_rd", label.na = FALSE), # , value.na = "grey40"
              fill.legend = tm_legend("# NTMs", position = c("left", "bottom"), frame = FALSE, bg.color = NA,
                                      text.size = 1.5, title.size = 2),
              lwd = 2) +
  tm_layout(frame = FALSE) 
dev.off()


