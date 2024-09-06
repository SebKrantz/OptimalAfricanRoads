library(fastverse)
fastverse_extend(qs, dggridR)

wld50kmhex <- dgconstruct(res = 9)

# Example on how model_calibration_data is computed from higher-resolution geospatial data

model_calibration_data <- qread("data/other_inputs/imputed_wealth_GDP_10km_hex.qs") |> 
  fmutate(cell = dgGEO_to_SEQNUM(wld50kmhex, lon, lat)$seqnum) |> 
  fgroup_by(cell) |> 
  fsummarise(across(c(IWI, RWI, GDP_per_capita_PPP, log_GDP_per_capita_PPP), fmean), 
             GDP_PPP = fsum(GDP_PPP),
             log_GDP_PPP = log(fsum(exp(log_GDP_PPP)))) 

# Comparison
fread("data/QSE/model_calibration_data.csv") %>%
  get_vars(intersect(names(.), names(model_calibration_data))) %>%
  join(model_calibration_data, on = "cell", how = "inner", 
       suffix = c("_loaded", "_computed")) %>% # colorderv() %>% View() 
  {pwcor(gvr(., "_loaded"), gvr(., "_computed"))}


# Example of aggregating a raster layer into a discrete global grid

# Global Motorized Friction Surface from https://data.malariaatlas.org
# Direct download link (might not work): https://data.malariaatlas.org/e56db3d3-8cdd-46f0-be11-ab55e462ddb7

# Friction 
MAP_MFS_19 <- terra::rast("/Users/sebastiankrantz/Documents/Data/MAPFrictionSurface2019/202001_Global_Motorized_Friction_Surface_2019/202001_Global_Motorized_Friction_Surface_2019.tif")
# Bounding box 
africa_bbox <- terra::ext(c(xmin = -27, xmax = 59, ymin = -36, ymax = 38))
# Clip raster
MAP_MFS_19_Africa <- terra::crop(MAP_MFS_19, africa_bbox)
plot(MAP_MFS_19_Africa)
# To data frame
MAP_MFS_19_Africa <- as.data.frame(MAP_MFS_19_Africa, xy = TRUE) |> qDT() |> setNames(c("lon", "lat", "rugg"))
# Aggregation
MAP_MFS_19_Africa <- MAP_MFS_19_Africa |> 
  fgroup_by(cell = dgGEO_to_SEQNUM(wld50kmhex, lon, lat)$seqnum) |> 
  fselect(rugg) |> 
  fsummarise(N = GRPN(), rugg = fmean(rugg))
# Descriptive Stats
MAP_MFS_19_Africa |> descr()
# Plot
MAP_MFS_19_Africa |> 
  ftransform(dgSEQNUM_to_GEO(wld50kmhex, cell)) |>
  ggplot(aes(x = lon_deg, y = lat_deg, colour = rugg)) + geom_point() +
  scale_colour_viridis_c(option = "H")
