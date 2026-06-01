#####################################################################
# Transport Network: Analyze PIDA Projects
#####################################################################

library(fastverse)
set_collapse(mask = c("manip", "helper", "special"))
fastverse_extend(qs2, sf, units, mapview, tmap, install = TRUE)
source("code/helpers/helpers.R")
fastverse_conflicts()

# Load essential network data
load("data/transport_network/trans_africa_network_param.RData")
edges_real <- qs_read("data/transport_network/edges_real_simplified.qs2") |> 
  select(from, to) |> rmapshaper::ms_simplify(keep = 0.06) |> st_make_valid()

# Load PIDA Data
st_layers("data/PIDA/PIDA_Transport_Projects_Harmonized.gpkg")
PIDA_lines <- st_read("data/PIDA/PIDA_Transport_Projects_Harmonized.gpkg", layer = "Routes") |> 
  subset(st_geometry_type(geom) %ilike% "LINESTRING" & 
         (name %ilike% "road|highway|freeway|RN06" | pida_code %ilike% "T.21.03.02|T.05.03.03") & 
         !pida_code %ilike% "T.RD.26.")

PIDA_OSBP <- st_read("data/PIDA/PIDA_Transport_Projects_Harmonized.gpkg", layer = "Border_Posts_OSBP")
PIDA_Other <- st_read("data/PIDA/PIDA_Transport_Projects_Harmonized.gpkg", layer = "Other_Points")

# Complement as needed
PIDA_lines_orig <- geojsonsf::geojson_sf("data/PIDA/AID_Africa_Transport_Projects_Routes.geojson") |> 
  subset(st_geometry_type(geometry) %ilike% "LINESTRING" & (name %ilike% "road" | name %like% "TAH"))
# mapview(PIDA_lines) + mapview(PIDA_lines_orig, color = "red")
# View(subset(PIDA_lines_orig, name %ilike% "osbp"))

PIDA_points_orig <- geojsonsf::geojson_sf("data/PIDA/AID_Africa_Transport_Projects_Points.geojson")
# mapview(PIDA_OSBP) + mapview(subset(PIDA_points_orig, name %ilike% "osbp"), color = "red")

# -> Nothing added in either case

# Now examining alternative downloads
PIDA_alt <- jsonlite::fromJSON("data/PIDA/PIDA_Map_Export_2026-04-18.json")$projects |>
  subset(!is.na(geometry_type))
PIDA_alt$location_map <- lapply(PIDA_alt$location_map, function(x) tryCatch(st_make_valid(st_as_sfc(x, crs = 4326)), error = function(e) NULL))
PIDA_alt %<>% subset(vlengths(location_map) > 0) %>% 
  mutate(location_map = do.call(c, location_map)) %>% 
  st_sf()

# mapview(subset(PIDA_lines, pida_code != "")) + mapview(subset(PIDA_alt,subsector %ilike% "road" & geometry_type %ilike% "linestring"), color = "red")
# -> Also nothing new

# View(fread("data/PIDA/pida_projects_export.csv"))
# View(fread("data/PIDA/pida-africa-projects-2026-04-18.csv"))

PIDA_roads <- subset(PIDA_alt,subsector %ilike% "road" & geometry_type %ilike% "linestring")
PIDA_border <- subset(PIDA_alt,subsector %ilike% "border" & geometry_type %ilike% "point")

# Now matching segments and borders ----------------------

edges_all_param <- rowbind(edges_param, 
                        rename(add_links_param, duration_65kmh = duration,
                               duration_100kmh = duration_imp, 
                               total_time_65kmh = total_time,
                               total_time_100kmh = total_time_imp) |> 
                        select(-id, -total_dist), fill = TRUE)
edges_all_param$geometry[seq_row(edges_real)] <- edges_real$geometry

PIDA_roads_segments <- nngeo::st_segments(PIDA_roads)
angle <- abs(stplanr::line_bearing(edges_all_param, bidirectional = TRUE) -
             stplanr::line_bearing(PIDA_roads_segments[st_nearest_feature(edges_all_param, PIDA_roads_segments), ], bidirectional = TRUE))

m <- unique(st_nearest_feature(st_centroid(PIDA_roads_segments), edges_all_param))
edges_all_param$PIDA <- "No"
edges_all_param$PIDA[m] <- "Yes"
mapview(edges_all_param, zcol = "PIDA")

mapview(edges_all_param) + mapview(PIDA_roads, color = "red")
PIDA_ind <- c(451, 450, 449, 448, 457, 471, 472, 452, 454, 509, 476, 432, 433, # Algeria <> Niger
              978, 970, 968, 971, 1007, 746, 745, 808, 928, 1005, 1018, # Lybia <> Niger
              6, 1071, 1292, 2621, 1408, 1462, 1572, 1713, 1912, 1958, # Chad <> Sudan
              1070, 1093, 1114, 1134, 1150, 1151, 1169, 1102, 1129, 1174, # CAR and Cameroon
              1039, 1038, 1024, 994, 993, 1002, 1008, 1010, 976, 920, 918, 
              919, 927, 1015, 820, 806, 791, 1109, 1095, 1054, 1026, 1027,  # More CEMAC
              989, 956, 957, 942, 910, 883, 825, 1194, 1208, 1213, 1238,
              1250, 1253, 1294, 1319, 1338, 1349, 1371, 1398, 1399, 1401,
              1492, 1519, 1612, 1400, 1455, 1554, 1555, 1602, 1637, 1636, 1664, # Through Congo (DRC) and Burundi
              1663, 1675, 2699, 1721, 1722, 1723, 1720, 1792, 1905, 1925, # Rwanda and Tanzania West
              2142, 2185, 1902, 1943, 1962, 1990, 2000, 2750, 2039, 2096, 
              2135, 2163, 2180, 2021, 2022, 2749, # 1923, 
              1984, 1938, 1894, # Kenya and Tanzania East
              1857, 1829, 2775, 2228, 2306, 2301, 2283, 2253, # Ethiopia
              1908, 1896, 1866, 1849, 1821, 1800, 1725, 1536, # Zambia
              1846, 1771, 1709, 1560, 1761, 1729, 1728, 1741, 1740, 1698, # Zimbabwe
              1614, 1580, 1173, 1160, 1147, 1113, # Angola
              466, 430, 406, 381, 363, 341, 334, 309, 300, 292, 2333, 249, 222, # Lagos <> Abidjan
              2549, # Bangui <> Sangha new link
              808, 779, 732, 688, 654, 585, 505, 481, 419, 391, 389, 336, 219, 216, 203, 179, 170, 141, # Niger <> Mali (Eyeballed from paper)
              1828, 2726, 1834, 1897 # Up north from Juba
)

edges_all_param$PIDA[PIDA_ind] <- "Yes"
edges_all_param$PIDA[-PIDA_ind] <- "No"

mapview(edges_all_param[PIDA_ind, ]) + mapview(PIDA_roads, color = "red") + 
  mapview(PIDA_OSBP[-c(18, 20), ]) # -c(18, 20) aare wrong. 

borders <- lapply(PIDA_border$location_map, st_make_valid) |> st_as_sfc(crs = 4326)

bdist <- st_distance(PIDA_OSBP[-c(18, 20), ], borders)
missing <- fmin(bdist) > as_units(100, "m")

mapview(PIDA_OSBP[-c(18, 20), ]) + mapview(PIDA_border[missing, ], color = "red")

names(PIDA_OSBP)
names(PIDA_border)

PIDA_OSBP_all <- # Merge PIDA_OSBP[-c(18, 20), ] and PIDA_border[missing, ] into one file. 


# pdf(sprintf("figures/transport_network/GE_dual/total/trans_africa_network_GE_%s_perc_ug.pdf", res_name), width = 8, height = 8)
tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges_all_param) +
  tm_lines(col = "PIDA", 
           col.scale = tm_scale_categorical(values = c("grey", "orange")),
           col.legend = tm_legend("PIDA Link", 
                                  position = c("left", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  # tm_shape(PIDA_roads) + tm_lines(col = "red") +
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey30") +
  tm_layout(frame = FALSE) 
# dev.off()


