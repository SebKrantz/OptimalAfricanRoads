#####################################################################
# PIDA: manual edge flag, merged OSBP + border points, map export
# Run from project root. Uses collapse (set_collapse mask); no dplyr.
#####################################################################

library(fastverse)
set_collapse(mask = c("manip", "helper", "special"))
fastverse_extend(qs2, sf, units, tmap, install = TRUE)
source("code/helpers/helpers.R")
fastverse_conflicts()

# -------------------------------------------------------------------
# Load network and PIDA inputs
# -------------------------------------------------------------------
load("data/transport_network/trans_africa_network_param.RData")
edges_real <- qs_read("data/transport_network/edges_real_simplified.qs2") |>
  select(from, to) |>
  rmapshaper::ms_simplify(keep = 0.06) |>
  st_make_valid()

PIDA_OSBP <- st_read("data/PIDA/PIDA_Transport_Projects_Harmonized.gpkg", layer = "Border_Posts_OSBP")

PIDA_alt <- jsonlite::fromJSON("data/PIDA/PIDA_Map_Export_2026-04-18.json")$projects |>
  subset(!is.na(geometry_type))
PIDA_alt$location_map <- lapply(PIDA_alt$location_map, function(x) {
  tryCatch(st_make_valid(st_as_sfc(x, crs = 4326)), error = function(e) NULL)
})
PIDA_alt <- subset(PIDA_alt, vlengths(location_map) > 0)
PIDA_alt <- ftransform(PIDA_alt, location_map = do.call(c, location_map))
PIDA_alt <- st_sf(PIDA_alt)

PIDA_border <- subset(PIDA_alt, subsector %ilike% "border" & geometry_type %ilike% "point")

# -------------------------------------------------------------------
# Edge table with simplified real geometries
# -------------------------------------------------------------------
edges_all_param <- rowbind(
  edges_param,
  rename(
    add_links_param,
    duration_65kmh = duration,
    duration_100kmh = duration_imp,
    total_time_65kmh = total_time,
    total_time_100kmh = total_time_imp
  ) |>
    select(-id, -total_dist),
  fill = TRUE
)
edges_all_param$geometry[seq_row(edges_real)] <- edges_real$geometry

# -------------------------------------------------------------------
# PIDA column from manual edge index list (authoritative)
# -------------------------------------------------------------------
PIDA_ind <- c(
  451, 450, 449, 448, 457, 471, 472, 452, 454, 509, 476, 432, 433, # Algeria <> Niger
  978, 970, 968, 971, 1007, 746, 745, 808, 928, 1005, 1018, # Lybia <> Niger
  6, 1071, 1292, 2621, 1408, 1462, 1572, 1713, 1912, 1958, # Chad <> Sudan
  1070, 1093, 1114, 1134, 1150, 1151, 1169, 1102, 1129, 1174, # CAR and Cameroon
  1039, 1038, 1024, 994, 993, 1002, 1008, 1010, 976, 920, 918,
  919, 927, 1015, 820, 806, 791, 1109, 1095, 1054, 1026, 1027, # More CEMAC
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

edges_all_param$PIDA <- "No"
edges_all_param$PIDA[PIDA_ind] <- "Yes"

# -------------------------------------------------------------------
# Merge harmonized OSBP (excluding wrong rows) with distant map-export borders
# -------------------------------------------------------------------
OSBP_DROP <- c(18L, 20L)
PIDA_OSBP_kept <- PIDA_OSBP[-OSBP_DROP, ]

borders <- lapply(PIDA_border$location_map, st_make_valid) |> st_as_sfc(crs = 4326)
bdist <- st_distance(PIDA_OSBP_kept, borders)
min_d <- apply(as.matrix(drop_units(bdist)), 1L, min)
missing <- min_d > 100

PIDA_border_extra <- PIDA_border[missing, ]
if (!identical(sf::st_crs(PIDA_OSBP_kept), sf::st_crs(PIDA_border_extra))) {
  sf::st_crs(PIDA_border_extra) <- sf::st_crs(PIDA_OSBP_kept)
}
# Harmonized OSBP (gpkg) vs map-export border points (JSON -> sf):
#   same: name, pida_code, sector, status, countries
#   raw_budget (numeric) -> capex_usd_millions
#   subsector -> mode ("Border Post" in both samples)
#   no route flags in export -> has_route = FALSE, route_segments = 0L
#   active geometry column location_map -> geom (collapse select(geom = ...) drops sf_column)
PIDA_border_osbp_aligned <- select(
  PIDA_border_extra,
  name,
  pida_code,
  sector,
  status,
  capex_usd_millions = raw_budget,
  countries,
  mode = subsector,
  location_map
)
PIDA_border_osbp_aligned <- ftransform(
  PIDA_border_osbp_aligned,
  has_route = FALSE,
  route_segments = 0L
)
PIDA_border_osbp_aligned <- rename(PIDA_border_osbp_aligned, geom = location_map)
attr(PIDA_border_osbp_aligned, "sf_column") <- "geom"
PIDA_border_osbp_aligned <- select(
  PIDA_border_osbp_aligned,
  name,
  pida_code,
  sector,
  status,
  capex_usd_millions,
  countries,
  has_route,
  route_segments,
  mode,
  geom
)
PIDA_OSBP_all <- rowbind(PIDA_OSBP_kept, PIDA_border_osbp_aligned)

merged_path <- "data/PIDA/PIDA_OSBP_border_merged.qs2"
qs_save(PIDA_OSBP_all, merged_path)
message("Wrote merged OSBP + border points: ", merged_path)

# -------------------------------------------------------------------
# Map: PIDA edges and nodes
# -------------------------------------------------------------------
dir.create("figures/transport_network/PIDA", recursive = TRUE, showWarnings = FALSE)
fig_path <- "figures/transport_network/PIDA/trans_africa_edges_PIDA_manual.pdf"

PIDA_OSBP_map <- ftransform(PIDA_OSBP_all, pid_osbp_leg = "OSBP")

pid_map <- tm_basemap("Esri.WorldGrayCanvas", zoom = 4) +
  tm_shape(edges_all_param) +
  tm_lines(
    col = "PIDA",
    col.scale = tm_scale_categorical(values = c("grey", "orange")),
    col.legend = tm_legend(
      "PIDA Link",
      position = c("left", "bottom"),
      frame = FALSE,
      text.size = 1.5,
      title.size = 2
    ),
    lwd = 2
  ) +
  tm_shape(subset(nodes, population > 0)) +
  tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) +
  tm_dots(size = 0.1, fill = "grey30") +
  tm_shape(PIDA_OSBP_map) +
  tm_dots(
    size = 0.2,
    fill = "pid_osbp_leg",
    col = "red",
    fill.scale = tm_scale_categorical(values = "red"),
    fill.legend = tm_legend(
      "Border posts",
      position = c("left", "bottom"),
      frame = FALSE,
      text.size = 1.5,
      title.size = 2
    )
  ) +
  tm_layout(frame = FALSE)

tmap_save(pid_map, fig_path, width = 8, height = 8, units = "in")
message("Wrote map: ", fig_path)
