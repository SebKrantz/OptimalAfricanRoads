###########################################
# Obtaining Distance and Travel Time Matrix
###########################################

# OPEN SOURCE ROUTING MACHINE: http://project-osrm.org/
# High-Performance Routing Service based on OSM
# https://map.project-osrm.org
# in R: install.packages("osrm")

library(fastverse)
fastverse_extend(osrm, qs, install = TRUE)

# ------------------------------------------------------------------
# Option 1: Run your own OSRM Server (Fastest, but need 50-60GB RAM)
# ------------------------------------------------------------------

# OSRM docker container
# Probably a memory issue with Africa: https://stackoverflow.com/questions/53157542/osrm-extract-silently-fails
# https://github.com/Project-OSRM/osrm-backend#using-docker

# Executed on Ubuntu server with 64Gb RAM:  ------------------------------

# See avialable RAM: 
# free -m

# (1) Install Docker --------
# sudo apt-get update
# 
# sudo apt-get install \
# ca-certificates \
# curl \
# gnupg \
# lsb-release
# 
# curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
#
# echo \
# "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu \
#   $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null 
#
# sudo apt-get update
# sudo apt-get install docker-ce docker-ce-cli containerd.io
# sudo systemctl status docker
#
# Add your user to the docker group
# sudo usermod -aG docker ${USER}
# Logout nd back in and verify: 
# docker ps

# (2) Start Docker, get Africa OSM, start OSRM Server
# sudo systemctl start docker.socket
# sudo systemctl start docker

# Get Data 
# sudo chown -R krantz:krantz /srv/shiny-server
# wget https://download.geofabrik.de/africa-latest.osm.pbf

# docker run -t -v "${PWD}:/data" ghcr.io/project-osrm/osrm-backend osrm-extract -p /opt/car.lua /data/africa-latest.osm.pbf || echo "osrm-extract failed"

# https://github.com/Project-OSRM/osrm-backend/issues/5022 # increase table size
# https://github.com/Project-OSRM/osrm-backend/issues/5500 # CH is better than MLD for large table queries
# see options with osrm-routed --help

# Setup with MLD
# docker run -t -v "${PWD}:/data" ghcr.io/project-osrm/osrm-backend osrm-partition /data/africa-latest.osrm || echo "osrm-partition failed"
# docker run -t -v "${PWD}:/data" ghcr.io/project-osrm/osrm-backend osrm-customize /data/africa-latest.osrm || echo "osrm-customize failed"
# docker run -t -i -p 5000:5000 -v "${PWD}:/data" ghcr.io/project-osrm/osrm-backend osrm-routed --algorithm mld --max-table-size 15000 /data/africa-latest.osrm

# # Setup with CH
# # docker run -t -v "${PWD}:/data" ghcr.io/project-osrm/osrm-backend osrm-contract /data/africa-latest.osrm || echo "osrm-contract failed"
# # docker run -t -i -p 5000:5000 -v "${PWD}:/data" ghcr.io/project-osrm/osrm-backend osrm-routed --algorithm ch --max-table-size 15000 /data/africa-latest.osrm

# Stop server with ctrl + C

# Set maximum header size to a large number !!
# export CURL_MAX_HTTP_HEADER=10240000000

# (3) Open a new terminal, run R, and set 
options(osrm.server = "http://0.0.0.0:5000/",
        osrm.profile = "car")

getOption("osrm.server")

# Testing
cape_town_to_cairo = rbind(data.frame(lon = 18.38030253247056, lat = -34.03336853793626), 
                           data.frame(lon = 31.425366291541163, lat = 30.187212885598868))

madagascar_to_mainland = rbind(data.frame(lon = 48.53008716640623, lat = -18.026478958576668), 
                               data.frame(lon = 38.246318779991775, lat = -15.049170402887114))

nowhere_to_cape_town = rbind(data.frame(lon = -6.5, lat = 23),
                             data.frame(lon = 18.38030253247056, lat = -34.03336853793626))

osrmTable(loc = cape_town_to_cairo, measure = c('duration', 'distance'))
osrmTable(loc = madagascar_to_mainland, measure = c('duration', 'distance'))
osrmTable(loc = nowhere_to_cape_town, measure = c('duration', 'distance'))

# Load Grid centroids
calib_data <- qDF(fread("data/QSE/QSE_model_calibration_data.csv"))
rownames(calib_data) <- calib_data$cell
source("code/helpers.R")

# Should be only a few seconds with your own server...
result <- split_large_dist_matrix(fselect(calib_data, lon, lat), chunk_size = 3000, verbose = TRUE)

result$centroids <- fselect(calib_data, cell, ISO3, lon, lat, pop_gpw4, pop_wpop)
qsave(result, "data/africa_full_distance_matrix_r9.qs")


# ----------------------------------------------------------------------------------
# Option 2: Query Public OSRM Server ~15000 times (slightly abusive, but far easier)
# ----------------------------------------------------------------------------------

# Load Grid centroids
calib_data <- qDF(fread("data/QSE/QSE_model_calibration_data.csv"))
rownames(calib_data) <- calib_data$cell
source("code/helpers.R")

# Probably takes at least 1 hour...
result <- split_large_dist_matrix(fselect(calib_data, lon, lat), chunk_size = 100, verbose = TRUE)

result$centroids <- fselect(calib_data, cell, ISO3, lon, lat, pop_gpw4, pop_wpop)
qsave(result, "data/africa_full_distance_matrix_r9.qs")


# ----------------------------------------------------------------------------------
# Option 3: Simply Download my Matrix from Google Drive
# ----------------------------------------------------------------------------------

# Link: https://drive.google.com/file/d/1oE_9i3SdqvYKcdl880uS9q774dD9KXmS/view?usp=sharing
