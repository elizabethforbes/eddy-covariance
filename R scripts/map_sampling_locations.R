# ============================================================
# Map: Two-panel figure — Organic & Conventional sites
# ============================================================
# Inputs (all in working directory):
#   - Shape_File.shp                 field polygons (EPSG:32115)
#   - 2020_soilcoring_locations.csv  soil cores (EPSG:32115)
#   - EC_tower_locations.xlsx        EC towers (WGS84 lat/lon)
#   - C-5.shp                        all 5 conventional collars (WGS84)
#   - O-4.shp                        all 4 organic collars (WGS84)
#
# NOTE: collar .shp files are missing companion files (.shx/.dbf/.prj).
#   SHAPE_RESTORE_SHX handles the missing index; CRS is assigned
#   manually as WGS84 (EPSG:4326) — confirmed from coordinates.
#
# NOTE on field IDs: the field_ID column in the CSV has EF01/EF02
#   labels that are spatially swapped relative to the EC tower xlsx.
#   The management column is authoritative; both datasets are grouped
#   by spatial location (organic site ≈ 42.166°N; conventional ≈ 42.116°N).
#
# The conventional field polygon is NOT in Shape_File.shp, so the
#   conventional panel shows only the point layers on a blank background.
# ============================================================

library(sf)
library(ggplot2)
library(dplyr)
library(readxl)
library(ggspatial)
library(patchwork)   # side-by-side panels

library(rnaturalearth)
library(rnaturalearthdata)


Sys.setenv(SHAPE_RESTORE_SHX = "YES")


################################################################################
# load and re-project all the different layers: towers, collars, and cores
################################################################################

setwd('/Users/elizabethforbes/Documents/Hudson Carbon/eddy covariance')

################################################################################
# cores
################################################################################

# List of core location shapefile paths (adjust paths as needed)
shapefile_paths <- c(
  "/Users/elizabethforbes/Documents/Hudson Carbon/eddy covariance/eddy_covariance_fluxdata/fall 2020 core sampling sites/conv_Nov2019_1.shp",
  "/Users/elizabethforbes/Documents/Hudson Carbon/eddy covariance/eddy_covariance_fluxdata/fall 2020 core sampling sites/conv_Nov2019_2.shp",
  "/Users/elizabethforbes/Documents/Hudson Carbon/eddy covariance/eddy_covariance_fluxdata/fall 2020 core sampling sites/conv_Nov2019_3.shp",
  "/Users/elizabethforbes/Documents/Hudson Carbon/eddy covariance/eddy_covariance_fluxdata/fall 2020 core sampling sites/org_Nov2019_1.shp",
  "/Users/elizabethforbes/Documents/Hudson Carbon/eddy covariance/eddy_covariance_fluxdata/fall 2020 core sampling sites/org_Nov2019_2.shp",
  "/Users/elizabethforbes/Documents/Hudson Carbon/eddy covariance/eddy_covariance_fluxdata/fall 2020 core sampling sites/org_Nov2019_3.shp",
  "/Users/elizabethforbes/Documents/Hudson Carbon/eddy covariance/eddy_covariance_fluxdata/fall 2020 core sampling sites/conv_Oct2018_1.shp",
  "/Users/elizabethforbes/Documents/Hudson Carbon/eddy covariance/eddy_covariance_fluxdata/fall 2020 core sampling sites/conv_Oct2018_2.shp",
  "/Users/elizabethforbes/Documents/Hudson Carbon/eddy covariance/eddy_covariance_fluxdata/fall 2020 core sampling sites/conv_Oct2018_3.shp",
  "/Users/elizabethforbes/Documents/Hudson Carbon/eddy covariance/eddy_covariance_fluxdata/fall 2020 core sampling sites/org_Oct2018_1.shp",
  "/Users/elizabethforbes/Documents/Hudson Carbon/eddy covariance/eddy_covariance_fluxdata/fall 2020 core sampling sites/org_Oct2018_2.shp",
  "/Users/elizabethforbes/Documents/Hudson Carbon/eddy covariance/eddy_covariance_fluxdata/fall 2020 core sampling sites/org_Oct2018_3.shp"
)

# Read all shapefiles into a list
sf_list <- lapply(shapefile_paths, st_read)

# combine all core locations into a list
cores_sf <- do.call(rbind, sf_list)

# add management, year columns (tho reminder these were all collected in Nov. 2020):
cores_sf <- cores_sf %>% 
  mutate(year = case_when(
    grepl("2019", layer) ~ "2019",
    TRUE ~ "2018"
  )) %>% 
  mutate(management = case_when(
    st_coordinates(st_centroid(geometry))[, "Y"] < 42.15 ~ "conventional",
    st_coordinates(st_centroid(geometry))[, "Y"] > 42.15 ~ "organic"
  ))

cores_sf_nyeast <- st_transform(cores_sf, crs = 4326)  # Transform to NAD83 / New York East

################################################################################
# EC towers
################################################################################

towers_sf <- read_excel("eddy_covariance_fluxdata/EC_tower_locations.xlsx")
towers_sf <- towers_sf %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

################################################################################
# collars
################################################################################

# Collars — read all files per treatment, combine, and deduplicate.
# The files are overlapping (each was saved with prior points still
# included), so deduplication by rounded coordinate is necessary.
# O-2 + O-3 together give all 5 unique organic collars;
# C-5 alone has all 5 conventional, but we read all for safety.

read_collars <- function(files, mgmt) {
  lapply(files, function(f) {
    st_read(f, quiet = TRUE) %>% st_set_crs(4326)
  }) %>%
    bind_rows() %>%
    mutate(lon_r = round(st_coordinates(.)[, 1], 5),
           lat_r = round(st_coordinates(.)[, 2], 5)) %>%
    distinct(lon_r, lat_r, .keep_all = TRUE) %>%
    select(-lon_r, -lat_r) %>%
    st_transform(target_crs) %>%
    mutate(management = mgmt)
}

conv_collars <- read_collars(
  c("additional_data/GHG Flux Rings/C-1.shp",
    "additional_data/GHG Flux Rings/C-2.shp",
    "additional_data/GHG Flux Rings/C-3.shp",
    "additional_data/GHG Flux Rings/C-4.shp",
    "additional_data/GHG Flux Rings/C-5.shp"),
  "conventional"
)

org_collars <- read_collars(
  c("additional_data/GHG Flux Rings/O-2.shp",
    "additional_data/GHG Flux Rings/O-3.shp",
    "additional_data/GHG Flux Rings/O-4.shp",
    "additional_data/GHG Flux Rings/O-5.shp"),
  "organic"
)

conv_collars_sf <- st_transform(conv_collars, crs = 4326) 
org_collars_sf <- st_transform(org_collars, crs = 4326)

################################################################################
# create shared visual settings
################################################################################

field_fill   <- "#b2df8a"   # soft green for field polygon
collar_col   <- "#7570b3"   # purple  — collars
core_col     <- "#d95f02"   # orange  — soil cores
tower_col    <- "#1b9e77"   # teal    — EC tower

# ── make map panels, then stitch together ───────────────────────────────

# function to reproject points relative to their tower
make_relative <- function(pts_sf, tower_sf) {
  # get tower coordinates as the origin
  tower_coords <- st_coordinates(tower_sf)
  tower_lon    <- tower_coords[1]
  tower_lat    <- tower_coords[2]
  
  # create a local azimuthal equidistant projection centred on the tower
  # this gives metres N/S and E/W from the tower
  local_crs <- sprintf(
    "+proj=aeqd +lat_0=%f +lon_0=%f +units=m",
    tower_lat, tower_lon
  )
  
  st_transform(pts_sf, crs = local_crs)
}

# reproject everything to local coordinates for each field
org_tower    <- towers_sf %>% filter(management == "organic")
conv_tower   <- towers_sf %>% filter(management == "conventional")

org_cores_rel    <- make_relative(cores_sf_nyeast %>% filter(management == "organic"),   org_tower)
org_collars_rel  <- make_relative(org_collars_sf,   org_tower)
org_tower_rel    <- make_relative(org_tower,         org_tower)

conv_cores_rel   <- make_relative(cores_sf_nyeast %>% filter(management == "conventional"), conv_tower)
conv_collars_rel <- make_relative(conv_collars_sf,   conv_tower)
conv_tower_rel   <- make_relative(conv_tower,         conv_tower)

make_panel <- function(panel_title,
                       cores_data,
                       collars_data,
                       tower_data,
                       x_range_m,   # shared geographic window in metres
                       y_range_m,
                       show_legend = FALSE) {
  
  # centre on tower (which is at 0,0 in local CRS)
  xlim <- c(-x_range_m / 2, x_range_m / 2)
  ylim <- c(-y_range_m / 2, y_range_m / 2)
  
  aspect_ratio <- diff(xlim) / diff(ylim)
  
  p <- ggplot() +
    
    geom_sf(data = cores_data,
            colour = core_col, fill = core_col,
            shape = 21, size = 1) +
    
    geom_sf(data = collars_data,
            colour = collar_col,
            shape = 1, size = 1, stroke = 1.0) +
    
    geom_sf(data = tower_data,
            colour = tower_col,
            shape = 8, size = 3, stroke = 1.0) +
    
    # axis labels in metres from tower
    scale_x_continuous(labels = function(x) paste0(x, " m")) +
    scale_y_continuous(labels = function(x) paste0(x, " m")) +
    
    annotate("text",
             x = xlim[1] + 0.05 * diff(xlim),
             y = ylim[2] - 0.05 * diff(ylim),
             label = panel_title,
             hjust = 0, vjust = 1,
             size = 4, fontface = "bold") +
    
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    annotation_scale(location = "br", width_hint = 0.3, text_cex = 0.7) +
    
    theme_bw(base_size = 11) +
    theme(
      axis.text        = element_text(size = 7),
      axis.title       = element_blank(),
      panel.grid.major = element_line(colour = "grey85", linewidth = 0.25),
      legend.position  = if (show_legend) "bottom" else "none",
      aspect.ratio     = aspect_ratio
    )
  
  p
}

# calculate shared window size from whichever field is more spread out
all_org_coords  <- st_coordinates(bind_rows(org_cores_rel, 
                                            org_collars_rel, 
                                            org_tower_rel))
all_conv_coords <- st_coordinates(bind_rows(conv_cores_rel,
                                            conv_collars_rel,
                                            conv_tower_rel))

x_range_m <- max(diff(range(all_org_coords[,1])),
                 diff(range(all_conv_coords[,1]))) * 1.5  # 1.5x for padding
y_range_m <- max(diff(range(all_org_coords[,2])),
                 diff(range(all_conv_coords[,2]))) * 1.5

# build panels
org <- make_panel(
  panel_title  = "Organic field",
  cores_data   = org_cores_rel,
  collars_data = org_collars_rel,
  tower_data   = org_tower_rel,
  x_range_m    = x_range_m,
  y_range_m    = y_range_m
) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

conv <- make_panel(
  panel_title  = "Conventional field",
  cores_data   = conv_cores_rel,
  collars_data = conv_collars_rel,
  tower_data   = conv_tower_rel,
  x_range_m    = x_range_m,
  y_range_m    = y_range_m
) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

combined_plot <- org + conv
print(combined_plot)


