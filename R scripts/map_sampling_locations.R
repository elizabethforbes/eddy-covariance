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

Sys.setenv(SHAPE_RESTORE_SHX = "YES")


# ── 1. Load & reproject all layers ──────────────────────────
setwd('/Users/elizabethforbes/Documents/Hudson Carbon/eddy covariance')

fields <- st_read("eddy_covariance_fluxdata/shape files_fields/Shape_File.shp", quiet = TRUE)   # EPSG:32115
target_crs <- st_crs(fields)

# Soil cores — fix the swapped field_ID by relying on spatial location
cores_raw <- read.csv("eddy_covariance_fluxdata/2020_soilcoring_locations.csv")
cores_raw <- cores_raw[, !grepl("Unnamed", names(cores_raw))]

cores_sf <- st_as_sf(cores_raw, coords = c("X", "Y"), crs = target_crs) %>%
  # Re-assign management by latitude (organic ~42.166°N, conventional ~42.116°N)
  mutate(
    lat_approx = st_coordinates(st_transform(., 4326))[, 2],
    management = if_else(lat_approx > 42.14, "organic", "conventional")
  )

# EC towers
towers_sf <- read_xlsx("eddy_covariance_fluxdata/EC_tower_locations.xlsx") %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(target_crs)

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

message("Conventional collars: ", nrow(conv_collars), " unique points")
message("Organic collars: ",      nrow(org_collars),  " unique points")


# ── 2. Organic field polygon (HF 12, Csite = EF01) ──────────

organic_field <- fields %>% filter(Csite == "EF01")

# ── 3. Shared visual settings ───────────────────────────────

BUFFER  <- 80    # metres of padding around each panel extent

field_fill   <- "#b2df8a"   # soft green for field polygon
collar_col   <- "#7570b3"   # purple  — collars
core_col     <- "#d95f02"   # orange  — soil cores
tower_col    <- "#1b9e77"   # teal    — EC tower

# ── 4. Helper: build a panel ────────────────────────────────

make_panel <- function(panel_title,
                       field_poly,       # sf polygon or NULL
                       cores_data,
                       collars_data,
                       tower_data,
                       buf = BUFFER,
                       show_legend = FALSE) {

  # Bounding box from all points + optional polygon
  all_pts <- bind_rows(cores_data, collars_data, tower_data)
  if (!is.null(field_poly) && nrow(field_poly) > 0) {
    bbox <- st_bbox(c(st_bbox(all_pts), st_bbox(field_poly)))
  } else {
    bbox <- st_bbox(all_pts)
  }

  xlim <- c(bbox["xmin"] - buf, bbox["xmax"] + buf)
  ylim <- c(bbox["ymin"] - buf, bbox["ymax"] + buf)

  p <- ggplot()

  # Field polygon (organic only — conventional has none in shapefile)
  # if (!is.null(field_poly) && nrow(field_poly) > 0) {
  #   p <- p +
  #     geom_sf(data = field_poly, fill = field_fill,
  #             colour = "grey40", linewidth = 0.6, alpha = 0.5)
  # }

  # Soil cores
  p <- p +
    geom_sf(data = cores_data,
            colour = core_col, fill = core_col,
            shape = 21, size = 3) +

    # Collars
    geom_sf(data = collars_data,
            colour = collar_col,
            shape = 1, size = 2, stroke = 1.0) +

    # EC tower
    geom_sf(data = tower_data,
            colour = tower_col,
            shape =13, size = 3, stroke = 1.0) +

    # Panel title as annotation
    annotate("text",
             x = xlim[1] + 0.05 * diff(xlim),
             y = ylim[2] - 0.05 * diff(ylim),
             label = panel_title,
             hjust = 0, vjust = 1,
             size = 4, fontface = "bold") +

    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    annotation_scale(location = "br", width_hint = 0.3,
                     text_cex = 0.7) +

    theme_bw(base_size = 11) +
    theme(
      axis.text        = element_text(size = 7),
      axis.title       = element_blank(),
      panel.grid.major = element_line(colour = "grey85", linewidth = 0.25),
      legend.position  = if (show_legend) "bottom" else "none"
    )

  p
}

# ── 5. Build each panel ─────────────────────────────────────

p_organic <- make_panel(
  panel_title  = "Organic field",
  field_poly   = organic_field,
  cores_data   = cores_sf %>% filter(management == "organic"),
  collars_data = org_collars,
  tower_data   = towers_sf %>% filter(management == "organic")
)

p_conv <- make_panel(
  panel_title  = "Conventional field",
  field_poly   = NULL,   # not in shapefile
  cores_data   = cores_sf %>% filter(management == "conventional"),
  collars_data = conv_collars,
  tower_data   = towers_sf %>% filter(management == "conventional")
)


# ── 6. Shared legend via a dummy plot ───────────────────────

legend_data <- data.frame(
  x    = 1:3,
  y    = 1,
  type = c("Soil core", "Gas exchange collar", "EC tower")
)
shapes <- c("Soil core" = 21, "Gas exchange collar" = 0, "EC tower" = 8)
cols   <- c("Soil core" = core_col,
            "Gas exchange collar" = collar_col,
            "EC tower" = tower_col)

legend_plot <- ggplot(legend_data, aes(x, y, shape = type, colour = type)) +
  geom_point(size = 3, stroke = 1.2) +
  scale_shape_manual(values = shapes, name = NULL) +
  scale_colour_manual(values = cols,  name = NULL) +
  theme_void() +
  theme(legend.position  = "bottom",
        legend.direction = "horizontal",
        legend.text      = element_text(size = 10),
        legend.key.size  = unit(1, "lines"))

shared_legend <- cowplot::get_legend(legend_plot)


# ── 7. Assemble with patchwork ──────────────────────────────

library(patchwork)
library(cowplot)   # for get_legend

combined <- (p_organic | p_conv) /
  wrap_elements(shared_legend) +
  plot_layout(heights = c(10, 1)) +
  plot_annotation(
    # title   = "Sampling locations",
    caption = "CRS: NAD83 / New York West (EPSG:32115)",
    theme   = theme(plot.title = element_text(face = "bold", size = 13))
  )

print(combined)


# ── 8. Save ─────────────────────────────────────────────────

ggsave("map_sampling_locations.png",
       plot   = combined,
       width  = 10,
       height = 6,
       dpi    = 300)

message("Saved → map_sampling_locations.png")
