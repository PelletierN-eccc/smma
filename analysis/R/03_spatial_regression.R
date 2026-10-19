# Distance-decay regressions of factor contributions against marine
# infrastructure (Section 3.2, Supplementary Section S15).
#
# Download the two Fisheries and Oceans Canada layers (Open Government
# Licence - Canada) and unzip them under analysis/data/ancillary/:
#   - Canadian Anchorages and Anchorage Areas
#     https://open.canada.ca/data/en/dataset/622a7f72-4a00-4f9e-b04f-af6551c77db3
#   - Small Craft Harbours Locations and Information
#     https://open.canada.ca/data/en/dataset/262451e7-6416-47b5-8453-f31d212ea657
# and set the file names below.
source("analysis/R/00_config.R")
library(sf)

layer_files <- list(
  anchorage     = file.path(cfg$data_dir, "ancillary", "ACHARE_A.shp"),
  small_harbour = file.path(cfg$data_dir, "ancillary", "small_craft_harbour.shp")
)
missing <- !file.exists(unlist(layer_files))
if (any(missing)) stop("Missing layer(s): ", paste(unlist(layer_files)[missing], collapse = ", "),
                       "\nSee the header of this script for download links.")

contrib <- read_factor_contributions(cfg$ex("bc_factor_contributions.csv"))
sites   <- read.csv(file.path(cfg$data_dir, "sites.csv"))
d <- merge(contrib, sites, by = "sample")

pts <- st_as_sf(d, coords = c("longitude", "latitude"), crs = 4326)
pts <- st_transform(pts, 3005)                         # BC Albers (metres)
layers <- lapply(layer_files, st_read, quiet = TRUE)

d <- cbind(d, distance_to_features(pts, layers))
xy <- st_coordinates(pts); d$X <- xy[, 1]; d$Y <- xy[, 2]

factors <- grep("^Factor", names(contrib), value = TRUE)
d <- relative_contributions(d, factors)
res <- fit_distance_decay(d,
                          responses = c(factors, paste0(factors, "_logit")),
                          distances = grep("^dist_", names(d), value = TRUE))
res$response_type <- ifelse(grepl("_logit$", res$response), "relative_contribution", "magnitude")
res$significance <- cut(res$p_value, c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", ""))
out(res, "distance_decay_gls.csv")
print(res[, c("response", "distance", "estimate", "p_value", "significance")])
