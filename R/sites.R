#' Default population tiers for site classification
#'
#' Each settlement is assigned a tier by population, and each tier has an
#' "urban" radius (site within the settlement's built-up/traffic
#' footprint) and a wider "mixed" radius (partial influence). Tiers
#' without an urban radius can only produce "Mixed" classifications.
#' These are the values used for coastal British Columbia in the
#' accompanying study; adapt them to your region.
#'
#' @return Data frame with `tier`, `min_pop`, `urban_km`, `mixed_km`.
#' @export
default_population_tiers <- function() {
  data.frame(
    tier = c("large", "medium", "small_town", "village", "hamlet"),
    min_pop = c(90000, 10000, 2000, 300, 0),
    urban_km = c(15, 8, 4, NA, NA),
    mixed_km = c(30, 20, 12, 6, 3),
    stringsAsFactors = FALSE
  )
}

#' Great-circle distance (km)
#'
#' @param lat1,lon1,lat2,lon2 Coordinates in decimal degrees (vectorized).
#' @return Distance in kilometres (haversine, Earth radius 6371 km).
#' @export
haversine_km <- function(lat1, lon1, lat2, lon2) {
  rad <- pi / 180
  dlat <- (lat2 - lat1) * rad; dlon <- (lon2 - lon1) * rad
  a <- sin(dlat / 2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon / 2)^2
  2 * 6371 * asin(pmin(1, sqrt(a)))
}

#' Classify sites as Urban, Mixed or Remote from population centres
#'
#' A site is "Urban" if it lies within the urban radius of any settlement,
#' "Mixed" if within any settlement's mixed radius but no urban radius,
#' and "Remote" otherwise.
#'
#' @param sites Data frame with a site identifier and coordinates.
#' @param centres Data frame of settlements with name, coordinates and population.
#' @param site_col,site_lat,site_lon Column names in `sites`.
#' @param centre_name,centre_lat,centre_lon,centre_pop Column names in `centres`.
#' @param tiers Tier table, see [default_population_tiers()].
#' @return `sites` with added columns `site_class`, `nearest_centre` and
#'   `nearest_km` (distance to the nearest settlement).
#' @examples
#' sites <- data.frame(id = c("a", "b"), lat = c(49.28, 52.0), lon = c(-123.12, -128.0))
#' towns <- data.frame(name = "Vancouver", latitude = 49.2827, longitude = -123.1207,
#'                     population = 662248)
#' classify_sites(sites, towns, site_col = "id", site_lat = "lat", site_lon = "lon")
#' @export
classify_sites <- function(sites, centres, site_col = "site", site_lat = "latitude",
                           site_lon = "longitude", centre_name = "name",
                           centre_lat = "latitude", centre_lon = "longitude",
                           centre_pop = "population", tiers = default_population_tiers()) {
  check_cols(sites, c(site_col, site_lat, site_lon), "sites")
  check_cols(centres, c(centre_name, centre_lat, centre_lon, centre_pop), "centres")
  tiers <- tiers[order(-tiers$min_pop), ]
  pop <- as.numeric(centres[[centre_pop]])
  ti <- vapply(pop, function(p) which(p >= tiers$min_pop)[1], integer(1))
  urb <- tiers$urban_km[ti]; mix <- tiers$mixed_km[ti]
  cls <- character(nrow(sites)); near <- character(nrow(sites)); nkm <- numeric(nrow(sites))
  for (i in seq_len(nrow(sites))) {
    d <- haversine_km(sites[[site_lat]][i], sites[[site_lon]][i],
                      centres[[centre_lat]], centres[[centre_lon]])
    cls[i] <- if (any(!is.na(urb) & d <= urb, na.rm = TRUE)) "Urban" else
      if (any(d <= mix, na.rm = TRUE)) "Mixed" else "Remote"
    j <- which.min(d)
    near[i] <- as.character(centres[[centre_name]][j]); nkm[i] <- round(d[j], 1)
  }
  sites$site_class <- factor(cls, levels = c("Urban", "Mixed", "Remote"))
  sites$nearest_centre <- near
  sites$nearest_km <- nkm
  sites
}
