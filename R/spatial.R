#' Distance from points to the nearest feature of each layer
#'
#' @param points An `sf` object of sample locations.
#' @param layers Named list of `sf` objects (e.g. anchorages, harbours).
#'   Polygons are converted to their boundaries so that points inside a
#'   polygon get the distance to its edge rather than zero; set
#'   `polygon_as_boundary = FALSE` to keep zero distances inside polygons.
#' @param crs Projected CRS (EPSG code or `crs` object) used for
#'   distances; defaults to the CRS of `points`. Use a metric projection.
#' @param polygon_as_boundary See `layers`.
#' @return Data frame with one `dist_<layer>` column (metres) per layer,
#'   in the row order of `points`.
#' @export
distance_to_features <- function(points, layers, crs = NULL, polygon_as_boundary = TRUE) {
  need_pkg("sf", "for distance calculations")
  if (!is.null(crs)) points <- sf::st_transform(points, crs)
  target <- sf::st_crs(points)
  if (isTRUE(sf::st_is_longlat(points))) {
    warning("points are in geographic coordinates; distances use ellipsoidal ",
            "calculations. Supply a projected 'crs' for faster, planar distances.", call. = FALSE)
  }
  out <- lapply(names(layers), function(nm) {
    ly <- sf::st_transform(layers[[nm]], target)
    if (polygon_as_boundary && any(sf::st_geometry_type(ly) %in% c("POLYGON", "MULTIPOLYGON"))) {
      ly <- sf::st_boundary(ly)
    }
    d <- sf::st_distance(points, ly)
    apply(matrix(as.numeric(d), nrow = nrow(points)), 1, min)
  })
  out <- as.data.frame(out)
  names(out) <- paste0("dist_", names(layers))
  out
}

#' Relative factor contributions and their logits
#'
#' @param contributions Data frame of factor contributions.
#' @param factors Factor columns.
#' @param eps Bounds used before the logit transform.
#' @return `contributions` with `<factor>_pct` and `<factor>_logit`
#'   columns added (negative contributions set to 0 first).
#' @export
relative_contributions <- function(contributions, factors, eps = 0.001) {
  m <- pmax(as.matrix(as_numeric_cols(contributions, factors)[, factors, drop = FALSE]), 0)
  tot <- rowSums(m, na.rm = TRUE)
  for (f in factors) {
    p <- ifelse(tot > 0, m[, f] / tot, NA_real_)
    contributions[[paste0(f, "_pct")]] <- 100 * p
    contributions[[paste0(f, "_logit")]] <- stats::qlogis(pmin(pmax(p, eps), 1 - eps))
  }
  contributions
}

#' Spatial distance-decay regressions
#'
#' Fits, for every response x distance pair, a generalized least-squares
#' model `response ~ log1p(distance)` with a spatial correlation structure
#' on the sample coordinates (\pkg{nlme}), and returns the slope of the
#' distance term. A negative slope means the response increases toward
#' the features.
#'
#' Duplicate coordinates are averaged first (`average_duplicates = TRUE`),
#' because coincident points make spatial correlation structures singular.
#'
#' @param data Data frame with responses, distances and coordinates.
#' @param responses Response columns (e.g. factor contributions and
#'   their `_logit` relative contributions).
#' @param distances Distance columns (e.g. from [distance_to_features()]).
#' @param x,y Coordinate columns (projected units).
#' @param correlation One of `"exp"`, `"gaus"`, `"spher"`.
#' @param log_distance Use `log1p(distance)` as the predictor.
#' @param average_duplicates See Details.
#' @return Data frame with `response`, `distance`, `estimate`,
#'   `std_error`, `t_value`, `p_value`, `n`; failed fits are returned
#'   with `NA` and a `note`.
#' @export
fit_distance_decay <- function(data, responses, distances, x = "X", y = "Y",
                               correlation = c("exp", "gaus", "spher"),
                               log_distance = TRUE, average_duplicates = TRUE) {
  need_pkg("nlme", "for spatial GLS models")
  correlation <- match.arg(correlation)
  check_cols(data, c(responses, distances, x, y), "data")
  d <- data[, c(responses, distances, x, y)]
  if (average_duplicates) {
    d <- stats::aggregate(d[, c(responses, distances)], by = list(X = d[[x]], Y = d[[y]]),
                          FUN = function(v) mean(v, na.rm = TRUE))
    x <- "X"; y <- "Y"
  }
  grid <- expand.grid(response = responses, distance = distances, stringsAsFactors = FALSE)
  rows <- lapply(seq_len(nrow(grid)), function(k) {
    rs <- grid$response[k]; ds <- grid$distance[k]
    pred <- if (log_distance) log1p(d[[ds]]) else d[[ds]]
    dd <- data.frame(resp = d[[rs]], pred = pred, X = d[[x]], Y = d[[y]])
    dd <- dd[stats::complete.cases(dd) & is.finite(dd$resp) & is.finite(dd$pred), ]
    cs <- switch(correlation,
                 exp = nlme::corExp(form = ~ X + Y, nugget = TRUE),
                 gaus = nlme::corGaus(form = ~ X + Y, nugget = TRUE),
                 spher = nlme::corSpher(form = ~ X + Y, nugget = TRUE))
    fit <- if (nrow(dd) >= 10) try(nlme::gls(resp ~ pred, data = dd, correlation = cs,
                                             method = "REML",
                                             control = nlme::glsControl(msMaxIter = 200, opt = "optim")),
                                   silent = TRUE) else structure("too few observations", class = "try-error")
    if (inherits(fit, "try-error")) {
      return(data.frame(response = rs, distance = ds, estimate = NA, std_error = NA,
                        t_value = NA, p_value = NA, n = nrow(dd),
                        note = trimws(as.character(fit)), stringsAsFactors = FALSE))
    }
    tt <- summary(fit)$tTable
    data.frame(response = rs, distance = ds, estimate = tt["pred", 1], std_error = tt["pred", 2],
               t_value = tt["pred", 3], p_value = tt["pred", 4], n = nrow(dd), note = "",
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
