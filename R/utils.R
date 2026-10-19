# Internal helpers --------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# Coerce selected columns of a data frame to numeric (silently turning
# non-numeric entries such as "<MDL" or "" into NA).
as_numeric_cols <- function(df, cols) {
  for (cl in cols) df[[cl]] <- suppressWarnings(as.numeric(df[[cl]]))
  df
}

# Row-normalize a numeric matrix so each row sums to 1 (NA ignored).
# Rows that are all NA or sum to zero become all NA.
normalize_rows <- function(m) {
  s <- rowSums(m, na.rm = TRUE)
  bad <- !is.finite(s) | s <= 0 | rowSums(!is.na(m)) == 0
  out <- m / s
  out[bad, ] <- NA_real_
  out
}

# Weighted empirical quantile: smallest x at which the cumulative
# (normalized) weight of the sorted values reaches each probability.
weighted_quantile <- function(x, w, probs = c(0.5, 0.025, 0.975)) {
  ok <- is.finite(x) & is.finite(w)
  x <- x[ok]; w <- w[ok]
  if (!length(x) || sum(w) <= 0) return(rep(NA_real_, length(probs)))
  o <- order(x)
  x <- x[o]; w <- w[o] / sum(w)
  cw <- cumsum(w)
  vapply(probs, function(p) x[which(cw >= p)[1]], numeric(1))
}

weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  x <- x[ok]; w <- w[ok]
  if (!length(x) || sum(w) == 0) return(NA_real_)
  sum(x * w) / sum(w)
}

# Relative RMSE (optionally weighted) between residuals and observations.
rel_rmse <- function(resid, obs, w = NULL) {
  if (is.null(w)) return(sqrt(mean(resid^2)) / sqrt(mean(obs^2)))
  ok <- is.finite(resid) & is.finite(obs) & is.finite(w) & w > 0
  resid <- resid[ok]; obs <- obs[ok]; w <- w[ok]
  if (!length(resid) || sum(w) <= 0) return(NA_real_)
  w <- w / sum(w)
  sqrt(sum(w * resid^2)) / sqrt(sum(w * obs^2))
}

check_cols <- function(df, cols, what = "data") {
  miss <- setdiff(cols, names(df))
  if (length(miss)) {
    stop(sprintf("Column(s) not found in %s: %s", what, paste(miss, collapse = ", ")),
         call. = FALSE)
  }
  invisible(TRUE)
}

need_pkg <- function(pkg, why) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required %s. Install it with install.packages('%s').",
                 pkg, why, pkg), call. = FALSE)
  }
}
