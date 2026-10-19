#' Per-analyte partitioning (Kow) weights
#'
#' Builds one weight per analyte from octanol-water partition
#' coefficients so that particle-affinitive compounds weigh more in a
#' sediment matrix (or less, for water). Weights are scaled to a mean of 1
#' over the analytes that have a usable value; analytes without one get 1.
#'
#' @param analytes Character vector of analyte names.
#' @param kow Named numeric vector (see [read_kow()]); `NULL` gives all 1.
#' @param input_scale Scale of the values in `kow`: `"log10"` if they are
#'   already log10 Kow (the usual tabulated form), `"linear"` if they are
#'   Kow. Weights are proportional to log10 Kow in both cases.
#' @param transform Optional function applied to the log10 Kow values
#'   before scaling (e.g. `log10` reproduces analyses in which tabulated
#'   log Kow values were log-transformed a second time).
#' @param matrix `"sediment"` (weight increases with Kow) or `"water"`
#'   (reversed, so more soluble compounds weigh more).
#' @return Named numeric vector of weights.
#' @examples
#' kow <- read_kow(system.file("extdata", "kow.csv", package = "smma"))
#' head(kow_weights(names(kow), kow))
#' @export
kow_weights <- function(analytes, kow = NULL, input_scale = c("log10", "linear"),
                        transform = NULL, matrix = c("sediment", "water")) {
  input_scale <- match.arg(input_scale)
  matrix <- match.arg(matrix)
  w <- stats::setNames(rep(1, length(analytes)), analytes)
  if (is.null(kow)) return(w)
  k <- suppressWarnings(as.numeric(kow[analytes]))
  k <- if (input_scale == "linear") ifelse(is.finite(k) & k > 0, log10(k), NA_real_) else k
  if (!is.null(transform)) k <- suppressWarnings(transform(k))
  k[!is.finite(k)] <- NA_real_
  if (matrix == "water" && any(is.finite(k))) {
    rng <- range(k, na.rm = TRUE)
    k <- sum(rng) - k
  }
  ok <- is.finite(k) & k > 0
  if (any(ok)) w[ok] <- k[ok] / mean(k[ok])
  w
}

#' Factor-specific signal-to-noise weights
#'
#' Extracts the S/N row of one factor as a weight vector over `analytes`.
#' Analytes absent from the table get weight 1.
#'
#' @param sn Data frame of S/N values (see [read_sn()]) or `NULL`.
#' @param factor_id Factor identifier.
#' @param analytes Analyte names.
#' @param factor_col Factor column in `sn`.
#' @return Named numeric vector.
#' @export
sn_weights <- function(sn, factor_id, analytes, factor_col = "factor") {
  w <- stats::setNames(rep(1, length(analytes)), analytes)
  if (is.null(sn)) return(w)
  if (!factor_col %in% names(sn)) {
    alt <- names(sn)[tolower(names(sn)) == tolower(factor_col)]
    if (!length(alt)) return(w)
    factor_col <- alt[1]
  }
  row <- sn[as.character(sn[[factor_col]]) == factor_id, , drop = FALSE]
  if (!nrow(row)) return(w)
  common <- intersect(analytes, names(sn))
  w[common] <- suppressWarnings(as.numeric(unlist(row[1, common, drop = TRUE])))
  w
}

# Combine weights; non-finite or non-positive values are replaced by
# `fallback` (1 by default, i.e. neutral weight).
combine_weights <- function(..., fallback = 1) {
  w <- Reduce(`*`, list(...))
  w[!is.finite(w) | w <= 0] <- fallback
  w
}
