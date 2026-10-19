#' Example data: British Columbia intertidal sediments
#'
#' Loads the example files shipped in `inst/extdata`: the four-factor
#' EPA PMF 5.0 solution for 211 intertidal sediment samples from coastal
#' British Columbia (profiles, contributions and S/N), the 371-profile
#' PAC/biomarker source library, and log10 Kow values. Field data are
#' published in the ECCC Data Catalogue ("ESTS Shoreline surveys data").
#'
#' @return A list with `profiles`, `contributions`, `sn`, `library`, `kow`.
#' @examples
#' ex <- smma_example()
#' names(ex)
#' @export
smma_example <- function() {
  f <- function(x) system.file("extdata", x, package = "smma", mustWork = TRUE)
  list(
    profiles = read_factor_profiles(f("bc_factor_profiles.csv")),
    contributions = read_factor_contributions(f("bc_factor_contributions.csv")),
    sn = read_sn(f("bc_factor_sn.csv")),
    library = read_source_library(f("source_library.csv")),
    kow = read_kow(f("kow.csv"))
  )
}

#' Simulate factor profiles from known library mixtures
#'
#' Builds synthetic "factors" as known mixtures of library profiles, with
#' optional multiplicative log-normal noise, to check how well SMMA
#' recovers a known answer with a given library and analyte suite.
#'
#' @param source_library Source library.
#' @param mixtures Named list; each element is a named numeric vector of
#'   source weights (names are source identifiers).
#' @param noise_sd Standard deviation of log-normal noise (0 = none).
#' @param analytes Analytes to keep (default: all matching `analyte_pattern`).
#' @param analyte_pattern Pattern for analyte columns.
#' @param source_col Source identifier column.
#' @param seed Optional seed.
#' @return Data frame of factor profiles (column `factor` + analytes) with
#'   the true mixtures in attribute `"truth"`.
#' @examples
#' ex <- smma_example()
#' src <- ex$library$source[1:2]
#' sim <- smma_simulate(ex$library, list(F1 = stats::setNames(c(0.7, 0.3), src)))
#' attr(sim, "truth")
#' @export
smma_simulate <- function(source_library, mixtures, noise_sd = 0, analytes = NULL,
                          analyte_pattern = "^(pah_|tph_|biom_)", source_col = "source",
                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(analytes)) analytes <- grep(analyte_pattern, names(source_library), value = TRUE)
  lib <- as_numeric_cols(source_library, analytes)
  M <- as.matrix(lib[, analytes, drop = FALSE])
  M[!is.finite(M)] <- 0
  M <- M / pmax(rowSums(M), .Machine$double.eps)
  rownames(M) <- lib[[source_col]]
  rows <- lapply(names(mixtures), function(nm) {
    w <- mixtures[[nm]]
    miss <- setdiff(names(w), rownames(M))
    if (length(miss)) stop("Unknown sources in mixture ", nm, ": ", paste(miss, collapse = ", "))
    v <- colSums(M[names(w), , drop = FALSE] * (w / sum(w)))
    if (noise_sd > 0) v <- v * exp(stats::rnorm(length(v), 0, noise_sd))
    v
  })
  out <- data.frame(factor = names(mixtures), do.call(rbind, rows), check.names = FALSE,
                    stringsAsFactors = FALSE)
  attr(out, "truth") <- mixtures
  out
}
