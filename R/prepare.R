#' Prepare inputs for Source Mixture Matching Analysis
#'
#' Aligns factor profiles and a source-profile library on their shared
#' analytes and normalizes every profile to unit sum over those analytes.
#'
#' @param factor_profiles Data frame with one row per factor (see
#'   [read_factor_profiles()]), or a sample table when SMMA is applied
#'   directly to samples (see [smma_samples()]).
#' @param source_library Data frame with one row per source profile (see
#'   [read_source_library()]).
#' @param factor_col Column in `factor_profiles` identifying factors.
#' @param source_col,category_col Columns in `source_library` holding the
#'   source identifier and its category.
#' @param analytes Optional character vector of analyte columns to use.
#'   By default, all columns matching `analyte_pattern` that occur in both
#'   tables.
#' @param analyte_pattern Regular expression used when `analytes` is `NULL`.
#' @param min_shared Minimum number of shared analytes required.
#' @return An object of class `smma_input`: a list with
#'   \describe{
#'     \item{F}{Factor matrix (factors x analytes), rows sum to 1.}
#'     \item{L}{Library matrix (analytes x sources), columns normalized,
#'       missing values set to 0.}
#'     \item{L_missing}{Logical matrix flagging values that were missing in
#'       the library before zero-filling.}
#'     \item{categories}{Data frame `source`, `category`.}
#'     \item{analytes}{The analytes used.}
#'   }
#' @examples
#' ex <- smma_example()
#' inp <- smma_prepare(ex$profiles, ex$library)
#' inp
#' @export
smma_prepare <- function(factor_profiles, source_library,
                         factor_col = "factor", source_col = "source",
                         category_col = "category", analytes = NULL,
                         analyte_pattern = "^(pah_|tph_|biom_)", min_shared = 12) {
  check_cols(factor_profiles, factor_col, "factor_profiles")
  check_cols(source_library, c(source_col, category_col), "source_library")
  if (is.null(analytes)) {
    analytes <- intersect(grep(analyte_pattern, names(factor_profiles), value = TRUE),
                          grep(analyte_pattern, names(source_library), value = TRUE))
  } else {
    check_cols(factor_profiles, analytes, "factor_profiles")
    check_cols(source_library, analytes, "source_library")
  }
  if (length(analytes) < min_shared) {
    stop(sprintf("Only %d shared analytes between factors and library (min_shared = %d).",
                 length(analytes), min_shared), call. = FALSE)
  }
  fp <- as_numeric_cols(factor_profiles, analytes)
  sl <- as_numeric_cols(source_library, analytes)

  Fm <- normalize_rows(as.matrix(fp[, analytes, drop = FALSE]))
  rownames(Fm) <- as.character(fp[[factor_col]])
  colnames(Fm) <- analytes

  Ls <- normalize_rows(as.matrix(sl[, analytes, drop = FALSE]))
  L <- t(Ls)
  rownames(L) <- analytes
  colnames(L) <- as.character(sl[[source_col]])
  L_missing <- !is.finite(L)
  L[L_missing] <- 0

  structure(list(
    F = Fm, L = L, L_missing = L_missing,
    categories = data.frame(source = as.character(sl[[source_col]]),
                            category = as.character(sl[[category_col]]),
                            stringsAsFactors = FALSE),
    analytes = analytes, factor_col = factor_col
  ), class = "smma_input")
}

#' Print an smma_input object
#'
#' @param x An `smma_input` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.smma_input <- function(x, ...) {
  cat("<smma_input>\n")
  cat(sprintf("  %d factors x %d analytes; library of %d sources in %d categories\n",
              nrow(x$F), length(x$analytes), ncol(x$L), length(unique(x$categories$category))))
  cat(sprintf("  library completeness: %.0f%% of source x analyte values present\n",
              100 * mean(!x$L_missing)))
  invisible(x)
}
