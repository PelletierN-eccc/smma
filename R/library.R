#' Fill missing analytes in a source library (sensitivity analysis)
#'
#' Literature source profiles often report only part of the analyte
#' suite (for example, no petroleum biomarkers). By default SMMA treats a
#' missing value as zero, which can structurally favour profiles that do
#' report those analytes. This function fills missing values with
#' user-supplied bounding values so that the effect of that assumption
#' can be tested by refitting (see the vignette).
#'
#' @param source_library Source library data frame.
#' @param fill Data frame with columns `category`, `analyte` and `value`,
#'   and optionally `ref_analyte`. If `ref_analyte` is given (non-NA), the
#'   filled value is `value * library[[ref_analyte]]` (ratio basis);
#'   otherwise `value` is used as is (absolute basis).
#' @param category_col Category column of the library.
#' @param only_missing Only replace missing (NA) values; existing values
#'   are never overwritten when `TRUE`.
#' @return The library with filled values, plus a logical attribute
#'   `"filled"` marking the modified cells.
#' @export
fill_library_gaps <- function(source_library, fill, category_col = "category",
                              only_missing = TRUE) {
  check_cols(fill, c("category", "analyte", "value"), "fill")
  lib <- source_library
  marks <- list()
  for (k in seq_len(nrow(fill))) {
    a <- fill$analyte[k]
    if (!a %in% names(lib)) lib[[a]] <- NA_real_
    lib[[a]] <- suppressWarnings(as.numeric(lib[[a]]))
    rows <- which(lib[[category_col]] == fill$category[k])
    if (only_missing) rows <- rows[is.na(lib[[a]][rows])]
    if (!length(rows)) next
    ref <- if ("ref_analyte" %in% names(fill)) fill$ref_analyte[k] else NA
    val <- if (!is.na(ref) && nzchar(ref)) {
      fill$value[k] * suppressWarnings(as.numeric(lib[[ref]][rows]))
    } else rep(fill$value[k], length(rows))
    lib[[a]][rows] <- val
    marks[[length(marks) + 1]] <- data.frame(row = rows, analyte = a)
  }
  attr(lib, "filled") <- if (length(marks)) do.call(rbind, marks) else NULL
  lib
}

#' Replace library values for selected sources
#'
#' Overwrites analyte values of the sources listed in `values` (one row
#' per source, one column per analyte to replace). Typical use: build
#' library variants with bounding values for analytes that literature
#' profiles do not report, then refit to test sensitivity.
#'
#' @param source_library Source library.
#' @param values Data frame with the source identifier column and the
#'   analyte columns to overwrite (other columns are ignored).
#' @param source_col Source identifier column (in both tables).
#' @param analytes Analyte columns to take from `values` (default: those
#'   matching `analyte_pattern`).
#' @param analyte_pattern Pattern used when `analytes` is `NULL`.
#' @return The modified library.
#' @export
update_library <- function(source_library, values, source_col = "source", analytes = NULL,
                           analyte_pattern = "^(pah_|tph_|biom_)") {
  check_cols(values, source_col, "values")
  if (is.null(analytes)) analytes <- grep(analyte_pattern, names(values), value = TRUE)
  lib <- source_library
  idx <- match(trimws(values[[source_col]]), trimws(lib[[source_col]]))
  if (anyNA(idx)) {
    warning(sum(is.na(idx)), " source(s) in 'values' not found in the library.", call. = FALSE)
  }
  for (a in analytes) {
    if (!a %in% names(lib)) lib[[a]] <- NA_real_
    lib[[a]] <- suppressWarnings(as.numeric(lib[[a]]))
    ok <- !is.na(idx)
    lib[[a]][idx[ok]] <- suppressWarnings(as.numeric(values[[a]][ok]))
  }
  lib
}

#' Apply SMMA directly to samples
#'
#' Fits each sample's own normalized composition against the library,
#' bypassing the receptor model. Useful for testing whether a specific
#' source is present in particular samples, ideally with a short,
#' well-justified candidate list (`sources`). Without the cross-sample
#' pooling of a receptor model, low-concentration samples often give
#' poorly constrained mixtures.
#'
#' @param samples Data frame of sample concentrations, one row per sample.
#' @param source_library Source library.
#' @param sample_col Sample identifier column.
#' @param sources Optional character vector restricting the candidate sources.
#' @param B Bootstrap draws per sample (0 for point estimates only).
#' @param ... Passed to [smma_prepare()] and [smma_bootstrap()] (e.g.
#'   `analytes`, `kow`, `seed`).
#' @return An `smma_bootstrap` object (B > 0) or the output of
#'   [smma_fit()] (B = 0), with samples in place of factors.
#' @export
smma_samples <- function(samples, source_library, sample_col = "sample", sources = NULL,
                         B = 300, ...) {
  dots <- list(...)
  lib <- if (is.null(sources)) source_library else
    source_library[source_library[[dots$source_col %||% "source"]] %in% sources, , drop = FALSE]
  prep_args <- dots[intersect(names(dots), names(formals(smma_prepare)))]
  inp <- do.call(smma_prepare, c(list(factor_profiles = samples, source_library = lib,
                                      factor_col = sample_col), prep_args))
  fit_args <- dots[intersect(names(dots), names(formals(smma_bootstrap)))]
  if (B > 0) do.call(smma_bootstrap, c(list(input = inp, B = B), fit_args)) else
    do.call(smma_fit, c(list(input = inp), dots[intersect(names(dots), names(formals(smma_fit)))]))
}
