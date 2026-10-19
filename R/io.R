#' Read a table from CSV or Excel
#'
#' Small convenience reader used by the other `read_*()` functions. CSV
#' files are read with `utils::read.csv()` (a UTF-8 byte-order mark is
#' handled); `.xlsx`/`.xls` files require the \pkg{readxl} package.
#'
#' @param path File path.
#' @param sheet Sheet name or number for Excel files (ignored for CSV).
#' @return A data frame.
#' @export
read_table_any <- function(path, sheet = NULL) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    need_pkg("readxl", "to read Excel files")
    as.data.frame(readxl::read_excel(path, sheet = sheet %||% 1))
  } else {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                    fileEncoding = "UTF-8-BOM")
  }
}

#' Read receptor-model factor profiles
#'
#' Reads factor profiles (concentration of each analyte in each factor)
#' and returns one row per factor. Accepts the "long" layout exported by
#' EPA PMF 5.0 (one row per species, one column per factor) or a "wide"
#' layout that already has one row per factor.
#'
#' @param x A file path (CSV or Excel) or a data frame.
#' @param sheet Excel sheet (default `"Profiles_conc"` when reading an
#'   Excel file; ignored for CSV and data frames).
#' @param species_col Name of the species column in the long layout. If
#'   this column is absent, `x` is assumed to be wide already.
#' @param factor_col Name of the factor-identifier column in the output
#'   (and in the input, for the wide layout).
#' @return A data frame with a `factor_col` column and one numeric column
#'   per analyte.
#' @examples
#' prof <- read_factor_profiles(system.file("extdata", "bc_factor_profiles.csv",
#'                                          package = "smma"))
#' prof[, 1:5]
#' @export
read_factor_profiles <- function(x, sheet = "Profiles_conc", species_col = "species",
                                 factor_col = "factor") {
  df <- if (is.character(x)) {
    read_table_any(x, sheet = if (grepl("\\.xlsx?$", x, ignore.case = TRUE)) sheet else NULL)
  } else as.data.frame(x)
  if (species_col %in% names(df)) {
    df <- df[!is.na(df[[species_col]]) & nzchar(trimws(df[[species_col]])), , drop = FALSE]
    species <- as.character(df[[species_col]])
    vals <- df[, setdiff(names(df), species_col), drop = FALSE]
    m <- t(as.matrix(as_numeric_cols(vals, names(vals))))
    out <- data.frame(names(vals), m, check.names = FALSE, stringsAsFactors = FALSE,
                      row.names = NULL)
    names(out) <- c(factor_col, species)
    out
  } else {
    check_cols(df, factor_col, "factor profiles")
    df[[factor_col]] <- as.character(df[[factor_col]])
    df
  }
}

#' Read receptor-model factor contributions (sample scores)
#'
#' @param x A file path (CSV or Excel) or a data frame with one row per
#'   sample and one column per factor.
#' @param sheet Excel sheet (default `"Contributions"`).
#' @param sample_col Name of the sample-identifier column.
#' @return A data frame with `sample_col` and numeric factor columns.
#' @export
read_factor_contributions <- function(x, sheet = "Contributions", sample_col = "sample") {
  df <- if (is.character(x)) {
    read_table_any(x, sheet = if (grepl("\\.xlsx?$", x, ignore.case = TRUE)) sheet else NULL)
  } else as.data.frame(x)
  check_cols(df, sample_col, "factor contributions")
  df[[sample_col]] <- as.character(df[[sample_col]])
  as_numeric_cols(df, setdiff(names(df), sample_col))
}

#' Read factor-specific signal-to-noise (S/N) values
#'
#' @param x File path or data frame with one row per factor and one
#'   column per analyte (as in the `SN` sheet written for EPA PMF output).
#' @param sheet Excel sheet (default `"SN"`).
#' @param factor_col Name of the factor-identifier column.
#' @return A data frame.
#' @export
read_sn <- function(x, sheet = "SN", factor_col = "factor") {
  df <- if (is.character(x)) {
    read_table_any(x, sheet = if (grepl("\\.xlsx?$", x, ignore.case = TRUE)) sheet else NULL)
  } else as.data.frame(x)
  if (!factor_col %in% names(df)) {
    alt <- names(df)[tolower(names(df)) == tolower(factor_col)]
    if (length(alt)) names(df)[names(df) == alt[1]] <- factor_col
  }
  check_cols(df, factor_col, "S/N table")
  df[[factor_col]] <- as.character(df[[factor_col]])
  as_numeric_cols(df, setdiff(names(df), factor_col))
}

#' Read octanol-water partition coefficients
#'
#' Accepts either a one-row "wide" table (one column per analyte) or a
#' two-column "long" table (`analyte`, `value`).
#'
#' @param x File path or data frame.
#' @param analyte_col,value_col Column names for the long layout.
#' @return A named numeric vector (names are analyte codes).
#' @export
read_kow <- function(x, analyte_col = "analyte", value_col = "kow") {
  df <- if (is.character(x)) read_table_any(x) else as.data.frame(x)
  if (all(c(analyte_col, value_col) %in% names(df))) {
    stats::setNames(suppressWarnings(as.numeric(df[[value_col]])), df[[analyte_col]])
  } else {
    v <- suppressWarnings(as.numeric(unlist(df[1, , drop = TRUE])))
    stats::setNames(v, names(df))
  }
}

#' Read a source-profile library
#'
#' @param x File path or data frame with one row per source profile, an
#'   identifier column, a category column and one column per analyte.
#' @param source_col,category_col Identifier and category column names.
#' @return A data frame (source identifiers are made unique if needed).
#' @export
read_source_library <- function(x, source_col = "source", category_col = "category") {
  df <- if (is.character(x)) read_table_any(x) else as.data.frame(x)
  check_cols(df, c(source_col, category_col), "source library")
  df[[source_col]] <- trimws(as.character(df[[source_col]]))
  if (anyDuplicated(df[[source_col]])) {
    warning("Duplicated source identifiers were made unique with make.unique().", call. = FALSE)
    df[[source_col]] <- make.unique(df[[source_col]])
  }
  df
}
