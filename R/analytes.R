#' Analyte naming scheme
#'
#' All functions that need to know *which* column is phenanthrene, which
#' columns are alkylated homologues, or which belong to the EPA 16
#' priority list take an analyte scheme. `ests_analyte_scheme()` returns
#' the column-naming convention used by the Emergency Science and
#' Technology Section (ESTS) of Environment and Climate Change Canada
#' (e.g. `pah_c0_p` = phenanthrene, `pah_c2_p` = C2-phenanthrenes,
#' `biom_c30_a_b_h` = C30 alpha-beta hopane). To use your own column
#' names, copy this list and edit the entries, or build one with
#' [analyte_scheme()].
#'
#' @return A list with elements:
#' \describe{
#'   \item{pac_pattern}{Regular expression identifying PAC columns.}
#'   \item{biomarker_pattern, tph_pattern}{Patterns for biomarker and TPH columns.}
#'   \item{alkylated_pattern}{Pattern identifying alkylated (C1-C4) homologues among PAC columns.}
#'   \item{exclude_pattern}{Pattern for summary columns (sums, indices) never treated as analytes.}
#'   \item{epa16}{Character vector of the 16 EPA priority PAH columns.}
#'   \item{ratio}{Named character vector mapping compounds used in
#'     diagnostic ratios (`ant`, `phe`, `fla`, `pyr`, `baa`, `chr`, `icdp`, `bghip`)
#'     to column names.}
#'   \item{homologues}{Named list of homologue series (parent first, then C1, C2, ...).}
#'   \item{labels}{Optional named character vector of display labels.}
#' }
#' @examples
#' sch <- ests_analyte_scheme()
#' sch$ratio
#' @export
ests_analyte_scheme <- function() {
  analyte_scheme(
    pac_pattern = "^pah_",
    biomarker_pattern = "^biom_",
    tph_pattern = "^tph_",
    alkylated_pattern = "_c[1-4]_",
    exclude_pattern = "_sum_|^pah_sum|^pah_total|pyrogenic_index",
    epa16 = c("pah_c0_n", "pah_acl", "pah_ace", "pah_c0_f", "pah_c0_p", "pah_an",
              "pah_fl", "pah_py", "pah_ba_a", "pah_c0_c", "pah_bb_f", "pah_bk_f",
              "pah_ba_p", "pah_ip", "pah_da", "pah_bg_p"),
    ratio = c(ant = "pah_an", phe = "pah_c0_p", fla = "pah_fl", pyr = "pah_py",
              baa = "pah_ba_a", chr = "pah_c0_c", icdp = "pah_ip", bghip = "pah_bg_p"),
    homologues = list(
      "Naphthalene"           = c("pah_c0_n", "pah_c1_n", "pah_c2_n", "pah_c3_n", "pah_c4_n"),
      "Phenanthrene"          = c("pah_c0_p", "pah_c1_p", "pah_c2_p", "pah_c3_p", "pah_c4_p"),
      "Dibenzothiophene"      = c("pah_c0_d", "pah_c1_d", "pah_c2_d", "pah_c3_d"),
      "Fluorene"              = c("pah_c0_f", "pah_c1_f", "pah_c2_f", "pah_c3_f"),
      "Fluoranthene/Pyrene"   = c("pah_c0_fl", "pah_c1_fl", "pah_c2_fl", "pah_c3_fl", "pah_c4_fl"),
      "Benzonaphthothiophene" = c("pah_c0_b", "pah_c1_b", "pah_c2_b", "pah_c3_b", "pah_c4_b"),
      "Chrysene"              = c("pah_c0_c", "pah_c1_c", "pah_c2_c", "pah_c3_c")
    )
  )
}

#' Build a custom analyte scheme
#'
#' @param pac_pattern,biomarker_pattern,tph_pattern Regular expressions
#'   identifying PAC, biomarker and TPH columns.
#' @param alkylated_pattern Regular expression identifying alkylated
#'   homologues among the PAC columns.
#' @param exclude_pattern Regular expression for summary columns to drop
#'   (use `NULL` for none).
#' @param epa16 Character vector of EPA 16 priority PAH column names.
#' @param ratio Named character vector; see [ests_analyte_scheme()].
#' @param homologues Named list of homologue series (parent first).
#' @param labels Optional named character vector of display labels.
#' @return A list of class `smma_scheme`.
#' @export
analyte_scheme <- function(pac_pattern, biomarker_pattern = NULL, tph_pattern = NULL,
                           alkylated_pattern = NULL, exclude_pattern = NULL,
                           epa16 = character(), ratio = character(),
                           homologues = list(), labels = NULL) {
  structure(list(pac_pattern = pac_pattern, biomarker_pattern = biomarker_pattern,
                 tph_pattern = tph_pattern, alkylated_pattern = alkylated_pattern,
                 exclude_pattern = exclude_pattern, epa16 = epa16, ratio = ratio,
                 homologues = homologues, labels = labels),
            class = "smma_scheme")
}

#' Define analyte groupings for source attribution
#'
#' Returns the four groupings used in the accompanying manuscript
#' (EPA16, Parent, Alkylated, Total) for a given set of analyte columns.
#' Any named list of character vectors can be used instead in
#' [attribute_sources()].
#'
#' @param analytes Character vector of available analyte column names
#'   (e.g. `names(profiles)`).
#' @param scheme Analyte scheme, see [ests_analyte_scheme()].
#' @return A named list of character vectors.
#' @examples
#' prof <- read_factor_profiles(system.file("extdata", "bc_factor_profiles.csv",
#'                                          package = "smma"))
#' str(analyte_groupings(names(prof)))
#' @export
analyte_groupings <- function(analytes, scheme = ests_analyte_scheme()) {
  pac <- grep(scheme$pac_pattern, analytes, value = TRUE)
  if (!is.null(scheme$exclude_pattern)) pac <- pac[!grepl(scheme$exclude_pattern, pac)]
  alk <- if (is.null(scheme$alkylated_pattern)) character() else
    pac[grepl(scheme$alkylated_pattern, pac)]
  parent <- setdiff(pac, alk)
  epa <- intersect(scheme$epa16, analytes)
  if (length(scheme$epa16) && length(epa) < length(scheme$epa16)) {
    warning("Some EPA16 columns are missing: ",
            paste(setdiff(scheme$epa16, analytes), collapse = ", "), call. = FALSE)
  }
  list(EPA16 = epa, Parent = parent, Alkylated = alk, Total = c(parent, alk))
}
