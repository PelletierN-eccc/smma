% Generated from roxygen comments by tools/roxy2rd.py
\name{ests_analyte_scheme}
\alias{ests_analyte_scheme}
\title{Analyte naming scheme}
\usage{
ests_analyte_scheme()
}
\description{
All functions that need to know *which* column is phenanthrene, which
columns are alkylated homologues, or which belong to the EPA 16
priority list take an analyte scheme. \code{ests_analyte_scheme()} returns
the column-naming convention used by the Emergency Science and
Technology Section (ESTS) of Environment and Climate Change Canada
(e.g. \code{pah_c0_p} = phenanthrene, \code{pah_c2_p} = C2-phenanthrenes,
\code{biom_c30_a_b_h} = C30 alpha-beta hopane). To use your own column
names, copy this list and edit the entries, or build one with
\code{\link{analyte_scheme}()}.
}
\value{
A list with elements:
\describe{
\item{pac_pattern}{Regular expression identifying PAC columns.}
\item{biomarker_pattern, tph_pattern}{Patterns for biomarker and TPH columns.}
\item{alkylated_pattern}{Pattern identifying alkylated (C1-C4) homologues among PAC columns.}
\item{exclude_pattern}{Pattern for summary columns (sums, indices) never treated as analytes.}
\item{epa16}{Character vector of the 16 EPA priority PAH columns.}
\item{ratio}{Named character vector mapping compounds used in
diagnostic ratios (\code{ant}, \code{phe}, \code{fla}, \code{pyr}, \code{baa}, \code{chr}, \code{icdp}, \code{bghip})
to column names.}
\item{homologues}{Named list of homologue series (parent first, then C1, C2, ...).}
\item{labels}{Optional named character vector of display labels.}
}
}
\examples{
sch <- ests_analyte_scheme()
sch$ratio
}
