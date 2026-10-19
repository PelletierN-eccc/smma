% Generated from roxygen comments by tools/roxy2rd.py
\name{smma_prepare}
\alias{smma_prepare}
\title{Prepare inputs for Source Mixture Matching Analysis}
\usage{
smma_prepare(factor_profiles, source_library, factor_col = "factor", source_col = "source", category_col = "category", analytes = NULL, analyte_pattern = "^(pah_|tph_|biom_)", min_shared = 12)
}
\arguments{
\item{factor_profiles}{Data frame with one row per factor (see
\code{\link{read_factor_profiles}()}), or a sample table when SMMA is applied
directly to samples (see \code{\link{smma_samples}()}).}
\item{source_library}{Data frame with one row per source profile (see
\code{\link{read_source_library}()}).}
\item{factor_col}{Column in \code{factor_profiles} identifying factors.}
\item{source_col,category_col}{Columns in \code{source_library} holding the
source identifier and its category.}
\item{analytes}{Optional character vector of analyte columns to use.
By default, all columns matching \code{analyte_pattern} that occur in both
tables.}
\item{analyte_pattern}{Regular expression used when \code{analytes} is \code{NULL}.}
\item{min_shared}{Minimum number of shared analytes required.}
}
\description{
Aligns factor profiles and a source-profile library on their shared
analytes and normalizes every profile to unit sum over those analytes.
}
\value{
An object of class \code{smma_input}: a list with
\describe{
\item{F}{Factor matrix (factors x analytes), rows sum to 1.}
\item{L}{Library matrix (analytes x sources), columns normalized,
missing values set to 0.}
\item{L_missing}{Logical matrix flagging values that were missing in
the library before zero-filling.}
\item{categories}{Data frame \code{source}, \code{category}.}
\item{analytes}{The analytes used.}
}
}
\examples{
ex <- smma_example()
inp <- smma_prepare(ex$profiles, ex$library)
inp
}
