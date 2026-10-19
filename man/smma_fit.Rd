% Generated from roxygen comments by tools/roxy2rd.py
\name{smma_fit}
\alias{smma_fit}
\title{Fit each factor as a mixture of library sources (single fit)}
\usage{
smma_fit(input, kow = NULL, sn = NULL, kow_args = list(), sn_factor_col = "factor", min_analytes = 12)
}
\arguments{
\item{input}{An \code{smma_input} object from \code{\link{smma_prepare}()}.}
\item{kow}{Optional named Kow vector (see \code{\link{read_kow}()}).}
\item{sn}{Optional factor-specific S/N table (see \code{\link{read_sn}()}).}
\item{kow_args}{List of further arguments to \code{\link{kow_weights}()}
(e.g. \code{list(input_scale = "log10", matrix = "sediment")}).}
\item{sn_factor_col}{Factor column name in \code{sn}.}
\item{min_analytes}{Minimum number of usable analytes per factor.}
}
\description{
Solves, for every factor, the weighted non-negative least-squares
problem \code{min || sqrt(W) (f - L s) ||^2, s >= 0} on the analytes shared by
the factor and the library, then renormalizes \code{s} to sum to 1. This is
the point estimate; use \code{\link{smma_bootstrap}()} for uncertainty.
}
\value{
Data frame with columns \code{factor}, \code{source}, \code{category},
\code{weight} and \code{rel_rmse} (relative RMSE of the reconstruction).
}
\examples{
ex <- smma_example()
fit <- smma_fit(smma_prepare(ex$profiles, ex$library))
head(fit[order(-fit$weight), ])
}
