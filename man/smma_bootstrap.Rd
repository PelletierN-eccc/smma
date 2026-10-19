% Generated from roxygen comments by tools/roxy2rd.py
\name{smma_bootstrap}
\alias{smma_bootstrap}
\title{Bootstrapped Source Mixture Matching Analysis}
\usage{
smma_bootstrap(input, B = 1000, seed = NULL, kow = NULL, sn = NULL, kow_args = list(), sn_factor_col = "factor", min_analytes = 12, keep_best = NULL, keep_min = 20, keep_draws = FALSE, progress = interactive())
}
\arguments{
\item{input}{An \code{smma_input} object from \code{\link{smma_prepare}()}.}
\item{B}{Number of bootstrap draws per factor.}
\item{seed}{Integer seed, or \code{NULL} to use the current RNG state.}
\item{kow}{Optional named Kow vector (see \code{\link{read_kow}()}).}
\item{sn}{Optional factor-specific S/N table (see \code{\link{read_sn}()}).}
\item{kow_args}{List of further arguments to \code{\link{kow_weights}()}
(e.g. \code{list(input_scale = "log10", matrix = "sediment")}).}
\item{sn_factor_col}{Factor column name in \code{sn}.}
\item{min_analytes}{Minimum number of usable analytes per factor.}
\item{keep_best}{Optional proportion in (0, 1]: keep only this share of
draws with the lowest relative RMSE (at least \code{keep_min} draws).}
\item{keep_min}{Minimum number of draws kept when \code{keep_best} is set.}
\item{keep_draws}{Store the full matrix of draws and their weights
(needed for \code{\link{smma_group_share}()} and distribution plots; can be large).}
\item{progress}{Show a progress bar.}
}
\description{
For each factor, repeatedly resamples the shared analytes with
replacement (same number of analytes, uniform probability),
renormalizes the factor and every library profile on the resampled
set, and refits the weighted NNLS mixture. Draws are weighted by the
inverse of their weighted relative RMSE, and mixing coefficients are
summarized per source and per category as weighted medians, 95\%
weighted percentile intervals, weighted means and selection rates
(weighted share of draws in which the source received a positive
coefficient).

Per-analyte weights are the product of \code{\link{kow_weights}()} and the
factor-specific S/N from \code{sn} (\code{\link{sn_weights}()}); non-finite or
non-positive products are set to 1 (neutral).

Random numbers: when \code{seed} is given, one seed per factor is derived
from it, so results for a factor do not depend on how many other
factors are fitted or on their order.
}
\value{
An object of class \code{smma_bootstrap}, a list with data frames
\code{per_source} and \code{per_category} (columns \code{factor}, \code{source} or
\code{category}, \code{median}, \code{p025}, \code{p975}, \code{mean}, \code{sel_rate}), \code{fit}
(relative RMSE summary per factor), optionally \code{draws}, and \code{settings}.
}
\examples{
ex <- smma_example()
inp <- smma_prepare(ex$profiles, ex$library)
bt <- smma_bootstrap(inp, B = 50, seed = 1)
head(bt$per_category)
}
