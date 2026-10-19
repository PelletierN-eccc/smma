\docType{package}
\name{smma-package}
\alias{smma}
\alias{smma-package}
\title{smma: Source Mixture Matching Analysis}
\description{
Interpret receptor-model factors (e.g. from EPA PMF 5.0) by fitting
each factor profile as a non-negative mixture of an external library
of measured source profiles, and propagate the result to samples.

Typical workflow:
\enumerate{
\item Read factor profiles, contributions and S/N (\code{\link{read_factor_profiles}()},
\code{\link{read_factor_contributions}()}, \code{\link{read_sn}()}) and a source library
(\code{\link{read_source_library}()}).
\item \code{\link{smma_prepare}()} aligns analytes and normalizes profiles.
\item \code{\link{smma_bootstrap}()} fits mixtures with analyte resampling;
\code{\link{smma_best_mixture}()} extracts a best-guess mixture.
\item \code{\link{attribute_sources}()} and \code{\link{summarise_attribution}()} attribute
sample concentrations to sources for any analyte grouping.
}
Supporting tools: \code{\link{profile_distances}()}, \code{\link{diagnostic_ratios}()},
\code{\link{alkylation_indices}()}, \code{\link{classify_sites}()}, \code{\link{distance_to_features}()},
\code{\link{fit_distance_decay}()}, \code{\link{fill_library_gaps}()}, \code{\link{smma_samples}()} and
\code{\link{smma_simulate}()}.
}
\keyword{internal}
