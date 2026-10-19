#' smma: Source Mixture Matching Analysis
#'
#' Interpret receptor-model factors (e.g. from EPA PMF 5.0) by fitting
#' each factor profile as a non-negative mixture of an external library
#' of measured source profiles, and propagate the result to samples.
#'
#' Typical workflow:
#' \enumerate{
#'   \item Read factor profiles, contributions and S/N ([read_factor_profiles()],
#'     [read_factor_contributions()], [read_sn()]) and a source library
#'     ([read_source_library()]).
#'   \item [smma_prepare()] aligns analytes and normalizes profiles.
#'   \item [smma_bootstrap()] fits mixtures with analyte resampling;
#'     [smma_best_mixture()] extracts a best-guess mixture.
#'   \item [attribute_sources()] and [summarise_attribution()] attribute
#'     sample concentrations to sources for any analyte grouping.
#' }
#' Supporting tools: [profile_distances()], [diagnostic_ratios()],
#' [alkylation_indices()], [classify_sites()], [distance_to_features()],
#' [fit_distance_decay()], [fill_library_gaps()], [smma_samples()] and
#' [smma_simulate()].
#'
#' @keywords internal
"_PACKAGE"

utils::globalVariables(".data")
