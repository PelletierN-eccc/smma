# smma 0.1.0

First public release, accompanying Pelletier et al. (submitted to *Environmental
Pollution*).

## Package

* Weighted non-negative least-squares fitting of receptor-model factor profiles as
  mixtures of library source profiles (`smma_prepare()`, `smma_fit()`,
  `smma_forward_select()`).
* Analyte bootstrap with fit-quality-weighted medians, percentile intervals and
  selection rates, per source and per category (`smma_bootstrap()`), best-guess
  mixtures (`smma_best_mixture()`) and group shares (`smma_group_share()`).
* Propagation of source weights to samples for any analyte grouping
  (`attribute_sources()`, `summarise_attribution()`, `attribution_shares()`).
* Profile distances, diagnostic ratios and alkylation indices
  (`profile_distances()`, `diagnostic_ratios()`, `alkylation_indices()`); analyte
  naming is set by an analyte scheme (`default_analyte_scheme()`, `analyte_scheme()`).
* Library sensitivity tools (`fill_library_gaps()`, `update_library()`,
  `smma_simulate()`), SMMA on raw samples (`smma_samples()`), site classification
  and spatial distance-decay regressions (`classify_sites()`,
  `distance_to_features()`, `fit_distance_decay()`).
* Example data: the four-factor EPA PMF 5.0 solution for 211 British Columbia
  intertidal sediments and a 371-profile PAC/biomarker source library.

## Analysis scripts (`analysis/`)

* Reproduce every SMMA-dependent result of the article. Per-analyte Kow weights are
  proportional to log10 Kow; `kow_double_log = TRUE` reproduces the weighting of the
  original submission draft, in which tabulated log10 Kow values were
  log10-transformed a second time.
