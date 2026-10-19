# smma: Source Mixture Matching Analysis

**smma** is an R package for interpreting receptor-model factors (for example from
EPA PMF 5.0) with an external library of measured source profiles. Each factor
profile is fitted as a non-negative mixture of library sources. Analyte resampling
gives uncertainty on that fit. The result can then be propagated back to every
sample for any grouping of analytes.

The method was developed at the Emergency Science and Technology Section (ESTS) of
Environment and Climate Change Canada for polycyclic aromatic compounds (PACs) and
petroleum biomarkers in shoreline sediments. Nothing in the package is tied to
those analytes. Column names, analyte groupings, source categories and weighting
choices are all arguments.

## What it does

| Step | Function |
|---|---|
| Read EPA PMF output, S/N, Kow and a source library (CSV or Excel) | `read_factor_profiles()`, `read_factor_contributions()`, `read_sn()`, `read_kow()`, `read_source_library()` |
| Align analytes and normalize profiles | `smma_prepare()` |
| Weighted NNLS mixture fit (point estimate) | `smma_fit()`, `smma_forward_select()` |
| Bootstrapped fit with fit-quality-weighted summaries | `smma_bootstrap()`, `smma_best_mixture()`, `smma_group_share()` |
| Source-resolved concentrations per sample and analyte grouping | `attribute_sources()`, `summarise_attribution()`, `attribution_shares()` |
| Other lines of evidence | `profile_distances()`, `best_matches()`, `diagnostic_ratios()`, `alkylation_indices()` |
| Library sensitivity tests | `fill_library_gaps()`, `update_library()`, `smma_simulate()` |
| SMMA on raw samples (no receptor model) | `smma_samples()` |
| Spatial context | `classify_sites()`, `distance_to_features()`, `relative_contributions()`, `fit_distance_decay()` |
| Plots | `plot_smma()`, `plot_reconstruction()` |

## Installation

```r
# install.packages("remotes")
remotes::install_git("https://gccode.ssc-spc.gc.ca/PelletierN-ECCC/smma")
```

Required: R ≥ 4.1 and the `nnls` package. Optional: `readxl` (Excel input),
`ggplot2` (plots), `sf` and `nlme` (spatial tools).

## Quick start

```r
library(smma)

profiles <- read_factor_profiles("my_pmf_output.xlsx", sheet = "Profiles_conc")
contrib  <- read_factor_contributions("my_pmf_output.xlsx", sheet = "Contributions")
sn       <- read_sn("my_pmf_output.xlsx", sheet = "SN")      # optional
kow      <- read_kow("kow.csv")                               # optional, log10 Kow
lib      <- read_source_library("source_library.csv")         # source, category, analytes...

inp  <- smma_prepare(profiles, lib)
boot <- smma_bootstrap(inp, B = 2000, seed = 1, kow = kow, sn = sn)
mix  <- smma_best_mixture(boot)
mix$categories                                                # what each factor is made of

att <- attribute_sources(contrib, profiles, mix$sources,
                         groupings = analyte_groupings(names(profiles)))
summarise_attribution(att)                                    # mean, 95% CI, % share
```

`vignette("smma")` walks through the example data shipped with the package (a
four-factor PMF solution for 211 British Columbia sediments and a 371-profile
PAC/biomarker library).

### Input conventions

- Analyte columns must have the same names in the factor profiles, the source library
  and the S/N and Kow tables. The default pattern is `^(pah_|tph_|biom_)`. Pass
  `analytes =` or `analyte_pattern =` for other naming.
- Library values may be missing. They are set to zero in the fit, and
  `print(smma_prepare(...))` reports library completeness. Test that assumption with
  `fill_library_gaps()` or `update_library()`.
- Diagnostic ratios, homologue series and the EPA16 list are defined in an analyte
  scheme. `default_analyte_scheme()` matches the example data; build your own with
  `analyte_scheme()`.

## Reproducing the BC intertidal sediment study

`analysis/` contains the scripts, settings and data used for:

> Pelletier N., Yang Z., Laforest S., Wojnarowicz P., Hollebone B.P., Beaulac V.
> Source apportionment of PACs and petroleum biomarkers in British Columbia
> intertidal sediments: an urban-to-remote baseline for spill response. Submitted to
> *Environmental Pollution*.

See [`analysis/README.md`](analysis/README.md). The field data are also published in
the ECCC Data Catalogue as
[ESTS Shoreline surveys data](https://data-donnees.az.ec.gc.ca/data/water/scientificknowledge/ests-shoreline-surveys-data?lang=en).

## Method summary

For factor profile **f** and library matrix **L** (analytes × sources), both
normalized to unit sum over their shared analytes, SMMA solves

    minimize || sqrt(W) (f − L s) ||²   subject to  s ≥ 0,   then  s ← s / Σ s

where **W** is diagonal. Its entries are the product of a log10 Kow weight (scaled to
mean 1) and the factor-specific S/N from the receptor model. In each of B bootstrap
draws, the shared analytes are resampled with replacement and all profiles are
renormalized on the resampled set before refitting. Draws are weighted by the
inverse of their weighted relative RMSE. Results are summarized as weighted medians,
95% weighted percentile intervals and selection rates. The best-guess mixture keeps
sources with selection rate ≥ 0.20 or median weight ≥ 0.01 and renormalizes their
weights.

## Citation

If you use smma, please cite the article above and this software (see
`CITATION.cff`).

## Licence

MIT. Copyright © His Majesty the King in Right of Canada, as represented by the
Minister of Environment and Climate Change, 2026.
