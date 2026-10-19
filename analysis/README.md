# Reproducing the British Columbia intertidal sediment analysis

These scripts reproduce the receptor-model interpretation in Pelletier et al.
(submitted to *Environmental Pollution*). They use the **smma** package, so install
it first (see the main README). Run them from the repository root:

```r
source("analysis/R/01_smma_factors.R")        # Sections 2.5–2.7 and 3.1; Figure 3; Tables 1, S14
source("analysis/R/02_source_attribution.R")  # Section 3.3; Table 2; Table S13; Figure S26
source("analysis/R/03_spatial_regression.R")  # Section 3.2; Section S15 (needs DFO layers)
source("analysis/R/04_biomarker_sensitivity.R") # Sections 2.8 and 3.1.6; Section S12 (slow)
source("analysis/R/05_smma_direct_samples.R") # Section 3.6; Section S11; Table S10 (after 02)
source("analysis/R/06_burden_and_texture.R")  # Sections 3.4-3.5; Tables S4, S6 (after 02)
source("analysis/R/08_figures.R")             # Figure 3; Figures S4-S15, S21, S24-S26
source("analysis/R/09_tables.R")              # Table 2; Tables S4, S6, S10, S12-S14 (CSV, xlsx)
```

Figures are written to `analysis/figures/` as 600-dpi LZW TIFF (17 cm wide) with PNG
previews; tables to `analysis/output/tables/`. All figures use one colour per
nine-category source family, in the same order.

`07_compare_kow_weighting.R` tabulates every weight-dependent number under both Kow
weightings (run scripts 01-06 with `options(smma.kow_double_log = TRUE)` and without it
first; it is run with `Rscript` from the repository root).

All settings are in `analysis/R/00_config.R`: B = 10,000 draws, seed = 42 and the
retention thresholds. Outputs go to `analysis/output/`; for a quick run, set
`options(smma.B = 1000)` before sourcing a script.

## Data

| File | Content |
|---|---|
| `inst/extdata/bc_factor_profiles.csv` | EPA PMF 5.0 four-factor profiles (concentration units, ng g⁻¹ dw) |
| `inst/extdata/bc_factor_contributions.csv` | Factor contributions for the 211 samples |
| `inst/extdata/bc_factor_sn.csv` | Factor-specific signal-to-noise per analyte |
| `inst/extdata/source_library.csv` | 371 source profiles (13 references; see Table S1 of the article) |
| `inst/extdata/kow.csv` | log10 Kow per analyte (placeholder values for TPH and biomarkers) |
| `analysis/data/pmf_input_samples.csv`, `pmf_input_uncertainties.csv` | PMF input concentration and uncertainty matrices |
| `analysis/data/sites.csv` | Sample coordinates and campaign |
| `analysis/data/population_centres_bc.csv` | Coastal BC population centres (2021 Census) |
| `analysis/data/source_category_mapping.csv` | Source → 9-category (main text) and 17-category (SI) schemes |
| `analysis/data/sample_texture.csv` | % fines (< 63 µm) and TOC per sample, built from the ECCC Data Catalogue release by `prep_sample_texture.R` |
| `analysis/data/biomarker_bounds_*.csv` | Low, central and high bounding biomarker values used in the sensitivity test (not measurements; see `biomarker_bounds_README.txt`) |

The complete field dataset (sediment chemistry, grain size, TOC and site metadata
for all ESTS pre-SCAT surveys) is published in the ECCC Data Catalogue:
<https://data-donnees.az.ec.gc.ca/data/water/scientificknowledge/ests-shoreline-surveys-data?lang=en>.

The two infrastructure layers for script 03 are Fisheries and Oceans Canada open
data. They are not redistributed here; the script header has download links.

## Notes on reproducing the published numbers

- **Kow weighting.** `kow.csv` tabulates log10 Kow. By default (`kow_double_log = FALSE`
  in `00_config.R`) the per-analyte Kow weights are proportional to log10 Kow, scaled
  to a mean of 1, as described in the Methods. In the original submission draft these
  values went through a second log10 transform, so the weights were proportional to
  log10(log10 Kow). Set `kow_double_log = TRUE` to reproduce the weighting of the
  original submission draft; outputs then go to `analysis/output_double_log/`.
- **Diagnostic ratios.** `diagnostic_ratios()` computes Ant/(Ant+Phe) with
  phenanthrene (`pah_c0_p`). The SI heat map (Figure S4) was produced with an earlier
  script that used perylene (`pah_pe`) in that denominator.
- **Random numbers.** The bootstrap reproduces the per-factor seeding of the original
  scripts, so results are identical for the same R version, seed, library and input
  order.
