# Factor identification: distances, diagnostic ratios, alkylation indices,
# bootstrapped SMMA and best-guess mixtures (Sections 2.5-2.7, 3.1).
source("analysis/R/00_config.R")

profiles <- read_factor_profiles(cfg$ex("bc_factor_profiles.csv"))
sn       <- read_sn(cfg$ex("bc_factor_sn.csv"))
srclib  <- read_source_library(cfg$ex("source_library.csv"))
kow      <- read_kow(cfg$ex("kow.csv"))
catmap   <- read.csv(file.path(cfg$data_dir, "source_category_mapping.csv"), check.names = FALSE)

# --- Line of evidence 2-3: diagnostic ratios and alkylation patterns -------
out(diagnostic_ratios(profiles), "factor_diagnostic_ratios.csv")
out(alkylation_indices(profiles), "factor_alkylation_indices.csv")

# --- Line of evidence 4: distance to library profiles ----------------------
wfun <- function(fid, a) kow_weights(a, kow, input_scale = "log10") *
  sn_weights(sn, fid, a)
dist <- profile_distances(profiles, srclib, min_shared = cfg$min_shared, weights = wfun)
out(dist, "factor_library_distances.csv")
out(best_matches(dist, n = 3), "factor_best_matches.csv")

# --- Line of evidence 5: SMMA -----------------------------------------------
inp <- smma_prepare(profiles, srclib, min_shared = cfg$min_shared)
print(inp)
boot <- smma_bootstrap(inp, B = cfg$B, seed = cfg$seed, kow = kow, sn = sn,
                       kow_args = cfg$kow_args, min_analytes = cfg$min_shared,
                       keep_draws = FALSE, progress = TRUE)
print(boot)
mix <- smma_best_mixture(boot, min_sel_rate = cfg$min_sel_rate, min_weight = cfg$min_weight)

# regroup library categories into the 9 (main text) and 17 (SI) categories
mix$sources$category_maintext <- catmap$category_maintext[match(mix$sources$source, catmap$source)]
mix$sources$category_detailed <- catmap$category_detailed[match(mix$sources$source, catmap$source)]
agg <- function(col) {
  a <- stats::aggregate(mix$sources$weight, list(factor = mix$sources$factor,
                                                 category = mix$sources[[col]]), sum)
  names(a)[3] <- "weight"; a[order(a$factor, -a$weight), ]
}
out(boot$per_source, "smma_bootstrap_per_source.csv")
out(boot$per_category, "smma_bootstrap_per_library_category.csv")
out(boot$fit, "smma_bootstrap_fit.csv")
out(mix$sources, "smma_best_mixture_sources.csv")
out(agg("category_maintext"), "smma_best_mixture_9cat.csv")     # Figure 3 (main-text scheme)
out(agg("category_detailed"), "smma_best_mixture_17cat.csv")    # Figure 3 / Table S14
saveRDS(list(input = inp, boot = boot, mix = mix), file.path(cfg$out_dir, "smma_factors.rds"))
