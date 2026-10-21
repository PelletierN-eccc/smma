# Biomarker-coverage sensitivity test (Sections 2.8, 3.1.6; Supplementary
# Section S12). Refits SMMA with literature-informed low / central / high
# biomarker values for the library profiles that report none, with and
# without three capped-ratio profiles. Slow: 8 variants x B draws.
source("analysis/R/00_config.R")

profiles <- read_factor_profiles(cfg$ex("bc_factor_profiles.csv"))
profiles <- profiles[profiles$factor %in% c("Factor 1", "Factor 2"), ]
sn       <- read_sn(cfg$ex("bc_factor_sn.csv"))
library0 <- read_source_library(cfg$ex("source_library.csv"))
kow      <- read_kow(cfg$ex("kow.csv"))

# profiles with no non-zero PAC data, removed from every variant
drop10 <- c("BAL.SP2", "TNB 5", "TNB 8", "TNB 7", "BAL.SKR1", "BAL.SKR2",
            "BAL.SKR3", "BAL.SKR5", "TNB 6", "LKH.SC3")
cap3   <- c("DEK 2", "TNB 2", "TNB 4")
petro_categories <- c("Diesel (fresh)", "Diesel (weathered)", "Fuel oil (fresh)",
                      "Fuel oil (weathered)", "Lubricating oil (fresh)",
                      "Lubricating oil (weathered)", "Asphalt")
is_petro <- function(category) ifelse(category %in% petro_categories, "Petrogenic", "Other")

bounds <- lapply(c(low = "low", central = "central", high = "high"), function(b)
  read.csv(file.path(cfg$data_dir, sprintf("biomarker_bounds_%s.csv", b)), check.names = FALSE))
bounds <- lapply(bounds, function(x) x[!trimws(x$source) %in% drop10, ])

variants <- expand.grid(magnitude = c("baseline", "low", "central", "high"),
                        drop_cap = c(FALSE, TRUE), stringsAsFactors = FALSE)
results <- lapply(seq_len(nrow(variants)), function(k) {
  v <- variants[k, ]
  lib <- library0[!trimws(library0$source) %in% drop10, ]
  if (v$magnitude != "baseline") lib <- update_library(lib, bounds[[v$magnitude]])
  if (v$drop_cap) lib <- lib[!trimws(lib$source) %in% cap3, ]
  inp <- smma_prepare(profiles, lib, min_shared = cfg$min_shared)
  bt <- smma_bootstrap(inp, B = cfg$B, seed = cfg$seed, kow = kow, sn = sn,
                       kow_args = cfg$kow_args, keep_draws = TRUE)
  sh <- smma_group_share(bt, is_petro, categories = inp$categories)
  ps <- bt$per_source[order(bt$per_source$factor, -bt$per_source$median), ]
  top <- do.call(rbind, lapply(split(ps, ps$factor), utils::head, 10))
  list(share = cbind(magnitude = v$magnitude, capped_profiles_excluded = v$drop_cap,
                     n_sources = ncol(inp$L), sh[sh$group == "Petrogenic", ]),
       top = cbind(magnitude = v$magnitude, capped_profiles_excluded = v$drop_cap, top))
})
res <- do.call(rbind, lapply(results, `[[`, "share"))
out(res, "biomarker_sensitivity_petrogenic_share.csv")   # Table S12
out(do.call(rbind, lapply(results, `[[`, "top")), "biomarker_sensitivity_top10_sources.csv")
print(res)
