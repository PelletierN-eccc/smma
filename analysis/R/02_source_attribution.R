# Sample-level source attribution by analyte basis and site class
# (Sections 2.9, 3.3; Table 2, Table S13, Figure S26).
source("analysis/R/00_config.R")

res      <- readRDS(file.path(cfg$out_dir, "smma_factors.rds"))
profiles <- read_factor_profiles(cfg$ex("bc_factor_profiles.csv"))
contrib  <- read_factor_contributions(cfg$ex("bc_factor_contributions.csv"))
catmap   <- read.csv(file.path(cfg$data_dir, "source_category_mapping.csv"), check.names = FALSE)

# --- Site classification (Supplementary Section S14) -----------------------
sites   <- read.csv(file.path(cfg$data_dir, "sites.csv"))
centres <- read.csv(file.path(cfg$data_dir, "population_centres_bc.csv"))
sites   <- classify_sites(sites, centres, site_col = "sample")
sites   <- sites[sites$sample %in% contrib$sample, ]
print(table(sites$site_class))
out(sites, "site_classification.csv")

# --- Attribution -------------------------------------------------------------
groupings <- analyte_groupings(names(profiles))
for (scheme in c("category_maintext", "category_detailed")) {
  att <- attribute_sources(contrib, profiles, res$mix$sources, groupings,
                           category_map = data.frame(source = catmap$source,
                                                     category = catmap[[scheme]]))
  tag <- if (scheme == "category_maintext") "9cat" else "17cat"
  out(att, sprintf("attribution_by_sample_%s.csv", tag))

  # Table 2: mean (bootstrap 95% CI) and share of mean burden by site class
  tab <- summarise_attribution(att, groups = sites[, c("sample", "site_class")],
                               group_col = "site_class", B = cfg$ci_B, seed = cfg$ci_seed)
  out(tab, sprintf("attribution_mean_ci_by_site_class_%s.csv", tag))

  # Table S13: median per-sample percentage share
  shr <- attribution_shares(att, groups = sites[, c("sample", "site_class")], group_col = "site_class")
  out(shr, sprintf("attribution_median_share_by_site_class_%s.csv", tag))
}
