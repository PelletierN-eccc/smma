# Source shares versus PAC burden and sediment texture (Sections 3.4, 3.5;
# Tables S4, S6). Needs the output of scripts 01-02 and
# analysis/data/sample_texture.csv (see prep_sample_texture.R).
source("analysis/R/00_config.R")

res      <- readRDS(file.path(cfg$out_dir, "smma_factors.rds"))
profiles <- read_factor_profiles(cfg$ex("bc_factor_profiles.csv"))
contrib  <- read_factor_contributions(cfg$ex("bc_factor_contributions.csv"))
samples  <- read.csv(file.path(cfg$data_dir, "pmf_input_samples.csv"), check.names = FALSE)
texture  <- read.csv(file.path(cfg$data_dir, "sample_texture.csv"))
groupings <- analyte_groupings(names(profiles))

# Broad categories of Tables S4 and S6 (library categories not listed are "Other").
recode_broad <- function(x) {
  out <- rep("Other", length(x))
  out[x %in% c("Diesel combustion", "Gasoline combustion")] <- "Diesel and gasoline combustion"
  out[x %in% c("Fuel oil (weathered)", "Lubricating oil (fresh)")] <- "Petrogenic"
  out[x %in% c("Street dust near unsealed lot", "Individual lot dust, sealed lot",
               "Coal-tar sealant and creosote")] <- "Street dust"
  out[x == "Industrial boiler emission"] <- "Industrial / coal"
  out[x == "Biomass combustion"] <- "Biomass combustion"
  out
}
cmap <- data.frame(source = res$input$categories$source,
                   category = recode_broad(res$input$categories$category))
att <- attribute_sources(contrib, profiles, res$mix$sources, groupings["EPA16"], category_map = cmap)

# per-sample percentage shares (EPA16 basis), samples x categories
wide <- stats::xtabs(value ~ sample + category, data = att)
pct <- 100 * unclass(wide) / rowSums(wide)
cats <- colnames(pct)
sum16 <- stats::setNames(rowSums(samples[, groupings$EPA16], na.rm = TRUE), samples$pit_code)
sum16 <- sum16[rownames(pct)]

# Section 3.4: Spearman correlation of each share with measured sum of EPA16
rho <- do.call(rbind, lapply(cats, function(k) {
  ct <- suppressWarnings(stats::cor.test(pct[, k], sum16, method = "spearman"))
  data.frame(category = k, n = length(sum16), rho = unname(ct$estimate), p_value = ct$p.value)
}))
out(rho, "share_vs_sum16_spearman.csv")

# Table S4: Tukey outliers of log10 sum16 and their two largest attributed sources
l16 <- log10(sum16)
thr <- stats::quantile(l16, 0.75) + 1.5 * stats::IQR(l16)
outl <- names(sum16)[l16 > thr]
s4 <- do.call(rbind, lapply(outl, function(s) {
  r <- sort(pct[s, ], decreasing = TRUE)
  data.frame(sample = s, sum16 = sum16[[s]], first = names(r)[1], first_pct = r[[1]],
             second = names(r)[2], second_pct = r[[2]])
}))
s4 <- s4[order(-s4$sum16), ]
out(s4, "tukey_outliers_table_s4.csv")

# top 5% of samples by sum16: largest single source
top <- names(sort(sum16, decreasing = TRUE))[seq_len(ceiling(0.05 * length(sum16)))]
lead <- data.frame(sample = top, sum16 = sum16[top],
                   largest = cats[apply(pct[top, , drop = FALSE], 1, which.max)])
out(lead, "top5pct_largest_source.csv")

# Section 3.5, Table S6: Pearson correlation of each share with % fines
fines <- texture$fines_pct[match(rownames(pct), texture$sample)]
ok <- is.finite(fines)
s6 <- do.call(rbind, lapply(cats, function(k) {
  ct <- stats::cor.test(pct[ok, k], fines[ok])
  data.frame(category = k, n = sum(ok), r = unname(ct$estimate), p_value = ct$p.value)
}))
out(s6, "share_vs_fines_table_s6.csv")

cat(sprintf("Tukey threshold: %.0f ng/g\n", 10^thr))
print(rho, digits = 3, row.names = FALSE)
print(s4, digits = 3, row.names = FALSE)
print(table(lead$largest))
print(s6, digits = 3, row.names = FALSE)
