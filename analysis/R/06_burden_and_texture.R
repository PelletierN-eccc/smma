# Source shares versus PAC burden and sediment texture (Sections 3.4, 3.5;
# Tables S4, S6). Needs the output of scripts 01-02 and
# analysis/data/sample_texture.csv (see prep_sample_texture.R).
source("analysis/R/00_config.R")

samples  <- read.csv(file.path(cfg$data_dir, "pmf_input_samples.csv"), check.names = FALSE)
texture  <- read.csv(file.path(cfg$data_dir, "sample_texture.csv"))
att      <- read.csv(file.path(cfg$out_dir, "attribution_by_sample_9cat.csv"))   # from script 02
att      <- att[att$grouping == "EPA16", ]
epa16    <- analyte_groupings(names(samples))$EPA16

# per-sample attributed concentrations and percentage shares, samples x categories,
# plus the three petroleum categories combined
conc <- unclass(stats::xtabs(value ~ sample + category, data = att))
petro <- grep("petroleum|Other oils", colnames(conc), value = TRUE)
pct <- 100 * conc / rowSums(conc)
pct <- cbind(pct, `All petroleum` = rowSums(pct[, petro, drop = FALSE]))
cats <- colnames(conc)
sum16 <- stats::setNames(rowSums(samples[, epa16], na.rm = TRUE), samples$pit_code)[rownames(pct)]

# Section 3.4: Spearman correlation of each share with the measured sum of EPA16
rho <- do.call(rbind, lapply(colnames(pct), function(k) {
  ct <- suppressWarnings(stats::cor.test(pct[, k], sum16, method = "spearman"))
  data.frame(category = k, n = length(sum16), rho = unname(ct$estimate), p_value = ct$p.value)
}))
out(rho, "share_vs_sum16_spearman.csv")

# Table S4: Tukey outliers of log10 sum16 and their two largest attributed categories
l16 <- log10(sum16)
thr <- stats::quantile(l16, 0.75) + 1.5 * stats::IQR(l16)
outl <- names(sum16)[l16 > thr]
s4 <- do.call(rbind, lapply(outl, function(s) {
  r <- sort(pct[s, cats], decreasing = TRUE)
  data.frame(sample = s, sum16 = sum16[[s]], first = names(r)[1], first_pct = r[[1]],
             second = names(r)[2], second_pct = r[[2]], petroleum_pct = pct[s, "All petroleum"])
}))
s4 <- s4[order(-s4$sum16), ]
out(s4, "tukey_outliers_table_s4.csv")

# elevation at the outliers: median attributed concentration of each category at
# the outliers over its median at the remaining samples
rest <- setdiff(rownames(conc), outl)
elev <- data.frame(category = cats,
                   ratio = apply(conc, 2, function(v) stats::median(v[outl]) / stats::median(v[rest])))
out(elev, "tukey_outliers_elevation.csv")

# top 5% of samples by sum16: largest single category
top <- names(sort(sum16, decreasing = TRUE))[seq_len(ceiling(0.05 * length(sum16)))]
lead <- data.frame(sample = top, sum16 = sum16[top],
                   largest = cats[apply(pct[top, cats, drop = FALSE], 1, which.max)])
out(lead, "top5pct_largest_source.csv")

# Section 3.5, Table S6: Pearson correlation of each share with % fines
fines <- texture$fines_pct[match(rownames(pct), texture$sample)]
ok <- is.finite(fines)
s6 <- do.call(rbind, lapply(colnames(pct), function(k) {
  ct <- stats::cor.test(pct[ok, k], fines[ok])
  data.frame(category = k, n = sum(ok), r = unname(ct$estimate), p_value = ct$p.value)
}))
out(s6, "share_vs_fines_table_s6.csv")

# Section 3.5: measured sum16 against texture (does not depend on SMMA)
toc <- texture$toc_pct[match(rownames(pct), texture$sample)]
ok2 <- is.finite(toc) & toc > 0
tex <- rbind(
  data.frame(variable = "% fines", n = sum(ok), r = stats::cor(l16[ok], fines[ok]),
             p_value = stats::cor.test(l16[ok], fines[ok])$p.value),
  data.frame(variable = "log10 TOC", n = sum(ok2), r = stats::cor(l16[ok2], log10(toc[ok2])),
             p_value = stats::cor.test(l16[ok2], log10(toc[ok2]))$p.value))
out(tex, "sum16_vs_texture.csv")

cat(sprintf("Tukey threshold: %.0f ng/g\n", 10^thr))
print(rho, digits = 3, row.names = FALSE)
print(s4, digits = 3, row.names = FALSE)
print(elev, digits = 3, row.names = FALSE)
print(table(lead$largest))
print(s6, digits = 3, row.names = FALSE)
print(tex, digits = 3, row.names = FALSE)
