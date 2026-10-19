ex <- smma_example()

test_that("example data load", {
  expect_equal(nrow(ex$profiles), 4)
  expect_equal(nrow(ex$contributions), 211)
  expect_true(all(c("source", "category") %in% names(ex$library)))
})

test_that("SMMA recovers a known noise-free mixture", {
  lib <- ex$library
  src <- c(lib$source[lib$category == "Fuel oil (fresh)"][1],
           lib$source[lib$category == "Biomass combustion"][1])
  sim <- smma_simulate(lib, list(F1 = stats::setNames(c(0.7, 0.3), src)))
  fit <- smma_fit(smma_prepare(sim, lib))
  w <- stats::setNames(fit$weight, fit$source)
  expect_gt(w[src[1]] + w[src[2]], 0.9)
  expect_equal(unname(w[src[1]] / (w[src[1]] + w[src[2]])), 0.7, tolerance = 0.05)
})

test_that("bootstrap returns consistent summaries", {
  inp <- smma_prepare(ex$profiles, ex$library)
  bt <- smma_bootstrap(inp, B = 20, seed = 1, kow = ex$kow, sn = ex$sn, keep_draws = TRUE)
  expect_s3_class(bt, "smma_bootstrap")
  expect_setequal(unique(bt$per_source$factor), ex$profiles$factor)
  expect_true(all(bt$per_source$sel_rate >= 0 & bt$per_source$sel_rate <= 1))
  # draws sum to one
  expect_true(all(abs(rowSums(bt$draws[[1]]) - 1) < 1e-8))
  # same seed, same answer
  bt2 <- smma_bootstrap(inp, B = 20, seed = 1, kow = ex$kow, sn = ex$sn)
  expect_equal(bt$per_source$median, bt2$per_source$median)
  mix <- smma_best_mixture(bt)
  s <- tapply(mix$sources$weight, mix$sources$factor, sum)
  expect_true(all(abs(s - 1) < 1e-8))
  sh <- smma_group_share(bt, function(cat) ifelse(grepl("oil|Diesel \\(", cat), "Petro", "Other"),
                         categories = inp$categories)
  expect_true(all(c("Petro", "Other") %in% sh$group))
})

test_that("attribution conserves the factor-reconstructed burden", {
  inp <- smma_prepare(ex$profiles, ex$library)
  mix <- smma_best_mixture(smma_bootstrap(inp, B = 10, seed = 2))
  g <- analyte_groupings(names(ex$profiles))
  att <- attribute_sources(ex$contributions, ex$profiles, mix$sources, g["EPA16"])
  tot <- c(tapply(att$value, att$sample, sum))  # c() drops the 1-d array dim
  prof <- ex$profiles
  fsum <- rowSums(as.matrix(prof[, g$EPA16]))
  names(fsum) <- prof$factor
  G <- as.matrix(ex$contributions[, names(fsum)]); G[G < 0] <- 0
  expect_equal(unname(tot[ex$contributions$sample]), unname(as.vector(G %*% fsum)), tolerance = 1e-8)
  # every sample x category present
  expect_equal(nrow(att), length(unique(att$sample)) * length(unique(att$category)))
  s <- summarise_attribution(att, B = 50, seed = 1)
  expect_equal(sum(s$share_pct), 100)
})

test_that("bootstrap handles a single category and a single draw", {
  lib <- ex$library
  lib$category <- "Only"
  bt <- smma_bootstrap(smma_prepare(ex$profiles[1, ], lib), B = 5, seed = 1)
  expect_s3_class(bt$per_category, "data.frame")
  expect_equal(bt$per_category$category, "Only")
  expect_equal(bt$per_category$median, 1)
  bt1 <- smma_bootstrap(smma_prepare(ex$profiles, ex$library), B = 1, seed = 1)
  expect_s3_class(bt1$per_source, "data.frame")
  expect_equal(nrow(bt1$per_source), 4 * nrow(ex$library))
})

test_that("category summaries equal the summed source draws", {
  inp <- smma_prepare(ex$profiles, ex$library)
  bt <- smma_bootstrap(inp, B = 30, seed = 3, keep_draws = TRUE)
  D <- bt$draws[["Factor 2"]]
  w <- bt$draw_weights[["Factor 2"]]
  cats <- inp$categories$category[match(colnames(D), inp$categories$source)]
  pc <- bt$per_category[bt$per_category$factor == "Factor 2", ]
  man <- vapply(pc$category, function(cc)
    smma:::weighted_quantile(rowSums(D[, cats == cc, drop = FALSE]), w)[1], numeric(1))
  expect_equal(unname(man), pc$median)
  expect_equal(rownames(bt$per_source), as.character(seq_len(nrow(bt$per_source))))
})

test_that("attribution works with one factor and one grouping", {
  inp <- smma_prepare(ex$profiles, ex$library)
  mix <- smma_best_mixture(smma_bootstrap(inp, B = 10, seed = 2))
  g <- analyte_groupings(names(ex$profiles))
  m3 <- mix$sources[mix$sources$factor == "Factor 3", ]
  att <- attribute_sources(ex$contributions, ex$profiles, m3, g["EPA16"])
  G3 <- pmax(ex$contributions[["Factor 3"]], 0)
  f3 <- sum(ex$profiles[ex$profiles$factor == "Factor 3", g$EPA16])
  tot <- c(tapply(att$value, att$sample, sum))
  expect_equal(unname(tot[ex$contributions$sample]), G3 * f3, tolerance = 1e-8)
  sh <- attribution_shares(att)
  expect_true(all(c("n", "median_pct") %in% names(sh)))
  expect_type(sh$median_pct, "double")
})
