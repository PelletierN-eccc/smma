test_that("Ant/(Ant+Phe) uses phenanthrene", {
  x <- data.frame(factor = "a", pah_an = 1, pah_c0_p = 3, pah_pe = 100, pah_fl = 1, pah_py = 1,
                  pah_ba_a = 1, pah_c0_c = 1, pah_ip = 1, pah_bg_p = 1)
  r <- diagnostic_ratios(x)
  expect_equal(r$Ant_AntPhe, 0.25)
  expect_equal(r$Fla_FlaPyr, 0.5)
})

test_that("alkylation indices", {
  x <- data.frame(factor = "a", pah_c0_n = 1, pah_c1_n = 1, pah_c2_n = 1, pah_c3_n = 1, pah_c4_n = 0)
  a <- alkylation_indices(x)
  expect_equal(a$parent_fraction, 0.25)
  expect_equal(a$alkyl_parent_ratio, 3)
  expect_equal(a$c2plus_fraction, 0.5)
})

test_that("site classification", {
  towns <- data.frame(name = c("Big", "Hamlet"), latitude = c(49, 50), longitude = c(-123, -125),
                      population = c(100000, 100))
  s <- data.frame(site = c("u", "m", "r", "h"), latitude = c(49.05, 49.2, 51, 50.01),
                  longitude = c(-123, -123, -127, -125))
  out <- classify_sites(s, towns)
  expect_equal(as.character(out$site_class), c("Urban", "Mixed", "Remote", "Mixed"))
  expect_equal(haversine_km(0, 0, 0, 1), 111.19, tolerance = 0.01)
})

test_that("profile distances are sensible", {
  ex <- smma_example()
  d <- profile_distances(ex$profiles[1, ], ex$library[1:20, ])
  expect_true(all(d$cosine[is.finite(d$cosine)] <= 1 + 1e-12))
  b <- best_matches(d, n = 2)
  expect_true(all(b$rank %in% 1:2))
})
