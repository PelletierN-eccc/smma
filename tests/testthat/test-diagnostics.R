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

test_that("distance_to_features and fit_distance_decay work on projected data", {
  skip_if_not_installed("sf")
  skip_if_not_installed("nlme")
  set.seed(1)
  n <- 40
  xy <- cbind(stats::runif(n, 0, 1e4), stats::runif(n, 0, 1e4))
  pts <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]), coords = c("x", "y"), crs = 3005)
  hub <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(5000, 5000)), crs = 3005))
  sq <- sf::st_polygon(list(rbind(c(0, 0), c(1000, 0), c(1000, 1000), c(0, 1000), c(0, 0))))
  box <- sf::st_sf(geometry = sf::st_sfc(sq, crs = 3005))
  d <- distance_to_features(pts, list(hub = hub, box = box))
  expect_equal(names(d), c("dist_hub", "dist_box"))
  expect_equal(d$dist_hub, sqrt((xy[, 1] - 5000)^2 + (xy[, 2] - 5000)^2))
  df <- data.frame(F1 = 10 * exp(-d$dist_hub / 3000) + stats::rnorm(n, 0, 0.5), d,
                   X = xy[, 1], Y = xy[, 2])
  df <- rbind(df, df[1:3, ])
  r <- fit_distance_decay(df, "F1", "dist_hub")
  expect_equal(r$n, n)
  expect_lt(r$estimate, 0)
  expect_lt(r$p_value, 0.001)
})
