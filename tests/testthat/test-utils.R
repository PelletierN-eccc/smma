test_that("weighted_quantile matches unweighted quantile type 1", {
  x <- c(5, 1, 3, 2, 4)
  expect_equal(smma:::weighted_quantile(x, rep(1, 5), 0.5), 3)
  expect_equal(smma:::weighted_quantile(c(1, 10), c(0.9, 0.1), 0.5), 1)
})

test_that("normalize_rows gives unit sums and NA rows for empty input", {
  m <- rbind(c(1, 3, NA), c(NA, NA, NA), c(0, 0, 0))
  n <- smma:::normalize_rows(m)
  expect_equal(sum(n[1, ], na.rm = TRUE), 1)
  expect_true(all(is.na(n[2, ])))
  expect_true(all(is.na(n[3, ])))
})

test_that("kow_weights have mean 1 and reverse for water", {
  k <- c(a = 3, b = 5, c = 7)
  w <- kow_weights(names(k), k)
  expect_equal(mean(w), 1)
  expect_true(w["c"] > w["a"])
  ww <- kow_weights(names(k), k, matrix = "water")
  expect_true(ww["a"] > ww["c"])
  expect_equal(unname(kow_weights("z", k)), 1)
})
