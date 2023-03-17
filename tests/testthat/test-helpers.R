test_that("rmssd estimates correctly", {
  expect_equal(rmssd(1:4), 1)
})

test_that("rmssd by ID estimates correctly", {
  expect_equal(
    rmssd_id(c(1:3, 1, 3, 5), rep(1:2, each = 3), FALSE),
    array(c(1, 2)), ignore_attr = TRUE)
})

test_that("sd by ID estimates correctly", {
  expect_equal(
    sd_id(c(1:3, 1, 3, 5), rep(1:2, each = 3), FALSE),
    array(c(1, 2)), ignore_attr = TRUE)
})

test_that("rolling_diff estimates correctly", {
  expect_equal(
    rolling_diff(1:4, 4),
    3)
})

test_that("rolling_diff by ID estimates correctly", {
  expect_equal(
    rolling_diff_id(c(1:4, 1, 3, 5, 7), rep(1:2, each = 4), FALSE, 4),
    array(c(3, 6)), ignore_attr = TRUE)
})
