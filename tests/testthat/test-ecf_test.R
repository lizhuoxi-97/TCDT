test_that("ecf_test works with basic inputs", {
  set.seed(123)
  x1 <- rnorm(50)
  x2 <- rnorm(50)
  y1 <- rnorm(50)
  y2 <- rnorm(50)

  # Test with default parameters
  result <- ecf_test(x1, x2, y1, y2, B = 10) # Using a small B for speed

  # Check output structure
  expect_type(result, "list")
  expect_named(result, c("h", "stat", "pvalue", "Tn_b"))

  # Check output values
  expect_type(result$stat, "double")
  expect_type(result$pvalue, "double")
  expect_gte(result$pvalue, 0)
  expect_lte(result$pvalue, 1)
  expect_type(result$Tn_b, "double")
  expect_length(result$Tn_b, 10)
})
