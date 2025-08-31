test_that("conf_test works with basic inputs", {
  set.seed(123)
  X1 <- matrix(rnorm(100), ncol = 2)
  X2 <- matrix(rnorm(100), ncol = 2)
  Y1 <- matrix(rnorm(50), ncol = 1)
  Y2 <- matrix(rnorm(50), ncol = 1)

  # Test with default parameters
  result <- conf_test(X1, X2, Y1, Y2)

  # Check output structure
  expect_type(result, "list")
  expect_named(result, c("U", "var", "z", "pvalue"))

  # Check output values
  expect_type(result$U, "double")
  expect_type(result$var, "double")
  expect_type(result$z, "double")
  expect_type(result$pvalue, "double")
  expect_gte(result$pvalue, 0)
  expect_lte(result$pvalue, 1)
})
