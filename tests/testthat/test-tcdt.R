test_that("tcdt works for global tests", {
  set.seed(123)
  X1 <- matrix(rnorm(100), ncol = 2)
  X2 <- matrix(rnorm(100), ncol = 2)
  Y1 <- matrix(rnorm(50), ncol = 1)
  Y2 <- matrix(rnorm(50), ncol = 1)

  # Test global CED
  result_ced <- tcdt(X1, X2, Y1, Y2, stat = "ced", B = 10)
  expect_type(result_ced, "list")
  expect_named(result_ced, c("problem", "B", "kern_smooth", "stat", "h", "Tn", "h_b", "Tn_b", "pvalue"))
  expect_equal(result_ced$problem, "global")
  expect_gte(result_ced$pvalue, 0)
  expect_lte(result_ced$pvalue, 1)

  # Test global CMMD
  result_cmmd <- tcdt(X1, X2, Y1, Y2, stat = "cmmd", B = 10)
  expect_type(result_cmmd, "list")
  expect_named(result_cmmd, c("problem", "B", "kern_smooth", "stat", "kern_mmd", "h", "Tn", "h_b", "Tn_b", "pvalue"))
  expect_equal(result_cmmd$problem, "global")
  expect_gte(result_cmmd$pvalue, 0)
  expect_lte(result_cmmd$pvalue, 1)
})

test_that("tcdt works for local tests", {
  set.seed(123)
  X1 <- matrix(rnorm(100), ncol = 2)
  X2 <- matrix(rnorm(100), ncol = 2)
  Y1 <- matrix(rnorm(50), ncol = 1)
  Y2 <- matrix(rnorm(50), ncol = 1)
  x0 <- c(0, 0)

  # Test local CED
  result_ced <- tcdt(X1, X2, Y1, Y2, x0 = x0, stat = "ced", B = 10)
  expect_type(result_ced, "list")
  expect_named(result_ced, c("problem", "B", "kern_smooth", "stat", "h", "Tn", "h_b", "Tn_b", "pvalue"))
  expect_equal(result_ced$problem, "local")
  expect_gte(result_ced$pvalue, 0)
  expect_lte(result_ced$pvalue, 1)

  # Test local CMMD
  result_cmmd <- tcdt(X1, X2, Y1, Y2, x0 = x0, stat = "cmmd", B = 10)
  expect_type(result_cmmd, "list")
  expect_named(result_cmmd, c("problem", "B", "kern_smooth", "stat", "kern_mmd", "h", "Tn", "h_b", "Tn_b", "pvalue"))
  expect_equal(result_cmmd$problem, "local")
  expect_gte(result_cmmd$pvalue, 0)
  expect_lte(result_cmmd$pvalue, 1)
})
