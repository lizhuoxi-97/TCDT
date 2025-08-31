#' @title Test for the equality of two regression curves via empirical characteristic functions
#' @description Implements the test for the equality of two regression curves
#' proposed by Pardo-Fernandez et al. (2015).
#' The test is based on the difference between the empirical characteristic
#' functions of the residuals in the two groups.
#' @param x1 A numeric vector of the covariate for the first sample.
#' @param x2 A numeric vector of the covariate for the second sample.
#' @param y1 A numeric vector of the outcome for the first sample.
#' @param y2 A numeric vector of the outcome for the second sample.
#' @param B The number of bootstrap replications to compute the p-value. Default is 299.
#' @param h The bandwidth parameter for the kernel smoothing.
#' If NULL (the default), it is selected automatically.
#' @param C A constant to multiply the automatically selected bandwidth. Default is 1.
#' @param rescale A logical indicating whether to rescale the covariate `x` to
#' the interval `[0, 1]`.
#' Default is FALSE.
#' @param kern_smooth The kernel function to be used for smoothing. One of
#' "epanechnikov", "gaussian", or "uniform".
#' Default is "epanechnikov".
#' @param reg_method The regression method to be used.
#' One of "llinear" (local linear) or "lconst" (local constant, Nadaraya-Watson).
#' Default is "llinear".
#' @return A list with the following components:
#' \item{stat}{The value of the test statistic.}
#' \item{pvalue}{The bootstrap p-value.}
#' \item{h}{The bandwidth used.}
#' \item{Tn_b}{The bootstrap values of the test statistic.}
#' @references Pardo‐Fernández, J. C., Jiménez‐Gamero, M. D., & Ghouch, A. E. (2015).
#' A Non‐parametric anova‐type Test for Regression Curves Based on Characteristic Functions.
#' Scandinavian Journal of Statistics, 42(1), 197-213.
#' @export
ecf_test <- function(x1, x2, y1, y2,
                     B = 299, h = NULL, C = NULL, rescale = FALSE,
                     kern_smooth = c("epanechnikov", "gaussian", "uniform"),
                     reg_method = c("llinear", "lconst")) {
  list_return <- list()

  kern_smooth <- match.arg(kern_smooth)
  reg_method <- match.arg(reg_method)

  # Coerce and check
  for (name_var in c("x1", "x2", "y1", "y2")) {
    obj <- eval(rlang::parse_expr(name_var))
    if (!is.data.frame(obj) && !is.matrix(obj) && !is.vector(obj)) {
      stop(paste(name_var, "must be a vector, matrix or data frame."))
    }

    if (is.data.frame(obj) || is.matrix(obj)) {
      if (ncol(obj) != 1) {
        stop(paste(name_var, "must have exactly 1 column."))
      } else if (ncol(obj) == 1) {
        assign(name_var, as.vector(obj[, 1]))
      }
    }
  }

  xpool <- c(x1, x2)
  ypool <- c(y1, y2)

  n1 <- length(y1)
  n2 <- length(y2)
  n <- length(ypool)

  if (rescale) {
    xpool <- scales::rescale(xpool)
    x1 <- xpool[1:n1]
    x2 <- xpool[(n1 + 1):n]
  }

  K <- kernel_for_smooth(kern_smooth, 2)

  if (is.null(C)) {
    C <- 1
  }

  if (is.null(h)) {
    h <- C * n^(- 3 / 8)
    # h <- 1.06 * n^(-1 / (1 / 2 + 2)) * sd(xpool) * 0.9
  }
  list_return$h <- h

  # Local constant weight
  Kpool <- K(outer(xpool, xpool, "-"), h)
  Kpool_1 <-  Kpool[, 1:n1]
  Kpool_2 <-  Kpool[, -(1:n1)]

  fpool <- rowSums(Kpool) / n
  f1 <- rowSums(Kpool_1) / n1
  f2 <- rowSums(Kpool_2) / n2

  wpool_nw <- sweep(Kpool, 1, fpool * n, "/")
  w1_nw <- sweep(Kpool_1, 1, f1 * n1, "/")
  w2_nw <- sweep(Kpool_2, 1, f2 * n2, "/")

  # Fall back to nearest neighbor
  for (idx in which(fpool == 0)) {
    wpool_nw[idx, ] <- as.numeric(seq_along(Kpool[idx, ]) == which.max(Kpool[idx, ]))
  }
  for (idx in which(f1 == 0)) {
    w1_nw[idx, ] <- as.numeric(seq_along(Kpool_1[idx, ]) == which.max(Kpool_1[idx, ]))
  }
  for (idx in which(f2 == 0)) {
    w2_nw[idx, ] <- as.numeric(seq_along(Kpool_2[idx, ]) == which.max(Kpool_2[idx, ]))
  }

  # Local linear
  if (reg_method == "llinear") {
    wpool <- matrix(nrow = n, ncol = n)
    w1 <- matrix(nrow = n, ncol = n1)
    w2 <- matrix(nrow = n, ncol = n2)

    for (i in 1:n) {
      x0 <- xpool[i]

      Kpool_0 <- K(xpool - x0, h)
      K1_0 <- Kpool_0[1:n1]
      K2_0 <- Kpool_0[-(1:n1)]

      xpoolaug <- cbind(1, xpool - x0)
      x1aug <- xpoolaug[1:n1, ]
      x2aug <- xpoolaug[-(1:n1), ]

      wpool[i, ] <- tryCatch(
        (solve(t(xpoolaug) %*% diag(Kpool_0) %*% xpoolaug) %*% t(xpoolaug) %*% diag(Kpool_0))[1, ],
        error = function(e) as.numeric(seq_along(Kpool_0) == which.max(Kpool_0)) # Fall back to nearest neighbor
      )
      w1[i, ] <- tryCatch(
        (solve(t(x1aug) %*% diag(K1_0) %*% x1aug) %*% t(x1aug) %*% diag(K1_0))[1, ],
        error = function(e) as.numeric(seq_along(K1_0) == which.max(K1_0)) # Fall back to nearest neighbor
      )
      w2[i, ] <- tryCatch(
        (solve(t(x2aug) %*% diag(K2_0) %*% x2aug) %*% t(x2aug) %*% diag(K2_0))[1, ],
        error = function(e) as.numeric(seq_along(K2_0) == which.max(K2_0)) # Fall back to nearest neighbor
      )
    }
  } else if (reg_method == "lconst") {
    wpool <- wpool_nw
    w1 <- w1_nw
    w2 <- w2_nw
  }

  m1 <- drop(w1 %*% y1)
  m2 <- drop(w2 %*% y2)
  mpool <- (f1  * m1 * n1 + f2 * m2 * n2) / fpool / n

  ypool_squared <- ypool^2

  # Estimate of conditional variance
  sigma2_1 <- pmax(
    drop(w1_nw %*% ypool_squared[1:n1]) - drop(w1_nw %*% ypool[1:n1])^2,
    1e-4
  )
  sigma2_2 <- pmax(
    drop(w2_nw %*% ypool_squared[-(1:n1)]) - drop(w2_nw %*% ypool[-(1:n1)])^2,
    1e-4
  )

  # Standardized residuals
  eps_1 <- (y1 - m1[1:n1]) / sqrt(sigma2_1[1:n1])
  eps_2 <- (y2 - m2[-(1:n1)]) / sqrt(sigma2_2[-(1:n1)])
  eps0_1 <- (y1 - mpool[1:n1]) / sqrt(sigma2_1[1:n1])
  eps0_2 <- (y2 - mpool[-(1:n1)]) / sqrt(sigma2_2[-(1:n1)])

  # Test statistic
  Tn <- ecf_reg_test_stat(eps_1, eps0_1, eps_2, eps0_2)
  list_return$stat <- Tn

  # Bootstrap
  if (B > 0) {
    a1 <- 2 * n1 ^ (- 3 / 10)
    a2 <- 2 * n2 ^ (- 3 / 10)
    Tn_b <- numeric(B)
    for (b in 1:B) {
      eps_1_star <- ifelse(runif(n1) <= a1, rnorm(n1), sample(eps_1, n1, replace = TRUE))
      eps_2_star <- ifelse(runif(n2) <= a2, rnorm(n2), sample(eps_2, n2, replace = TRUE))

      y1_b <- mpool[1:n1] + sqrt(sigma2_1[1:n1]) * eps_1_star
      y2_b <- mpool[-(1:n1)] + sqrt(sigma2_2[-(1:n1)]) * eps_2_star

      ypool_b <- c(y1_b, y2_b)

      m1_b <- drop(w1 %*% y1_b)
      m2_b <- drop(w2 %*% y2_b)
      mpool_b <- (f1  * m1_b * n1 + f2 * m2_b * n2) / fpool / n

      ypool_squared_b <- ypool_b^2

      # Estimate of conditional variance
      sigma2_1_b <- pmax(
        drop(w1_nw %*% ypool_squared_b[1:n1]) - drop(w1_nw %*% ypool_b[1:n1])^2,
        1e-4
      )
      sigma2_2_b <- pmax(
        drop(w2_nw %*% ypool_squared_b[-(1:n1)]) - drop(w2_nw %*% ypool_b[-(1:n1)])^2,
        1e-4
      )

      # Standardized residuals
      eps_1_b <- (y1_b - m1_b[1:n1]) / sqrt(sigma2_1_b[1:n1])
      eps_2_b <- (y2_b - m2_b[-(1:n1)]) / sqrt(sigma2_2_b[-(1:n1)])
      eps0_1_b <- (y1_b - mpool_b[1:n1]) / sqrt(sigma2_1_b[1:n1])
      eps0_2_b <- (y2_b - mpool_b[-(1:n1)]) / sqrt(sigma2_2_b[-(1:n1)])

      # Test statistic
      Tn_b[b] <- ecf_reg_test_stat(eps_1_b, eps0_1_b, eps_2_b, eps0_2_b)
    }

    pvalue <- (sum(Tn_b >= Tn) + 1) / (B + 1)
    list_return$pvalue <- pvalue
    list_return$Tn_b <- Tn_b
  }

  list_return
}
