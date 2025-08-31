#' @title Two-sample conditional distribution test
#' @description
#' Performs a two-sample conditional distribution test proposed by Yan, Li and Zhang (2025),
#' using either conditional energy distance (CED) or conditional maximum mean discrepancy (CMMD).
#' The function can conduct both global tests for the equality of conditional
#' distributions over the entire support, and local tests at a specific point `x0`.
#'
#' @details
#' This function implements the testing framework proposed by Yan, Li and Zhang (2025).
#' When `x0` is `NULL`, a global test is performed using the integrated measure of
#' discrepancy between the conditional distributions.
#' When `x0` is specified, a local test is performed, focusing on the discrepancy
#' at that specific point.
#'
#' Bandwidth selection for the smoothing kernel on `X` is crucial.
#' The default method `h = "cv"` uses leave-one-out cross-validation (`select_bw_cv`)
#' to find an optimal bandwidth.
#'
#' P-values are computed using a local bootstrap procedure which generates
#' bootstrap samples under the null hypothesis while preserving the marginal
#' distribution of `X`.
#' 
#' @references Yan, J., Li, Z., & Zhang, X. (2025). Distance and Kernel-Based 
#' Measures for Global and Local Two-Sample Conditional Distribution Testing
#' arXiv preprint arXiv:2210.08149.
#'
#' @param X1,X2 Matrices or data frames of the predictor variables for the two
#'   samples. Each row is an observation, and columns are variables.
#' @param Y1,Y2 Matrices or data frames of the response variables for the two
#'   samples.
#' @param x0 An optional numeric vector specifying the point for a local test.
#'   If `NULL` (the default), a global test is performed.
#' @param B The number of bootstrap replications for calculating the p-value.
#'   Default is 299. Set to 0 to compute only the test statistic.
#' @param h A character string or a function specifying the bandwidth selection
#'   method for the smoothing kernel on `X`.
#'   \itemize{
#'     \item `"cv"`: (Default) Uses leave-one-out cross-validation.
#'     \item `"undersmoothing"`: Uses a rule-of-thumb bandwidth that satisfies
#'       the undersmoothing condition required by the theory.
#'     \item A `function`: A user-provided function that takes `(n, p, nu, C)`
#'       and returns a bandwidth vector.
#'   }
#' @param adj_bw A small numeric value used to adjust the exponent of the sample
#'   size for the `"undersmoothing"` bandwidth to ensure the theoretical
#'   undersmoothing condition is met. Default is 0.1.
#' @param C A numeric multiplier for the bandwidth selected by the
#'   `"undersmoothing"` method or passed to a user-defined `h` function.
#' @param kern_smooth The type of smoothing kernel to use for the predictors `X`.
#'   One of `"gaussian"`, `"epanechnikov"`, or `"uniform"`. Default is `"gaussian"`.
#' @param nu The order of the smoothing kernel. Default is 2.
#' @param stat The type of test statistic. One of `"cmmd"` (conditional maximum
#'   mean discrepancy) or `"ced"` (conditional energy distance). Default is `"ced"`.
#' @param kern_mmd The type of RKHS kernel to use for the response `Y` when
#'   `stat = "cmmd"`. One of `"gaussian"` or `"laplacian"`. Default is `"gaussian"`.
#'
#' @return A list containing the following components:
#' \item{problem}{A string indicating the type of test performed: `"global"` or `"local"`.}
#' \item{B}{The number of bootstrap replications used.}
#' \item{kern_smooth}{The smoothing kernel function used for `X`.}
#' \item{stat}{The type of test statistic used (`"ced"` or `"cmmd"`).}
#' \item{kern_mmd}{The RKHS kernel function used for `Y` (if `stat="cmmd"`).}
#' \item{h}{A matrix of the final bandwidths used for each sample and dimension.}
#' \item{Tn}{The calculated value of the test statistic.}
#' \item{h_b}{The bandwidths used for the bootstrap procedure (if `B > 0`).}
#' \item{Tn_b}{A vector of the bootstrap test statistics (if `B > 0`).}
#' \item{pvalue}{The computed p-value of the test (if `B > 0`).}
#' @export
tcdt <- function(X1, X2, Y1, Y2, x0 = NULL, B = 299,
                 h = "cv", adj_bw = 0.1, C = 1,
                 kern_smooth = c("gaussian", "epanechnikov", "uniform"), nu = 2,
                 stat = c("ced", "cmmd"), kern_mmd = c("gaussian", "laplacian")) {
  # --- 1. Argument Matching and Data Preparation ---
  stat <- match.arg(stat)
  kern_mmd <- match.arg(kern_mmd)
  kern_smooth <- match.arg(kern_smooth)

  # Coerce inputs to matrix
  for (name_var in c("X1", "X2", "Y1", "Y2")) {
    obj <- get(name_var)
    if (is.numeric(obj) && is.null(dim(obj))) {
      assign(name_var, matrix(obj, ncol = 1))
    } else {
      assign(name_var, as.matrix(obj))
    }
  }
  if (ncol(X1) != ncol(X2)) stop("`X1` and `X2` must have the same number of columns.")
  if (ncol(Y1) != ncol(Y2)) stop("`Y1` and `Y2` must have the same number of columns.")

  n1 <- nrow(X1)
  n2 <- nrow(X2)
  p <- ncol(X1)

  n <- n1 + n2
  idx1 <- 1:n1
  idx2 <- (n1 + 1):n

  if (is.null(x0)) {
    problem <- "global"
  } else {
    problem <- "local"
    x0 <- as.numeric(x0)
    if (length(x0) != p) stop("Dimension of `x0` must match dimension of `X`.")
  }

  list_return <- list(
    problem = problem, B = B, kern_smooth = kern_smooth, stat = stat
  )
  if (stat == "cmmd") list_return$kern_mmd <- kern_mmd

  # --- 2. Construct Kernel Matrix for Y ---
  Y_pool <- rbind(Y1, Y2)
  if (stat == "ced" || (stat == "cmmd" && kern_mmd == "gaussian")) {
    dY_dist <- dist(Y_pool)
  } else if (stat == "cmmd" && kern_mmd == "laplacian") {
    dY_dist <- dist(Y_pool, method = "manhattan")
  }
  dY <- as.matrix(dY_dist)

  if (stat == "ced") {
    # Conditional Energy Distance (CED)
    kY <- -dY
  } else if (stat == "cmmd") {
    # Conditional MMD
    gamma <- median(dY_dist)
    k <- kernel_rpdc(kern_mmd)
    kY <- k(dY, gamma)
  }
  kY11 <- kY[idx1, idx1]
  kY22 <- kY[idx2, idx2]
  kY12 <- kY[idx1, idx2]

  # --- 3. Select Smoothing Bandwidths for X ---
  G <- kernel_for_smooth(kern_smooth, nu)
  nu <- attr(G, "nu")

  A1 <- pmin(apply(X1, 2, IQR) / 1.34, apply(X1, 2, sd))
  A2 <- pmin(apply(X2, 2, IQR) / 1.34, apply(X2, 2, sd))

  # Handle zero variance
  A1[A1 == 0] <- 1
  A2[A2 == 0] <- 1

  if (h == "cv") {
    kappa <- 1 / (p + 2 * nu)

    h1 <- select_bw_cv(
        X = X1, kY = kY11, x0 = x0, kappa = kappa, nu = nu, kern_smooth = kern_smooth
    )$h_vec_optim
    h2 <- select_bw_cv(
        X = X2, kY = kY22, x0 = x0, kappa = kappa, nu = nu, kern_smooth = kern_smooth
    )$h_vec_optim

  } else if (h == "undersmoothing") {
    if (problem == "global") {
      h1 <- C * n1^(-1 / (p / 2 + nu - adj_bw)) * A1
      h2 <- C * n2^(-1 / (p / 2 + nu - adj_bw)) * A2
    } else { # local
      h1 <- C * n1^(-1 / (p + nu - adj_bw)) * A1
      h2 <- C * n2^(-1 / (p + nu - adj_bw)) * A2
    }
  } else if (is.function(h)) {
    h1 <- h(n1, p, nu, C)
    h2 <- h(n2, p, nu, C)
  }
  list_return$h <- cbind(h1, h2)

  # --- 4. Calculate Test Statistic ---
  if (problem == "global") {
    sign_G1_11 <- matrix(1, nrow = n1, ncol = n1)
    sign_G1_12 <- matrix(1, nrow = n1, ncol = n2)
    sign_G2_21 <- matrix(1, nrow = n2, ncol = n1)
    sign_G2_22 <- matrix(1, nrow = n2, ncol = n2)
    log_abs_G1_11 <- matrix(0, nrow = n1, ncol = n1)
    log_abs_G1_12 <- matrix(0, nrow = n1, ncol = n2)
    log_abs_G2_21 <- matrix(0, nrow = n2, ncol = n1)
    log_abs_G2_22 <- matrix(0, nrow = n2, ncol = n2)

    for (j in 1:p) {
      G_vals_11 <- G(outer(X1[, j], X1[, j], "-"), h1[j])
      G_vals_12 <- G(outer(X1[, j], X2[, j], "-"), h1[j])
      G_vals_21 <- G(outer(X2[, j], X1[, j], "-"), h2[j])
      G_vals_22 <- G(outer(X2[, j], X2[, j], "-"), h2[j])
      
      sign_G1_11 <- sign_G1_11 * sign(G_vals_11)
      sign_G1_12 <- sign_G1_12 * sign(G_vals_12)
      sign_G2_21 <- sign_G2_21 * sign(G_vals_21)
      sign_G2_22 <- sign_G2_22 * sign(G_vals_22)
      
      abs_G11 <- abs(G_vals_11)
      abs_G12 <- abs(G_vals_12)
      abs_G21 <- abs(G_vals_21)
      abs_G22 <- abs(G_vals_22)

      log_abs_G1_11 <- log_abs_G1_11 + ifelse(abs_G11 > 0, log(abs_G11), -Inf)
      log_abs_G1_12 <- log_abs_G1_12 + ifelse(abs_G12 > 0, log(abs_G12), -Inf)
      log_abs_G2_21 <- log_abs_G2_21 + ifelse(abs_G21 > 0, log(abs_G21), -Inf)
      log_abs_G2_22 <- log_abs_G2_22 + ifelse(abs_G22 > 0, log(abs_G22), -Inf)
    }
    G1_11 <- sign_G1_11 * exp(log_abs_G1_11)
    G1_12 <- sign_G1_12 * exp(log_abs_G1_12)
    G2_21 <- sign_G2_21 * exp(log_abs_G2_21)
    G2_22 <- sign_G2_22 * exp(log_abs_G2_22)

    S1_1 <- colSums(G2_21)
    S1_2 <- colSums(G1_12)
    S2_1 <- S1_1^2 - colSums(G2_21^2)
    S2_2 <- S1_2^2 - colSums(G1_12^2)
    Tn <- teststatg(n1, n2, kY11, kY22, kY12, G1_11, G1_12, G2_22, G2_21, S1_1, S1_2, S2_1, S2_2)
  } else { # local
    sign_G1X <- rep(1, n1)
    sign_G2X <- rep(1, n2)
    log_abs_G1X <- rep(0, n1)
    log_abs_G2X <- rep(0, n2)

    for (j in 1:p) {
      G_vals_1 <- G(X1[, j] - x0[j], h1[j])
      G_vals_2 <- G(X2[, j] - x0[j], h2[j])
      
      sign_G1X <- sign_G1X * sign(G_vals_1)
      sign_G2X <- sign_G2X * sign(G_vals_2)
      
      abs_G1 <- abs(G_vals_1)
      abs_G2 <- abs(G_vals_2)
      
      log_abs_G1X <- log_abs_G1X + ifelse(abs_G1 > 0, log(abs_G1), -Inf)
      log_abs_G2X <- log_abs_G2X + ifelse(abs_G2 > 0, log(abs_G2), -Inf)
    }
    G1X <- sign_G1X * exp(log_abs_G1X)
    G2X <- sign_G2X * exp(log_abs_G2X)

    S1_1 <- sum(G1X)
    S1_2 <- sum(G2X)
    S2_1 <- S1_1^2 - sum(G1X^2)
    S2_2 <- S1_2^2 - sum(G2X^2)

    Tn <- teststatl(n1, n2, kY11, kY22, kY12, G1X, G2X, S1_1, S1_2, S2_1, S2_2)
  }
  list_return$Tn <- Tn

  # --- 5. Local Bootstrap for p-value ---
  if (B > 0) {
    G_b <- kernel_for_smooth("gaussian", nu = 2)

    # Rule of thumb
    h1_b <- n1^(-1 / (p + 4)) * A1
    h2_b <- n2^(-1 / (p + 4)) * A2
    list_return$h_b <- cbind(h1_b, h2_b)

    log_G1_11_b <- matrix(0, nrow = n1, ncol = n1)
    log_G1_12_b <- matrix(0, nrow = n1, ncol = n2)
    log_G2_21_b <- matrix(0, nrow = n2, ncol = n1)
    log_G2_22_b <- matrix(0, nrow = n2, ncol = n2)
    for (j in 1:p) {
      G_vals_11_b <- G_b(outer(X1[, j], X1[, j], "-"), h1_b[j])
      G_vals_12_b <- G_b(outer(X1[, j], X2[, j], "-"), h1_b[j])
      G_vals_21_b <- G_b(outer(X2[, j], X1[, j], "-"), h2_b[j])
      G_vals_22_b <- G_b(outer(X2[, j], X2[, j], "-"), h2_b[j])

      log_G1_11_b <- log_G1_11_b + ifelse(G_vals_11_b > 0, log(G_vals_11_b), -Inf)
      log_G1_12_b <- log_G1_12_b + ifelse(G_vals_12_b > 0, log(G_vals_12_b), -Inf)
      log_G2_21_b <- log_G2_21_b + ifelse(G_vals_21_b > 0, log(G_vals_21_b), -Inf)
      log_G2_22_b <- log_G2_22_b + ifelse(G_vals_22_b > 0, log(G_vals_22_b), -Inf)
    }
    G1_11_b <- exp(log_G1_11_b)
    G1_12_b <- exp(log_G1_12_b)
    G2_21_b <- exp(log_G2_21_b)
    G2_22_b <- exp(log_G2_22_b)

    Gp <- rbind(cbind(G1_11_b, G1_12_b), cbind(G2_21_b, G2_22_b))
    w <- sweep(Gp, 2, colSums(Gp), "/", check.margin = FALSE)

    Tn_b <- numeric(B)
    for (b in 1:B) {
      idx_b <- numeric(n)

      # Local bootstrap
      for (i in 1:n) {
        idx_b[i] <- sample(1:n, 1, prob = w[, i])
      }
      idx1_b <- idx_b[idx1]
      idx2_b <- idx_b[idx2]

      kYb11 <- kY[idx1_b, idx1_b]
      kYb22 <- kY[idx2_b, idx2_b]
      kYb12 <- kY[idx1_b, idx2_b]

      if (problem == "global") {
        Tn_b[b] <- teststatg(n1, n2, kYb11, kYb22, kYb12, G1_11, G1_12, G2_22, G2_21, S1_1, S1_2, S2_1, S2_2)
      } else if (problem == "local") {
        Tn_b[b] <- teststatl(n1, n2, kYb11, kYb22, kYb12, G1X, G2X, S1_1, S1_2, S2_1, S2_2)
      }
    }
    pvalue <- (1 + sum(Tn_b >= Tn)) / (1 + B)

    list_return$Tn_b <- Tn_b
    list_return$pvalue <- pvalue
  }

  list_return
}

select_bw_cv <- function(X, kY,
                         x0 = NULL,
                         kappa = NULL, nu = 2,
                         kern_smooth = c("gaussian", "epanechnikov", "uniform"),
                         C_init = NULL, lower_C = 1e-4,
                         optim_control = list(maxit = 100),
                         verbose = FALSE) {
  # Initial Checks and Setup
  kern_smooth <- match.arg(kern_smooth)
  n <- nrow(X)
  p <- ncol(X)
  if (nrow(kY) != n || ncol(kY) != n) stop("kY matrix dimensions must match X nrow.")
  if (n <= 1) stop("Need at least 2 data points for LOO-CV.")

  if (is.null(kappa)) {
    kappa <- 1 / (p + 2 * nu)
  }
  kern_map <- list(gaussian = 1, epanechnikov = 2, uniform = 3)
  kern_type_int <- kern_map[[kern_smooth]]

  local_weights <- rep(1.0, n)

  # Optimization
  if (is.null(C_init)) C_init <- rep(1.0, p)
  if (verbose) cat("Starting optimization for sample...\n")

  objective_for_optim <- function(C) {
    cv_objective(
      C = C, X = X, kY = kY, local_weights = local_weights,
      kappa = kappa, nu = nu, kern_type = kern_type_int, lower_C = lower_C
    )
  }

  opt_result <- optim(
    par = C_init, fn = objective_for_optim, method = "L-BFGS-B",
    lower = rep(lower_C, p), control = optim_control
  )

  if (opt_result$convergence != 0) {
    warning("Optimization did not converge. Result may be unreliable. Code: ", opt_result$convergence, " Message: ", opt_result$message)
  } else {
    if (verbose) cat("Optimization converged.\n")
  }

  C_optim <- opt_result$par
  h_vec_optim <- C_optim * (n^(-kappa))

  list(
    C_optim = C_optim, h_vec_optim = h_vec_optim,
    optim_result = opt_result, cv_value = opt_result$value
  )
}
