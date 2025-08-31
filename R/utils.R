# Scaled kernel for local smoothing
kernel_for_smooth <- function(type = c("gaussian", "epanechnikov", "uniform"),
                              nu = 2, param = NULL) {
  type <- match.arg(type)

  if (type == "gaussian") {
    if (nu == 2) {
      kernel <- function(x, h) {
        exp(-(x / h)^2 / 2) / sqrt(2 * pi) / h
      }
    } else if (nu == 4) {
      kernel <- function(x, h) {
        u_sq <- (x / h)^2
        (3 - u_sq) / 2 * exp(-u_sq / 2) / sqrt(2 * pi) / h
      }
    } else if (nu == 6) {
      kernel <- function(x, h) {
        u_sq <- (x / h)^2
        (15 - 10 * u_sq + u_sq^2) / 8 * exp(-u_sq / 2) / sqrt(2 * pi) / h
      }
    } else {
      stop("For gaussian smoothing kernel, only support `nu` = 2, 4, or 6 currently.")
    }
  } else if (type == "epanechnikov") {
    if (nu == 2) {
      kernel <- function(x, h) {
        ifelse(
          abs(x / h) <= 1,
          (1 - (x / h)^2) * 3 / 4 / h,
          0
        )
      }
    } else {
      stop("For epanechnikov smoothing kernel, only support `nu` = 2 currently.")
    }
  } else if (type == "uniform") {
    if (nu == 2) {
      if (is.null(param)) {
        boundary <- 1
      } else {
        boundary <- param
      }
      kernel <- function(x, h) {
        ifelse(
          (x / h) >= -boundary & (x / h) <= boundary,
          1 / (2 * boundary * h),
          0
        )
      }
    } else {
      stop("For uniform smoothing kernel, only support `nu` = 2 currently.")
    }
  }

  attr(kernel, "nu") <- nu
  kernel
}

# Reproducing kernel
kernel_rpdc <- function(type = c("gaussian", "laplacian")) {
  type <- match.arg(type)

  if (type == "gaussian") {
    kernel <- function(x, gamma) {
      exp(-x^2 / gamma^2 / 2)
    }
  } else if (type == "laplacian") {
    kernel <- function(x, gamma) {
      exp(-x / gamma)
    }
  }

  kernel
}

