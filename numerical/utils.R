pkg_all <- c("mixtools", "mvtnorm", "sn", "ggplot2", "ggthemes")
for (pkg in pkg_all) {
  suppressWarnings(suppressMessages(library(pkg, character.only = TRUE)))
}

#' Data generating process with regression framework
#' @param n1,n2 Sample sizes of the first and the second sample, respectively.
#' @param p,q Number of dimensions of the covariate and the response, respectively.
#' @param dist_x1,param1_x1,param2_x1 Specify the marginal distribution of
#' covariate of the first sample.
#' See [dgp_margin()] for details.
#' #' @param dist_x2,param1_x2,param2_x2 Specify the marginal distribution of
#' covariate of the second sample.
#' See [dgp_margin()] for details.
#' @param form_m1,form_m2,signal_m Specify the regression function for two
#' samples, and the (potential) signal on the second sample.
#' - `signal_m` is a parameter determining the signal level of the discrepancy
#' between the conditional mean functions for the two samples.
#' - `form_m1` and `form_m2` should be function of x.
#' - `form_m2` can provide a potential argument `signal`.
#' When `form_m2` does not have the argument `signal`, the `signal_m` must be
#' specified as `NULL`.
#' @param form_v1,form_v2,signal_v Specify the conditional variance for
#' tow samples, and the (potential) on the second sample.
dgp_reg <- function(n1, n2, p, q,
                    dist_x1, param1_x1 = NULL, param2_x1 = NULL,
                    dist_x2, param1_x2 = NULL, param2_x2 = NULL,
                    dist_eps1, param_eps1 = NULL,
                    dist_eps2, param_eps2 = NULL,
                    form_m1, form_m2, signal_m = NULL,
                    form_v1, form_v2, signal_v = NULL) {

  # Generate x
  X1 <- component(n1, p, dist_x1, param1_x1, param2_x1)
  X2 <- component(n2, p, dist_x2, param1_x2, param2_x2)

  # Generate conditional mean
  m1 <- component_cond(X1, form_m1)
  m2 <- component_cond(X2, form_m2, signal_m)

  # Generate error term
  epsilon1 <- component(n1, q, dist_eps1, NULL, param_eps1)
  epsilon2 <- component(n2, q, dist_eps2, NULL, param_eps2)

  # Generate conditional variance
  v1 <- component_cond(X1, form_v1)
  v2 <- component_cond(X2, form_v2, signal_v)

  # Obtain y (also work for multivariate case)
  Y1 <- m1 + sqrt(v1) * epsilon1
  Y2 <- m2 + sqrt(v2) * epsilon2

  list(X1 = X1, X2 = X2, Y1 = Y1, Y2 = Y2)
}

#' Wrapper function to generate random sample from certain distributions
#' @param n Sample size
#' @param d Number of dimensions
#' @param dist,param1,param2 Specify the distribution.
#' - `dist` can be "normal", "uniform" or "beta"
#' - For `dist = "normal"`, `param1` is the mean, and `param2` is the variance
#' (for `d = 1`) or covariance matrix (for `d > 1`)
component <- function(n, d,
                      dist = c("normal", "t", "uniform", "beta", "gaussian_mixture"),
                      param1 = NULL, param2 = NULL) {
  dist <- match.arg(dist)

  n <- as.integer(n)
  if (!(n >= 1)) {
    stop()
  }

  d <- as.integer(d)
  if (!(d >= 1)) {
    stop()
  }

  if (dist == "normal") {
    # mean vector (mean for d = 1)
    mu <- param1
    if (is.null(mu)) {
      mu <- rep(0, d)
    } else if (!is.null(mu)) {
      mu <- drop(as.matrix(mu))
      if (!(inherits(mu, "numeric") && length(mu) == d)) {
        stop()
      }
    }

    # covariance matrix (or sd for d = 1)
    sigma <- param2
    if (is.null(sigma)) {
      if (d == 1) {
        sigma <- 1
      } else if (d >= 2) {
        sigma <- diag(1, d)
      }
    } else if (!is.null(sigma)) {
      sigma <- drop(as.matrix(sigma))
      if (d == 1 && !(inherits(sigma, "numeric") && length(sigma) == 1)) {
        stop()
      } else if (d >= 2 && !(inherits(sigma, "matrix") && all(dim(sigma) == c(d, d)))) {
        stop()
      }
    }

    component <- matrix(rnorm(n * d), nrow = n)
    if (d == 1) {
      component <- sigma * component + mu
    } else if (d >= 2) {
      component <- sweep(component %*% chol(sigma), 2, mu, "+")
    }
  } else if (dist == "t") {
    df <- param1
    if (is.null(df)) {
      df <- 5
    }

    sd <- param2

    component <- matrix(rt(n * d, df = df), nrow = n)

    # If a desired standard deviation is specified, scale the sample
    if (!is.null(sd)) {
      # The standard deviation of a standard t-distribution is sqrt(df / (df - 2)) for df > 2
      if (df <= 2) {
        stop("Standard deviation is undefined for degrees of freedom less than or equal to 2.")
      }
      # Scale the sample to have a standard deviation of 1
      component <- component / sqrt(df / (df - 2))
      # Adjust the scaled sample to the desired standard deviation
      component <- component * sd
    }


  } else if (dist == "uniform") {
    param_min <- param1
    if (is.null(param_min)){
      param_min <- 0
    }

    param_max <- param2
    if (is.null(param_max)) {
      param_max <- 1
    }
    component <- drop(replicate(d, runif(n, min = param_min, max = param_max)))
  } else if (dist == "beta") {
    shape1 <- param1
    shape2 <- param2

    component <- rbeta(n, shape1 = shape1, shape2 = shape2)
  } else if (dist == "gaussian_mixture") {
    mu2 <- param1
    sigma2 <- param2

    lambda <- c(0.5, 0.5)
    mu <- c(0, mu2)
    sigma <- c(1, sigma2)

    component <- rnormmix(n, lambda, mu, sigma)
  }

  component
}

component_cond <- function(X, form, signal = NULL) {
  if (is.null(signal)) {
    cond_part <- form(X)
  } else {
    cond_part <- form(X, signal)
  }

  cond_part
}

extract_result <- function(obj_extract, params_varying, params_fixed) {
  params_varying <- lapply(params_varying, as.character)
  params_fixed <- lapply(params_fixed, as.character)

  for (name_param in names(params_fixed)) {
    assign(name_param, params_fixed[[name_param]])
    params_fixed[[name_param]] <- NULL
  }
  
  dir_result <- file.path(
    "output", "results", "simulations", problem, setting_shift, setting_dgp
  )
  
  if (
    problem == "global" && 
    setting_shift == "covariate_shift" && 
    setting_dgp == "different_mean"
  ) {
    str_glue <- paste(
      "m_{setting_m}",
      "{signal_m}",
      "n_{n1}_{n2}",
      "p_{p}",
      "dist_x_{dist_x}",
      "shift_x_{shift_x}",
      "dist_eps_{dist_eps}",
      "num_simu_{num_simu}.RData",
      sep = "__"
    )
  } else if (
    problem == "global" && 
    setting_shift == "covariate_shift" && 
    setting_dgp == "different_var"
  ) {
    str_glue <- paste(
      "v_{setting_v}",
      "{signal_v}",
      "m_{setting_m}",
      "n_{n1}_{n2}",
      "p_{p}",
      "dist_x_{dist_x}",
      "shift_x_{shift_x}",
      "dist_eps_{dist_eps}",
      "num_simu_{num_simu}.RData",
      sep = "__"
    )
  } else if (
    problem == "global" && 
    setting_shift == "label_shift" && 
    setting_dgp == "different_mean"
  ) {
    str_glue <- paste(
      "m_{setting_m}",
      "{signal_m}",
      "n_{n1}_{n2}",
      "p_{p}",
      "shift_y_{shift_y}",
      "num_simu_{num_simu}.RData",
      sep = "__"
    )
  }

  # Force "method" be the 2nd dimension (i.e., the column)
  dim_array <- append(params_varying, list(method = methods), after = 1)

  result_array <- array(
    dim = lapply(dim_array, length),
    dimnames = dim_array
  )

  params_grid <- expand.grid(params_varying, stringsAsFactors = FALSE)
  for (ind_param in seq_len(nrow(params_grid))) {
    param <- params_grid[ind_param, ]
    for (name_param in names(param)) {
      assign(name_param, param[, name_param])
      if (name_param == "n") {
        n1 <- n2 <- n
      }
    }
    
    path_result <- file.path(dir_result, glue(str_glue))

    obj_result <- env()
    load(path_result, obj_result)

    if (obj_extract == "power_array") {
      result <- colMeans(obj_result$reject)
    } else if (obj_extract == "time_array") {
      result <- colMeans(obj_result$timecost)
    }
    result_array <- do.call( # Equivalent to `result_array[methods, ...] <- result[methods]`
      "[<-",
      c(
        list(object = result_array),
        append(param, list(method = methods), after = 1), # adaptive to the param
        list(value = result[methods])
      )
    )
  }

  result_array
}

my_theme_base <- function() {
  theme_foundation() +
    theme(
      line = element_line(colour = "black", lineend = "round", linetype = "solid"),
      rect = element_rect(fill = "white", colour = "black", linetype = "solid"),
      text = element_text(colour = "black", face = "plain", family = "", size = 16, vjust = 0.5, hjust = 0.5, lineheight = 1),
      panel.grid = element_line(colour = "gray95"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.background = element_blank(),
      strip.background = element_rect(colour = NA, fill = NA),
      legend.key = element_rect(colour = NA),
      title = element_text(size = rel(1)),
      plot.title = element_text(hjust = 0.5),
      strip.text = element_text(),
      axis.ticks.length = unit(0.5, "lines"),
      axis.text.x=element_text(angle = 60, hjust = 1),
      legend.position = "bottom"
    )
}
