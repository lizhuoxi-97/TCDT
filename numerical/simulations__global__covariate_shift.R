pkg_all <- c("purrr", "glue", "digest", "TCDT")
for (pkg in pkg_all) {
  suppressWarnings(suppressMessages(library(pkg, character.only = TRUE)))
}

source("numerical/utils.R")

# Overall settings ----
B <- 299 # number of bootstrap replications

# Parameterized settings ----
args <- commandArgs(TRUE)
if (length(args) > 0) {
  type_args <- "shell"
} else if (length(args) == 0) {
  type_args <- "Rscript"
}

if (type_args == "shell") {
  setting_dgp <- args[1]
  n1 <- as.numeric(args[2])
  n2 <- as.numeric(args[3])
  p <- as.numeric(args[4])
} else if (type_args == "Rscript") {
  setting_dgp <- "different_mean"
  n1 <- 50
  n2 <- 50
  # p <- 1
  p <- 4
  dist_x <- "norm"
  shift_x <- "true"
  dist_eps <- "norm"
  # setting_m <- "x2+x(x+2)(x-2)"
  setting_m <- "sindex_quad"
  signal <- 0
  num_simu <- 3
}

settings_dgp_reg <- c("different_mean", "different_var")
if (setting_dgp %in% settings_dgp_reg) {
  dgp_framework <- "reg"
} else {
  dgp_framework <- "custom"
}

if (dgp_framework == "reg") {
  # Regression DGP framework ----
  if (type_args == "shell") {
    dist_x <- args[5]
    shift_x <- args[6]
  }

  ## Marginal distribution of X ----
  if (dist_x == "norm") {
    ### normal ----
    dist_x1 <- dist_x2 <- "normal"

    if (p == 1) {
      #### p = 1 ----
      param1_x1 <- 0 # mean of normal
      param2_x1 <- 1 # sd of normal

      if (shift_x == "false") {
        param1_x2 <- param1_x1
        param2_x2 <- param2_x1
      } else if (shift_x == "true") {
        param1_x2 <- param1_x1 + 1
        param2_x2 <- param2_x1
      }
    } else if (p == 4) {
      #### p = 4 ----
      param1_x1 <- rep(0, p) # mean vector of normal
      param2_x1 <- diag(1, p) # covariance matrix of normal

      if (shift_x == "false") {
        param1_x2 <- param1_x1
        param2_x2 <- param2_x1
      } else if (shift_x == "true") {
        param1_x2 <- c(1, 1, -1, 0)
        param2_x2 <- param2_x1
      }
    }
  } else if (dist_x == "gaumix") {
    ### gaussian mixture ----
    dist_x1 <- dist_x2 <- "gaussian_mixture"

    # 0.5 N(0, 1) + 0.5 N(param1, param2^2)

    param1_x1 <- 1 # mean of mix2 for X1
    param2_x1 <- 1 # sd of mix2 for X1

    if (shift_x == "false") {
      param1_x2 <- param1_x1
      param2_x2 <- param2_x1
    } else if (shift_x == "true") {
      param1_x2 <- 0 # mean of mix2 for X2
      param2_x2 <- sqrt(1.5) # sd of mix2 for X2
    }
  } else if (dist_x == "normt") {
    ### gaussian and t ----
    dist_x1 <- "norm"
    dist_x2 <- "t"

    param1_x1 <- 0 # mean of normal
    param2_x1 <- 1 # sd of normal

    if (shift_x == "false") {
      stop("dist_x = 'normt' only allows shift_x = 'true'")
    } else if (shift_x == "true") {
      param1_x2 <- 5 # degree of freedom of t
      param2_x2 <- NULL # sd of t (take default, that is t itself)
    }
  }

  ## Conditional distribution settings, and methods to be compared ----
  if (setting_dgp == "different_mean") {
    ### global, covariate shift, different mean ----
    # params:
    # 1. setting_dgp
    # 2. n1
    # 3. n2
    # 4. p
    # 5. dist_x
    # 6. shift_x
    # 7. dist_eps
    # 8. setting_m
    # 9. signal
    # 10. num_simu

    if (type_args == "shell") {
      dist_eps <- args[7]
      setting_m <- args[8]
      signal <- as.numeric(args[9])
    }

    name_simulation <- glue(
      "m_{setting_m}",
      "{signal}",
      "n_{n1}_{n2}",
      "p_{p}",
      "dist_x_{dist_x}",
      "shift_x_{shift_x}",
      "dist_eps_{dist_eps}",
      .sep = "__"
    )

    q <- 1

    if (p == 1) {
      #### p = 1 ----
      nu <- 2

      methods <- c("CED", "CMMD", "CONF", "ECF")

      if (dist_eps == "norm") {
        dist_eps1 <- dist_eps2 <- "normal"
        param_eps1 <- param_eps2 <- NULL
      } else if (dist_eps == "t") {
        dist_eps1 <- dist_eps2 <- "t"
        param_eps1 <- param_eps2 <- NULL
      }

      signal_m <- signal
      if (setting_m == "x2+x(x+2)(x-2)") {
        ##### x2+x(x+2)(x-2) ----
        form_m1 <- function(x) 1 + x^2
        form_m2 <- function(x, signal) 1 + x^2 + signal * x*(x+2)*(x-2)/5

        form_v1 <- form_v2 <- function(x) 0.5
        signal_v <- NULL
      } else if (setting_m == "exp+sin") {
        ##### exp+sin ----
        form_m1 <- function(x) 1 + exp(x)
        form_m2 <- function(x, signal) 1 + exp(x) + signal * sin(2 * pi * x)

        form_v1 <- form_v2 <- function(x) 0.5
        signal_v <- NULL
      } else {
        stop("`setting_m` not matched.")
      }
    } else if (p == 4) {
      #### p = 4 ----
      nu <- 4

      methods <- c("CED", "CMMD", "CONF")

      if (dist_eps == "norm") {
        dist_eps1 <- dist_eps2 <- "normal"
        param_eps1 <- param_eps2 <- 1
      } else if (dist_eps == "t") {
        dist_eps1 <- dist_eps2 <- "t"
        param_eps1 <- param_eps2 <- NULL
      }

      signal_m <- signal
      if (setting_m == "sindex_quad") {
        ##### single index cubic ----
        form_m1 <- function(X) {
          theta <- X %*% c(1, 1, 1, 1) + 1
          drop(theta^2)
        }
        form_m2 <- function(X, signal) {
          theta <- X %*% c(1, 1, 1, 1) + 1
          drop(theta^2 + signal * theta)
        }
        form_v1 <- form_v2 <- function(x) 0.5
        signal_v <- NULL
      } else {
        stop("`setting_m` not matched.")
      }
    }
  } else if (setting_dgp == "different_var") {
    ### global, covariate shift, different variance ----
    # params:
    # 1. setting_dgp
    # 2. n1
    # 3. n2
    # 4. p
    # 5. dist_x
    # 6. shift_x
    # 7. dist_eps
    # 8. setting_m
    # 9. setting_v
    # 10. signal
    # 11. num_simu

    if (type_args == "shell") {
      dist_eps <- args[7]
      setting_m <- args[8]
      setting_v <- args[9]
      signal <- as.numeric(args[10])
    }

    name_simulation <- glue(
      "v_{setting_v}",
      "{signal}",
      "m_{setting_m}",
      "n_{n1}_{n2}",
      "p_{p}",
      "dist_x_{dist_x}",
      "shift_x_{shift_x}",
      "dist_eps_{dist_eps}",
      .sep = "__"
    )

    q <- 1

    nu <- 2

    methods <- c("CED", "CMMD", "CONF", "ECF")

    if (dist_eps == "norm") {
      dist_eps1 <- dist_eps2 <- "normal"
      param_eps1 <- param_eps2 <- 1
    } else if (dist_eps == "t") {
      dist_eps1 <- dist_eps2 <- "t"
      param_eps1 <- param_eps2 <- NULL
    }

    if (setting_m == "sin") {
      form_m1 <- form_m2 <- function(x) sin(2 * pi * x)
      signal_m <- NULL
    } else if (setting_m == "quad") {
      form_m1 <- form_m2 <- function(x) 1 + x^2
      signal_m <- NULL
    } else {
      stop("`setting_m` not matched.")
    }

    signal_v <- signal
    if (setting_v == "homo") {
      #### homo ----
      form_v1 <- function(x) 0.25
      form_v2 <- function(x, signal) 0.25 + signal
    } else if (setting_v == "hetero") {
      #### hetero ----
      form_v1 <- function(x) 0.25^2 * exp(2 * x)
      form_v2 <- function(x, signal) (0.25 * (1 + 2 * signal))^2 * exp(2 * x)
    } else {
      stop("`setting_v` not matched.")
    }
  }

  dgp <- partial(
    dgp_reg,
    n1 = n1, n2 = n2, p = p, q = q,
    dist_x1 = dist_x1, param1_x1 = param1_x1, param2_x1 = param2_x1,
    dist_x2 = dist_x2, param1_x2 = param1_x2, param2_x2 = param2_x2,
    dist_eps1 = dist_eps1, param_eps1 = param_eps1,
    dist_eps2 = dist_eps2, param_eps2 = param_eps2,
    form_m1 = form_m1, form_m2 = form_m2, signal_m = signal_m,
    form_v1 = form_v1, form_v2 = form_v2, signal_v = signal_v
  )

} else if (dgp_framework == "custom") {
  cat("> Using custom DGP framework.")

  # Specify your custom dgp function as `dgp()` here:

}

# Path and filename for storing results ----
dir <- file.path(
  "output", "results", "simulations", "global", "covariate_shift", setting_dgp
)
if (type_args == "shell") {
  num_simu <- as.integer(tail(args, 1))
}
name_simulation <- glue(
  name_simulation,
  "num_simu_{num_simu}",
  .sep = "__"
)
simu_all <- 1:num_simu

if (!dir.exists(dir)) {
  dir.create(dir, recursive = TRUE)
}

# Simulation ----
reject <- timecost <- matrix(
  nrow = length(simu_all), ncol = length(methods),
  dimnames = list(simu_all, methods)
)

for (ind_simu in seq_along(simu_all)) {
  simu <- simu_all[ind_simu]

  # Fixed seed for dgp
  set.seed(strtoi(substr(digest(simu, algo = "xxhash32"), 1, 7), 16))

  obj_dgp <- dgp()
  X1 <- obj_dgp$X1
  X2 <- obj_dgp$X2
  Y1 <- obj_dgp$Y1
  Y2 <- obj_dgp$Y2

  for (method in methods) {
    # Fixed seed for each method
    set.seed(strtoi(substr(digest(simu, algo = "xxhash32"), 1, 7), 16))

    t1 <- proc.time()
    if (method == "CONF") {
      obj_test <- conf_test(X1, X2, Y1, Y2)
    } else if (method == "ECF") {
      obj_test <- ecf_test(X1, X2, Y1, Y2, B = B)
    } else if (method == "CED") {
      obj_test <- tcdt(X1, X2, Y1, Y2, stat = "ced", h = "cv", B = B, nu = nu)
    } else if (method == "CMMD") {
      obj_test <- tcdt(X1, X2, Y1, Y2, stat = "cmmd", h = "cv", B = B, nu = nu)
    } else {
      stop("`method` not matched.")
    }
    t2 <- proc.time()

    reject[ind_simu, method] <- 1 * (obj_test$pvalue <= 0.05)
    timecost[ind_simu, method] <- (t2 - t1)[3]
  }

  cat(glue("\r> Simulation {simu} finished."))
}

save.image(file.path(dir, glue("{name_simulation}.RData")))
