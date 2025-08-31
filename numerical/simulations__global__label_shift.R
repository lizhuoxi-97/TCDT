pkg_all <- c("TCDT", "digest", "purrr", "glue", "mixtools", "mvtnorm")
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
  p <- as.numeric(args[4]) # dimension of x (y in original dgp)
} else if (type_args == "Rscript") {
  setting_dgp <- "different_mean"
  n1 <- 50
  n2 <- 50
  p <- 5
  setting_m <- "nonlinear"
  shift_y <- "true"
  signal <- 1
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
    shift_y <- args[5]
  }

  ## Marginal distribution of Y ----

  q <- 1 # dimension of y (x in original dgp)
  nu <- 2

  dist_y1 <- dist_y2 <- "normal"
  param1_y1 <- 0
  param2_y1 <- 1

  if (shift_y == "false") {
    param1_y2 <- param1_y1
    param2_y2 <- param2_y1
  } else if (shift_y == "true") {
    param1_y2 <- 0
    param2_y2 <- sqrt(1.5)
  }

  ## Conditional distribution settings, and methods to be compared ----
  if (setting_dgp == "different_mean") {
    ### global, label shift, different mean ----
    # params:
    # 1. setting_dgp
    # 2. n1
    # 3. n2
    # 4. p
    # 5. shift_y
    # 6. setting_m
    # 7. signal
    # 8. num_simu

    if (type_args == "shell") {
      setting_m <- args[6]
      signal <- as.numeric(args[7])
    }

    name_simulation <- glue(
      "m_{setting_m}",
      "{signal}",
      "n_{n1}_{n2}",
      "p_{p}",
      "shift_y_{shift_y}",
      .sep = "__"
    )

    methods <- c("CED", "CMMD", "CONF")

    g0_mix <- list(
      d = function(x) {
        p <- length(x)
        0.5 * dmvnorm(x, rep(0, p), diag(1, p)) + 0.5 * dmvnorm(x, rep(1, p), diag(1, p))
      },
      r = function(m, p) {
        rmvnormmix(
          m,
          c(0.5, 0.5),
          matrix(
            c(rep(0, p), rep(1, p)),
            nrow = 2,
            byrow = TRUE
          ),
          matrix(
            c(rep(1, p), rep(1, p)),
            nrow = 2,
            byrow = TRUE
          )
        )
      }
    )

    dist_eps1 <- dist_eps2 <- "normal"
    rho <- 0.5
    sigma <- rho^abs(outer(1:p, 1:p, "-"))
    param_eps1 <- param_eps2 <- sigma

    form_v1 <- form_v2 <- function(y) {
      1
    }
    signal_v <- NULL

    signal_m <- signal
    if (setting_m == "nonlinear") {
      #### nonlinear ----
      form_m1 <- function(y) {
        if (p == 5) {
          mu_base <- c(-0.5, -0.5, 0.5, 0.5, 1)
        } else if (p == 20) {
          mu_base <- c(rep(-0.5, 3), rep(0.5, 2), rep(1, 5), rep(0, 10))
        }
        theta <- drop(y) %o% mu_base
        (theta^2 + 3 * theta + 2) / (theta^2 + 1)
      }
      form_m2 <- function(y, signal) {
        if (p == 5) {
          mu_base <- (1 + signal) * c(-0.5, -0.5, 0.5, 0.5, 1)
        } else if (p == 20) {
          mu_base <- (1 + signal) * c(rep(-0.5, 3), rep(0.5, 2), rep(1, 5), rep(0, 10))
        }
        theta <- drop(y) %o% mu_base
        (theta^2 + 3 * theta + 2) / (theta^2 + 1)
      }
    } else {
      stop("`setting_m` not matched.")
    }
  } else if (setting_dgp == "different_var") {
    ### global, label shift, different var ----
    # params:
    # 1. setting_dgp
    # 2. n1
    # 3. n2
    # 4. p
    # 5. shift_y
    # 6. setting_m
    # 7. setting_v
    # 8. signal
    # 9. num_simu

    if (type_args == "shell") {
      setting_m <- args[6]
      setting_v <- args[7]
      signal <- as.numeric(args[8])
    }

    name_simulation <- glue(
      "v_{setting_v}",
      "{signal}",
      "m_{setting_m}",
      "n_{n1}_{n2}",
      "p_{p}",
      "shift_y_{shift_y}",
      .sep = "__"
    )

    methods <- c("CED", "CMMD", "CONF")

    dist_eps1 <- dist_eps2 <- "normal"
    param_eps1 <- param_eps2 <- diag(1, p)

    signal_m <- NULL
    if (setting_m == "y") {
      form_m1 <- form_m2 <- function(y) {
        if (p == 5) {
          mu_base <- c(-0.5, -0.5, 0.5, 0.5, 0)
        } else if (p == 20) {
          mu_base <- c(rep(-0.5, 5), rep(0.5, 5), rep(0, 10))
        }

        drop(y) %o% mu_base
      }
    }

    signal_v <- signal
    if (setting_v == "hetero") {
      #### hetero ----
      form_v1 <- function(y) {
        # Notice we generate p-vector of cond. var from 1-dim y
        if (p == 5) {
          factor_base <- c(1.5, 1.5, 0.5, 0.5, 1)
        } else if (p == 20) {
          factor_base <- c(rep(1.5, 5), rep(0.5, 5), rep(1, 10))
        }

        (0.1 / (1 + drop(y)^2)) %o% factor_base
      }
      form_v2 <- function(y, signal) {
        # Notice we generate p-vector of cond. var from 1-dim y
        if (p == 5) {
          factor_base <- c(1.5, 1.5, 0.5, 0.5, 1)
        } else if (p == 20) {
          factor_base <- c(rep(1.5, 5), rep(0.5, 5), rep(1, 10))
        }

        ((1 + 3 * signal) * 0.1 / (1 + drop(y)^2)) %o% factor_base
      }
    }
  }

  dgp <- partial(
    dgp_reg,
    n1 = n1, n2 = n2, p = q, q = p,
    dist_x1 = dist_y1, param1_x1 = param1_y1, param2_x1 = param2_y1,
    dist_x2 = dist_y2, param1_x2 = param1_y2, param2_x2 = param2_y2,
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
  "output", "results", "simulations", "global", "label_shift", setting_dgp
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
  set.seed(strtoi(substr(digest(simu, algo = "xxhash32"), 1, 7), 16)) # Fixed seed for dgp

  obj_dgp <- dgp()
  X1 <- obj_dgp$Y1
  X2 <- obj_dgp$Y2
  Y1 <- obj_dgp$X1
  Y2 <- obj_dgp$X2

  for (method in methods) {
    set.seed(strtoi(substr(digest(simu, algo = "xxhash32"), 1, 7), 16)) # Fixed seed for dgp

    t1 <- proc.time()
    if (method == "CONF") {
      obj_test <- conf_test(Y1, Y2, X1, X2)
    } else if (method == "CED") {
      obj_test <- tcdt(
        Y1, Y2, X1, X2, stat = "ced", h = "cv", B = B, nu = nu
      )
    } else if (method == "CMMD") {
      obj_test <- tcdt(
        Y1, Y2, X1, X2, stat = "cmmd",h = "cv", B = B, nu = nu
      )
    } else {
      stop("`method` not matched.")
    }
    t2 <- proc.time()

    reject[ind_simu, method] <- 1 * (obj_test$pvalue <= 0.05)
    timecost[ind_simu, method] <- (t2 - t1)[3]
  }

  cat(glue("\r> Simulation {simu} finished."))
}

rm(obj_dgp)
rm(X1)
rm(X2)

save.image(file.path(dir, glue("{name_simulation}.RData")))
