pkg_all <- c("TCDT", "purrr", "glue", "mvtnorm", "sn", "digest")
for (pkg in pkg_all) {
  suppressWarnings(suppressMessages(library(pkg, character.only = TRUE)))
}

source("numerical/utils.R")

# Fixed settings ----
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
  setting_m <- args[5]
  dist_x <- args[6]
  shift_x <- args[7]
  mode_local <- args[8]
} else if (type_args == "Rscript") {
  setting_dgp <- "different_mean"
  n1 <- 150
  n2 <- 150
  p <- 1
  setting_m <- "x"
  dist_x <- "normt"
  shift_x <- "true"
  num_simu <- 3
  mode_local <- "parallel"
  x0 <- 1
}

name_simulation <- glue(
  "m_{setting_m}",
  "n_{n1}_{n2}",
  "p_{p}",
  "dist_x_{dist_x}",
  "shift_x_{shift_x}",
  "mode_local_{mode_local}",
  .sep = "__"
)

if (mode_local == "parallel") {
  if (type_args == "shell") {
    x0 <- as.numeric(args[9])
  }

  name_simulation <- glue(
    name_simulation,
    "x0_{x0}",
    .sep = "__"
  )
}


if (setting_dgp %in% c("different_mean")) {
  dgp_framework <- "reg"
} else {
  dgp_framework <- "custom"
}

if (dgp_framework == "reg") {
  # Regression DGP dgp_framework ----
  ## Marginal distribution of X ----
  if (dist_x == "unif") {
    ### unifx ----
    dist_x1 <- dist_x2 <- "uniform"
    param1_x1 <- 0
    param2_x1 <- 1

    if (shift_x == "false") {
      param1_x2 <- param1_x1
      param2_x2 <- param2_x1
    } else if (shift_x == "true") {
      param1_x2 <- param1_x1 + 1
      param2_x2 <- param2_x1 + 1
    }
  } else if (dist_x == "norm") {
    ### normx ----
    dist_x1 <- dist_x2 <- "normal"

    if (p == 1) {
      #### p = 1 ----
      param1_x1 <- -0.5
      param2_x1 <- 1

      if (shift_x == "false") {
        param1_x2 <- param1_x1
        param2_x2 <- param2_x1
      } else if (shift_x == "true") {
        param1_x2 <- 0.5
        param2_x2 <- param2_x1
      }
    } else if (p >= 2) {
      #### p >= 2 ----
      param1_x1 <- rep(0, p) # mean vector
      param2_x1 <- diag(1, p) # Covariance matrix

      if (shift_x == "false") {
        param1_x2 <- param1_x1
        param2_x2 <- param2_x1
      } else if (shift_x == "true") {
        param1_x2 <- param1_x1 + 1
        param2_x2 <- param2_x1
      }
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

    q <- 1

    ### local, different mean ----
    methods <- c("CED", "CMMD")
    nu <- 2

    dist_eps1 <- dist_eps2 <- "normal"

    if (setting_m == "x") {
      ##### x ----
      form_m1 <- function(x) 0.5 * x
      form_m2 <- function(x) -0.5 * x
      signal_m <- NULL

      form_v1 <- form_v2 <- function(x) 0.5
      signal_v <- NULL

      if (mode_local == "single") {
        x0_all <- seq(-1.3, 1.3, 0.1)
      } else if (mode_local == "parallel") {
        x0_all <- x0
      }
    } else {
      stop("`setting_m` not matched.")
    }
  }

  dgp <- partial(
    dgp_reg,
    n1 = n1, n2 = n2, p = p, q = q,
    dist_x1 = dist_x1, param1_x1 = param1_x1, param2_x1 = param2_x1,
    dist_x2 = dist_x2, param1_x2 = param1_x2, param2_x2 = param2_x2,
    dist_eps1 = dist_eps1, dist_eps2 = dist_eps2,
    form_m1 = form_m1, form_m2 = form_m2, signal_m = signal_m,
    form_v1 = form_v1, form_v2 = form_v2, signal_v = signal_v
  )
} else {
  cat("> Using custom DGP framework.")

  # Specify your custom dgp function as `dgp()` here:

}

# Path and filename for storing results ----
dir <- file.path(
  "output", "results", "simulations", "local", "univariate", setting_dgp
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

reject <- array(
  dim = c(length(simu_all), length(methods), length(x0_all)),
  dimnames = list(simu_all, methods, format(x0_all, trim = TRUE, nsmall = 1))
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

  for (x0 in x0_all) {
    x0_chr <- format(x0, trim = TRUE, nsmall = 1)

    for (method in methods) {
      # Fixed seed for each method
      set.seed(strtoi(substr(digest(simu, algo = "xxhash32"), 1, 7), 16))

      t1 <- proc.time()
      if (method == "CED") {
        obj_test <- tcdt(X1, X2, Y1, Y2, x0, stat = "ced", h = "cv", B = B, nu = nu)
      } else if (method == "CMMD") {
        obj_test <- tcdt(X1, X2, Y1, Y2, x0, stat = "cmmd", h = "cv", B = B, nu = nu)
      } else {
        stop("`method` not matched.")
      }

      reject[ind_simu, method, x0_chr] <- 1 * (obj_test$pvalue <= 0.05)
    }
  }

  cat(glue("\r> Simulation {simu} finished."))
}

save.image(file.path(dir, glue("{name_simulation}.RData")))
