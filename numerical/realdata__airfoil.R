pkg_all <- c("purrr", "dplyr", "readr", "glue", "TCDT", "tictoc", "digest")
for (pkg in pkg_all) {
  suppressWarnings(suppressMessages(library(pkg, character.only = TRUE)))
}

data <- read_table(
  "data/airfoil_self_noise.dat",
  col_names = c("freq", "angle", "chord", "velocity", "thickness", "sound")
)

data <- data %>%
  mutate(logfreq = log(freq), logthick = log(thickness)) %>%
  select(logfreq, angle, chord, velocity, logthick, sound)
X <- as.matrix(select(data, -sound))
y <- as.matrix(select(data, sound))

N <- nrow(data)
p <- ncol(data)

args <- commandArgs(TRUE)
if (length(args) > 0) {
  setting <- args[1]
} else if (length(args) == 0) {
  setting <- "cov_shift_null"
}

B <- 299

name <- glue("realdata__airfoil__{setting}")

dgp_realdata2 <- function(X, y,
                          setting = c("cov_shift_null", "prior_shift_null",
                                      "cov_shift_alt", "prior_shift_alt")) {
  setting <- match.arg(setting)

  if (setting == "cov_shift_null") {
    n1 <- round(N * 0.2)
    n2 <- N - n1

    idx <- sample(1:N, n1, replace = FALSE)
    X1 <- X[idx, ]
    X2 <- X[-idx, ]
    Y1 <- y[idx, , drop = FALSE]
    Y2 <- y[-idx, , drop = FALSE]

    weight <- exp(X2[, c(1, 5)] %*% c(-1, 1))
    idx2 <- sample(1:n2, size = ceiling(n2 * 0.25), prob = weight, replace = FALSE)
    X2 <- X2[idx2, ]
    Y2 <- Y2[idx2, , drop = FALSE]
  } else if (setting == "prior_shift_null") {
    n1 <- round(N * 0.2)
    n2 <- N - n1

    idx <- sample(1:N, n1, replace = FALSE)
    X1 <- y[idx, , drop = FALSE]  # Consider testing P_{X|Y} later, so here swap X and y for convenience
    X2 <- y[-idx, , drop = FALSE]
    Y1 <- X[idx, ]
    Y2 <- X[-idx, ]

    weight <- exp(X2)
    idx2 <- sample(1:n2, size = ceiling(n2 * 0.25), prob = weight, replace = FALSE)
    X2 <- X2[idx2, , drop = FALSE]
    Y2 <- Y2[idx2, ]
  } else if (setting == "cov_shift_alt") {
    n1 <- round(N / 2)
    n2 <- N - n1

    idx <- order(y)[1:n1]
    X1 <- X[idx, ]
    X2 <- X[-idx, ]
    Y1 <- y[idx, , drop = FALSE]
    Y2 <- y[-idx, , drop = FALSE]

    n1 <- nrow(X1)
    n2 <- nrow(X2)

    flipsize <- round(n1 * 0.05)
    idx1 <- sample(1:n1, flipsize, replace = FALSE)
    idx2 <- sample(1:n2, flipsize, replace = FALSE)
    X2_tmp <- X2[idx2, ]
    Y2_tmp <- Y2[idx2, ]
    X1_tmp <- X1[idx1, ]
    Y1_tmp <- Y1[idx1, ]

    X1[idx1, ] <- X2_tmp
    Y1[idx1, ] <- Y2_tmp
    X2[idx2, ] <- X1_tmp
    Y2[idx2, ] <- Y1_tmp
  } else if (setting == "prior_shift_alt") {
    n1 <- round(N / 2)
    n2 <- N - n1

    idx <- order(X[, "chord"])[1:n1]
    X1 <- y[idx, , drop = FALSE] # Consider testing P_{X|Y} later, so here swap X and y for convenience
    X2 <- y[-idx, , drop = FALSE]
    Y1 <- X[idx, ]
    Y2 <- X[-idx, ]

    n1 <- nrow(X1)
    n2 <- nrow(X2)

    flipsize <- round(n1 * 0.05)
    idx1 <- sample(1:n1, flipsize, replace = FALSE)
    idx2 <- sample(1:n2, flipsize, replace = FALSE)
    X2_tmp <- X2[idx2, ]
    Y2_tmp <- Y2[idx2, ]
    X1_tmp <- X1[idx1, ]
    Y1_tmp <- Y1[idx1, ]

    X1[idx1, ] <- X2_tmp
    Y1[idx1, ] <- Y2_tmp
    X2[idx2, ] <- X1_tmp
    Y2[idx2, ] <- Y1_tmp
  }

  list(
    X1 = X1, X2 = X2, Y1 = Y1, Y2 = Y2
  )
}

if (setting %in% c("cov_shift_null", "prior_shift_null")) {
  num_simu <- 500
} else if (setting %in% c("cov_shift_alt", "prior_shift_alt")) {
  num_simu <- 1
}

local_bootstrap <- "sep"
if (setting %in% c("cov_shift_null", "cov_shift_alt")) {
  nu <- 4
} else if (setting %in% c("prior_shift_null", "prior_shift_alt")) {
  nu <- 2
}

methods <- c("CED", "CMMD", "CONF")

pvalue <- reject <- matrix(nrow = num_simu, ncol = length(methods), dimnames = list(NULL, methods))

for (simu in 1:num_simu) {
  set.seed(strtoi(substr(digest(simu, algo="xxhash32"), 1, 7), 16)) # Fixed seed for dgp
  cat(glue("> Simulation {simu}"), "\n")

  obj_data <- dgp_realdata2(X, y, setting)
  X1 <- obj_data$X1
  X2 <- obj_data$X2
  Y1 <- obj_data$Y1
  Y2 <- obj_data$Y2

  tic()
  for (method in methods) {
    if (method == "CONF") {
      obj_test <- conf_test(X1, X2, Y1, Y2)
    } else if (method == "CED") {
      obj_test <- tcdt(
        X1, X2, Y1, Y2, stat = "ced", h = "cv", B = B, nu = nu
      )
    } else if (method == "CMMD") {
      obj_test <- tcdt(
        X1, X2, Y1, Y2, stat = "cmmd", h = "cv", B = B, nu = nu
      )
    } else {
      stop("`method` not matched.")
    }
    pvalue[simu, method] <- obj_test$pvalue
    reject[simu, method] <- 1 * (pvalue[simu, method] <= 0.05)
  }
  toc()
}

path <- file.path("output", "results", "realdata")
if (!dir.exists(path)) {
  dir.create(path, recursive = TRUE)
}
save.image(file.path(path, glue("{name}.RData")))
