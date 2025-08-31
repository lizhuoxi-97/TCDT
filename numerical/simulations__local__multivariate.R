pkg_all <- c(
  "ggplot2", "dplyr", "glue", "TruncatedNormal", "gridExtra", "TCDT", "doFuture"
)
for (pkg in pkg_all) {
  suppressWarnings(suppressMessages(library(pkg, character.only = TRUE)))
}

plan(multisession(workers = 120))

source("numerical/utils.R")

seed <- 1L
set.seed(seed, "L'Ecuyer-CMRG")

f1 <- function(x, p = 2) {
  x <- matrix(x, ncol = p)
  rowSums(x^2) / 2
}
f2 <- function(x, p = 2) {
  x <- matrix(x, ncol = p)
  mu <- rowSums(x^2) / 2
  ifelse(x[, 1] * x[, 2] >= 0, mu, -mu)
}

grid_size <- 21
grid_range <- c(-2, 2)
Xgrid <- expand.grid(
  seq(grid_range[1], grid_range[2], length.out=grid_size),
  seq(grid_range[1], grid_range[2], length.out=grid_size)
)

colnames(Xgrid) <- c("X1", "X2")
mu1 <- f1(as.matrix(Xgrid))
mu2 <- f2(as.matrix(Xgrid))
class <- ifelse(mu1 != mu2, 1, 0)

B <- 299
nu <- 2

n1 <- n2 <- 500
p <- 2
signal <- 1
xmarg <- "TN_U"

if (xmarg == "TN_U") {
  X1 <- matrix(rtnorm(n1 * p, lb = grid_range[1], ub = grid_range[2]), ncol = p)
  X2 <- matrix(runif(n2 * p, grid_range[1], grid_range[2]), ncol = p)
} else {
  stop("`xmarg` not matched.")
}

Y1 <- signal * f1(X1) + rnorm(n1)
Y2 <- signal * f2(X2) + rnorm(n2)


methods <- c("CED", "CMMD")
pvalue <- reject <- matrix(
  nrow = nrow(Xgrid), ncol = length(methods), dimnames = list(NULL, methods)
)
for (method in methods) {
  set.seed(seed, "L'Ecuyer-CMRG")
  cat(glue("\r> Using `foreach` to calculate method: {method}."))
  if (method == "CED") {
    obj_foreach <- foreach(ind = seq_len(nrow(Xgrid)), .options.future = list(seed = TRUE)) %dofuture% {
      x0 <- Xgrid[ind, ]
      tcdt(X1, X2, Y1, Y2, x0, stat = "ced", h = "cv", B = B, nu = nu)$pvalue
    }
  } else if (method == "CMMD") {
    obj_foreach <- foreach(ind = seq_len(nrow(Xgrid)), .options.future = list(seed = TRUE)) %dofuture% {
      x0 <- Xgrid[ind, ]
      tcdt(X1, X2, Y1, Y2, x0, stat = "cmmd", h = "cv", B = B, nu = nu)$pvalue
    }
  } else {
    stop("`method` not matched.")
  }

  pvalue[, method] <- unlist(obj_foreach)
  reject[, method] <- ifelse(pvalue[, method] <= 0.05, 1, NA)
}

name_file <- glue(
  "xmarg_{xmarg}",
  "n_{n1}",
  "signal_{signal}",
  "grid_{grid_size}",
  .sep = "__"
)

dir_result <- file.path("output", "results", "simulations", "local", "multivariate")
if (!dir.exists(dir_result)) {
  dir.create(dir_result, recursive = TRUE)
}
save.image(
  file.path(dir_result, glue("{name_file}.RData"))
)

