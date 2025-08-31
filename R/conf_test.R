#' @title Conformal test for conditional distributions
#' @description This function implements the two-sample conditional distribution 
#' test using conformal prediction, proposed by Hu and Lei (2024).
#' @param X1 A matrix or data frame of covariates for the first sample.
#' @param X2 A matrix or data frame of covariates for the second sample.
#' @param Y1 A matrix or data frame of outcomes for the first sample.
#' @param Y2 A matrix or data frame of outcomes for the second sample.
#' @param prop The proportion of data to be used for calculating the test statistic. 
#' Default is 0.5.
#' @param folds The number of folds for cross-validation. Default is NULL.
#' @param bandwidth The bandwidth for kernel density estimation. Default is NULL.
#' @param est_method The method for estimating the density ratio. 
#' One of "KLR", "LL", "QL", "NN", "KLIEP". Default is "KLR".
#' @return A list containing the test statistic, variance, z-score, and p-value.
#' @references Hu, X., & Lei, J. (2024). A two-sample conditional distribution 
#' test using conformal prediction and weighted rank sum. Journal of the American 
#' Statistical Association, 119(546), 1136-1154.
#' @export
conf_test <- function(X1, X2, Y1, Y2, prop = 0.5, folds = NULL, bandwidth = NULL,
                      est_method = c("KLR", "LL", "QL", "NN", "KLIEP")) {
  list_return <- list()

  # Match the method argument
  est_method <- match.arg(est_method)

  # Coerce and check
  for (name_var in c("X1", "X2", "Y1", "Y2")) {
    obj <- eval(rlang::parse_expr(name_var))
    if (inherits(obj, "data.frame")) {
      assign(name_var, as.matrix(obj))
    } else if (inherits(obj, "numeric")) {
      assign(name_var, matrix(obj, nrow = length(obj)))
    } else if (!inherits(obj, "matrix")) {
      stop(paste(name_var, "should be a vector, matrix or data frame."))
    }
  }
  if (ncol(X1) != ncol(X2)) {
    stop("`X1` and `X2` should have the same number of columns.")
  }
  if (ncol(Y1) != ncol(Y2)) {
    stop("`Y1` and `Y2` should have the same number of columns.")
  }

  # Sample sizes
  n1 <- nrow(X1)
  n2 <- nrow(X2)

  # Data shuffling
  perm1 <- sample(1:n1)
  perm2 <- sample(1:n2)
  X1 <- X1[perm1, , drop = FALSE]
  X2 <- X2[perm2, , drop = FALSE]
  Y1 <- Y1[perm1, , drop = FALSE]
  Y2 <- Y2[perm2, , drop = FALSE]
  
  # Data splitting
  n1_test <- ceiling(n1 * prop)
  n2_test <- ceiling(n2 * prop)
  n1_estr <- n1 - n1_test
  n2_estr <- n2 - n2_test
  
  X1_estr <- X1[1:n1_estr, , drop = FALSE]
  X1_test <- X1[-(1:n1_estr), , drop = FALSE]
  Y1_estr <- Y1[1:n1_estr, , drop = FALSE]
  Y1_test <- Y1[-(1:n1_estr), , drop = FALSE]
  
  X2_estr <- X2[1:n2_estr, , drop = FALSE]
  X2_test <- X2[-(1:n2_estr), , drop = FALSE]
  Y2_estr <- Y2[1:n2_estr, , drop = FALSE]
  Y2_test <- Y2[-(1:n2_estr), , drop = FALSE]
  
  # Estimate density ratios
  ratios <- estimate_r(
    X1_estr, X1_test, X2_estr, X2_test,
    Y1_estr, Y1_test, Y2_estr, Y2_test, est_method
  )
  g1 <- ratios$g1
  g2 <- ratios$g2
  v1 <- ratios$v1
  v2 <- ratios$v2
  
  # Calculate test statistic
  ecdf1 <- ecdf(v2)(v1)
  ecdf2 <- sapply(v1, function(x) mean(x == v2))
  ecdf_est1 <- (ecdf1 + (ecdf1 - ecdf2)) / 2
  zeta <- runif(n2_test)
  inner_sum <- sapply(v1, function(x) mean(zeta * (x == v2)))
  U <- mean(g1 * ((1 - ecdf1) + inner_sum)) / mean(g1)
  
  # Calculate the asymptotic variance under 1st sample
  ss <- var(g1 * (1 - ecdf_est1))
  var1 <- (ss / n1_test + 1 / 12 / n2_test) +
    (var(g1) / n1_test) / 4 -
    cov(g1 * (1 - ecdf_est1), g1) / n1_test
  
  # Calculate the asymptotic variance under 2nd sample
  ecdf3 <- ecdf(v2)(v2)
  ecdf4 <- sapply(v2, function(x) mean(x == v2))
  ecdf_est2 <- (ecdf3 + (ecdf3 - ecdf4)) / 2
  ss <- pmax(mean(g2 * (1 - ecdf_est2)^2) - 1 / 4, 0)
  var2 <- (ss / n1_test + 1 / 12 / n2_test) +
    pmax((mean(g2) - 1) / n1_test, 0) / 4 -
    (mean(g2 * (1 - ecdf3)) + mean(g2 * ecdf4) / 2 - 1 / 2) / n1_test
  
  var_hm <- 2 / (1 / var1 + 1 / var2) # harmonic mean
  
  z <- (U - 0.5) / sqrt(var_hm)
  pvalue <- pnorm(z)
  
  list_return$U <- U
  list_return$var <- var_hm
  list_return$z <- z
  list_return$pvalue <- pvalue

  list_return
}

estimate_r <- function(x1_train, x1_test, x2_train, x2_test,
                       y1_train, y1_test, y2_train, y2_test,
                       est_method = c("KLR", "LL", "QL", "NN", "KLIEP")) {
  est_method <- match.arg(est_method)

  n1_train <- nrow(x1_train)
  n1_test <- nrow(x1_test)
  n2_train <- nrow(x2_train)
  n2_test <- nrow(x2_test)
  label_fit <- factor(c(rep(0, n1_train), rep(1, n2_train)))

  if (est_method == "LL") {
    x_fit <- rbind(x1_train, x2_train)
    fit_marginal <- glm(label_fit ~ ., data = as.data.frame(x_fit), family = binomial())
    new_data <- rbind(x1_test, x2_test)
    prob_marginal <- predict(fit_marginal, newdata = as.data.frame(new_data), type = "response")

    xy_fit <- cbind(rbind(x1_train, x2_train), rbind(y1_train, y2_train))
    fit_joint <- glm(label_fit ~ ., data = as.data.frame(xy_fit), family = binomial())
    new_data <- cbind(new_data, rbind(y1_test, y2_test))
    prob_joint <- predict(fit_joint, newdata = as.data.frame(new_data), type = "response")
  } else if (est_method == "QL") {
    x_fit <- poly(rbind(x1_train, x2_train), degree = 2, raw = TRUE)
    fit_marginal <- glm(label_fit ~ ., data = as.data.frame(x_fit), family = binomial())
    new_data <- poly(rbind(x1_test, x2_test), degree = 2, raw = TRUE)
    prob_marginal <- predict(fit_marginal, newdata = as.data.frame(new_data), type = "response")

    xy_fit <- poly(cbind(rbind(x1_train, x2_train), rbind(y1_train, y2_train)), degree = 2, raw = TRUE)
    fit_joint <- glm(label_fit ~ ., data = as.data.frame(xy_fit), family = binomial())
    new_data <- poly(cbind(rbind(x1_test, x2_test), rbind(y1_test, y2_test)), degree = 2, raw = TRUE)
    prob_joint <- predict(fit_joint, newdata = as.data.frame(new_data), type = "response")
  } else if (est_method == "KLR") {
    klrlearner <- CVST::constructKlogRegLearner()
    sigma <- 0.005
    lambda_seq <- 10^seq(log(1, 10), log(1e-3, 10), length.out = 10)
    params_CV <- CVST::constructParams(
      kernel = "rbfdot",
      sigma = sigma,
      lambda = lambda_seq / (n1_test + n2_test),
      tol = 10e-6,
      maxiter = 500
    )

    x_fit <- rbind(x1_train, x2_train)
    data_fit <- CVST::constructData(x_fit, label_fit)
    newlabel <- as.factor(c(rep(0, n1_test), rep(1, n2_test)))
    newdata <- rbind(x1_test, x2_test)
    data_validation <- CVST::constructData(newdata, newlabel)
    params <- cv_KLR(data_fit, klrlearner, params_CV, fold = FALSE, data_validation = data_validation, verbose = FALSE)[[1]]
    fit_marginal <- klrlearner$learn(data_fit, params)
    K <- kernlab::kernelMult(fit_marginal$kernel, newdata, fit_marginal$data, fit_marginal$alpha)
    prob_marginal <- 1 / (1 + exp(-as.vector(K)))

    xy_fit <- cbind(rbind(x1_train, x2_train), rbind(y1_train, y2_train))
    data_fit <- CVST::constructData(xy_fit, label_fit)
    newlabel <- as.factor(c(rep(0, n1_test), rep(1, n2_test)))
    newdata <- cbind(rbind(x1_test, x2_test), rbind(y1_test, y2_test))
    data_validation <- CVST::constructData(newdata, newlabel)
    params <- cv_KLR(data_fit, klrlearner, params_CV, fold = FALSE, data_validation = data_validation, verbose = FALSE)[[1]]
    fit_joint <- klrlearner$learn(data_fit, params)
    K <- kernlab::kernelMult(fit_joint$kernel, newdata, fit_joint$data, fit_joint$alpha)
    prob_joint <- 1 / (1 + exp(-as.vector(K)))
  } else if (est_method == "NN") {
    hidden_layers <- c(10, 10)
    learn_rates <- 0.001
    n_epochs <- 500

    x_fit <- rbind(x1_train, x2_train)
    newdata <- rbind(x1_test, x2_test)
    prob_marginal <- NNfun(x_fit, label_fit, newdata,
                           nnrep = 5, hidden_layers = hidden_layers,
                           n_epochs = n_epochs, learn_rates = learn_rates
    )

    xy_fit <- cbind(rbind(x1_train, x2_train), rbind(y1_train, y2_train))
    newdata <- cbind(rbind(x1_test, x2_test), rbind(y1_test, y2_test))
    prob_joint <- NNfun(xy_fit, label_fit, newdata,
                        nnrep = 5, hidden_layers = hidden_layers,
                        n_epochs = n_epochs, learn_rates = learn_rates
    )
  } else if (est_method == "KLIEP") {
    fit_marginal <- densratio::densratio(x2_train, x1_train, method = "KLIEP", verbose = FALSE)

    xy1 <- cbind(x1_train, y1_train)
    xy2 <- cbind(x2_train, y2_train)
    fit_joint <- densratio::densratio(xy2, xy1, method = "KLIEP", verbose = FALSE)
  }

  if (est_method == "KLIEP") {
    g1 <- fit_marginal$compute_density_ratio(x1_test)
    g2 <- fit_marginal$compute_density_ratio(x2_test)
    v1 <- fit_joint$compute_density_ratio(cbind(x1_test, y1_test)) * g1
    v2 <- fit_joint$compute_density_ratio(cbind(x2_test, y2_test)) * g2
  } else {
    prob_marginal[prob_marginal < 0.01] <- 0.01
    prob_marginal[prob_marginal > 0.99] <- 0.99
    prob_joint[prob_joint < 0.01] <- 0.01
    prob_joint[prob_joint > 0.99] <- 0.99

    g1 <- prob_marginal[1:n1_test] / (1 - prob_marginal[1:n1_test]) * n1_train / n2_train
    g2 <- prob_marginal[(n1_test + 1):(n1_test + n2_test)] / (1 - prob_marginal[(n1_test + 1):(n1_test + n2_test)]) * n1_train / n2_train
    v1 <- (1 - prob_joint[1:n1_test]) / prob_joint[1:n1_test] * g1
    v2 <- (1 - prob_joint[(n1_test + 1):(n1_test + n2_test)]) / prob_joint[(n1_test + 1):(n1_test + n2_test)] * g2
  }

  list(g1 = g1, g2 = g2, v1 = v1, v2 = v2)
}

cv_KLR <- function(data, learner, params, fold = 5, data_validation, verbose = TRUE) {
  stopifnot(inherits(learner, "CVST.learner") &&
    inherits(data, "CVST.data") &&
    inherits(params, "CVST.params"))
  nParams <- length(params)

  data <- CVST::shuffleData(data)

  if (fold) {
    dimnames <- list(as.character(1:fold), names(params))
    results <- matrix(0, fold, nParams, dimnames = dimnames)
    size <- CVST::getN(data) / fold
  } else {
    results <- rep(0, nParams)
    names(results) <- names(params)
  }

  for (ind in 1:nParams) {
    p <- params[[ind]]
    if (fold) {
      for (f in 1:fold) {
        validationIndex <- seq((f - 1) * size + 1, f * size)
        curTrain <- CVST::getSubset(data, -validationIndex)
        curTest <- CVST::getSubset(data, validationIndex)
        # either mean squared error or mean classification error
        results[f, ind] <- mean(eval_KLR(curTrain, curTest, learner, p))
      }
      if (verbose) {
        cat(names(params)[ind], "(", mean(results[, ind]), ")\n")
      }
      winner <- which.min(apply(results, 2, mean))
    } else {
      results[ind] <- mean(eval_KLR(data, data_validation, learner, p))
      winner <- which.min(results)
    }
  }

  if (length(winner) == 0) {
    return(NULL)
  } else {
    return(params[winner])
  }
}

eval_KLR <- function(train, test, learner, param) {
  stopifnot(inherits(learner, "CVST.learner") && inherits(train, "CVST.data") && inherits(test, "CVST.data"))
  model <- try(learner$learn(train, param))
  if (inherits(model, "try-error")) {
    prob_pre <- rep(NA, length(test$y))
  } else {
    K <- try(kernlab::kernelMult(model$kernel, test$x, model$data, model$alpha))
    if (inherits(K, "try-error")) {
      prob_pre <- rep(NA, length(test$y))
    } else {
      prob_pre <- 1 / (1 + exp(-as.vector(K)))
      prob_pre[prob_pre < 0.01] <- 0.01
      prob_pre[prob_pre > 0.99] <- 0.99
    }
  }

  y <- 1 * (test$y != levels(test$y)[1])
  res <- mean(-y * log(prob_pre) - (1 - y) * log(1 - prob_pre))
  res
}

NNfun <- function(x, z, xpre, nnrep = 10, hidden_layers = NA,
                  acfun = "sigmoid", optim_type = "sgd", n_epochs = 500,
                  learn_rates = 0.001, L1 = 0) {
  n <- nrow(xpre)
  prob <- matrix(0, nrow = n, ncol = nnrep)

  for (i in 1:nnrep) {
    fit_nn <- ANN2::neuralnetwork(x, z,
      hidden.layers = hidden_layers, optim.type = optim_type,
      val.prop = 0, learn.rates = learn_rates, L1 = L1,
      n.epochs = n_epochs, activ.functions = acfun, verbose = 0
    )

    prob[, i] <- predict(fit_nn, xpre)$probabilities[, 2]
  }

  prob_fit <- rowMeans(prob)
  prob_fit[prob_fit < 0.01] <- 0.01
  prob_fit[prob_fit > 0.99] <- 0.99

  prob_fit
}

