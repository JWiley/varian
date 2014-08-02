#' Use a Variability Model to predict another outcome
#'
#' This function uses a linear mixed effects model that assumes the level 1 residual
#' variance varies by Level 2 units.  That is rather than assuming a homogenous residual
#' variance, it assumes the residual standard deviations come from a Gamma distribution.
#' In the first stage of this model, each Level 2's residual standard deviation is
#' estimated, and in the second stage, these standard deviations are used to predict
#' another Level 2 outcome.  The interface uses an intuitive formula interface, but
#' the underlying model is implemented in Stan, with minimally informative priors for all
#' parameters.
#'
#' @param formula A formula describing the model to estimate
#' @param var The Level 1 outcome from which to estimate each Level 2 units residual
#'   standard deviation
#' @param data A long data frame containing an both the Level 2 and Level 1 outcomes,
#'   as well as all covariates and an ID variable.
#' @param totaliter The total number of iterations to be used (not including the
#'   warmup iterations), these are distributed equally across multiple independent
#'   chains.
#' @param warmup The number of warmup iterations.  Each independent chain
#'   has the same number of warmup iterations, before it starts the iterations
#'   that will be used for inference.
#' @param chains The number of independent chains to run (default to 1).
#' @param useU A logical value whether the latent intercept estimated in Stage 1 should
#'   also be used as a predictor.  Defaults to \code{TRUE}.
#' @param inits Initial values passed on to \code{stan}.  If \code{NULL}, the default,
#'   initial values are estimated means, standard deviations, and coefficients from a
#'   single level linear regression.
#' @param opts A list giving options.  Currently only \code{SD_Tol} which controls
#'   the tolerance for how small a variables standard deviation may be without
#'   stopping estimation (this ensures that duplicate variables, or variables without
#'   any variability are included as predictors).
#' @param \dots Additional arguments passed to \code{stan}.
#' @return A named list containing the model \code{Output},
#'   the \code{variableNames}, the \code{data}, and the initial function
#'   \code{.call}.
#' @author Joshua F. Wiley <josh@@elkhartgroup.com>
#' @export
#' @keywords models
#' @examples
#' \dontrun{
#'   sim.data <- with(simulate_gvm(4, 60, 0, 1, 3, 2, 94367), {
#'     set.seed(265393)
#'     x2 <- MASS::mvrnorm(k, c(0, 0), matrix(c(1, .3, .3, 1), 2))
#'     y2 <- rnorm(k, cbind(Int = 1, x2) %*% matrix(c(3, .5, .7)) + sigma, sd = 3)
#'     data.frame(
#'       y = Data$y,
#'       y2 = y2[Data$ID2],
#'       x1 = x2[Data$ID2, 1],
#'       x2 = x2[Data$ID2, 2],
#'       ID = Data$ID2)
#'   })
#'   m <- vm_predict(y2 ~ x1 + x2 + y + ID, var = "y", data = sim.data,
#'     totaliter = 10000, warmup = 1500, thin = 10, chains = 4, verbose=TRUE)
#'
#'   # check diagnostics
#'   vm_diagnostics(m)
#'
#'   sim.data2 <- with(simulate_gvm(21, 250, 0, 1, 3, 2, 94367), {
#'     set.seed(265393)
#'     x2 <- MASS::mvrnorm(k, c(0, 0), matrix(c(1, .3, .3, 1), 2))
#'     y2 <- rnorm(k, cbind(Int = 1, x2) %*% matrix(c(3, .5, .7)) + sigma, sd = 3)
#'     data.frame(
#'       y = Data$y,
#'       y2 = y2[Data$ID2],
#'       x1 = x2[Data$ID2, 1],
#'       x2 = x2[Data$ID2, 2],
#'       ID = Data$ID2)
#'   })
#'   # warning: may take several minutes
#'   m2 <- vm_predict(y2 ~ x1 + x2 + y + ID, var = "y", data = sim.data2,
#'     totaliter = 10000, warmup = 1500, thin = 10, chains = 4, verbose=TRUE)
#'   # check diagnostics
#'   vm_diagnostics(m2)
#' }
vm_predict <- function(formula, var, data, totaliter = 2000, warmup = 1000, chains = 1,
  useU = TRUE, inits = NULL, opts = list(SD_Tol = .01), ...) {
  stopifnot(is.data.frame(data))

  storedCall <- match.call()

  # drop any missing levels to avoid redudant dummy codes in model matrix
  data <- droplevels(data)

  key <- list()

  IDvar <- as.character(formula[[3]])
  IDvar <- IDvar[length(IDvar)]

  # make sure ID is a numeric/integer or a factor
  stopifnot(class(data[, IDvar]) %in% c("numeric", "integer", "factor"))

  test.IDVar <- sd_id(data[, var], data[, IDvar], long=FALSE)

  if (!all(test.IDVar != 0, na.rm=TRUE)) {
    stop(sprintf("The following IDs have no variability in the first stage outcome:\n%s\nTry using\n%s\nto remove these from the data.",
      paste(names(test.IDVar)[which(test.IDVar == 0)], collapse = ', '),
      paste0("subset(your_data, sd_id(", var, ", ", IDvar, ") != 0)")))
  }

  key$OriginalID <- data[, IDvar]

  data[, IDvar] <- as.integer(data[, IDvar])

  key$IntegerID <- data[, IDvar]

  data <- data[order(data[, IDvar]), ]
  formula2 <- terms.formula(formula, data = data, keep.order = TRUE)
  mf <- model.frame(formula2, data = data, na.action = na.omit)

  key$MFOriginalID <- mf[, IDvar]

  mf[, IDvar] <- as.integer(factor(mf[, IDvar]))

  key$MFIntegerID <- mf[, IDvar]

  key <- with(key, {
    tmp1 <- data.frame(OriginalID, IntegerID)[!duplicated(IntegerID), ]
    tmp2 <- data.frame(MFOriginalID, MFIntegerID)[!duplicated(MFIntegerID), ]
    data.frame(OriginalID = tmp1[match(tmp2[, 1], tmp1[, 2]), 1],
               InternalID = tmp2[, 2])
  })

  mm <- model.matrix(formula2, data = mf)

  i <- which(colnames(mm) == var)
  if((ncol(mm) - 1) > i) {
    iX1 <- c(1, (i + 1):(ncol(mm) - 1))
  } else {
    iX1 <- 1
  }
  varNames <- list(
    Y2 = colnames(mf)[1],
    Y1 = var,
    X2 = colnames(mm)[1:(i - 1)],
    X1 = colnames(mm)[iX1],
    ID1 = colnames(mm)[ncol(mm)])

  Y1 <- mf[, varNames$Y1]
  X1 <- mm[, varNames$X1, drop = FALSE]
  ID1 <- mf[, varNames$ID1]
  keepObs <- !duplicated(ID1)
  Y2 <- mf[keepObs, varNames$Y2]
  X2 <- mm[keepObs, varNames$X2, drop = FALSE]

  # Code to check the scaling of variables
  v.sds <- c(sd(Y1, na.rm=TRUE), apply(X1, 2, sd, na.rm=TRUE),
             sd(Y2, na.rm=TRUE), apply(X2, 2, sd, na.rm=TRUE))
  names(v.sds) <- c(varNames$Y1, varNames$X1, varNames$Y2, varNames$X2)

  v.sds <- v.sds[names(v.sds) != "(Intercept)"]

  v.sds.index <- is.na(v.sds) | v.sds < opts$SD_Tol | v.sds > 50
  if (any(v.sds.index)) {
    stop(sprintf("The follow variables SDs are either too small or too large.\n  Remove or rescale variables before modelling.\n  Variables: %s",
                 paste(names(v.sds)[v.sds.index], collapse = ", ")))
  }

  class.tests <- c(
    Y1 = is.numeric(Y1) & is.vector(Y1),
    X1 = is.null(X1) || (is.matrix(X1) & is.numeric(X1)),
    Y2 = is.numeric(Y1) & is.vector(Y1),
    X2 = is.null(X2) || (is.matrix(X2) & is.numeric(X2)))

  if (!all(class.tests)) {
    stop(c("Y1 must be a numeric vector", "X1 must be NULL or a numeric matrix",
      "Y2 must be a numeric vector", "X2 must be NULL or a numeric matrix")[class.tests])
  }
  n <- length(Y1)
  k <- length(unique(ID1))

  if (!all(is.null(X1) || identical(nrow(X1), n), identical(length(ID1), n))) {
    stop("The length and rows of Y1, X1, and ID1 must all be equal")
  }

  if (!all(is.null(X2) || identical(nrow(X2), k), identical(length(Y2), k))) {
    stop("The length and rows of Y2, X2, and the unique IDs in ID1 must all be equal")
  }

  p1 <- ncol(X1)
  p2 <- ncol(X2)

  if (useU) {
    p3 <- 2L
    Y2_model <- "Y2_hat[i] <- (X2[i] * B2) + alpha[1] * Sigma_Y1[i] + alpha[2] * U[i];"
  } else {
    p3 <- 1L
    Y2_model <- "Y2_hat[i] <- (X2[i] * B2) + alpha[1] * Sigma_Y1[i];"
  }

  stan.data <- list(
    Y1 = Y1, X1 = X1, ID1 = ID1,
    Y2 = Y2, X2 = X2,
    n = n, k = k,
    p1 = p1, p2 = p2, p3 = p3)


  model <- gsub("Y2_model", Y2_model, "
  // upper case letters indicate vectors/matrices
  // lower case letters indicate scalars
  data {
    int<lower=1> n;
    int<lower=1> k;
    int<lower=1> p1;
    int<lower=1> p2;
    int<lower=1> p3;
    real Y1[n];
    real Y2[k];

    matrix[n, p1] X1;
    matrix[k, p2] X2;

    int ID1[n];
  }
  parameters {
    vector[p1] B1;
    real U[k];
    real<lower=0> sigma_U;
    real<lower=0> shape;
    real<lower=0> rate;
    real<lower=0> Sigma_Y1[k];

    vector[p2] B2;
    vector[p3] alpha;
    real<lower=0> sigma_Y2;
  }
  transformed parameters {
    real Y1_hat[n];
    real Sigma_Y1_hat[n];
    real Y2_hat[k];

    for (i in 1:n) {
      Y1_hat[i] <- (X1[i] * B1) + U[ID1[i]];
      Sigma_Y1_hat[i] <- Sigma_Y1[ID1[i]];
    }
    for (i in 1:k) {
      Y2_model
    }
  }
  model {
    // Priors for Stage 1 Location
    B1 ~ normal(0, 1000);
    U ~ normal(0, sigma_U);
    // Priors for Stage 1 scale of random location
    sigma_U ~ cauchy(0, 5);
    // Priors for Stage 1 Scale
    shape ~ cauchy(0, 5);
    rate ~ cauchy(0, 5);
    // Model for Stage 1 Scale
    Sigma_Y1 ~ gamma(shape, rate);

    // Priors for Stage 2 location
    B2 ~ normal(0, 1000);
    alpha ~ normal(0, 1000);
    // Priors for Stage 2 scale
    sigma_Y2 ~ cauchy(0, 5);

    // Likelihood for Stage 1
    Y1 ~ normal(Y1_hat, Sigma_Y1_hat);
    // Likelihood for Stage 2
    Y2 ~ normal(Y2_hat, sigma_Y2);
  }")

  stan_inits <- function() {
    rg <- res_gamma(Y1, ID1)

    out <- with(stan.data, list(
      B1 = as.array(coef(lm.fit(X1, Y1))),
      U = as.array(by_id(Y1, ID1, mean, FALSE, na.rm=TRUE)),
      shape = rg$alpha,
      rate = rg$beta,
      Sigma_Y1 = as.array(sd_id(Y1, ID1, FALSE)),
      sigma_Y2 = sd(Y2, na.rm=TRUE)
    ))

    out$sigma_U <- sd(out$U, na.rm=TRUE)

    if (useU) {
      b2 <- with(stan.data, coef(lm.fit(cbind(X2, Res = out$Sigma_Y1, U = out$U), Y2)))
    } else {
      b2 <- with(stan.data, coef(lm.fit(cbind(X2, Res = out$Sigma_Y1), Y2)))
    }

    nc <- ncol(stan.data$X2)

    out$B2 <- as.array(b2[1:nc])
    out$alpha <- as.array(b2[(nc + 1):length(b2)])

    return(out)
  }

  if (is.null(inits)) inits <- stan_inits

  res <- parallel_stan(model_code = model, standata = stan.data,
    totaliter = totaliter, warmup = warmup, chains = chains,
    pars = c("B1", "U", "sigma_U", "shape", "rate", "Sigma_Y1",
             "B2", "alpha", "sigma_Y2"), init = inits, ...)

  out <- list(
    Output = res,
    variableNames = varNames,
    data = c(stan.data, list(IDkey = key)),
    .call = storedCall
  )

  class(out) <- c("vmp", "list")

  return(out)
}