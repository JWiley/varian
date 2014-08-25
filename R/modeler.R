## data {
##     int<lower=1> N;
##     real rrt[N];                    //outcome
##     real so[N];                     //predictor
##     int<lower=1> I;                 //number of subjects
##     int<lower=1> K;                 //number of items
##     int<lower=1, upper=I> subj[N];  //subject id
##     int<lower=1, upper=K> item[N];  //item id
##     vector[2] mu_prior;             //vector of zeros passed in from R
## }
## parameters {
##     vector[2] beta;                 // intercept and slope
##     vector[2] u[I];                 // random intercept and slope
##     real w[K];                      // random intercept item
##     real<lower = 0> sigma_e;        // residual sd
##     vector<lower=0>[2] sigma_u;     // subj sd
##     real<lower=0> sigma_w;          // item sd
##     corr_matrix[2] Omega;           // correlation matrix for random intercepts and slopes
## }
## transformed parameters {
##     matrix[2,2] D;
##     D <- diag_matrix(sigma_u);
## }
## model {
##     matrix[2,2] L;
##     matrix[2,2] DL;
##     // priors
##     beta ~ normal(0,5);
##     sigma_e ~ cauchy(0,2);
##     sigma_u ~ cauchy(0,2);
##     sigma_w ~ cauchy(0,2);
##     Omega ~ lkj_corr(2.0);
##     L <- cholesky_decompose(Omega);
##     DL <- D * L;
##     for (i in 1:I)                             // loop for subj random effects
##         u[i] ~ multi_normal_cholesky(mu_prior, DL);
##     for (k in 1:K)                             // loop for item random effects
##         w[k] ~ normal(0,sigma_w);
##     // likelihood
##     for (n in 1:N) {
##         rrt[n] ~ normal(beta[1] + beta[2]*so[n] + u[subj[n], 1] + u[subj[n], 2]*so[n],
##                         sigma_e);
##     }
## }
## generated quantities {
##     cov_matrix[2] Sigma;
##     Sigma <- D * Omega * D;
## }



#' Create a Stan class VM model
#'
#' Internal function to create and compile a Stan model.
#'
#' @param type A character string indicating the type of model,
#'   currently only \dQuote{vm_predict} and \dQuote{vm_predict_mediation}.
#' @param useU A logical value whether the latent intercept estimated in Stage 1 should
#'   also be used as a predictor.  Defaults to \code{TRUE}.
#' @param \dots Additional arguments passed to \code{stan_model}.
#' @return A compiled Stan model.
#' @author Joshua F. Wiley <josh@@elkhartgroup.com>
#' @keywords models
#' @examples
#' # Make Me!
#' \dontrun{
#'   test1 <- vm_stan("vm_predict", useU=TRUE)
#'   test2 <- vm_stan("vm_predict", useU=FALSE)
#'   test3 <- vm_stan("vm_predict_mediation", useU=TRUE)
#'   test4 <- vm_stan("vm_predict_mediation", useU=FALSE)
#' }
vm_stan <- function(type = c("vm_predict", "vm_predict_mediation"), useU=TRUE, ...) {

  type <- match.arg(type)

  model.core.V <- list(
    data = "
      int<lower=1> n;
      int<lower=1> k;
      int VID[n];

      // Data related to V (variability)
      int<lower=1> pVX;
      real V[n];
      matrix[n, pVX] VX;
    ",
    parameters = "
      // Parameters related to V
      vector[pVX] VB;
      real U[k];
      real<lower=0> sigma_U;
      real<lower=0> shape;
      real<lower=0> rate;
      real<lower=0> Sigma_V[k];
    ",
    tparameters.declarations = "
      // Params related to V
      real V_hat[n];
      real Sigma_V_hat[n];
    ",
    tparameters.statements = "
      // Params related to V
      for (i in 1:n) {
        V_hat[i] <- (VX[i] * VB) + U[VID[i]];
        Sigma_V_hat[i] <- Sigma_V[VID[i]];
      }
    ",
    model.declarations = "
      // Priors for V Location
      VB ~ normal(0, 1000);
      U ~ normal(0, sigma_U);
      // Priors for Stage 1 scale of random location
      sigma_U ~ cauchy(0, 5);
      // Priors for Stage 1 Scale
      shape ~ cauchy(0, 5);
      rate ~ cauchy(0, 5);
      // Model for Stage 1 Scale
      Sigma_V ~ gamma(shape, rate);
    ",
    model.statements = "
      // Likelihood for V
      V ~ normal(V_hat, Sigma_V_hat);
    ")

  model.core.Y <- list(
    data = "
      // Data related to Y (outcome)
      int<lower=1> pYX;
      int<lower=1> pYX2;
      real Y[k];
      matrix[k, pYX] YX;
    ",
    parameters = "
      // Parameters related to Y
      vector[pYX] YB;
      vector[pYX2] Yalpha;
      real<lower=0> sigma_Y;
    ",
    tparameters.declarations = "
      // Params related to Y
      real Y_hat[k];
    ",
    tparameters.statements = "
      // Params related to Y
      for (i in 1:k) {
        Y_hat[i] <- (YX[i] * YB) + Yalpha[1] * Sigma_V[i]YuseU;
      }
    ",
    model.delcarations = "
      // Priors for Y location
      YB ~ normal(0, 1000);
      Yalpha ~ normal(0, 1000);
      // Priors for Y scale
      sigma_Y ~ cauchy(0, 5);
    ",
    model.statements = "
      // Likelihood for Y
      Y ~ normal(Y_hat, sigma_Y);
    ")
  model.core.Y$tparameters.statements <- gsub("YuseU", ifelse(useU, " + Yalpha[2] * U[i]", ""), model.core.Y$tparameters.statements)

  model.core.M <- list(
    data = "
      // Data related to M (mediator)
      int<lower=1> pMX;
      int<lower=1> pMX2;
      real M[k];
      matrix[k, pMX] MX;
    ",
    parameters = "
      // Parameters related to M
      vector[pMX] MB;
      vector[pMX2] Malpha;
      real<lower=0> sigma_M;
    ",
    tparameters.declarations = "
      // Params related to M
      real M_hat[k];
    ",
    tparameters.statements = "
      // Params related to M
      for (i in 1:k) {
        M_hat[i] <- (MX[i] * MB) + Malpha[1] * Sigma_V[i]MuseU;
      }
    ",
    model.declarations = "
      // Priors for M location
      MB ~ normal(0, 1000);
      Malpha ~ normal(0, 1000);
      // Priors for M scale
      sigma_M ~ cauchy(0, 5);
    ",
    model.statements = "
      // Likelihood for M
      M ~ normal(M_hat, sigma_M);
    ")
  model.core.M$tparameters.statements <- gsub("MuseU", ifelse(useU, " + Malpha[2] * U[i]", ""), model.core.M$tparameters.statements)

  model_builder <- function(...) {
    pieces <- list(...)
    n <- names(pieces[[1]])

    combined <- lapply(n, function(i) {
      do.call(paste, c(lapply(pieces, function(x) x[[i]]), list(collapse = "\n")))
    })
    names(combined) <- n

    with(combined, sprintf("
    // upper case letters indicate vectors/matrices
    // lower case letters indicate scalars
    data {
    %s
    }
    parameters {
    %s
    }
    transformed parameters {
    %s
    %s
    }
    model {
    %s
    %s
    }
    ", data, parameters,
       tparameters.declarations, tparameters.statements,
       model.declarations, model.statements))
  }

  model <- switch(type,
    vm_predict =  model_builder(model.core.V, model.core.Y),
    vm_predict_mediation = model_builder(model.core.V, model.core.Y, model.core.M)
  )

  stan_model(model_code = model, save_dso=TRUE, ...)
}

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
#' @param y.formula A formula describing a model for the outcome.  At present,
#'   this must be a continuous, normally distributed variable.
#' @param v.formula A formula describing a model for the variability. Note
#'   this must end with \code{ | ID}, where \code{ID} is the name of the
#'   ID variable in the dataset.  At present, this must be a continuous,
#'   normally distributed variable.
#' @param m.formula An optional formula decribing a model for a mediatior variable.
#'   At present, this must be a continuous normally distributed variable.
#' @param data A long data frame containing an both the Level 2 and Level 1 outcomes,
#'   as well as all covariates and an ID variable.
#' @param useU A logical value whether the latent intercept estimated in Stage 1 should
#'   also be used as a predictor.  Defaults to \code{TRUE}.  Note if there is a
#'   mediator as well as main outcome, the latent intercepts will be used as a predictor
#'   for both.
#' @param totaliter The total number of iterations to be used (not including the
#'   warmup iterations), these are distributed equally across multiple independent
#'   chains.
#' @param warmup The number of warmup iterations.  Each independent chain
#'   has the same number of warmup iterations, before it starts the iterations
#'   that will be used for inference.
#' @param chains The number of independent chains to run (default to 1).
#' @param inits Initial values passed on to \code{stan}.  If \code{NULL}, the default,
#'   initial values are estimated means, standard deviations, and coefficients from a
#'   single level linear regression.
#' @param modelfit A compiled Stan model (e.g., from a previous run).
#' @param opts A list giving options.  Currently only \code{SD_Tol} which controls
#'   the tolerance for how small a variables standard deviation may be without
#'   stopping estimation (this ensures that duplicate variables, or variables without
#'   any variability are included as predictors).
#' @param \dots Additional arguments passed to \code{stan}.
#' @return A named list containing the model \code{results}, the \code{model},
#'   the \code{variable.names}, the \code{data}, the random \code{seeds},
#'   and the initial function \code{.call}.
#' @author Joshua F. Wiley <josh@@elkhartgroup.com>
#' @import Formula
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
#'   m <- vm_predict(y2 ~ x1 + x2, y ~ 1 | ID, data = sim.data,
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
vm_predict <- function(y.formula, v.formula, m.formula, data, useU = TRUE,
                       totaliter = 2000, warmup = 1000, chains = 1,
                       inits = NULL, modelfit, opts = list(SD_Tol = .01), ...) {
  stopifnot(is.data.frame(data))
  stopifnot(!missing(y.formula))
  stopifnot(!missing(v.formula))
  stopifnot(is.logical(useU))

  storedCall <- match.call()

  # logical flag for mediation
  med <- !missing(m.formula)

  # drop any missing levels to avoid redudant dummy codes in model matrix
  data <- droplevels(data)

  var.names <- list(
    Y = all.vars(terms(as.Formula(y.formula), lhs = 1, rhs = -1)),
    V = all.vars(terms(as.Formula(v.formula), lhs = 1, rhs = -c(1, 2))),
    VID = all.vars(terms(as.Formula(v.formula), lhs = -1, rhs = 2))
  )

  if (med) {
    all.formula <- as.Formula(y.formula, m.formula, v.formula)
    var.names$M <- all.vars(terms(as.Formula(m.formula), lhs = 1, rhs = -c(1, 2)))
  } else {
    all.formula <- as.Formula(y.formula, v.formula)
  }

  # make sure ID is a numeric/integer or a factor
  stopifnot(class(data[, var.names$VID]) %in% c("numeric", "integer", "factor"))

  test.VID <- sd_id(data[, var.names$V], data[, var.names$VID], long=FALSE)

  if (!all(test.VID != 0, na.rm=TRUE)) {
    stop(sprintf("The following IDs have no variability in the first stage outcome:\n%s\nTry using\n%s\nto remove these from the data.",
      paste(names(test.VID)[which(test.VID == 0)], collapse = ', '),
      paste0("subset(your_data, sd_id(", var.names$V, ", ", var.names$VID, ") != 0)")))
  }

  key <- list(OriginalID = data[, var.names$VID])

  data[, var.names$VID] <- as.integer(data[, var.names$VID])

  key$IntegerID <- data[, var.names$VID]

  data <- data[order(data[, var.names$VID]), ]

  mf <- model.frame(all.formula, data = data, na.action = na.omit)

  key$MFOriginalID <- mf[, var.names$VID]

  mf[, var.names$VID] <- as.integer(factor(mf[, var.names$VID]))

  key$MFIntegerID <- mf[, var.names$VID]

  key <- with(key, {
    tmp1 <- data.frame(OriginalID, IntegerID)[!duplicated(IntegerID), ]
    tmp2 <- data.frame(MFOriginalID, MFIntegerID)[!duplicated(MFIntegerID), ]
    data.frame(OriginalID = tmp1[match(tmp2[, 1], tmp1[, 2]), 1],
               InternalID = tmp2[, 2])
  })

  y.mm <- model.matrix(y.formula, data = mf)
  v.mm <- model.matrix(as.Formula(v.formula), data = mf, rhs = 1)
  if (med) {
    m.mm <- model.matrix(m.formula, data = mf)
  }

  V <- mf[, var.names$V]
  VX <- v.mm
  VID <- mf[, var.names$VID]

  keep.obs <- !duplicated(VID)

  Y <- mf[keep.obs, var.names$Y]
  YX <- y.mm[keep.obs, , drop=FALSE]

  var.names$YX <- colnames(y.mm)
  var.names$VX <- colnames(v.mm)

  if (med) {
    M <- mf[keep.obs, var.names$M]
    MX <- m.mm[keep.obs, , drop = FALSE]
    var.names$MX <- colnames(m.mm)
  }

  # Code to check the scaling of variables
  if (med) {
    v.sds <- c(sd(Y, na.rm=TRUE), sd(V, na.rm=TRUE), sd(M, na.rm=TRUE),
             apply(YX, 2, sd, na.rm=TRUE), apply(VX, 2, sd, na.rm=TRUE), apply(MX, 2, sd, na.rm=TRUE))
    names(v.sds) <- c(var.names$Y, var.names$V, var.names$M,
                      var.names$YX, var.names$VX, var.names$MX)
  } else {
    v.sds <- c(sd(Y, na.rm=TRUE), sd(V, na.rm=TRUE),
             apply(YX, 2, sd, na.rm=TRUE), apply(VX, 2, sd, na.rm=TRUE))
    names(v.sds) <- c(var.names$Y, var.names$V,
                      var.names$YX, var.names$VX)
  }

  v.sds <- v.sds[names(v.sds) != "(Intercept)"]

  v.sds.index <- is.na(v.sds) | v.sds < opts$SD_Tol | v.sds > 50
  if (any(v.sds.index)) {
    stop(sprintf("The follow variables SDs are either too small or too large.\n  Remove or rescale variables before modelling.\n  Variables: %s",
                 paste(names(v.sds)[v.sds.index], collapse = ", ")))
  }

  class.tests <- c(
    V = is.numeric(V) & is.vector(V),
    VX = (is.matrix(VX) & is.numeric(VX)),
    Y = is.numeric(Y) & is.vector(Y),
    YX = is.matrix(YX) & is.numeric(YX))

  if (med) {
    class.tests <- c(class.tests,
      M = is.numeric(M) & is.vector(M),
      MX = is.matrix(MX) & is.numeric(MX))
  } else {
    class.tests <- c(class.tests, M = TRUE, MX = TRUE)
  }

  if (!all(class.tests)) {
    stop(c("V must be a numeric vector", "VX must be a numeric matrix",
           "Y must be a numeric vector", "YX must be a numeric matrix",
           "M must be a numeric vector", "MX must be a numeric matrix")[!class.tests])
  }

  n <- length(V)
  k <- length(unique(VID))

  if (!all(identical(nrow(VX), n), identical(length(VID), n))) {
    stop("The length and rows of VX, V, and VID must all be equal")
  }

  if (!all(identical(nrow(YX), k), identical(length(Y), k))) {
    stop("The length and rows of Y, YX, and the unique IDs must all be equal")
  }

  if (med) {
    if (!all(identical(nrow(MX), k), identical(length(M), k))) {
      stop("The length and rows of M, MX, and the unique IDs must all be equal")
    }
  }

  pVX <- ncol(VX)
  pYX <- ncol(YX)
  pYX2 <- 1L + useU

  if (med) {
    pMX <- ncol(MX)
    pMX2 <- 1L + useU

    stan.data <- list(
      V = V, VX = VX, VID = VID,
      Y = Y, YX = YX,
      M = M, MX = MX,
      n = n, k = k,
      pVX = pVX, pYX = pYX, pYX2 = pYX2, pMX = pMX, pMX2 = pMX2)

  } else {
    stan.data <- list(
      V = V, VX = VX, VID = VID,
      Y = Y, YX = YX,
      n = n, k = k,
      pVX = pVX, pYX = pYX, pYX2 = pYX2)
  }


  stan_inits <- function() {
    rg <- res_gamma(V, VID)

    out <- with(stan.data, list(
      VB = as.array(coef(lm.fit(VX, V))),
      U = as.array(by_id(V, VID, mean, FALSE, na.rm=TRUE)-mean(V, na.rm=TRUE)),
      shape = rg$alpha,
      rate = rg$beta,
      Sigma_V = as.array(sd_id(V, VID, FALSE)),
      sigma_Y = sd(Y, na.rm=TRUE)
    ))

    out$sigma_U <- sd(out$U, na.rm=TRUE)

    if (useU) {
      yb2 <- with(stan.data, coef(lm.fit(cbind(YX, Res = out$Sigma_V, U = out$U), Y)))
    } else {
      yb2 <- with(stan.data, coef(lm.fit(cbind(YX, Res = out$Sigma_V), Y)))
    }

    ync <- ncol(stan.data$YX)

    out$YB <- as.array(yb2[1:ync])
    out$Yalpha <- as.array(yb2[(ync + 1):length(yb2)])

    if (med) {
      out$sigma_M <- sd(stan.data$M, na.rm=TRUE)
      if (useU) {
        mb2 <- with(stan.data, coef(lm.fit(cbind(MX, Res = out$Sigma_V, U = out$U), M)))
      } else {
        mb2 <- with(stan.data, coef(lm.fit(cbind(MX, Res = out$Sigma_V), M)))
      }

      mnc <- ncol(stan.data$MX)

      out$MB <- as.array(mb2[1:mnc])
      out$Malpha <- as.array(mb2[(mnc + 1):length(mb2)])
    }

    return(out)
  }

  if (is.null(inits)) inits <- stan_inits

  if (!missing(modelfit)) {
    model <- modelfit
  } else {
    if (med) {
      model <- vm_stan("vm_predict_mediation", useU=useU)
    } else {
      model <- vm_stan("vm_predict", useU=useU)
    }
  }

  if (med) {
    res <- parallel_stan(modelfit = model, standata = stan.data,
      totaliter = totaliter, warmup = warmup, chains = chains,
      pars = c("VB", "U", "sigma_U", "shape", "rate", "Sigma_V",
               "YB", "Yalpha", "sigma_Y",
               "MB", "Malpha", "sigma_M"), init = inits, ...)
  } else {
    res <- parallel_stan(modelfit = model, standata = stan.data,
      totaliter = totaliter, warmup = warmup, chains = chains,
      pars = c("VB", "U", "sigma_U", "shape", "rate", "Sigma_V",
               "YB", "Yalpha", "sigma_Y"), init = inits, ...)
  }

  out <- list(
    results = res$results,
    model = res$model,
    variable.names = var.names,
    data = c(stan.data, list(IDkey = key)),
    seeds = res$seeds,
    .call = storedCall
  )

  if (med) {
    class(out) <- c("vmp_med", "list")
  } else {
    class(out) <- c("vmp", "list")
  }

  return(out)
}
