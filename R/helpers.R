#' @name Variability_Measures
#' @title Variability Measures
#' @rdname VarMeasures
#' @aliases by_id
#' @aliases sd_id
#' @aliases rmssd
#' @aliases rmssid_id
#' @aliases rolling_diff
#' @aliases rolling_diff_id
#'
#' @note These are a set of functions designed to calculate various
#' measures of variability either on a single data vector, or
#' calculate them by an ID.
#'
#' @param x A data vector to operate on.  Should be a numeric or
#'   integer vector, or coercible to such (e.g., logical).
#' @param ID an ID variable indicating how to split up the \code{x}
#'   vector.  Should be the same length as \code{x}.
#' @param fun The function to calculate by ID
#' @param long A logical indicating whether to return results in
#'   \dQuote{long} form (the default) or wide (if \code{FALSE}).
#' @param \dots Additional arguments passed on to \code{fun}
#' @author Joshua F. Wiley <josh@@elkhartgroup.com>
NULL


#' Variability Measures
#'
#' \code{by_id} - Internal function to allow a simple statistic (e.g., SD)
#' to be calculated individually by an ID variable and returned
#' either as per ID (i.e., wide form) or for every observation of an
#' ID (i.e., long form).
#' @return \code{by_id} - A vector the same length as \code{x}
#'   if \code{long=TRUE}, or the length of unique \code{ID}s if
#'   \code{long=FALSE}.
#' @rdname VarMeasures
by_id <- function(x, ID, fun, long=TRUE, ...) {
  if (long) {
    ave(x, ID, FUN = function(x) fun(x, ...))
  } else {
    tapply(x, ID, FUN = function(x) fun(x, ...))
  }
}

#' Variability Measures
#'
#' \code{sd_id} - Calculates the standard deviation of observations by \code{ID}.
#'
#' @return \code{sd_id} - A vector of the standard deviations by ID
#' @keywords utilities
#' @export
#' @rdname VarMeasures
#' @examples
#'
#' sd_id(mtcars$mpg, mtcars$cyl, long=TRUE)
#' sd_id(mtcars$mpg, mtcars$cyl, long=FALSE)
sd_id <- function(x, ID, long=TRUE) {
  by_id(x, ID, fun = sd, long = long, na.rm=TRUE)
}

#' Variability Measures
#'
#' \code{rmssd} - Calculates the root mean square of successive differences (RMSSD).
#'   Note that missing values are removed.
#'
#' @return \code{rmssd} - The RMSSD for the data.
#' @export
#' @rdname VarMeasures
#' @examples
#' rmssd(1:4)
#' rmssd(c(1, 3, 2, 4))
rmssd <- function(x) {
  x <- na.omit(diff(x))
  as.vector(sqrt(mean(x^2)))
}

#' Variability Measures
#'
#' \code{rmssd_id} - Calculates the RMSSD by ID.
#'
#' @return \code{rmssd_id} - A vector of the RMSSDs by ID
#' @export
#' @rdname VarMeasures
#' @examples
#' rmssd_id(mtcars$mpg, mtcars$cyl)
#' rmssd_id(mtcars$mpg, mtcars$cyl, long=FALSE)
rmssd_id <- function(x, ID, long=TRUE) {
  by_id(x, ID, fun = rmssd, long = long)
}

#' Variability Measures
#'
#' \code{rolling_diff} - Calculates the average rolling difference of the data.
#'   Within each window, the difference between the maximum and minimum value is
#'   computed and these are averaged across all windows.  The equation is:
#' \deqn{\frac{\sum_{t = 1}^{N - k} max(x_{t}, \ldots, x_{t + k}) - min(x_{t}, \ldots, x_{t + k})}{N - k}}
#'
#' @param window An integer indicating the size of the rolling window.
#'   Must be at least the length of \code{x}.
#' @return \code{rolling_diff} - The average of the rolling differences between maximum and minimum.
#' @export
#' @rdname VarMeasures
#' @examples
#' rolling_diff(1:7, window = 4)
#' rolling_diff(c(1, 4, 3, 4, 5))
rolling_diff <- function(x, window = 4) {
  stopifnot(length(x) >= window)

  index <- 1:(length(x) + 1 - window)

  mean(sapply(index, function(i) {
    x <- na.omit(x[i:(i + window - 1)])
    if (length(x) < 2) {
      NA
    } else {
      diff(range(x))
    }
  }), na.rm=TRUE)
}

#' Variability Measures
#'
#' \code{rolling_diff_id} - Calculates the average rolling difference by ID
#'
#' @return \code{rolling_diff_id} - A vector of the average rolling differences by ID
#' @export
#' @rdname VarMeasures
#' @examples
#' rolling_diff_id(mtcars$mpg, mtcars$cyl, window = 3)
rolling_diff_id <- function(x, ID, long=TRUE, window = 4) {
  by_id(x, ID, fun = rolling_diff, long = long, window = window)
}


#' Estimate the parameters for a Gamma distribution
#'
#' This is a simple function to estimate what the parameters for a Gamma
#' distribution would be from a data vector.  It is used internally to
#' generate start values.
#'
#' @param x a data vector to operate on
#' @return a list of the shape (alpha) and rate (beta) parameters
#'   and the mean and variance
#' @author Joshua F. Wiley <josh@@elkhartgroup.com>
#' @keywords utilities
gamma_params <- function(x) {
  m <- mean(x, na.rm=TRUE)
  v <- var(x, na.rm=TRUE)
  beta <- m / v
  alpha <- m * beta
  list(alpha = alpha, beta = beta, mean = m, variance = v)
}

#' Estimates the parameters of a Gamma distribution from SDs
#'
#' This function calcualtes the parameters of a Gamma distribution
#' from the residuals from an individuals' own mean.
#' That is, the distribution of (standard) deviations from individuals'
#' own mean are calculated and then an estimate of the parameters of a
#' Gamma distribution are calculated.
#'
#' @param x A data vector to operate on
#' @param ID an ID variable of the same length as \code{x}
#' @return a list of the shape (alpha) and rate (beta) parameters
#'   and the mean and variance
#' @author Joshua F. Wiley <josh@@elkhartgroup.com>
#' @export
#' @keywords utilities
#' @examples
#'
#' set.seed(1234)
#' y <- rgamma(100, 3, 2)
#' x <- rnorm(100 * 10, mean = 0, sd = rep(y, each = 10))
#' ID <- rep(1:100, each = 10)
#' res_gamma(x, ID)
res_gamma <- function(x, ID) {
  gamma_params(sd_id(x, ID, long = FALSE))
}

#' Calculates an empirical p-value based on the data
#'
#' This function takes a vector of statistics and calculates
#' the empirical p-value, that is, how many fall on the other
#' side of zero.  It calculates a two-tailed p-value.
#'
#' @param x a data vector to operate on
#' @param na.rm Logical whether to remove NA values. Defaults to \code{TRUE}
#' @return a named vector with the number of values falling at
#'   or below zero, above zero, and the empirical p-value.
#' @author Joshua F. Wiley <josh@@elkhartgroup.com>
#' @export
#' @keywords utilities
#' @examples
#'
#' empirical_pvalue(rnorm(100))
empirical_pvalue <- function(x, na.rm = TRUE) {
  x <- as.integer(x <= 0)
  tmp <- table(factor(x, levels = 1:0, labels = c("<= 0", "> 0")))
  m <- mean(x, na.rm = na.rm)
  pval2 <- 2 * min(m, 1 - m)
  out <- c(as.vector(tmp), pval2)
  names(out) <- c(names(tmp), "p-value")

  out
}


#' nice formatting for p-values
#'
#' @param p a numeric pvalue
#' @param d the digits less than which should be displayed as less than
#' @param sd scientific digits for round
#' @author Joshua F. Wiley <josh@@elkhartgroup.com>
#' @keywords utilities
#' @examples
#' varian:::pval_smartformat(c(1, .15346, .085463, .05673, .04837, .015353462,
#'   .0089, .00164, .0006589, .0000000053326), 3, 5)
pval_smartformat <- function(p, d = 3, sd = 5) {
  p.out <- ifelse(p < 1/(10^d),
         paste0("< .", paste(rep(0, d - 1), collapse = ""), "1"),
                  format(round(p, digits = d), digits = d, nsmall = d, scientific = sd))
  gsub("0\\.", ".", p.out)
}


#' Calculates summaries for a parameter
#'
#' This function takes a vector of statistics and calculates
#' several summaries: mean, median, 95% CI, and
#' the empirical p-value, that is, how many fall on the other
#' side of zero.
#'
#' @param x a data vector to operate on
#' @param \dots Additional arguments passed to \code{pval_smartformat}
#'   to control p-value printing.
#' @param na.rm Logical whether to remove NA values. Defaults to \code{TRUE}
#' @return A data frame of summary statistics
#' @author Joshua F. Wiley <josh@@elkhartgroup.com>
#' @export
#' @keywords utilities
#' @examples
#'
#' param_summary(rnorm(100))
param_summary <- function(x, ..., na.rm = TRUE) {
  data.frame(
    Mean = mean(x, na.rm = na.rm),
    Median = median(x, na.rm = na.rm),
    SE = sd(x, na.rm = na.rm),
    LL2.5 = as.vector(quantile(x, probs = .025, na.rm = na.rm)),
    UL97.5 = as.vector(quantile(x, probs = .975, na.rm = na.rm)),
    pvalue = pval_smartformat(empirical_pvalue(x)[["p-value"]], ...))
}

#' Format a data frame of summary statistics
#'
#' This functions nicely formats a data frame of parameter summary
#' statistics and is designed to be used with the param_summary()
#' function.
#'
#' @param d A data frame of the parameter summary statistics
#' @param digits Number of digits to round to for printing
#' @param pretty Logical value whether prettified values should be returned.
#'   Defaults to \code{FALSE}.
#' @return A formatted data frame of summary statistics or a formated
#' vector (if \code{pretty = TRUE}).
#' @author Joshua F. Wiley <josh@@elkhartgroup.com>
#' @export
#' @keywords utilities
#' @examples
#' set.seed(1234)
#' xsum <- do.call(rbind, apply(matrix(rnorm(100*10), ncol = 10),
#'   2, param_summary))
#' rownames(xsum) <- letters[1:10]
#' param_summary_format(xsum)
#' param_summary_format(xsum, pretty = TRUE)
#'
#' rm(xsum)
param_summary_format <- function(d, digits = getOption("digits"), pretty = FALSE) {
  if (pretty) {
    string <- sprintf("%%01.%df [%%01.%df, %%01.%df], %%s", digits, digits, digits)

    out <- sprintf(string,
                   d$Mean,
                   d$LL2.5,
                   d$UL97.5,
                   ifelse(grepl("<", d$pvalue),
                          paste0("p ", d$pvalue),
                          paste0("p = ", d$pvalue)))
    names(out) <- rownames(d)
  } else {
    out <- as.data.frame(lapply(d, function(x) {
      if (is.numeric(x)) {
        sprintf(sprintf("%%01.%df", digits), x)
      } else {
        x
      }
    }), stringsAsFactors = FALSE)
    rownames(out) <- rownames(d)
    }

  return(out)
}

#' Summary method for variability model objects
#'
#' @param object An object of class \code{vm}.
#' @param digits The number of digits to use for rounding.
#' @param \dots Additional arguments. Not currently used.
#' @export
#' @method summary vm
summary.vm <- function(object, digits = getOption("digits"), ...) {
  tmp <- extract(object$results, permute = TRUE)

  out <- with(tmp, {
       vbsum <- do.call(rbind.data.frame, apply(VB, 2, param_summary))
       rownames(vbsum) <- paste0("V<-", object$variable.names[[4]])

       sigmaUsum <- param_summary(sigma_U)
       rownames(sigmaUsum) <- "U<->U (residual sd)"
       shapesum <- param_summary(shape)
       rownames(shapesum) <- "Gamma (shape)"
       ratesum <- param_summary(rate)
       rownames(ratesum) <- "Gamma (rate)"

       if (grepl("Y", object$design)) {
         ybsum <- do.call(rbind.data.frame, apply(YB, 2, param_summary))
         rownames(ybsum) <- paste0("Y<-", object$variable.names$YX)

         yalphasum <- do.call(rbind.data.frame, apply(Yalpha, 2, param_summary))
         rownames(yalphasum) <- paste0("Y<-", c("sigma", "mu", "sigma2", "mu2"))[c(TRUE, object$useU, object$IIVQ, object$UQ)]

         sigmaYsum <- param_summary(sigma_Y)
         rownames(sigmaYsum) <- "Y<->Y (residual sd)"
       }

       if (grepl("M", object$design)) {
         mbsum <- do.call(rbind.data.frame, apply(MB, 2, param_summary))
         rownames(mbsum) <- paste0("M<-", object$variable.names$MX)

         malphasum <- do.call(rbind.data.frame, apply(Malpha, 2, param_summary))
         rownames(malphasum) <- paste0("M<-", c("sigma", "mu", "sigma2", "mu2"))[c(TRUE, object$useU, object$IIVQ, object$UQ)]

         sigmaMsum <- param_summary(sigma_M)
         rownames(sigmaMsum) <- "M<->M (residual sd)"
       }

       if (grepl("Y", object$design) & grepl("M", object$design)) {
         do.call(rbind, list(ybsum, yalphasum, sigmaYsum,
                             mbsum, malphasum, sigmaMsum,
                             vbsum, sigmaUsum, shapesum, ratesum))
       } else if (grepl("Y", object$design)) {
         do.call(rbind, list(ybsum, yalphasum, sigmaYsum,
                             vbsum, sigmaUsum, shapesum, ratesum))
       } else {
         do.call(rbind, list(vbsum, sigmaUsum, shapesum, ratesum))
       }
     })

  param_summary_format(out, digits = digits, ...)
}

#' Simulate a Gamma Variability Model
#'
#' This function facilitates simulation of a Gamma Variability Model
#' and allows the number of units and repeated measures to be varied
#' as well as the degree of variability.
#'
#' @param n The number of repeated measures on each unit
#' @param k The number of units
#' @param mu The grand mean of the variable
#' @param mu.sigma The standard deviation of the random mean of the variable
#' @param sigma.shape the shape (alpha) parameter of the Gamma distribution
#'   controlling the residual variability
#' @param sigma.rate the rate (beta) parameter of the Gamma distribution
#'   controlling the residual variability
#' @param seed the random seed, used to make simulations reproductible.
#'   Defaults to 5346 (arbitrarily).
#' @return a list of the data, IDs, and the parameters used for the simulation
#' @author Joshua F. Wiley <josh@@elkhartgroup.com>
#' @export
#' @import MASS
#' @keywords utilities
#' @examples
#' raw.sim <- simulate_gvm(12, 140, 0, 1, 4, .1, 94367)
#' sim.data <- with(raw.sim, {
#'   set.seed(265393)
#'   x2 <- MASS::mvrnorm(k, c(0, 0), matrix(c(1, .3, .3, 1), 2))
#'   y2 <- rnorm(k, cbind(Int = 1, x2) %*% matrix(c(3, .5, .7)) + sigma, sd = 3)
#'   data.frame(
#'     y = Data$y,
#'     y2 = y2[Data$ID2],
#'     x1 = x2[Data$ID2, 1],
#'     x2 = x2[Data$ID2, 2],
#'     ID = Data$ID2)
#' })
simulate_gvm <- function(n, k, mu, mu.sigma, sigma.shape, sigma.rate, seed = 5346) {
  set.seed(seed)
  m <- rnorm(k, mu, mu.sigma)
  sigma <- rgamma(k, sigma.shape, sigma.rate)
  y <- rnorm(n * k, mean = rep(m, each = n),
    sd = rep(sigma, each = n))
  y2 <- rnorm(k, sigma/sqrt(sigma.shape/(sigma.rate^2)), sd = 3)

  list(
    Data = data.frame(y = y, y2 = y2,
    ID1 = 1:(n * k), ID2 = rep(1:k, each = n)),
    n = n, k = k, mu = mu, mu.sigma,
    sigma.shape, sigma.rate, sigma = sigma, seed)
}
