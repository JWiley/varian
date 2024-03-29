% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modeler.R
\name{stan_inits}
\alias{stan_inits}
\title{Calculate Initial Values for Stan VM Model}
\usage{
stan_inits(
  stan.data,
  design = c("V -> Y", "V -> M -> Y", "V", "X -> V", "X -> V -> Y", "X -> M -> V"),
  useU = TRUE,
  UQ = FALSE,
  IIVQ = FALSE,
  centerU = 0,
  centerIIV = 0,
  ...
)
}
\arguments{
\item{stan.data}{A list containing the data to be passed to Stan}

\item{design}{A character string indicating the type of model to be run.  One of
\dQuote{V -> Y} for variability predicting an outcome,
\dQuote{V -> M -> Y} for mediation of variability on an outcome,
\dQuote{V} to take posterior samples of individual variability estimates alone.}

\item{useU}{A logical value whether the latent intercept estimated in Stage 1 should
also be used as a predictor.  Defaults to \code{TRUE}.  Note if there is a
mediator as well as main outcome, the latent intercepts will be used as a predictor
for both.}

\item{UQ}{A logical value whether the latent intercept estimated in Stage 1 should
also be used as a predictor with a qudratice effect.  Defaults to \code{FALSE}.
Note if there is a mediator as well as main outcome, the latent intercepts
will be used as a predictor for both.}

\item{IIVQ}{A logical value whether the latent variabilities estimated in Stage 1 should
also be used as a predictor with a qudratice effect.  Defaults to \code{FALSE}.
Note if there is a mediator as well as main outcome, the latent intercepts
will be used as a predictor for both.}

\item{centerU}{A numeric vector of length one (scalar) that is used to center the
latent intercept estimates before using as a predictor of the outcome / mediator.
Uses the formula: (U - centerU). Particularly useful when including quadratic terms.}

\item{centerIIV}{A numeric vector of length one (scalar) that is used to center the
latent variability estimates before using as a predictor of the outcome / mediator.
Uses the formula: (IIV - centerIIV). Particularly useful when including quadratic terms.}

\item{...}{Additional arguments passed to \code{sampling}.}
}
\value{
A named list containing the initial values for Stan.
}
\description{
Internal function used to get rough starting values for a
variability model in Stan.  Uses inidivudal standard deviations, means,
and linear regressions.
}
\examples{
# make me!
}
\author{
Joshua F. Wiley <josh@elkhartgroup.com>
}
\keyword{models}
