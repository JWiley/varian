% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{simulate_gvm}
\alias{simulate_gvm}
\title{Simulate a Gamma Variability Model}
\usage{
simulate_gvm(n, k, mu, mu.sigma, sigma.shape, sigma.rate, seed = 5346)
}
\arguments{
\item{n}{The number of repeated measures on each unit}

\item{k}{The number of units}

\item{mu}{The grand mean of the variable}

\item{mu.sigma}{The standard deviation of the random mean of the variable}

\item{sigma.shape}{the shape (alpha) parameter of the Gamma distribution
controlling the residual variability}

\item{sigma.rate}{the rate (beta) parameter of the Gamma distribution
controlling the residual variability}

\item{seed}{the random seed, used to make simulations reproductible.
Defaults to 5346 (arbitrarily).}
}
\value{
a list of the data, IDs, and the parameters used for the simulation
}
\description{
This function facilitates simulation of a Gamma Variability Model
and allows the number of units and repeated measures to be varied
as well as the degree of variability.
}
\examples{
raw.sim <- simulate_gvm(12, 140, 0, 1, 4, .1, 94367)
sim.data <- with(raw.sim, {
  set.seed(265393)
  x2 <- MASS::mvrnorm(k, c(0, 0), matrix(c(1, .3, .3, 1), 2))
  y2 <- rnorm(k, cbind(Int = 1, x2) \%*\% matrix(c(3, .5, .7)) + sigma, sd = 3)
  data.frame(
    y = Data$y,
    y2 = y2[Data$ID2],
    x1 = x2[Data$ID2, 1],
    x2 = x2[Data$ID2, 2],
    ID = Data$ID2)
})
}
\author{
Joshua F. Wiley <jwiley.psych@gmail.com>
}
\keyword{utilities}
