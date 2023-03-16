# varian 0.3.0

## Major Changes
* `parallel_stan()` has been removed as the `rstan` sampling function
  can run in parallel now.  To use multiple cores now, follow
  the `rstan` approach of:
  `rstan_options(auto_write = TRUE)`
  `options(mc.cores = 4)` if you wanted 4 cores, for example, 

* `varian()` now only requires a single seed to be set, as this is
  now controlled by `rstan` rather than the removed `parallel_stan()`
  function.

## New Features
* `varian()` can now include **quadratic** effects of latent means and
  intraindividual variabilities using the new arguments, `UQ = TRUE`
  and `IIVQ = TRUE`.

* `summary.vm()` method now added for a convenient summary.

* `shinystan` package added as a suggested package. This implements
  interactive and high quality model diagnostics. This will likely
  replace the `vm_diagnostics()` function in the near future.

# varian 0.2.0

## Major Changes

* `vm_predict()` renamed to `varian()` reflecting a unification of
  separate functions into a more general purpose, variability analysis
  function.

## New Features

* `varian()` now allows different model `design`s including "V" to
  estimate intra-individual variability alone (without using it as a
  predictor) and "V -> M -> Y" to estimate a simple mediation model.

* Many back end changes including more pre-modeling data checks and
  better estimates for start values.

# varian 0.1.0

## Functions

* `vm_predict()` calculates the intraindividual variability and uses
  this to predict an outcome
