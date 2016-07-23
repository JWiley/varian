# varian 0.3.0

## New Features
* `varian()` can now include quadratic effects of latent means and
  intraindividual variabilities.

* `summary.vm()` method now added for a convenient summary.

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
