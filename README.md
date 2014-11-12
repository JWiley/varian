<img src="https://elkhartgroup.com/images/varian_logo.png"
	 alt="varian logo"
	 style = "width: 650px; height: 160px;"/>

`varian` is a free, open source `R` package for *vari*ability
*an*alysis using Bayesian inference.

It provides a formula interface in `R` to define a multilevel model
for a variable with repeated measures, in addition to standard
multilevel features such as fixed effects and a random
intercepts/slopes, the residual variance is a random variable by
ID.  The latent residual variance estimates are then used as
predictors in a regression of another outcome.

# Getting Started

At its core, `varian` uses [`Stan`](http://mc-stan.org) to estimate the
models.  Because the `rstan` package is not yet on CRAN, `varian` must
be installed from github.

## Prerequisites

### R

`R` version 3.1.1 or later is required.  You can download the latest
version of `R` here:

<http://www.r-project.org/>

### Tools

`Stan` and the `varian` package rely on being able to compile `C++`
code.

* For Windows, download and install Rtools at
  <http://cran.r-project.org/bin/windows/Rtools/>
* For Mac, make sure you have the latest version of Xcode installed.
  The installation instructions for `rstan` provide more details.
* For Linux, make sure you have a recent version of ge++ or clang++.

### RStan

Install the latest version of the
[`rstan`](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)
package by following the directions from the creators of `Stan`, and
ensure that `Stan` is working correctly.

That page also lists more detailed directions for getting the
necessary tools installed.

## Install varian

You can install:

* the latest stable release (0.2.0) from github

    ```R
	install.packages("devtools")
	devtools::install_github("varian", "ElkhartGroup", ref = "v0.2.0")
	```

* the latest development version from github

    ```R
	install.packages("devtools")
	devtools::install_github("varian", "ElkhartGroup")
	```

## Learn about the theory

See our open access pre-print on the arXiv <http://arxiv.org/abs/1411.2961>.


# More

* If you find a bug or have a feature request, please let us know on
  the issue tracker: <https://github.com/ElkhartGroup/varian/issues>.

* Have questions or need support?  Email us at info@elkhartgroup.com
  and mention "varian" in the subject.
