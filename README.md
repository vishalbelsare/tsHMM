
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tsHMM

<!-- badges: start -->
<!-- badges: end -->

The goal of tsHMM is to …

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Akai01/tsHMM")
```

## Example

This is a basic example which shows you how to solve a common problem:

We will use `earthquakes` data and try to forecast probability of number
of earthquakes in a year.

``` r
library(tsHMM)

data("earthquakes", package = "tsHMM")

m = 2 # two states

h = 4 # years ahead forecast

lambda <- c(15,25)

gamma <- matrix(c(0.9,0.1,0.1,0.9), m, m,byrow=TRUE)

model <- tsHMM$new(x = earthquakes, h = h, m = 2, lambda = lambda, gamma = gamma)

model$fit()

forecasts <- model$forecast(xf = 0:50)

model$plot(h = 1:4)
```

<img src="man/figures/README-example-1.png" width="100%" />
