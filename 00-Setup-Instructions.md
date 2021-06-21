Installation Instructions
================
Eric Novik, Generable Inc.
18 Jun 2021

  - [Introduction](#introduction)
  - [Installing RStan](#installing-rstan)
  - [Installing RStanArm Survival
    Branch](#installing-rstanarm-survival-branch)
  - [Getting Help and Other Resouces](#getting-help-and-other-resouces)
  - [References](#references)

## Introduction

This markdown document explains how to install required software for the
Bayesian Survival and Joint Models using Rstanarm for the [NCB 2021
confernece](https://community.amstat.org/biop/events/ncb/index) class
taught by [Jacqueline
Buros](https://scholar.google.com/citations?user=O_AWJTwAAAAJ&hl=en).

Please, exexecute all steps and verify that the software works before
the class. We will not have time to troubleshoot installation problems
during the class but if you experience any problems please send us an
email to [**ncb2021@generable.com**](mailto:ncb2021@generable.com).

[**Stan**](https://mc-stan.org/) (Carpenter et al. 2017) is a
programming language with its own compiler and in this class, we will be
accessing Stan from R directly via the
[RStan](https://mc-stan.org/users/interfaces/rstan) package and
indirectly from an
[RStanArm](https://mc-stan.org/users/interfaces/rstanarm) package which
provides a set of pre-compiled Stan models including Survival models.

## Installing RStan

For this class, you should be using R version 4.0 or later and we
recommend using RStudio version 1.4 or later. To install RStan, execute
the following commands.

``` r
install.packages("rstan", dependencies = TRUE)
```

On most systems, this should work without a problem, but if you are
experiencing installation issues, take a look a the [RStan Getting
Started](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)
document.

Execute the following code to test your installation.

``` r
# In general, you should write your Stan code in a separate file
# but this will work for testing:
stanmodelcode <- "
data {
  int<lower=0> N;
  real y[N];
} 

parameters {
  real mu;
} 

model {
  target += normal_lpdf(mu | 0, 10);
  target += normal_lpdf(y  | mu, 1);
} 
"
library(rstan)
y <- rnorm(20) 
dat <- list(N = 20, y = y); 

# This will compile the model and start the sampling process
# It may take 10-30 seconds depending on your machine
fit <- stan(model_code = stanmodelcode, 
            model_name = "example", 
            data = dat)

print(fit)
```

The last command should produce something like the following.

    Inference for Stan model: example.
    4 chains, each with iter=2000; warmup=1000; thin=1; 
    post-warmup draws per chain=1000, total post-warmup draws=4000.
    
           mean se_mean   sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
    mu    -0.08    0.01 0.23  -0.54  -0.24  -0.07   0.08   0.38  1236    1
    lp__ -30.48    0.02 0.78 -32.72 -30.63 -30.18 -29.99 -29.94  1246    1

In addition, install the following R packages.

``` r
p <- c('survival', 'bayesplot', 'RISCA', 'tidyverse', 'remotes')
install.packages(p)

# helper functions
remotes::install_github('generable/rgeco')
```

## Installing RStanArm Survival Branch

The other major package we will use is RStanArm. The advantage of
RStanArm is that it comes with a large set of pre-compiled models that
have been tested for computational efficiency. As of this writing, the
Survival functions are not available on CRAN and so you have to install
the package from a GitHub branch as follows.

``` r
# This will take a while -- more than a few minutes as lots of Stan
# models are being compiled. Go have a cup of coffee if that's your thing.
remotes::install_github("stan-dev/rstanarm", ref = "feature/survival", build_vignettes = FALSE)
```

To test if the installation worked execute the following code.

``` r
library(rstanarm)
library(survival)
library(dplyr)

fit <- stan_surv(formula = Surv(time, status) ~ x,
                 data = leukemia, # from Survival package
                 basehaz = "exp")

print(fit)
```

After executing the print command, you should see something like this.

    stan_surv
     baseline hazard: exponential
     formula:         Surv(time, status) ~ x
     observations:    23
     events:          18 (78.3%)
     right censored:  5 (21.7%)
     delayed entry:   no
    ------
                   Median MAD_SD exp(Median)
    (Intercept)    -4.1    0.4     NA       
    xNonmaintained  0.9    0.5    2.6       
    
    ------
    * For help interpreting the printed output see ?print.stanreg
    * For info on the priors used see ?prior_summary.stanreg

## Getting Help and Other Resouces

For informaiton on the Bayesian Survival analysis see the paper
[*Bayesian Survival Analysis Using the rstanarm
R*](https://arxiv.org/abs/2002.09633) (Brilleman et al. 2020). For a
description of the Joint Models in Bayesian context see the *Joint
longitudinal and time-to-event models for multilevel hierarchical data*
paper (Brilleman et al. 2019).

For the introduction to the RStanArm pacakge see [*How to Use the
rstanarm Package*](https://mc-stan.org/rstanarm/articles/rstanarm.html)
vignette. For more information on Joint Models in RStanArm see
[this](https://cran.r-project.org/web/packages/rstanarm/vignettes/jm.html)
and [this](http://mc-stan.org/rstanarm/reference/stan_jm.html).

Stan is extrememy well documented. Stan User Manual and other Stan
documentation are available
[here](https://mc-stan.org/users/documentation/).

If you have a question about Stan inlcluding the language and models,
check out a user-friendly [community
forum](https://discourse.mc-stan.org/).

We hope you engoy the workshop. If there anything we can do to improve,
please send us feebback to
[**ncb2021@generable.com**](mailto:ncb2021@generable.com)**.**

## References

<div id="refs" class="references hanging-indent">

<div id="ref-brilleman2019">

Brilleman, Samuel L, Michael J Crowther, Margarita Moreno-Betancur,
Jacqueline Buros Novik, James Dunyak, Nidal Al-Huniti, Robert Fox, Jeff
Hammerbacher, and Rory Wolfe. 2019. “Joint Longitudinal and
Time-to-Event Models for Multilevel Hierarchical Data.” *Statistical
Methods in Medical Research* 28 (12): 3502–15.
<https://doi.org/10.1177/0962280218808821>.

</div>

<div id="ref-brilleman2020">

Brilleman, Samuel L., Eren M. Elci, Jacqueline Buros Novik, and Rory
Wolfe. 2020. “Bayesian Survival Analysis Using the Rstanarm R Package.”
*arXiv:2002.09633 \[Stat\]*, February.
<http://arxiv.org/abs/2002.09633>.

</div>

<div id="ref-carpenter2017">

Carpenter, Bob, Andrew Gelman, Matthew D. Hoffman, Daniel Lee, Ben
Goodrich, Michael Betancourt, Marcus Brubaker, Jiqiang Guo, Peter Li,
and Allen Riddell. 2017. “Stan: A Probabilistic Programming Language.”
*Journal of Statistical Software* 76 (1): 1–32.
<https://doi.org/10.18637/jss.v076.i01>.

</div>

</div>
