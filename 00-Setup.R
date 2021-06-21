
install.packages('tidyverse')
install.packages('remotes')

## install feature/survival branch of rstanarm
## NOTE: This is a WIP branch. Please report any issues you 
## find so that we can improve it.
remotes::install_github('stan-dev/rstanarm', ref = 'feature/survival', build_vignettes = FALSE)


# mainly for accessing the Generable API (which we will not use).
# This contains some convenience functions which we will use.
remotes::install_github('generable/rgeco')

## install other prerequisites
packages <- c('RISCA', 'survival', 'tidybayes', 'bayesplot', 'loo', 'cowplot', 'bayesplot', 'simsurv')
install.packages(pkgs = packages)



