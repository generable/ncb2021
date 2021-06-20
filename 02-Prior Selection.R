
# Our goal is to fit a survival model to these data

# Recall that a standard CoxPH survival assumes that the instantaneous rate of failure (hazard rate) 
# for an individual i at any time t 
# 
#    h_i(t) = h_0(t)[exp(f_i(t))]
# 
# This function has two parts:
#
# h_0(t): baseline hazard at time t. This hazard is shared among all individuals i
# exp(f_i(t)): relative hazard for individual i at time t.
# 
# This model is a proportional hazards (the PH part of CoxPH) when the relative hazard for individual i is constant for all times t.
# That is, the relative hazards for two individuals are proportional to one another over the follow-up time.

# Aside: A common forms of f_i(t) is X_i*B

# Model selection will thus proceed in two parts:
# 1. Specifying the baseline hazard (functional form + priors)
# 2. Modeling the relative hazard Describing covariate effects on the hazard

# In this dataset, we have time-varying covariates so our f_i(t) will look like: X_it*B

# ---- Installing rstanarm: feature/survival branch ----

remotes::install_github('stan-dev/rstanarm', ref='feature/survival')
# After installing, restart you R session

# In Rstudio: Use the menu item Session > Restart R 
#   or the associated keyboard shortcut Ctrl+Shift+F10 (Windows and Linux) or Command+Shift+F10 (Mac OS).
# Otherwise: `q()`, then start a new R session
if (TRUE == FALSE) {
  .rs.restartR()
}

# ---- Setup ----
# This code uses a few helper-functions
source(file.path('functions', 'helper_functions.R'))

library(rstanarm)
library(survival)
library(bayesplot)
library(tidybayes)
library(cowplot)
options(mc.cores = parallel::detectCores())

# load data
d <- read_data(data_path = '~/Downloads/CMAPSSData/train_FD001.txt')
events <- read_data_as_events(data_path = '~/Downloads/CMAPSSData/train_FD001.txt')

# declare desired parameters for Stan
CHAINS <- 4
CORES <- 2
ITER <- 2000
SEED <- 42

# ---- Baseline hazard (exp) ----

# draw from the prior predictive distribution of the stan_surv survival model
prior_exp_hazard <- stan_surv(
  formula = Surv(time = end, event = event) ~ 1,
  data = events,
  basehaz = "exp",
  prior_PD = TRUE,
  chains = CHAINS,
  cores = CORES,
  iter = ITER,
  seed = SEED)

## What are the priors for this model?
prior_summary(prior_exp_hazard)

## Check prior predicted RMST vs observed RMST
rmst_check_plot(prior_exp_hazard, tau = 300)

## Plot prior hazard over time
pphaz <- posterior_survfit(prior_exp_hazard, type = 'haz', newdata = events, prob = 0.5)
ggplot(pphaz, aes(x = time, y = median, ymin = ci_lb, ymax = ci_ub, group = id)) + 
  tidybayes::geom_lineribbon(fill = 'grey') +
  scale_y_continuous('Prior predictive hazard rate') +
  ggtitle('Prior predicted hazard over time (events / machine-cycle') +
  labs(caption = stringr::str_c('Model using default priors for exp baseline hazard',
                                'Showing median + 50% CrI',
                                'y-axis is truncated at 0.1',
                                sep = '\n')) +
  coord_cartesian(ylim = c(0, 0.1))

# just how large are the CI around the hazard?
pphaz %>% dplyr::arrange(time) %>% head()

# ---- ... Updating priors (v2) ---- 

# let's reduce the prior scale on the intercept term
prior_exp_hazard2 <- update(prior_exp_hazard,
                           prior_intercept = normal(0, 1),
                           prior = normal(0, 0.5))

## Check prior predicted RMST vs observed RMST
## NOTE: in practice, you wouldn't check this against the data, but against your expectations
## What is a realistic survival time for a turbofan?
rmst_check_plot(prior_exp_hazard2, tau = 300)

## Plot prior hazard over time

pphaz2 <- posterior_survfit(prior_exp_hazard2, type = 'haz', newdata = events, prob = 0.5)
ggplot(pphaz2, aes(x = time, y = median,
                   ymin = ci_lb, ymax = ci_ub,
                   group = id, fill = 'prior')) + 
  ggdist::geom_lineribbon(alpha = 0.04) + 
  scale_y_continuous('Median hazard rate') +
  ggtitle('Prior predicted hazard over time') +
  labs(caption = stringr::str_c('Exp model using scale of 1 on prior intercept',
                                'Showing median + 50% CrI',
                                sep = '\n'))
  
# This is a more reasonable range for our priors

## You'll note, by the way, that the data structure returned by `posterior_survfit`
## by default is a little strange.
  
## We are getting back an estimate of the hazard for each row of our input data
## For this reason, I recommend you always pass in `newdata`, so that the predicted 
## inputs are clear.
pphaz2 %>% dplyr::left_join(events, by = 'id') %>% head()

# NOTE: to get standardized predicted hazards, you would do something like this:
# (not run)
if (TRUE == FALSE) {
  # here, newdata aren't necessary since the hazard at each time is 
  # averaged over the input args
  pphaz2_standardized <- seq(from = 0, to = 300, by = 10) %>%
    purrr::set_names() %>%
    purrr::map_dfr(~ posterior_survfit(prior_exp_hazard2, type = 'haz', times = .x,
                                       standardise = TRUE, prob = 0.5),
                   .id = 'time'
    ) %>%
    dplyr::mutate(time = as.integer(time))
  
  # plot standardized hazard estimates
  ggplot(pphaz2_standardized,
         aes(x = time, y = median,
             ymin = ci_lb, ymax = ci_ub, 
             fill = 'prior')) + 
    ggdist::geom_lineribbon(alpha = 0.4) + 
    scale_y_continuous('Median population-averaged hazard rate') +
    ggtitle('Prior predicted hazard over time') +
    labs(caption = stringr::str_c('Exp model using scale of 1 on prior intercept',
                                  'Showing median + 50% CrI',
                                  sep = '\n'))
  
  ## it's a bit of a curious thing, why the population-averaged rate varies over 
  ## time, vs the previous summary which appeared perfectly constant?
  
  # The reason is that we are sampling draws. 
  # Use the draws parameter to include all draws in the estimate
  pphaz2_standardized2 <- seq(from = 0, to = 300, by = 10) %>%
    purrr::set_names() %>%
    purrr::map_dfr(~ posterior_survfit(prior_exp_hazard2, type = 'haz', times = .x,
                                       standardise = TRUE, prob = 0.5, draws = 4000),
                   .id = 'time'
    ) %>%
    dplyr::mutate(time = as.integer(time))
  ggplot(pphaz2_standardized2,
         aes(x = time, y = median,
             ymin = ci_lb, ymax = ci_ub, 
             fill = 'prior')) + 
    ggdist::geom_lineribbon(alpha = 0.4) + 
    scale_y_continuous('Median population-averaged hazard rate') +
    ggtitle('Prior predicted hazard over time') +
    labs(caption = stringr::str_c('Exp model using scale of 1 on prior intercept',
                                  'Showing median + 50% CrI over all draws',
                                  sep = '\n'))
}


# It can be helpful to transform this to a more meaningful scale
# for example, the number of failures per 100 machine-cycles
ggplot(pphaz2, aes(x = time, y = median*100,
                   ymin = ci_lb*100, ymax = ci_ub*100,
                   group = id, fill = 'prior')) + 
  ggdist::geom_lineribbon(alpha = 0.01) + 
  scale_y_continuous('Prior predicted hazard rate\n(per 100 machine-cycles)') +
  ggtitle('Prior predicted hazard over time') +
  labs(caption = stringr::str_c('Exp model using scale of 1 on prior intercept',
                                'Showing median + 50% CrI',
                                sep = '\n'))


# In a survival model, we consider the cumulative hazard as
# well as a baseline hazard
## Since this is a constant-hazard model, this is (perhaps unsurprisingly) linear over time.

ppcumhaz2 <- posterior_survfit(prior_exp_hazard2, type = 'cumhaz', newdata = events, prob = 0.5)
ggplot(ppcumhaz2, aes(x = time, y = median,
                   ymin = ci_lb, ymax = ci_ub,
                   group = id, fill = 'prior')) + 
  ggdist::geom_lineribbon(alpha = 0.04) + 
  scale_y_continuous('Prior predicted cumulative hazard') +
  ggtitle('Prior predicted cumulative hazard over time') +
  labs(caption = stringr::str_c('Exp model using scale of 1 on prior intercept',
                                'Showing median + 50% CrI',
                                sep = '\n'))

# From this, we can estimate the cumulative failure (or 1-failure = survival) 

ppsurv2 <- posterior_survfit(prior_exp_hazard2, type = 'surv', newdata = events, prob = 0.5)
ggplot(ppsurv2, aes(x = time, y = median,
                      ymin = ci_lb, ymax = ci_ub,
                      group = id, fill = 'prior')) + 
  ggdist::geom_lineribbon(alpha = 0.04) + 
  scale_y_continuous('Prior predicted survival', labels = scales::percent) +
  ggtitle('Prior predicted survival over time') +
  labs(caption = stringr::str_c('Exp model using scale of 1 on prior intercept',
                                'Showing median + 50% CrI',
                                sep = '\n'))

# ---- Baseline hazard (m-spline) ----

## Now we repeat this process of testing how resonable our priors are, using the 
## spline model

# The m-splines have a nice property of being relatively easy to integrate analytically
# so they are both computationally tractable & flexible enough to fit a variety of baseline hazards

# draw from the prior predictive distribution of the stan_surv survival model
prior_ms_hazard_df_5 <- update(prior_exp_hazard2,
                               basehaz = "ms",
                               basehaz_ops = list(df = 5))

## Check prior predicted RMST vs observed RMST
rmst_check_plot(prior_ms_hazard_df_5, tau = 300)

# Note: priors need to be updated in order to 
# ---- ... Updating priors (v2) ----

prior_ms_hazard_df_5_v2 <- update(prior_ms_hazard_df_5,
                                  prior_intercept = normal(5, 1),
                                  prior = normal(0, 1))

rmst_check_plot(prior_ms_hazard_df_5_v2, tau = 300)

# NOTE: the appropriate priors vary depending on the model!

# ---- Baseline hazard for frailty model (exp) ----

events_per_time <- d %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(start = time-1,
                end = time,
                event = dplyr::if_else(time == max(time), TRUE, FALSE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(setting1_normalized = (setting1 - mean(setting1))/sd(setting1),
                setting2_normalized = (setting2 - mean(setting2))/sd(setting2))

prior_exp_hazard2_frailty <- stan_surv(
  formula = Surv(time = start, time2 = end, event = event) ~ 1 + (1|id),
  data = events_per_time,
  basehaz = "exp",
  prior_PD = TRUE,
  chains = CHAINS,
  cores = CORES,
  iter = ITER,
  seed = SEED,
  prior_intercept = normal(0, 1),
  prior = normal(0, 0.5)
)

rmst_check_plot(prior_exp_hazard2_frailty, tau = 300)
rmst_check_plot(prior_exp_hazard2, tau = 300)

cowplot::plot_grid(rmst_check_plot(prior_exp_hazard2) + coord_cartesian(xlim = c(0, 300)),
                   rmst_check_plot(prior_exp_hazard2_frailty) + coord_cartesian(xlim = c(0, 300)),
                   ncol = 1, labels = c('base model', 'frailty model')
                   )

# What are the priors for this model?
prior_summary(prior_exp_hazard2_frailty)


# ---- ... Updating priors (v2) ----

# Update the prior on the covariance (describes variation in hazard across subjects)
prior_exp_hazard2_frailty2 <- 
  update(prior_exp_hazard2_frailty, 
         prior_covariance = decov(regularization = 1, concentration = 1, shape = 1, scale = 0.2)
  )

# now, we've ended up with approximately the same prior predictive distribution
rmst_check_plot(prior_exp_hazard2_frailty2, tau = 300)
cowplot::plot_grid(rmst_check_plot(prior_exp_hazard2, tau = 300) + coord_cartesian(xlim = c(0, 300)),
                   rmst_check_plot(prior_exp_hazard2_frailty2, tau = 300) + coord_cartesian(xlim = c(0, 300)),
                   ncol = 1, labels = c('base model', 'frailty model')
)

# look at updated, prior predicted survival curves
prior_ppsurv_exp <- posterior_survfit(prior_exp_hazard2_frailty2, 
                                      newdata = events %>% dplyr::filter(id <= 10) %>% 
                                        dplyr::mutate(normalized_setting1 = 0, normalized_setting2 = 0))

# see ?posterior_survfit.stansurv for details

# plot results for first 10 ids
ggplot(prior_ppsurv_exp, aes(x = time, y = median, ymin = ci_lb, ymax = ci_ub, group = id, colour = id)) + 
  ggdist::geom_lineribbon(alpha = 0.2) +
  ggtitle('Prior predicted survival for IDs 1-10')


ggplot(prior_ppsurv_exp, aes(x = time, y = median, ymin = ci_lb, ymax = ci_ub, group = id, colour = id, fill = id)) + 
  ggdist::geom_lineribbon(alpha = 0.1) +
  ggtitle('Prior predicted survival for IDs 1-10',
          subtitle = 'Showing 95% CrI')


# ---- Baseline hazard for frailty model (m-spline) ----

prior_ms_hazard_df_5_frailty <- stan_surv(
  formula = Surv(time = start, time2 = end, event = event) ~ 1 + (1|id),
  data = events_per_time,
  basehaz = "ms",
  prior_PD = TRUE,
  chains = CHAINS,
  cores = CORES,
  iter = ITER,
  seed = SEED,
  basehaz_ops = list(df = 5),
  prior_intercept = normal(5, 1),
  prior = normal(0, 1),
  prior_covariance = decov(regularization = 1, concentration = 1, shape = 1, scale = 0.2)
)

rmst_check_plot(prior_ms_hazard_df_5_frailty, tau = 300)
cowplot::plot_grid(rmst_check_plot(prior_ms_hazard_df_5_v2, tau = 300) + coord_cartesian(xlim = c(0, 300)),
                   rmst_check_plot(prior_ms_hazard_df_5_frailty, tau = 300) + coord_cartesian(xlim = c(0, 300)),
                   ncol = 1, labels = c('base m-spline model', 'frailty m-spline model'), hjust = 0
)

# look at predicted survival curves
prior_ppsurv_ms <- posterior_survfit(prior_ms_hazard_df_5_frailty, 
                                  newdata = events %>% dplyr::filter(id <= 10) %>% 
                                    dplyr::mutate(normalized_setting1 = 0, normalized_setting2 = 0))

# see ?posterior_survfit.stansurv for details

# plot prior predicted survival for first 10 ids
ggplot(prior_ppsurv_ms, aes(x = time, y = median, ymin = ci_lb, ymax = ci_ub, group = id, colour = id)) + 
  ggdist::geom_lineribbon(alpha = 0.2) +
  ggtitle('Prior predicted survival for IDs 1-10', subtitle = 'M-spline baseline hazard with df = 5')

# Compare to priors using exp-based frailty model
ggplot(dplyr::bind_rows(ms = prior_ppsurv_ms,
                        exp = prior_ppsurv_exp,
                        .id = 'baseline_hazard'), 
       aes(x = time, y = median, ymin = ci_lb, ymax = ci_ub, group = id, colour = baseline_hazard)) + 
  ggdist::geom_lineribbon(alpha = 0.2) +
  ggtitle('Prior predicted survival for IDs 1-10', subtitle = 'M-spline baseline hazard with df = 5')


# ---- Comparison of priors ----

ggplot(dplyr::bind_rows(ms = prior_ppsurv_ms,
                        exp = prior_ppsurv_exp,
                        .id = 'baseline_hazard'), 
       aes(x = time, y = median, ymin = ci_lb, ymax = ci_ub, group = id, colour = baseline_hazard, fill = baseline_hazard)) + 
  ggdist::geom_lineribbon(alpha = 0.05) +
  scale_y_continuous('Prior Predicted Survival', labels = scales::percent) + 
  ggtitle('Prior predicted survival for IDs 1-10', subtitle = 'M-spline baseline hazard with df = 5')
