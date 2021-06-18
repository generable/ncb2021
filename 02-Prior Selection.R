
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

# We start by investigating various forms of the the baseline hazard
# In practice, we want the 
remotes::install_github('stan-dev/rstanarm', ref='feature/survival')
# After installing, restart you R session

# In Rstudio: Use the menu item Session > Restart R 
#   or the associated keyboard shortcut Ctrl+Shift+F10 (Windows and Linux) or Command+Shift+F10 (Mac OS).
# Otherwise: `q()`, then start a new R session
.rs.restartR()

# ---- Baseline hazard (exp) ----
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

# draw from the prior predictive distribution of the stan_surv survival model
prior_exp_hazard <- stan_surv(
  formula = Surv(end, event) ~ 1,
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
pphaz <- posterior_survfit(prior_exp_hazard, type = 'haz', prob = 0.5)
ggplot(pphaz, aes(x = time, y = median, ymin = ci_lb, ymax = ci_ub)) + 
  tidybayes::geom_lineribbon(fill = 'grey') +
  scale_y_continuous('Median hazard rate') +
  ggtitle('Prior predicted hazard over time (events / machine-cycle') +
  labs(caption = 'Model using default priors for exp baseline hazard') +
  coord_cartesian(ylim = c(0, 0.1))

ggplot(pphaz, aes(x = time, y = median*1000, ymin = ci_lb*1000, ymax = ci_ub*1000)) + 
  tidybayes::geom_lineribbon() + 
  scale_y_continuous('Median hazard rate \n per 1000 machine-cycles') +
  ggtitle('Prior predicted hazard over time') +
  labs(caption = 'Model using default priors for exp baseline hazard') +
  

# ---- ... Updating priors (v2) ---- 

prior_exp_hazard2 <- update(prior_exp_hazard,
                           prior_intercept = normal(0, 1),
                           prior = normal(0, 0.5))

## Check prior predicted RMST vs observed RMST
## NOTE: in practice, you wouldn't check this against the data, but against your expectations
## What is a realistic survival time for a turbofan?
rmst_check_plot(prior_exp_hazard2, tau = 300)

## Plot prior hazard over time
pphaz2 <- posterior_survfit(prior_exp_hazard2, type = 'haz')
ggplot(pphaz2, aes(x = time, y = median)) + 
  geom_line() + 
  scale_y_continuous('Median hazard rate') +
  ggtitle('Prior predicted hazard over time (events / machine-cycle') +
  labs(caption = 'Model using default priors for exp baseline hazard')

ggplot(pphaz2, aes(x = time, y = median*1000)) + 
  geom_line() + 
  scale_y_continuous('Median hazard rate \n per 1000 machine-cycles') +
  ggtitle('Prior predicted hazard over time') +
  labs(caption = 'Model using default priors for exp baseline hazard')

# ---- Baseline hazard (m-spline) ----

# draw from the prior predictive distribution of the stan_surv survival model
prior_ms_hazard_df_5 <- update(prior_exp_hazard2,
                               basehaz = "ms",
                               basehaz_ops = list(df = 5))

## Check prior predicted RMST vs observed RMST
rmst_check_plot(prior_ms_hazard_df_5, tau = 300)

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
