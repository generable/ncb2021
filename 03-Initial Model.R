
library(rstanarm)
library(tidybayes)
library(tidyverse)

# This code references data loaded in our EDA example, and helper-functions
source(here::here('functions', 'helper_functions.R'))
events <- read_data_as_events()

# declare desired parameters for Stan
CHAINS <- 4
CORES <- 2
ITER <- 2000
SEED <- 42

# ---- Baseline hazard (exp) ----
# prior draws (again)
prior_exp_hazard2 <- stan_surv(
  formula = Surv(time = end, event = event) ~ 1,
  data = events,
  basehaz = "exp",
  prior_PD = TRUE,
  chains = CHAINS,
  cores = CORES,
  iter = ITER,
  seed = SEED,
  prior_intercept = normal(0, 1),
  prior = normal(0, 0.5))

# fit posterior model
post_exp_hazard2 <- stan_surv(
  formula = Surv(time = end, event = event) ~ 1,
  data = events,
  basehaz = "exp",
  prior_PD = FALSE,
  chains = CHAINS,
  cores = CORES,
  iter = ITER,
  seed = SEED,
  prior_intercept = normal(0, 1),
  prior = normal(0, 0.5))

cowplot::plot_grid(rmst_check_plot(prior_exp_hazard2, tau = 300) + 
                     coord_cartesian(xlim = c(0, 300)) +
                     scale_x_continuous(''),
                   rmst_check_plot(post_exp_hazard2, tau = 300) + 
                     coord_cartesian(xlim = c(0, 300)) + 
                     scale_x_continuous('RMST with tau = 300') +
                     labs(caption = 'RMST from intercept only model with exp baseline hazard (model 2)'),
                   ncol = 1, labels = c('Priors', 'Posterior'), hjust = 0
)


# That's curious. Why is the observed RMST is different from the posterior predicted?

# try with mspline model

post_ms_hazard2 <- stan_surv(
  formula = Surv(time = end, event = event) ~ 1,
  data = events,
  basehaz = "ms",
  prior_PD = FALSE,
  chains = CHAINS,
  cores = CORES,
  iter = ITER,
  seed = SEED,
  prior_intercept = normal(5, 1),
  prior = normal(0, 1),
  basehaz_ops = list(df = 5))

cowplot::plot_grid(rmst_check_plot(post_exp_hazard2, tau = 300) + 
                     coord_cartesian(xlim = c(0, 300)) +
                     scale_x_continuous(''),
                   rmst_check_plot(post_ms_hazard2, tau = 300) + 
                     coord_cartesian(xlim = c(0, 300)) + 
                     scale_x_continuous('RMST with tau = 300') +
                     labs(caption = 'RMST from intercept only models'),
                   ncol = 1, labels = c('exp model', 'm-spline model'), hjust = 0
)


# ---- model comparison with loo-psis ----

loo_compare(loo(post_exp_hazard2), loo(post_ms_hazard2))

# ---- compare survival curves ----

ppsurv_exp <- posterior_survfit(post_exp_hazard2, newdata = events, type = 'surv', prob = 0.5)
ppsurv_ms <- posterior_survfit(post_ms_hazard2, newdata = events, type = 'surv', prob = 0.5)

ggplot(dplyr::bind_rows(exp = ppsurv_exp, ms = ppsurv_ms, .id = 'type'),
       aes(x = time, y = median,
                    ymin = ci_lb, ymax = ci_ub,
                    group = id, colour = type, fill = type)) + 
  ggdist::geom_lineribbon() + 
  scale_y_continuous('Prior predicted survival', labels = scales::percent) +
  ggtitle('Prior predicted survival over time') +
  labs(caption = stringr::str_c('Exp model using scale of 1 on prior intercept',
                                'Showing median + 50% CrI',
                                sep = '\n'))

