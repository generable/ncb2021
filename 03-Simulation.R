

library(simsurv)
library(rstanarm)
library(tidyverse)
library(RISCA)

ggplot2::theme_set(tidybayes::theme_tidybayes())

# ---- exponential model ----

# simulate data for an exponential model with a fixed treatment effect
set.seed(1234)
covs <- data.frame(id = 1:300, trt = stats::rbinom(300, 1L, 0.5))
s1 <- simsurv(lambdas = 0.1, betas = c(trt = -0.5),
              dist = 'exponential',
              x = covs, maxt = 50) %>%
  mutate(true_lambda = 0.1, true_beta = -0.5, dist = 'exponential')
sim_exp <- left_join(covs, s1, by = 'id')

# plot K-M curve
sim_km <- prep_km_data(sim_exp, Surv(time = eventtime, event = status) ~ trt)
sim_km %>%
  ggplot(., aes(x = time, y = estimate, group = trt, colour = trt)) +
  geom_step() +
  ggtitle('Simulated survival data (exponential hazard)')

# estimate RMST
sim_rmst <- RISCA::rmst(times = sim_km$time, surv.rates = sim_km$estimate, max.time = 20)

coxph(formula = Surv(time = eventtime, event = status) ~ trt,
      data = sim_exp)

# Fit stan_surv to these data
sim_fit <- stan_surv(
  formula = Surv(time = eventtime, event = status) ~ trt,
  data = sim_exp, 
  basehaz = 'exp',
  prior_intercept = normal(0, 1))

# print summary
sim_fit

# plot treatment effect
tidybayes::gather_draws(sim_fit, trt) %>%
  ggplot(., aes(x = .value)) +
  stat_dotsinterval() +
  geom_vline(data = sim_exp, aes(xintercept = true_beta)) +
  scale_y_continuous('density') +
  scale_x_continuous('Beta[trt]') +
  ggtitle('Posterior Estimates of Treatment Beta fit to Simulated Data') +
  labs(caption = glue::glue('Showing posterior draws for a single replicate truncated at 50',
                            '\nVertical line shows true value used for simulation ({scales::comma(-0.5, accuracy = 0.01)})'))
       


# ---- test sensitivity to sample size ----

simulate_data <- function(n, balance = 0.5, maxt = 50, seed = 4243) {
  set.seed(seed)
  covs <- data.frame(id = seq_len(n), trt = stats::rbinom(n, 1L, 0.5))
  s1 <- simsurv(lambdas = 0.1, betas = c(trt = -0.5), dist = 'exponential',
                x = covs, maxt = maxt) %>%
    mutate(true_lambda = 0.1, true_beta = -0.5, dist = 'exponential')
  sim_exp <- left_join(covs, s1, by = 'id')
}



scenarios <- expand_grid(
  n = seq(from = 50, to = 100, by = 50),
  maxt = seq(from = 40, to = 60, by = 10),
  seed = seq_len(3))

simulations <- purrr::pmap(scenarios, simulate_data)
all_sims <- bind_rows(simulations, .id = 'scenario')

fits <- simulations %>%
  purrr::map(
    ~ stan_surv(formula = Surv(time = eventtime, event = status) ~ trt,
                data = .x, basehaz = 'exp',
                prior_intercept = normal(0, 1))
  )

estimates <- fits %>%
  purrr::map_dfr(tidybayes::gather_draws, trt, .id = 'scenario') %>%
  dplyr::mutate(scenario = as.integer(scenario)) %>%
  dplyr::left_join(scenarios %>% dplyr::mutate(scenario = row_number()))

ggplot(estimates, aes(x = .value, y = factor(maxt),
                      group = factor(maxt), colour = factor(maxt))) +
  stat_dotsinterval(position = position_dodge(width = 0.5)) +
  facet_wrap( ~ stringr::str_c('Sample size = ', n), scale = 'free_y') +
  geom_vline(data = all_sims, aes(xintercept = true_beta), colour = 'gray10', linetype = 'dotted') +
  scale_y_discrete('Max follow-up time') + 
  scale_x_continuous('Treatment beta') +
  labs(caption = glue::glue('Showing posterior draws for three independent replicates truncated at different timepoints',
                            '\nVertical line shows true value used for simulation ({scales::comma(unique(all_sims$true_beta), accuracy = 0.01)})')
  ) +
  ggtitle('Simulation study of Posterior Beta')

ggplot(estimates, aes(x = .value, y = factor(maxt),
                      group = factor(maxt), fill = stat(x < 0))) +
  stat_dotsinterval(position = position_dodge(width = 0.5), quantiles = 50) +
  facet_wrap( ~ stringr::str_c('Sample size = ', n), scale = 'free_y') +
  geom_vline(data = all_sims, aes(xintercept = true_beta), colour = 'gray10', linetype = 'dotted') +
  scale_y_discrete('Max follow-up time') + 
  scale_x_continuous('Treatment beta') +
  labs(caption = glue::glue('Showing posterior draws for three independent replicates truncated at different timepoints',
                            '\nVertical line shows true value used for simulation ({scales::comma(unique(all_sims$true_beta), accuracy = 0.01)})')
  ) +
  ggtitle('Simulation study of Posterior Beta')

ggplot(estimates, aes(y = .value, x = factor(maxt), colour = factor(seed), group = factor(seed))) +
  stat_pointinterval(position = position_dodge(width = 0.2), .width = c(.5, .8)) +
  facet_wrap( ~ stringr::str_c('Sample size = ', n), scale = 'free_y', ncol = 1) +
  geom_hline(data = all_sims, aes(yintercept = true_beta), colour = 'gray10', linetype = 'dotted') +
  scale_x_discrete('Max followup time (maxt)') + 
  scale_y_continuous('Treatment beta') +
  labs(caption = glue::glue('Showing posterior median and {glue::glue_collapse(scales::percent(c(0.5, 0.8)), sep = "/")} CrI',
                            '\nHorizontal line shows true value used for simulation ({scales::comma(unique(all_sims$true_beta), accuracy = 0.01)})')
       ) +
  ggtitle('Summary of Posterior Beta estimate fit to simulated data')

# ---- check performance using loo ----

sim_loo <- loo(sim_fit)

# compare against a model with a different baseline hazard
sim_fit_alt <- stan_surv(
  formula = Surv(time = eventtime, event = status) ~ trt,
  data = sim_exp, 
  basehaz = 'ms',
  prior_intercept = normal(0, 1))

sim_loo_alt <- loo(sim_fit_alt)

# Compare models
loo_compare(sim_loo, sim_loo_alt)

