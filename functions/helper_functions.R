
library(RISCA)
library(tidyverse)

# Evaluate RMST for the fitted model
#
# @param fit The fitted stan_surv model.
# @param tau The time horizon for evaluating RMST.
rmst_check <- function(fit, tau = 300, newdata = NULL) {
  
  df_gk_ <- tibble(
    locations = c(
      0.991455371120813, 		
      0.949107912342759, 	 
      0.864864423359769, 	
      0.741531185599394, 	 
      0.586087235467691, 	
      0.405845151377397, 	 
      0.207784955007898, 	
      0.000000000000000
    ),
    weights = c(
      0.022935322010529,
      0.063092092629979,
      0.104790010322250,
      0.140653259715525,
      0.169004726639267,
      0.190350578064785,
      0.204432940075298,
      0.209482141084728
    )
  )
  
  df_gk <- df_gk_ %>% 
    dplyr::mutate(locations = -locations) %>% 
    bind_rows(df_gk_) %>% 
    distinct() %>% 
    arrange(locations)
  
  if (is.null(newdata)) {
    nd  <- data.frame(group = 'training', normalized_setting1 = 0, normalized_setting2 = 0, id = 1)
  } else {
    nd <- newdata
  }
  locs    <- df_gk$locations
  weights <- df_gk$weights
  qpts    <- tau/2 + tau/2 * locs
  # empty list to store predicted survival
  ps <- list() 
  
  # use all MCMC draws so that the sample draws
  # are used at each quadrature point
  ndraws <- nrow(as.matrix(fit))
  # evaluate predicted survival at each quadrature point
  for (q in 1:length(qpts)) {
    ps[[q]] <- posterior_survfit(object        = fit,
                                 newdata       = nd, 
                                 times         = qpts[q], 
                                 extrapolate   = FALSE, 
                                 condition     = FALSE,
                                 return_matrix = TRUE,
                                 draws         = ndraws)
  }
  
  # convert predicted survival to rmst using GK quadrature
  rmst <- 1:length(qpts) %>%
    map(~ ps[[.]][[1]] * weights[[.]] * tau/2) %>%
    Reduce('+', .)
  
  # return as data frame
  return(posterior::as_draws_df(rmst))
}

rmst_check_plot <- function(fit, tau = 300, newdata = NULL) {
  
  if (is.null(newdata)) {
    surv_data <- fit$data
  } else {
    surv_data <- newdata
  }
  surv_formula <- fit$formula$lhs_form
  km_data <- survival::survfit(data = surv_data, formula = surv_formula)
  km_rmst <- RISCA::rmst(times = km_data$time, surv.rates = km_data$surv, max.time = tau)
  
  df_rmst <- rmst_check(fit, tau = tau, newdata = NULL)
  
  df_rmst %>% 
    bayesplot::mcmc_hist(pars = c("1"), binwidth = 0.1) +
    bayesplot::vline_at(km_rmst, color='green') -> 
    p1
  
  p1
}

# Evaluate the approximate leave-one-out mean posterior survival function using Pareto Smoothed Importance Sampling
#
# @param fit The fitted stan_surv model.
# @param times Vector of times for which the survival functions should be evaluated.
get_loo_mean_surv <- function(fit, times) {
  
  ntimes <- length(times)
  
  # log likelihood for the fitted model
  log_lik <- log_lik(fit)
  
  # extract posterior draws
  post <- as.array(fit) 
  
  # evaluate relative efficiencies (effective samples / total samples)
  r_eff <- relative_eff(
    exp(log_lik), 
    chain_id = rep(1:dim(post)[2], each = dim(post)[1])) 
  
  # evaluate loo
  loo_object <- loo(log_lik, 
                    r_eff     = r_eff,
                    cores     = 2, 
                    save_psis = TRUE)
  # empty list to store predicted survival
  ps <- list() 
  
  # evaluate predicted survival at each time
  for (q in 1:ntimes) {
    ps[[q]] <- posterior_survfit(object        = fit,
                                 times         = times[q], 
                                 extrapolate   = FALSE, 
                                 condition     = FALSE,
                                 return_matrix = TRUE,
                                 draws         = nrow(as.matrix(fit)))
  }
  
  do.call(cbind,
          map(1:ntimes, 
              ~ E_loo(x           = ps[[.]][[1]], 
                      psis_object = loo_object$psis_object, 
                      type        = "mean", 
                      log_ratios  = -log_lik)$value)
  ) 
}

read_data <- function(
  data_dir = '~/Downloads/CMAPSSData',
  data_name = 'train_FD001.txt',
  data_path = file.path(DATA_DIR, data_name)
) {
  if (missing(data_name)) {
    data_name <- basename(data_path)
  }
  is_training <- stringr::str_detect(data_name, pattern = '^train')
  d <- readr::read_delim(data_path,
                         delim = ' ', 
                         col_names = c('id', 'time', 
                                       stringr::str_c('setting', seq_len(3)),
                                       stringr::str_c('sensor', seq_len(21)),
                                       'empty'),
                         trim_ws = TRUE
  ) %>%
    dplyr::select(-empty) %>% 
    dplyr::mutate_at(.vars = dplyr::vars(dplyr::starts_with('setting'), dplyr::starts_with('sensor')),
                     .funs = list(normalized = ~ (.x - mean(.x))/sd(.x))) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(start = time-1,
                  end = time,
                  event = dplyr::if_else(time == max(time), is_training, !is_training),
                  group = dplyr::if_else(is_training, 'training', 'test')) %>%
    dplyr::ungroup()
}

read_data_as_events <- function(
    data_dir = '~/Downloads/CMAPSSData',
    data_name = 'train_FD001.txt',
    data_path = file.path(DATA_DIR, data_name)
  ) {
  d <- read_data(data_path = data_path)
  events <- d %>%
    dplyr::group_by(id, group) %>%
    dplyr::summarize(start = min(time)-1,
                     end = max(time),
                     event = dplyr::last(event, order_by = time)) %>%
    dplyr::ungroup()
}