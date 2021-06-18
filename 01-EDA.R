
library(tidyverse)
ggplot2::theme_set(ggplot2::theme_minimal())

# install `rgeco` package
# we will use a few convenience functions from this package
remotes::install_github('generable/rgeco')
library(rgeco)

# data downloaded from: https://ti.arc.nasa.gov/tech/dash/groups/pcoe/prognostic-data-repository/#turbofan
DATA_DIR <- '~/Downloads/CMAPSSData'
DATA_PATH <- file.path(DATA_DIR, 'train_FD001.txt')
DATA_NAME <- basename(DATA_PATH)

# ---- Load data ----

# this is a dataset of simulated turbofan sensor readings over the lifetime of 100 turbofans
# In this training set, each machine is followed up to the failure time.
d <- readr::read_delim(DATA_PATH, delim = ' ', 
                       col_names = c('id', 'time', 
                                     stringr::str_c('setting', seq_len(3)),
                                     stringr::str_c('sensor', seq_len(21)),
                                     'empty'),
                       trim_ws = TRUE
                       ) %>%
  dplyr::select(-empty)

# ---- Format as time-to-event ----

## Let's encode the explicit time-to-failure for each machine.
## each machine will have an "event" code of 1 since the failure is observed in this dataset
events <- d %>%
  dplyr::group_by(id) %>%
  dplyr::summarize(start = min(time)-1,
                   end = max(time),
                   event = TRUE,
                   group = 'training') %>%
  dplyr::ungroup()

# we can use this `events` dataframe to plot the survival time over this set of machines
rgeco::prep_km_data(data = events, formula = survival::Surv(origin = start, time = end, event = event) ~ group) %>%
  ggplot(aes(x = time, y = estimate)) + 
  geom_step() +
  scale_y_continuous('Machines Running (%)', labels = scales::percent) + 
  scale_x_continuous('Time (cycles)') +
  ggtitle('Kaplan-Meier Estimates of Turbofan Lifetimes') +
  labs(caption = glue::glue('For the dataset: {DATA_NAME}'))

# ---- EDA of settings data ----

## This dataset records each of 3 operational settings, for each fan
## The operational settings are assumed to be exogenous, that is: they influence the failure rate for each machine, but are not influenced by it.

d %>% 
  tidyr::gather(setting_name, setting_value, dplyr::starts_with('setting')) %>%
  ggplot(., aes(x = time, y = setting_value, group = id)) + 
  geom_line(alpha = 0.1, size = 0.1) + 
  facet_wrap(~ setting_name, scale = 'free_y') +
  scale_x_continuous('Time (cycles') + 
  scale_y_continuous('Setting value') +
  ggtitle('Operational settings over time for each machine') +
  labs(caption = glue::glue('For the dataset: {DATA_NAME}'))

# We learn from this that the `setting3` value is 100 for all machines.
# Also, setting1 & setting2 are on different scales.

# let's plot as percent change from baseline
d %>% 
  tidyr::gather(setting_name, setting_value, dplyr::starts_with('setting')) %>%
  dplyr::group_by(id, setting_name) %>%
  dplyr::mutate(percent_change = (setting_value - dplyr::first(setting_value, order_by=time)) / dplyr::first(setting_value, order_by=time)) %>%
  dplyr::ungroup() %>%
  ggplot(., aes(x = time, y = percent_change, group = id)) + 
  geom_line(alpha = 0.1, size = 0.1) + 
  facet_wrap(~ setting_name, scale = 'free_y') +
  scale_x_continuous('Time (cycles') + 
  scale_y_continuous('Setting value (% change from baseline)', labels = scales::percent) +
  ggtitle('Operational settings over time for each machine') +
  labs(caption = glue::glue('For the dataset: {DATA_NAME}'))

# both setting1 & setting2 have a reasonable frequency of extreme values.
# also, It looks like setting2 has discrete values whereas setting1 is continuous.

# ---- EDA: Setting1 ----

d %>% 
  dplyr::filter(time == 1) %>%
  ggplot(aes(x = setting1)) + 
  geom_histogram() +
  scale_x_continuous('Setting1 value at baseline') +
  ggtitle('Distribution of Setting1 values at baseline') +
  labs(caption = glue::glue('For the dataset: {DATA_NAME}'))

# At baseline, this setting shows a skewed distribution

d %>% 
  ggplot(aes(x = setting1, group = time, fill = time)) + 
  geom_histogram(position = 'stack') +
  scale_x_continuous('Setting1 value') +
  ggtitle('Distribution of Setting1 values over all times') +
  labs(caption = glue::glue('For the dataset: {DATA_NAME}'))

# Over all time, this shows a distribution centered at 0
# with very little evidence of being skewed

d %>% 
  dplyr::mutate(setting1_normalized = (setting1 - mean(setting1))/sd(setting1)) %>%
  ggplot(aes(x = setting1_normalized, group = time, fill = time)) + 
  geom_histogram() +
  scale_x_continuous('[Setting1 - mean(Setting1)] / sd(Setting1)') +
  ggtitle('Distribution of Normalized Setting1 values') +
  labs(caption = glue::glue('For the dataset: {DATA_NAME}'))

# these data, like those of setting2, appear well-behaved when normalized

# ---- EDA: Setting2 ----

d %>% 
  dplyr::filter(time == 1) %>%
  ggplot(aes(x = setting2)) + 
  geom_histogram() +
  scale_x_continuous('Setting2 value at baseline') +
  ggtitle('Distribution of Setting2 values at baseline') +
  labs(caption = glue::glue('For the dataset: {DATA_NAME}'))

d %>% 
  ggplot(aes(x = setting2, group = time, fill = time)) + 
  geom_histogram() +
  scale_x_continuous('Setting2 value') +
  ggtitle('Distribution of Setting2 values') +
  labs(caption = glue::glue('For the dataset: {DATA_NAME}'))

d %>% 
  dplyr::mutate(setting2_normalized = (setting2 - mean(setting2))/sd(setting2)) %>%
  ggplot(aes(x = setting2_normalized, group = time, fill = time)) + 
  geom_histogram() +
  scale_x_continuous('[Setting2 - mean(Setting2)] / sd(Setting2)') +
  ggtitle('Distribution of Normalized Setting2 values') +
  labs(caption = glue::glue('For the dataset: {DATA_NAME}'))


# ---- EDA: correlation of settings ----

## Are setting1 & setting2 correlated?

d %>% 
  dplyr::mutate(setting1_normalized = (setting1 - mean(setting1))/sd(setting1),
                setting2_normalized = (setting2 - mean(setting2))/sd(setting2)) %>%
  ggplot(aes(x = setting1_normalized, y = setting2_normalized)) + 
  geom_jitter(alpha = 0.1) +
  scale_x_continuous('Setting1 (normalized)') +
  scale_y_continuous('Setting2 (normalized)') +
  ggtitle('Correlation of Setting1 and Setting2 values') +
  labs(caption = glue::glue('For the dataset: {DATA_NAME}'))

# There is almost no evidence of correlation between these two values

## Are setting1 & setting2 correlated with baseline values?

d %>% 
  dplyr::mutate(setting1_normalized = (setting1 - mean(setting1))/sd(setting1),
                setting2_normalized = (setting2 - mean(setting2))/sd(setting2)) %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(first_setting1 = dplyr::first(setting1_normalized, order_by = time)) %>%
  dplyr::ungroup() %>%
  ggplot(aes(x = first_setting1, y = setting2_normalized)) + 
  geom_jitter(alpha = 0.2) +
  scale_x_continuous('First Setting1 (normalized)') +
  scale_y_continuous('Setting2 (normalized)') +
  ggtitle('Correlation of First Setting1 and all Setting2 values') +
  labs(caption = glue::glue('For the dataset: {DATA_NAME}'))


d %>% 
  dplyr::mutate(setting1_normalized = (setting1 - mean(setting1))/sd(setting1),
                setting2_normalized = (setting2 - mean(setting2))/sd(setting2)) %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(first_setting1 = dplyr::first(setting1_normalized, order_by = time),
                first_setting2 = dplyr::first(setting2_normalized, order_by = time)) %>%
  dplyr::ungroup() %>%
  ggplot(aes(x = first_setting2, y = setting1_normalized)) + 
  geom_jitter(alpha = 0.2) +
  scale_x_continuous('First Setting2 (normalized)') +
  scale_y_continuous('Setting1 (normalized)') +
  ggtitle('Correlation of First Setting2 and all Setting1 values') +
  labs(caption = glue::glue('For the dataset: {DATA_NAME}'))

# It appears that setting1 & setting2 are set randomly for all machines
# irrespective of the initial setting or the other setting values

# In other words, we do not have a scenario where there are combinations of settings 
# that are more likely than others.

