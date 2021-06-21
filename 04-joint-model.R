
library(rstanarm)
data(pbcLong)
data(pbcSurv)
library(tidyverse)

# ---- explore long data ----

ids <- c(25,31:33,36,38:40)
pbcLong_subset <- pbcLong[pbcLong$id %in% ids, ]
pbcLong_subset <- merge(pbcLong_subset, pbcSurv)
pbcLong_subset$Died <- factor(pbcLong_subset$death,
                              labels = c("No", "Yes"))
patient_labels <- paste("Patient", 1:8)
names(patient_labels) <- ids
ggplot() + 
  geom_line(aes(y = logBili, x = year, lty = Died), 
            data = pbcLong_subset) + 
  facet_wrap(~ id, ncol = 4, labeller = labeller(id = patient_labels)) +
  theme_bw() +
  ylab("Log serum bilirubin") +
  xlab("Time (years)")

# ---- full JM ----

mod1 <- stan_jm(
  formulaLong = logBili ~ sex + trt + year + (year | id), 
  dataLong = pbcLong,
  formulaEvent = survival::Surv(futimeYears, death) ~ sex + trt, 
  dataEvent = pbcSurv,
  time_var = "year",
  chains = 1, refresh = 2000, seed = 12345)


# ---- Long part ----

long1 <- stan_glmer(
  formula = logBili ~ sex + trt + year + (year | id), 
  data = pbcLong,
  chains = 1, refresh = 2000, seed = 12345)

# ---- survival part ----

surv1 <- stan_surv(
  formula = survival::Surv(futimeYears, death) ~ sex + trt,
  data = pbcSurv,
  chains = 1, refresh = 2000, seed = 12345)
  
)

# ---- comparison of predicted values ----

pplong_jm <- posterior_traj(mod1, ids = ids)

# plot pplong 
ggplot(pplong_jm, 
       aes(x = year, y = yfit, ymin = ci_lb, ymax = ci_ub)) + 
  geom_ribbon(alpha = 0.2) + facet_wrap(~ id) +
  geom_line(aes(y = logBili, x = year, lty = Died), 
            data = pbcLong_subset, inherit.aes = FALSE) + 
  ggtitle('Posterior predictive checks for Joint Model')

# obtain standarized predictions
ppsurv_jm <- seq(from = 0, to = 10, by = 1) %>%
  purrr::map_dfr(~ posterior_survfit(mod1, standardise = TRUE, times = .x))

# plot ppsurv
ggplot(ppsurv_jm, 
       aes(x = year, y = median, ymin = ci_lb, ymax = ci_ub)) + 
  geom_ribbon(alpha = 0.2) + 
  ggtitle('Posterior predictive checks for Joint Model')


                                                                                                                        data = pbcLong_subset, inherit.aes = FALSE)
