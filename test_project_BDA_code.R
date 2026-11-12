library(tidyverse)
library(brms)
library(bayesplot)
library(loo)

# 1. Load and clean data
data_clean <- read.csv("https://vincentarelbundock.github.io/Rdatasets/csv/DAAG/nassCDS.csv") %>%
  drop_na(yearVeh) %>%
  mutate(
    y = if_else(dead == "dead", 1, 0),
    seatbelt = factor(seatbelt, levels = c("none", "belted")),
    airbag   = factor(airbag, levels = c("none", "airbag")),
    dvcat    = factor(dvcat, levels = c("1-9km/h", "10-24", "25-39", "40-54", "55+")),
    yearVeh  = as.factor(yearVeh)
  )

# 2. Define informative priors
priors_m1 <- c(
  prior(normal(-3.5, 1), class = "Intercept"),
  prior(normal(-0.2, 0.2), class = "b", coef = "seatbeltbelted"),
  prior(normal(-0.1, 0.3), class = "b", coef = "airbagairbag"),
  prior(normal(0.3, 1), class = "b", coef="dvcat10M24"),
  prior(normal(1.1, 0.8), class = "b", coef="dvcat25M39"),
  prior(normal(1.8, 0.7), class = "b", coef="dvcat40M54"),
  prior(normal(2.4, 0.3), class = "b", coef="dvcat55P")
)

priors_m2 <- c(
  priors_m1,
  prior(normal(0,1), class = "sd", group = "yearVeh")
)

# Model 1: Non-Hierarchical Logistic Model
m1_nonhier <- brm(
  formula = y ~ 1 + seatbelt + airbag + dvcat,
  data    = data_clean,
  family  = bernoulli(link = "logit"),
  prior   = priors_m1,
  chains  = 4,
  iter    = 3000,
  warmup  = 1000,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  cores   = 4,
  seed    = 42
)

# Model 2: Hierarchical Logistic Model (Varying Intercept by Vehicle Year)
m2_hier <- brm(
  formula = y ~ 1 + seatbelt + airbag + dvcat + (1 | yearVeh),
  data    = data_clean,
  family  = bernoulli(link = "logit"),
  prior   = priors_m2,
  chains  = 4,
  iter    = 3000,
  warmup  = 1000,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  cores   = 4,
  seed    = 42
)

# Summary checks for R-hat, ESS, and Divergences
summary(m1_nonhier)
summary(m2_hier)

# Extract MCMC diagnostic summary
np_m2 <- nuts_params(m2_hier$fit)
num_divergences <- sum(subset(np_m2, Parameter == "divergent__")$Value)
cat("Number of divergent transitions (Model 2):", num_divergences, "\n")

# Trace plot inspection for key parameters
mcmc_trace(m2_hier, pars = c("b_Intercept", "b_seatbeltbelted", "b_airbagairbag", "sd_yearVeh__Intercept"))

# Posterior predictive checks: Distribution of predicted proportions
pp_check(m1_nonhier, type = "dens_overlay", ndraws = 100) + ggtitle("PPC: Model 1 (Non-Hierarchical)")
pp_check(m2_hier, type = "dens_overlay", ndraws = 100) + ggtitle("PPC: Model 2 (Hierarchical)")

# Grouped PPC by vehicle year to inspect group-level fit
pp_check(m2_hier, type = "stat_grouped", group = "yearVeh", stat = "mean")

priors_weak <- c(
  prior(normal(0, 10), class = "Intercept"),
  prior(normal(0, 5), class = "b"),
  prior(student_t(3, 0, 2.5), class = "sd", group = "yearVeh")
)

m2_sensitivity <- brm(
  formula = y ~ 1 + seatbelt + airbag + dvcat + (1 | yearVeh),
  data    = data_clean,
  family  = bernoulli(link = "logit"),
  prior   = priors_weak,
  chains  = 4,
  iter    = 3000,
  warmup  = 1000,
  control = list(adapt_delta = 0.95),
  cores   = 4,
  seed    = 42
)

# Compare posterior estimates between informative and weak priors
posterior_comp <- bind_rows(
  as_draws_df(m2_hier) %>% mutate(Prior = "Informative"),
  as_draws_df(m2_sensitivity) %>% mutate(Prior = "Weakly Informative")
)

mcmc_areas(posterior_comp, pars = c("b_seatbeltbelted", "b_airbagairbag"), prob = 0.95)library(tidyverse)
library(brms)
library(bayesplot)
library(loo)

# 1. Load and clean data
data_clean <- read.csv("https://vincentarelbundock.github.io/Rdatasets/csv/DAAG/nassCDS.csv") %>%
  drop_na(yearVeh) %>%
  mutate(
    y = if_else(dead == "dead", 1, 0),
    seatbelt = factor(seatbelt, levels = c("none", "belted")),
    airbag   = factor(airbag, levels = c("none", "airbag")),
    dvcat    = factor(dvcat, levels = c("1-9km/h", "10-24", "25-39", "40-54", "55+")),
    yearVeh  = as.factor(yearVeh)
  )

# 2. Define informative priors
priors_m1 <- c(
  prior(normal(-3.5, 1), class = "Intercept"),
  prior(normal(-0.2, 0.2), class = "b", coef = "seatbeltbelted"),
  prior(normal(-0.1, 0.3), class = "b", coef = "airbagairbag"),
  prior(normal(0.3, 1), class = "b", coef="dvcat10M24"),
  prior(normal(1.1, 0.8), class = "b", coef="dvcat25M39"),
  prior(normal(1.8, 0.7), class = "b", coef="dvcat40M54"),
  prior(normal(2.4, 0.3), class = "b", coef="dvcat55P")
)

priors_m2 <- c(
  priors_m1,
  prior(normal(0,1), class = "sd", group = "yearVeh")
)

# Model 1: Non-Hierarchical Logistic Model
m1_nonhier <- brm(
  formula = y ~ 1 + seatbelt + airbag + dvcat,
  data    = data_clean,
  family  = bernoulli(link = "logit"),
  prior   = priors_m1,
  chains  = 4,
  iter    = 3000,
  warmup  = 1000,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  cores   = 4,
  seed    = 42
)

# Model 2: Hierarchical Logistic Model (Varying Intercept by Vehicle Year)
m2_hier <- brm(
  formula = y ~ 1 + seatbelt + airbag + dvcat + (1 | yearVeh),
  data    = data_clean,
  family  = bernoulli(link = "logit"),
  prior   = priors_m2,
  chains  = 4,
  iter    = 3000,
  warmup  = 1000,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  cores   = 4,
  seed    = 42
)

# Summary checks for R-hat, ESS, and Divergences
summary(m1_nonhier)
summary(m2_hier)

# Extract MCMC diagnostic summary
np_m2 <- nuts_params(m2_hier$fit)
num_divergences <- sum(subset(np_m2, Parameter == "divergent__")$Value)
cat("Number of divergent transitions (Model 2):", num_divergences, "\n")

# Trace plot inspection for key parameters
mcmc_trace(m2_hier, pars = c("b_Intercept", "b_seatbeltbelted", "b_airbagairbag", "sd_yearVeh__Intercept"))

# Posterior predictive checks: Distribution of predicted proportions
pp_check(m1_nonhier, type = "dens_overlay", ndraws = 100) + ggtitle("PPC: Model 1 (Non-Hierarchical)")
pp_check(m2_hier, type = "dens_overlay", ndraws = 100) + ggtitle("PPC: Model 2 (Hierarchical)")

# Grouped PPC by vehicle year to inspect group-level fit
pp_check(m2_hier, type = "stat_grouped", group = "yearVeh", stat = "mean")

priors_weak <- c(
  prior(normal(0, 10), class = "Intercept"),
  prior(normal(0, 5), class = "b"),
  prior(student_t(3, 0, 2.5), class = "sd", group = "yearVeh")
)

m2_sensitivity <- brm(
  formula = y ~ 1 + seatbelt + airbag + dvcat + (1 | yearVeh),
  data    = data_clean,
  family  = bernoulli(link = "logit"),
  prior   = priors_weak,
  chains  = 4,
  iter    = 3000,
  warmup  = 1000,
  control = list(adapt_delta = 0.95),
  cores   = 4,
  seed    = 42
)

# Compare posterior estimates between informative and weak priors
posterior_comp <- bind_rows(
  as_draws_df(m2_hier) %>% mutate(Prior = "Informative"),
  as_draws_df(m2_sensitivity) %>% mutate(Prior = "Weakly Informative")
)

mcmc_areas(posterior_comp, pars = c("b_seatbeltbelted", "b_airbagairbag"), prob = 0.95)
