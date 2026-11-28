# LIBRARIES
library(tidyverse)
library(brms)
library(bayesplot)
library(loo)
library(reliabilitydiag)
library(ggplot2)

# DATA LOADING AND CLEANING
data_clean <- read.csv("https://vincentarelbundock.github.io/Rdatasets/csv/DAAG/nassCDS.csv") %>%
  drop_na(yearVeh) %>%
  mutate(
    y = if_else(dead == "dead", 1, 0),
    seatbelt = factor(seatbelt, levels = c("none", "belted")),
    airbag   = factor(airbag, levels = c("none", "airbag")),
    dvcat    = factor(dvcat, levels = c("1-9km/h", "10-24", "25-39", "40-54", "55+")),
    yearVeh  = as.factor(yearVeh)
  )
  
#----------------------------------------------------------
# PRIORS
#----------------------------------------------------------
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

# MODELS
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

# Save brmsfit to RDS file
saveRDS(m1_nonhier, file = "notebooks/BDA Project/m1_nonhier.rds")
saveRDS(m2_hier, file = "notebooks/BDA Project/m2_hier.rds")
# Read brmsfit from RDS file
m1_nonhier <- readRDS("notebooks/BDA Project/m1_nonhier.rds")
m2_hier <- readRDS("notebooks/BDA Project/m2_hier.rds")

#----------------------------------------------------------
# CONVERGENCE DIAGNOSTICS
#-----------------------------------------------------------
# Summary checks for R-hat, ESS, and Divergences
summary(m1_nonhier)
summary(m2_hier)

# Extract MCMC diagnostic summary
np_m2 <- nuts_params(m2_hier$fit)
num_divergences <- sum(subset(np_m2, Parameter == "divergent__")$Value)
cat("Number of divergent transitions (Model 2):", num_divergences, "\n")

# Trace plot inspection for key parameters
mcmc_trace(m2_hier, pars = c("b_Intercept", "b_seatbeltbelted", "b_airbagairbag", "sd_yearVeh__Intercept"))

#-------------------------------------------------------
# POSTERIOR PREDICTIVE CHECKS
#-------------------------------------------------------

# Model 1: 
p_epred1 <- posterior_epred(
  m1_nonhier,
  newdata = data_clean,
  ndraws  = 1000         # moderate number of draws to avoid memory issues
)
p_hat1 <- colMeans(p_epred1)
y_obs <- data_clean$y
rd1 <- reliabilitydiag(
  EMOS = p_hat1,
  y    = y_obs
)

autoplot(rd1) +
  labs(
    x = "Predicted probability of fatality",
    y = "Conditional event probabilities"
  ) +
  ggtitle("Calibration (PAV) – Model 1: Non-hierarchical (subset)")

# Model 2:
p_epred2 <- posterior_epred(
  m2_hier,
  newdata = data_clean,
  ndraws  = 1000
)
p_hat2 <- colMeans(p_epred2)
rd2 <- reliabilitydiag(
  EMOS = p_hat2,
  y    = y_obs
)
autoplot(rd2) +
  labs(
    x = "Predicted probability of fatality",
    y = "Conditional event probabilities"
  ) +
  ggtitle("Calibration (PAV) – Model 2: Hierarchical (subset)")

#--------------------------------------------------------------------
# PRIOR AND LIKELIHOOD SENSITIVITY
#--------------------------------------------------------------------

main_vars_m1 <- c(
  "b_Intercept", "b_seatbeltbelted", "b_airbagairbag",
  "b_dvcat10M24", "b_dvcat25M39", "b_dvcat40M54", "b_dvcat55P"
)

main_vars_m2 <- c(main_vars_m1, "sd_yearVeh__Intercept")

set.seed(42)
sub_idx <- sample(seq_len(nrow(data_clean)), size = 4000)
data_ps <- data_clean[sub_idx, , drop = FALSE]

ps_m1 <- create_priorsense_data(m1_nonhier, newdata = data_ps)
ps_m2 <- create_priorsense_data(m2_hier,    newdata = data_ps)

cat("=== Model 1 Power-Scale Sensitivity (4,000 rows) ===\n")
powerscale_sensitivity(ps_m1, variable = main_vars_m1)

cat("\n=== Model 2 Power-Scale Sensitivity (4,000 rows) ===\n")
powerscale_sensitivity(ps_m2, variable = main_vars_m2)

powerscale_plot_ecdf(
  ps_m1,
  variable = c("b_Intercept", "b_seatbeltbelted", "b_airbagairbag")
)

powerscale_plot_ecdf(
  ps_m2,
  variable = c("b_Intercept", "b_seatbeltbelted", "b_airbagairbag", "sd_yearVeh__Intercept")
)

#--------------------------------------------------------------------
# MODEL COMPARISON
#--------------------------------------------------------------------

# PSIS-LOO-CV for both models
loo_m1 <- loo(m1_nonhier, save_psis = TRUE)
loo_m2 <- loo(m2_hier, save_psis = TRUE)
print(loo_m1)
print(loo_m2)
loo_compare(loo_m1, loo_m2)

# LOO predictive probabilities of fatality (out-of-sample)
p_loo1 <- loo_epred(m1_nonhier, psis_object = loo_m1$psis_object)
p_loo2 <- loo_epred(m2_hier, psis_object = loo_m2$psis_object)

# Accuracy and balanced accuracy at a given classification threshold
classify <- function(p, y, threshold) {
  pred <- as.numeric(p >= threshold)
  sensitivity <- mean(pred[y == 1] == 1)
  specificity <- mean(pred[y == 0] == 0)
  c(accuracy = mean(pred == y),
    balanced_accuracy = (sensitivity + specificity) / 2)
}

for (threshold in c(0.5, 0.05)) {
  cat("\nThreshold", threshold, "\n")
  print(rbind(non_hierarchical = classify(p_loo1, y_obs, threshold),
              hierarchical     = classify(p_loo2, y_obs, threshold)))
}

#--------------------------------------------------------------------
# MARGINAL POSTERIORS
#--------------------------------------------------------------------

fixed_pars <- c(
  "b_Intercept",
  "b_seatbeltbelted",
  "b_airbagairbag",
  "b_dvcat10M24",
  "b_dvcat25M39",
  "b_dvcat40M54",
  "b_dvcat55P"
)

summary_m1_logodds <- posterior_summary(m1_nonhier, variable = fixed_pars)
summary_m2_logodds <- posterior_summary(m2_hier, variable = c(fixed_pars, "sd_yearVeh__Intercept"))

print(round(summary_m1_logodds, 3))
print(round(summary_m2_logodds, 3))

or_m2 <- as_draws_df(m2_hier) %>%
  select(all_of(fixed_pars)) %>%
  mutate(across(everything(), exp)) %>%
  summarise_draws(
    mean = mean,
    median = median,
    sd = sd,
    ~quantile2(.x, probs = c(0.025, 0.975))
  )

print(or_m2, digits = 3)


p_marg_predictors <- mcmc_areas(
  m2_hier,
  pars = c(
    "b_seatbeltbelted",
    "b_airbagairbag",
    "b_dvcat10M24",
    "b_dvcat25M39",
    "b_dvcat40M54",
    "b_dvcat55P"
  ),
  prob = 0.80,
  prob_outer = 0.95
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Marginal Posteriors for Restraints and Crash Speed (m2_hier)",
    subtitle = "80% inner shading and 95% outer interval on log-odds scale",
    x = "Log-odds",
    y = NULL
  ) +
  theme_minimal()

print(p_marg_predictors)

p_marg_sd_year <- mcmc_areas(
  m2_hier,
  pars = "sd_yearVeh__Intercept",
  prob = 0.80,
  prob_outer = 0.95
) +
  labs(
    title = "Marginal Posterior of Vehicle Manufacturing Year SD",
    x = "sd_yearVeh__Intercept (log-odds)",
    y = NULL
  ) +
  theme_minimal()

print(p_marg_sd_year)
