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

# ==============================================================================
# 1. POSTERIOR PREDICTIVE CHECKS (Binary-appropriate checks for both models)
# ==============================================================================

pp_check(m1_nonhier, type = "dens_overlay", ndraws = 100) +
  ggtitle("PPC Bars: Model 1 (Non-Hierarchical)")

pp_check(m2_hier, type = "dens_overlay", ndraws = 100) +
  ggtitle("PPC Bars: Model 2 (Hierarchical)")


# ==============================================================================
# 2. SENSITIVITY ANALYSIS (Both Model 1 and Model 2)
# ==============================================================================

priors_weak_m1 <- c(
  prior(normal(0, 10), class = "Intercept"),
  prior(normal(0, 5), class = "b")
)

priors_weak_m2 <- c(
  priors_weak_m1,
  prior(student_t(3, 0, 2.5), class = "sd", group = "yearVeh")
)

m1_sensitivity <- brm(
  formula = y ~ 1 + seatbelt + airbag + dvcat,
  data    = data_clean,
  family  = bernoulli(link = "logit"),
  prior   = priors_weak_m1,
  chains  = 4,
  iter    = 3000,
  warmup  = 1000,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  cores   = 4,
  seed    = 42
)

m2_sensitivity <- brm(
  formula = y ~ 1 + seatbelt + airbag + dvcat + (1 | yearVeh),
  data    = data_clean,
  family  = bernoulli(link = "logit"),
  prior   = priors_weak_m2,
  chains  = 4,
  iter    = 3000,
  warmup  = 1000,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  cores   = 4,
  seed    = 42
)

# Numerical comparison of posterior summaries (Informative vs. Weakly Informative)
compare_posteriors <- function(fit_info, fit_weak, model_label) {
  s_info <- as_draws_df(fit_info) %>%
    select(starts_with("b_"), any_of("sd_yearVeh__Intercept")) %>%
    pivot_longer(everything(), names_to = "param", values_to = "val") %>%
    group_by(param) %>%
    summarise(info_mean = mean(val), info_sd = sd(val), .groups = "drop")
  
  s_weak <- as_draws_df(fit_weak) %>%
    select(starts_with("b_"), any_of("sd_yearVeh__Intercept")) %>%
    pivot_longer(everything(), names_to = "param", values_to = "val") %>%
    group_by(param) %>%
    summarise(weak_mean = mean(val), weak_sd = sd(val), .groups = "drop")
  
  left_join(s_info, s_weak, by = "param") %>%
    mutate(
      model = model_label,
      diff_mean = weak_mean - info_mean
    )
}

sens_table <- bind_rows(
  compare_posteriors(m1_nonhier, m1_sensitivity, "Model 1"),
  compare_posteriors(m2_hier, m2_sensitivity, "Model 2")
)
print(sens_table, n = Inf)

# Density overlay plot for key fixed effects (Model 2)
sens_draws_m2 <- bind_rows(
  as_draws_df(m2_hier) %>% select(b_seatbeltbelted, b_airbagairbag, b_dvcat55P) %>% mutate(Prior = "Informative"),
  as_draws_df(m2_sensitivity) %>% select(b_seatbeltbelted, b_airbagairbag, b_dvcat55P) %>% mutate(Prior = "Weakly Informative")
) %>%
  pivot_longer(-Prior, names_to = "Parameter", values_to = "Draw")

ggplot(sens_draws_m2, aes(x = Draw, fill = Prior)) +
  geom_density(alpha = 0.45) +
  facet_wrap(~ Parameter, scales = "free") +
  labs(title = "Model 2 Prior Sensitivity: Informative vs. Weakly Informative", x = "Log-Odds", y = "Density") +
  theme_minimal()


# ==============================================================================
# 3. MODEL COMPARISON (PSIS-LOO-CV + Practical Predictive Metrics)
# ==============================================================================

# Compute PSIS-LOO (use moment_match = TRUE if any Pareto k > 0.7)
loo_m1 <- loo(m1_nonhier, cores = 4)
loo_m2 <- loo(m2_hier, cores = 4)

print(loo_m1)
print(loo_m2)

comp_loo <- loo_compare(loo_m1, loo_m2)
print(comp_loo, simplify = FALSE)

# Practical classification metrics (e.g., threshold = 0.05 and 0.50)
eval_classification <- function(model, data, thresholds = c(0.05, 0.50)) {
  probs <- colMeans(posterior_epred(model, newdata = data))
  brier <- mean((probs - data$y)^2)
  cat(sprintf("Brier Score: %.4f\n", brier))
  
  for (th in thresholds) {
    pred_cls <- if_else(probs >= th, 1, 0)
    tp <- sum(pred_cls == 1 & data$y == 1)
    tn <- sum(pred_cls == 0 & data$y == 0)
    fp <- sum(pred_cls == 1 & data$y == 0)
    fn <- sum(pred_cls == 0 & data$y == 1)
    
    acc  <- (tp + tn) / length(data$y)
    sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    
    cat(sprintf("Threshold = %.2f | Acc: %.2f%% | Sens (Recall): %.2f%% | Spec: %.2f%% | Prec: %.2f%%\n",
                th, acc * 100, sens * 100, spec * 100, prec * 100))
  }
}

cat("--- Model 1 Practical Metrics ---\n")
eval_classification(m1_nonhier, data_clean)

cat("\n--- Model 2 Practical Metrics ---\n")
eval_classification(m2_hier, data_clean)


# ==============================================================================
# 4. HELPER FOR CONCLUSION / ODDS RATIOS
# ==============================================================================

# Posterior mean and 95% credible interval on Odds Ratio scale (Model 2)
or_summary_m2 <- as_draws_df(m2_hier) %>%
  select(starts_with("b_")) %>%
  mutate(across(everything(), exp)) %>%
  pivot_longer(everything(), names_to = "Parameter", values_to = "OR") %>%
  group_by(Parameter) %>%
  summarise(
    OR_mean = mean(OR),
    OR_l95  = quantile(OR, 0.025),
    OR_u95  = quantile(OR, 0.975),
    .groups = "drop"
  )
print(or_summary_m2)
