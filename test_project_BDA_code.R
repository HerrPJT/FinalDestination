library(tidyverse)
library(brms)
library(bayesplot)
library(loo)

# ==============================================================================
# 0. CHECKPOINTS AND OUTPUTS
# ==============================================================================
# fits/    : cached model fits, LOO objects and predictions (.rds). Large files;
#            not tracked in git. Delete a file to force that step to rerun.
# outputs/ : plots (.png), tables (.csv) and printed summaries (.txt) to inspect.
#
# brm(file = ..., file_refit = "on_change") loads the saved fit if it exists and
# only refits when the data, formula, priors or sampler settings change.

fit_dir <- "fits"
out_dir <- "outputs"
dir.create(fit_dir, showWarnings = FALSE)
dir.create(out_dir, showWarnings = FALSE)

fit_path <- function(name) file.path(fit_dir, name)
out_path <- function(name) file.path(out_dir, name)

# Run `expr` once and cache its result; later runs load it from disk.
cached <- function(name, expr) {
  path <- fit_path(paste0(name, ".rds"))
  if (file.exists(path)) {
    message("Loading cached ", path)
    return(readRDS(path))
  }
  res <- expr
  saveRDS(res, path)
  res
}

save_plot <- function(plot, name, width = 8, height = 5) {
  ggsave(out_path(name), plot, width = width, height = height, dpi = 150)
  print(plot)
}

save_text <- function(name, expr) {
  txt <- capture.output(expr)
  writeLines(txt, out_path(name))
  cat(txt, sep = "\n")
}

# ==============================================================================
# 1. LOAD AND CLEAN DATA
# ==============================================================================

data_src <- if (file.exists("nassCDS.csv")) "nassCDS.csv" else
  "https://vincentarelbundock.github.io/Rdatasets/csv/DAAG/nassCDS.csv"

data_raw <- read.csv(data_src)

# Drops 1 row (out of 26,217) with missing yearVeh. Applied to both models so
# loo_compare() compares them on identical observations.
data_clean <- data_raw %>%
  drop_na(yearVeh) %>%
  mutate(
    y = if_else(dead == "dead", 1, 0),
    seatbelt = factor(seatbelt, levels = c("none", "belted")),
    airbag   = factor(airbag, levels = c("none", "airbag")),
    dvcat    = factor(dvcat, levels = c("1-9km/h", "10-24", "25-39", "40-54", "55+")),
    yearVeh  = as.factor(yearVeh)
  )

save_text("01_data_summary.txt", {
  cat("Rows raw:", nrow(data_raw), " | Rows clean:", nrow(data_clean), "\n")
  cat("Death rate:", mean(data_clean$y), "\n\n")
  print(summary(data_clean[, c("y", "seatbelt", "airbag", "dvcat")]))
})

# Observed death rate by each predictor (useful for interpreting coefficients)
obs_rates <- bind_rows(
  data_clean %>% group_by(level = seatbelt) %>% summarise(n = n(), deaths = sum(y)) %>% mutate(variable = "seatbelt"),
  data_clean %>% group_by(level = airbag)   %>% summarise(n = n(), deaths = sum(y)) %>% mutate(variable = "airbag"),
  data_clean %>% group_by(level = dvcat)    %>% summarise(n = n(), deaths = sum(y)) %>% mutate(variable = "dvcat"),
  data_clean %>% group_by(level = yearVeh)  %>% summarise(n = n(), deaths = sum(y)) %>% mutate(variable = "yearVeh")
) %>%
  mutate(level = as.character(level), death_rate = deaths / n) %>%
  select(variable, level, n, deaths, death_rate)
write_csv(obs_rates, out_path("01_observed_death_rates.csv"))

# ==============================================================================
# 2. PRIORS AND MODEL FITS
# ==============================================================================

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
  seed    = 42,
  file    = fit_path("m1_nonhier"),
  file_refit = "on_change"
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
  seed    = 42,
  file    = fit_path("m2_hier"),
  file_refit = "on_change"
)

# ==============================================================================
# 3. MCMC DIAGNOSTICS
# ==============================================================================

save_text("03_summary_m1.txt", print(summary(m1_nonhier)))
save_text("03_summary_m2.txt", print(summary(m2_hier)))

# One row per model: worst R-hat, worst ESS ratio, divergences, treedepth hits
diagnostics_row <- function(fit, label, max_treedepth = 12) {
  np <- nuts_params(fit)
  tibble(
    model            = label,
    max_rhat         = max(rhat(fit), na.rm = TRUE),
    min_neff_ratio   = min(neff_ratio(fit), na.rm = TRUE),
    divergences      = sum(subset(np, Parameter == "divergent__")$Value),
    treedepth_hits   = sum(subset(np, Parameter == "treedepth__")$Value >= max_treedepth)
  )
}

diag_main <- bind_rows(
  diagnostics_row(m1_nonhier, "Model 1"),
  diagnostics_row(m2_hier, "Model 2")
)
write_csv(diag_main, out_path("03_diagnostics.csv"))
print(diag_main)

# Trace plot inspection for key parameters
p_trace <- mcmc_trace(m2_hier, pars = c("b_Intercept", "b_seatbeltbelted", "b_airbagairbag", "sd_yearVeh__Intercept"))
save_plot(p_trace, "03_trace_m2.png", width = 10, height = 6)

# Fixed-effect estimates for both models, in one table
fixef_table <- bind_rows(
  as_tibble(fixef(m1_nonhier), rownames = "param") %>% mutate(model = "Model 1"),
  as_tibble(fixef(m2_hier), rownames = "param") %>% mutate(model = "Model 2")
)
write_csv(fixef_table, out_path("03_fixed_effects.csv"))

# Vehicle-year intercepts (Model 2)
yearveh_effects <- as_tibble(ranef(m2_hier)$yearVeh[, , "Intercept"], rownames = "yearVeh")
write_csv(yearveh_effects, out_path("03_yearVeh_intercepts_m2.csv"))

# ==============================================================================
# 4. POSTERIOR PREDICTIVE CHECKS
# ==============================================================================

p_ppc_m1 <- pp_check(m1_nonhier, type = "dens_overlay", ndraws = 100) +
  ggtitle("PPC Bars: Model 1 (Non-Hierarchical)")
save_plot(p_ppc_m1, "04_ppc_m1.png")

p_ppc_m2 <- pp_check(m2_hier, type = "dens_overlay", ndraws = 100) +
  ggtitle("PPC Bars: Model 2 (Hierarchical)")
save_plot(p_ppc_m2, "04_ppc_m2.png")

# ==============================================================================
# 5. SENSITIVITY ANALYSIS (Both Model 1 and Model 2)
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
  seed    = 42,
  file    = fit_path("m1_sensitivity"),
  file_refit = "on_change"
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
  seed    = 42,
  file    = fit_path("m2_sensitivity"),
  file_refit = "on_change"
)

diag_sens <- bind_rows(
  diagnostics_row(m1_sensitivity, "Model 1 (weak priors)"),
  diagnostics_row(m2_sensitivity, "Model 2 (weak priors)")
)
write_csv(diag_sens, out_path("05_diagnostics_sensitivity.csv"))
print(diag_sens)

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
write_csv(sens_table, out_path("05_sensitivity_table.csv"))
print(sens_table, n = Inf)

# Density overlay plot for key fixed effects (Model 2)
sens_draws_m2 <- bind_rows(
  as_draws_df(m2_hier) %>% select(b_seatbeltbelted, b_airbagairbag, b_dvcat55P) %>% mutate(Prior = "Informative"),
  as_draws_df(m2_sensitivity) %>% select(b_seatbeltbelted, b_airbagairbag, b_dvcat55P) %>% mutate(Prior = "Weakly Informative")
) %>%
  pivot_longer(-Prior, names_to = "Parameter", values_to = "Draw")

p_sens <- ggplot(sens_draws_m2, aes(x = Draw, fill = Prior)) +
  geom_density(alpha = 0.45) +
  facet_wrap(~ Parameter, scales = "free") +
  labs(title = "Model 2 Prior Sensitivity: Informative vs. Weakly Informative", x = "Log-Odds", y = "Density") +
  theme_minimal()
save_plot(p_sens, "05_sensitivity_m2.png", width = 10, height = 4)

# ==============================================================================
# 6. MODEL COMPARISON (PSIS-LOO-CV + Practical Predictive Metrics)
# ==============================================================================

# Compute PSIS-LOO (cached: this is one of the slow steps)
# (moment_match = TRUE would also need save_pars = save_pars(all = TRUE) in brm)
loo_m1 <- cached("loo_m1", loo(m1_nonhier, cores = 4))
loo_m2 <- cached("loo_m2", loo(m2_hier, cores = 4))

comp_loo <- loo_compare(loo_m1, loo_m2)

save_text("06_loo.txt", {
  cat("=== Model 1 ===\n"); print(loo_m1)
  cat("\n=== Model 2 ===\n"); print(loo_m2)
  cat("\n=== Comparison ===\n"); print(comp_loo, simplify = FALSE)
})
write_csv(as_tibble(comp_loo, rownames = "model"), out_path("06_loo_compare.csv"))

# Posterior mean predicted probability per observation (cached: memory-heavy)
probs_m1 <- cached("probs_m1", colMeans(posterior_epred(m1_nonhier, newdata = data_clean)))
probs_m2 <- cached("probs_m2", colMeans(posterior_epred(m2_hier, newdata = data_clean)))

# Practical classification metrics (e.g., threshold = 0.05 and 0.50)
eval_classification <- function(probs, y, model_label, thresholds = c(0.05, 0.50)) {
  brier <- mean((probs - y)^2)

  map_dfr(thresholds, function(th) {
    pred_cls <- if_else(probs >= th, 1, 0)
    tp <- sum(pred_cls == 1 & y == 1)
    tn <- sum(pred_cls == 0 & y == 0)
    fp <- sum(pred_cls == 1 & y == 0)
    fn <- sum(pred_cls == 0 & y == 1)

    tibble(
      model     = model_label,
      threshold = th,
      brier     = brier,
      tp = tp, tn = tn, fp = fp, fn = fn,
      accuracy    = (tp + tn) / length(y),
      sensitivity = if ((tp + fn) > 0) tp / (tp + fn) else NA_real_,
      specificity = if ((tn + fp) > 0) tn / (tn + fp) else NA_real_,
      precision   = if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    )
  })
}

class_metrics <- bind_rows(
  eval_classification(probs_m1, data_clean$y, "Model 1"),
  eval_classification(probs_m2, data_clean$y, "Model 2")
)
write_csv(class_metrics, out_path("06_classification_metrics.csv"))
print(class_metrics, width = Inf)

# ==============================================================================
# 7. HELPER FOR CONCLUSION / ODDS RATIOS
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
write_csv(or_summary_m2, out_path("07_odds_ratios_m2.csv"))
print(or_summary_m2)

save_text("00_session_info.txt", print(sessionInfo()))
