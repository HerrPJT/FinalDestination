# Model comparison (PSIS-LOO-CV) for the BDA project, from the saved fits.
# Run from the folder that contains fits/:  Rscript model_comparison_fractran.R

library(brms)
library(loo)

m1_nonhier <- readRDS("fits/m1_nonhier.rds")
m2_hier    <- readRDS("fits/m2_hier.rds")
y_obs      <- m1_nonhier$data$y

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

