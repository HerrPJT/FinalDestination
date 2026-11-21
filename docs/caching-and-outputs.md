# Caching and outputs

The analysis script [`test_project_BDA_code.R`](../test_project_BDA_code.R) saves slow results to `fits/` and loads them on later runs instead of recomputing them. It also writes plots and tables for each section to `outputs/`.

## How to rerun something

**To rerun a step, delete its file in `fits/` and run the script again.** Anything not deleted is loaded from disk.

| To rerun... | Delete |
|---|---|
| Model 1 (informative priors) | `fits/m1_nonhier.rds` |
| Model 2 (informative priors) | `fits/m2_hier.rds` |
| Model 1 (weak priors, sensitivity) | `fits/m1_sensitivity.rds` |
| Model 2 (weak priors, sensitivity) | `fits/m2_sensitivity.rds` |
| LOO for Model 1 / Model 2 | `fits/loo_m1.rds` / `fits/loo_m2.rds` |
| Predicted probabilities for Model 1 / Model 2 (used by the classification metrics) | `fits/probs_m1.rds` / `fits/probs_m2.rds` |
| Everything | the whole `fits/` folder |

> **Important: after a main model refits, delete its LOO and probability files too.**
> The LOO and probability caches don't know when their model changed. If Model 1 refits (for example after a prior change), delete `fits/loo_m1.rds` and `fits/probs_m1.rds` as well, or the comparison and metrics will come from the old fit. The same goes for Model 2 with `loo_m2.rds` and `probs_m2.rds`.

From an R console or notebook cell:

```r
file.remove("fits/m2_hier.rds", "fits/loo_m2.rds", "fits/probs_m2.rds")  # rerun Model 2 and everything that depends on it
unlink("fits", recursive = TRUE)                                         # rerun everything
```

Plots and tables in `outputs/` are recreated on every run, so they never need deleting.

## When the models refit on their own

The four model fits use `brms`'s built-in caching (`brm(file = ..., file_refit = "on_change")`). A model refits automatically when its **data, formula, priors or sampler settings** (chains, iterations, `adapt_delta`, seed, ...) change. Otherwise the saved fit is loaded in seconds.

The LOO and probability files use a simpler helper, `cached()`, that only checks whether the file exists. That's why they need deleting by hand (see the warning above).

## Where files go

Both folders are created relative to R's **working directory**, so run the script from the `FinalDestination/` folder. In JupyterHub the working directory is the folder the notebook is in; check it with `getwd()`.

### `fits/`: cached results (not in git)

`.rds` files that can be tens of MB each. The folder is in `.gitignore`. To save teammates a long first run, share the folder through JupyterHub rather than git.

### `outputs/`: results to look at (in git)

Files are prefixed by section number:

| File | Contents |
|---|---|
| `00_session_info.txt` | R and package versions |
| `01_data_summary.txt` | Row counts, death rate, predictor summaries |
| `01_observed_death_rates.csv` | Raw death rate by each predictor level and vehicle year |
| `03_summary_m1.txt`, `03_summary_m2.txt` | `brms` model summaries |
| `03_diagnostics.csv` | Per model: worst R-hat, worst ESS ratio, divergences, treedepth hits |
| `03_trace_m2.png` | Trace plots for key Model 2 parameters |
| `03_fixed_effects.csv` | Coefficient estimates for both models |
| `03_yearVeh_intercepts_m2.csv` | Vehicle-year intercepts from Model 2 |
| `04_ppc_m1.png`, `04_ppc_m2.png` | Posterior predictive checks |
| `05_diagnostics_sensitivity.csv` | Diagnostics for the weak-prior models |
| `05_sensitivity_table.csv` | Posterior means and SDs, informative vs. weak priors |
| `05_sensitivity_m2.png` | Density comparison for key Model 2 coefficients |
| `06_loo.txt`, `06_loo_compare.csv` | PSIS-LOO results and model comparison |
| `06_classification_metrics.csv` | Brier score, accuracy, sensitivity, specificity, precision |
| `07_odds_ratios_m2.csv` | Model 2 odds ratios with 95% credible intervals |
