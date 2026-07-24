## ---------------------------------------------------------------
## Tobit (left-censored) regression of nitrogen on huc12-level
## covariates, clustered on huc10id
## ---------------------------------------------------------------

library(survival)
library(dplyr)
library(car)

df <- read.csv("waterdata.csv", stringsAsFactors = FALSE)

## ---------------------------------------------------------------
## 1. Build the censored outcome
## ---------------------------------------------------------------
## nitrogenmgl == 0 whenever a value was below the detection limit

DL <- 0.01   # reporting limit for NO2+NO3 (mg/L)

df <- df %>%
  mutate(
    n_censored = no2no3mgl == 0,
    n_value    = ifelse(n_censored, DL, no2no3mgl),
    n_status   = ifelse(n_censored, 0, 1),
    log_cafo_load = log(huc12_weighted_cafo_load + 4.33e-05 / 2)
  )


## sanity check
cat("\nNumber of censored observations:\n")
table(df$n_censored)

## ---------------------------------------------------------------
## 2. huc12-level predictors
## ---------------------------------------------------------------
## Final predictor set selected using VIF, cluster-robust inference,
## likelihood-ratio tests, AIC, and ecological interpretability.
##
## huc12_pct_cropland was removed because it did not improve model fit
## (LR p = 0.992; AIC improved after removal).
##
## log_cafo_load was retained despite a non-significant Wald test because
## likelihood-ratio testing showed that removing it significantly worsened
## model fit.
##
## CAFO density, wastewater density, landfill density, and population
## density remained stable and were retained in the final model.

predictors <- c(
  "huc12_upstream_cafo_density",
  "huc12_upstream_wastewater_densit",
  "huc12_upstream_landfill_density",
  "log_cafo_load",
  "huc12_pop_density"
)

## huc12s are nested inside huc10s,
## so clustering on the coarser unit captures correlation from both
## levels without discarding the identifying variation.
## check missingness before fitting

sapply(df[predictors], function(x) sum(is.na(x)))

model_df <- df %>%
  select(huc12id, huc10id, n_value, n_status, all_of(predictors)) %>%
  na.omit()

## ---------------------------------------------------------------
## 3. Check predictor collinearity (VIF)
## ---------------------------------------------------------------
## Rule of thumb: VIF > 5 worth scrutinizing, > 10 a real problem.

vif_check <- lm(as.formula(paste("n_value ~", paste(predictors, collapse = " + "))),
                data = model_df)
vif(vif_check)

## ---------------------------------------------------------------
## 4. Fit the tobit model
## ---------------------------------------------------------------
## dist = "lognormal" fits the model on the log scale internally,
## suited to the strong right-skew typical of nutrient concentrations,
## while Surv()/survreg() correctly handles the left-censoring on the
## original (untransformed) scale.

form <- as.formula(
  paste("Surv(n_value, n_status, type = 'left') ~",
        paste(predictors, collapse = " + "),
        "+ cluster(huc10id)")
)

fit <- survreg(form, data = model_df, dist = "lognormal")
summary(fit)

## ---------------------------------------------------------------
## 5. Interpreting coefficients
## ---------------------------------------------------------------
## Because dist = "lognormal", coefficients are on the log scale.
## exp(coef) gives the multiplicative change in the (latent) median
## nitrogen concentration per unit increase in the predictor.

exp(coef(fit))

## Scale parameter (sigma) is the log-scale residual SD -- printed
## automatically in summary(fit) as "Scale".

## ---------------------------------------------------------------
## 4c. Results table (exponentiated coefficients + 95% CI)
## ---------------------------------------------------------------
## Uses the robust (cluster-corrected) SEs -- summary(fit)$table's
## second column, not the naive one -- for the z/p already shown

tab <- summary(fit)$table
tab <- tab[rownames(tab) != "Log(scale)", , drop = FALSE]  # drop scale row; not a covariate effect

est <- tab[, 1]
se  <- tab[, 2]                       # robust (cluster) SE
zval <- tab[, ncol(tab) - 1]
pval <- tab[, ncol(tab)]

z975 <- qnorm(0.975)

results <- data.frame(
  term          = rownames(tab),
  estimate      = round(est, 4),
  std_error     = round(se, 4),
  z             = round(zval, 2),
  p_value       = signif(pval, 3),
  exp_estimate  = round(exp(est), 3),
  exp_ci_low    = round(exp(est - z975 * se), 3),
  exp_ci_high   = round(exp(est + z975 * se), 3),
  row.names     = NULL
)

results

## ---------------------------------------------------------------
## 4d. Rescaled version for reporting
## ---------------------------------------------------------------
## Raw per-1-unit effects are unusable for cafo_density/landfill_density/
## pop_density -- their "1 unit" is far outside the observed range (e.g.
## landfill_density ranges 0-0.10 with SD 0.0135, so exp(coef * 1) implies
## a jump ~7x larger than the entire observed range). Rescale each
## predictor to a step size that's actually meaningful to read:
##   cafo_density:     per +0.10  (range 0-0.97, SD 0.15)
##   landfill_density: per +0.01  (range 0-0.10, SD 0.014)
##   pop_density:        per +100  (persons/km^2; range 3-797, SD 143)
## intercept is left as-is (not a slope, rescaling doesn't apply).

scale_map <- c(
  "(Intercept)"                      = 1,
  "huc12_upstream_cafo_density"      = 0.10,
  "huc12_upstream_wastewater_densit" = 0.10,
  "huc12_upstream_landfill_density"  = 0.01,
  "log_cafo_load"                    = 1,
  "huc12_pop_density"                = 100
)
results_scaled <- results
results_scaled$unit_step   <- scale_map[results_scaled$term]
results_scaled$exp_estimate <- round(exp(est * results_scaled$unit_step), 3)
results_scaled$exp_ci_low   <- round(exp((est - z975 * se) * results_scaled$unit_step), 3)
results_scaled$exp_ci_high  <- round(exp((est + z975 * se) * results_scaled$unit_step), 3)

## drop the raw log-scale columns for the report version -- keep term,
## the step size so readers know what "the effect" is per, and the
## exponentiated effect + CI
results_scaled <- results_scaled[, c(
  "term",
  "unit_step",
  "p_value",
  "exp_estimate",
  "exp_ci_low",
  "exp_ci_high"
)]
results_scaled

## ---------------------------------------------------------------
## 4e. Cleaned-up version: clean labels, intercept dropped
## ---------------------------------------------------------------

label_map <- c(
  "huc12_upstream_cafo_density"     = "CAFO density (per +0.1 CAFOs/km2)",
  "huc12_upstream_wastewater_densit" = "Wastewater density (per +0.1 facilities/km2)",
  "huc12_upstream_landfill_density" = "Landfill density (per +0.01 landfills/km2)",
  "log_cafo_load"                    = "Log CAFO load (per +1 log unit)",
  "huc12_pop_density"                 = "Population density (per +100 people/km2)"
)

results_report <- results_scaled[results_scaled$term != "(Intercept)", ]
results_report$term <- label_map[results_report$term]
results_report$pct_change    <- round((results_report$exp_estimate - 1) * 100, 1)
results_report$pct_ci_low    <- round((results_report$exp_ci_low - 1) * 100, 1)
results_report$pct_ci_high   <- round((results_report$exp_ci_high - 1) * 100, 1)
results_report$significant  <-  ifelse(results_report$p_value < 0.05, "Yes", "No")
results_report <- results_report[, c("term", "pct_change", "pct_ci_low",
                                     "pct_ci_high", "p_value", "significant")]
names(results_report) <- c("Predictor", "% change in median nitrogen",
                           "95% CI low", "95% CI high", "p-value", "Significant (p<.05)")
results_report

## ---------------------------------------------------------------
## FIGURES
## ---------------------------------------------------------------

library(ggplot2)
library(scales)

##### PREDICTOR EFFECTS CIs #####
forest_df <- results_scaled %>%
  filter(term != "(Intercept)") %>%
  mutate(
    Predictor = factor(
      label_map[term],
      levels = rev(c(
        "CAFO density (per +0.1 facilities/km2)",
        "Wastewater density (per +0.1 facilities/km2)",
        "Landfill density (per +0.01 landfills/km2)",
        "Log CAFO load (per +1 log unit)",
        "Population density (per +100 people/km2)"
      ))
    ),
    pct_change  = (exp_estimate - 1) * 100,
    pct_ci_low  = (exp_ci_low - 1) * 100,
    pct_ci_high = (exp_ci_high - 1) * 100,
    Significant = p_value < 0.05,
    Label = sprintf("%.1f%%", pct_change)
  )

ggplot(forest_df,
       aes(x = pct_change,
           y = Predictor,
           colour = Significant)) +
  
  geom_vline(xintercept = 0,
             linewidth = 0.6,
             linetype = "dashed",
             colour = "grey50") +
  
  geom_errorbar(
    aes(xmin = pct_ci_low,
        xmax = pct_ci_high),
    orientation = "y",
    width = 0.18,
    linewidth = 0.8
  ) +
  
  geom_point(size = 3.8) +
  
  geom_text(
    aes(label = Label),
    nudge_x = 2,
    colour = "black",
    size = 3.8
  ) +
  
  scale_colour_manual(
    values = c(
      "TRUE" = "#1F78B4",
      "FALSE" = "grey60"
    ),
    guide = "none"
  ) +
  
  labs(
    title = "Effects of watershed characteristics on NO2 + NO3 concentrations",
    subtitle = "Points show estimated percent change in median concentration; bars indicate 95% confidence intervals",
    x = "Percent change in median NO2 + NO3 concentration",
    y = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(face = "bold")
  )

##### OBSERVED VS PREDICTED #####

model_df$predicted <- predict(
  fit,
  type = "response"
)

ggplot(model_df,
       aes(predicted,
           n_value)) +
  
  geom_point(alpha = 0.35) +
  
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed",
              color = "red") +
  
  labs(
    x = "Predicted nitrogen",
    y = "Observed nitrogen",
    title = "Observed versus predicted nitrogen"
  ) +
  
  theme_bw(base_size = 13)

##### PARTIAL EFFECTS #####
library(ggeffects)
library(patchwork)

p1 <- ggplot(ggpredict(fit,
                       terms = "huc12_upstream_cafo_density"),
             aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.2) +
  theme_bw() +
  labs(x="CAFO density",
       y="Predicted N")

p2 <- ggplot(ggpredict(fit,
                       terms = "huc12_upstream_wastewater_densit"),
             aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.2) +
  theme_bw() +
  labs(x="Wastewater density",
       y="")

p3 <- ggplot(ggpredict(fit,
                       terms = "huc12_upstream_landfill_density"),
             aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.2) +
  theme_bw() +
  labs(x="Landfill density",
       y="Predicted N")

p4 <- ggplot(ggpredict(fit,
                       terms = "huc12_pop_density"),
             aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.2) +
  theme_bw() +
  labs(x="Population density",
       y="")

p5 <- ggplot(ggpredict(fit,
                       terms = "log_cafo_load"),
             aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.2) +
  theme_bw() +
  labs(x="Log(CAFO load)",
       y="Predicted N")

(p1 | p2) /
  (p3 | p4) /
  (p5 | plot_spacer())
