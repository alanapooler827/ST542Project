## ---------------------------------------------------------------
## Tobit (left-censored) regression of phosphorus on huc12-level
## covariates, clustered on huc10id
## ---------------------------------------------------------------

library(survival)
library(dplyr)
library(car)

df <- read.csv("waterdata.csv", stringsAsFactors = FALSE)

## ---------------------------------------------------------------
## 1. Build the censored outcome
## ---------------------------------------------------------------
## phosphorusmgl == 0 whenever a value was below the detection limit

DL <- 0.01  # detection limit for phosphorus

df <- df %>%
  mutate(
    p_censored = phosphorusmgl == 0,                     # TRUE = below detection
    p_value    = ifelse(p_censored, DL, phosphorusmgl),  # value fed to Surv()
    p_status   = ifelse(p_censored, 0, 1)                 # 0 = left-censored, 1 = observed
  )


fit
table(df$p_censored)

## ---------------------------------------------------------------
## 2. huc12-level predictors
## ---------------------------------------------------------------
## Final set, after checking VIF, cluster-robust stability, and an
## AIC/LR comparison:
## - huc12_upstream_cafo_count dropped (VIF 7.9, unstable/counterintuitive
##   sign -- collinear with cafo_density and weighted_cafo_load)
## - huc12_weighted_cafo_load dropped (unstable cluster-robust SE driven
##   by a handful of high-leverage huc10 clusters, even after a log
##   transform its effect wasn't distinguishable from zero once
##   cafo_density was in the model, and dropping it didn't change any
##   other coefficient's story)
## - huc12_upstream_wastewater_densit dropped: weak correlation with
##   every other predictor (all |r| < 0.15, so unlikely to be doing
##   confounding-control work), cluster-robust p in the 0.6-0.9 range
##   across every spec tried, AIC negligibly prefers keeping it (Delta
##   ~1.1, below the usual "meaningfully different" threshold of 2),
##   and even the anti-conservative (non-cluster-corrected) LR test only
##   reaches p=0.08.
## - huc12_pct_cropland KEPT despite its own p-value being non-significant:
##   moderately correlated with both cafo_density (r=0.39) and
##   pop_density (r=-0.36), so it may still be doing confounding-control
##   work, and cropland/nutrient runoff is a substantively motivated
##   pathway worth reporting even as a null result.
## cafo_density and landfill_density are the two variables that stayed
## strong and stable across every specification tried.

predictors <- c(
  "huc12_upstream_cafo_density",
  "huc12_upstream_landfill_density",
  "huc12_pct_cropland",
  "huc12_pop_density"
)

## huc12s are nested inside huc10s,
## so clustering on the coarser unit captures correlation from both
## levels without discarding the identifying variation.
## check missingness before fitting

## check missingness before fitting
sapply(df[predictors], function(x) sum(is.na(x)))

model_df <- df %>%
  select(huc12id, huc10id, p_value, p_status, all_of(predictors)) %>%
  na.omit()

## ---------------------------------------------------------------
## 3. Check predictor collinearity (VIF)
## ---------------------------------------------------------------
## Rule of thumb: VIF > 5 worth scrutinizing, > 10 a real problem.

vif_check <- lm(as.formula(paste("p_value ~", paste(predictors, collapse = " + "))),
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
  paste("Surv(p_value, p_status, type = 'left') ~",
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
## phosphorus concentration per unit increase in the predictor.

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
##   pct_cropland:      per +10   (already in percentage points; range 0-60, SD 11)
##   pop_density:        per +100  (persons/km^2; range 3-797, SD 143)
## intercept is left as-is (not a slope, rescaling doesn't apply).

scale_map <- c(
  "(Intercept)"                      = 1,
  "huc12_upstream_cafo_density"      = 0.10,
  "huc12_upstream_landfill_density"  = 0.01,
  "huc12_pct_cropland"               = 10,
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
results_scaled <- results_scaled[, c("term", "unit_step", "p_value",
                                     "exp_estimate", "exp_ci_low", "exp_ci_high")]
results_scaled

## ---------------------------------------------------------------
## 4e. Cleaned-up version: clean labels, intercept dropped
## ---------------------------------------------------------------
## The intercept (predicted phosphorus when every predictor = 0) isn't
## a meaningful quantity to report here -- "0 CAFOs, 0 landfills, 0%
## cropland, 0 people/km^2 upstream" isn't a real scenario in this
## watershed, and it's not the kind of effect a results table is meant
## to communicate. Drop it for the report version; keep it in `results`
## above if a reviewer wants the full model output.

label_map <- c(
  "huc12_upstream_cafo_density"     = "CAFO density (per +0.1 CAFOs/km2)",
  "huc12_upstream_landfill_density" = "Landfill density (per +0.01 landfills/km2)",
  "huc12_pct_cropland"               = "Cropland (per +10 percentage points)",
  "huc12_pop_density"                 = "Population density (per +100 people/km2)"
)

results_report <- results_scaled[results_scaled$term != "(Intercept)", ]
results_report$term <- label_map[results_report$term]
results_report$pct_change    <- round((results_report$exp_estimate - 1) * 100, 1)
results_report$pct_ci_low    <- round((results_report$exp_ci_low - 1) * 100, 1)
results_report$pct_ci_high   <- round((results_report$exp_ci_high - 1) * 100, 1)
results_report$significant   <- ifelse(results_report$p_value < 0.05, "Yes", "No")

results_report <- results_report[, c("term", "pct_change", "pct_ci_low",
                                     "pct_ci_high", "p_value", "significant")]
names(results_report) <- c("Predictor", "% change in median phosphorus",
                           "95% CI low", "95% CI high", "p-value", "Significant (p<.05)")
results_report

## ---------------------------------------------------------------
## 6. Quick diagnostics
## ---------------------------------------------------------------
## % censored, by way of a sanity check on model fit
mean(model_df$p_status == 0)

## Predicted (latent, uncensored) values on the response scale
model_df$pred <- predict(fit, type = "response")

## ---------------------------------------------------------------
## Figures
## ---------------------------------------------------------------
## Nitrate/nitrite: same pattern, same DL logic
##   df$n_censored <- df$no2no3mgl == 0
##   df$n_value    <- ifelse(df$n_censored, 0.01, df$no2no3mgl)
##   df$n_status   <- ifelse(df$n_censored, 0, 1)
##
## Fecal coliform: censoring flagged via fecalcoliform_rmk instead of a
## fixed numeric pattern -- check which remark codes mean "below
## detection" (e.g. "U") vs. other qualifiers (e.g. "B1"-"B7" often mean
## blank-contamination flags, not censoring) before treating a code as
## left-censored. Then build f_value/f_status the same way.

## okay now figures

library(ggplot2)
library(dplyr)
library(ggeffects)
library(patchwork)

##### Forest plot #####
forest_df <- results_scaled %>%
  filter(term != "(Intercept)") %>%
  mutate(
    Predictor = factor(
      label_map[term],
      levels = rev(c(
        "CAFO density (per +0.1 CAFOs/km2)",
        "Landfill density (per +0.01 landfills/km2)",
        "Cropland (per +10 percentage points)",
        "Population density (per +100 people/km2)"
      ))
    ),
    
    pct_change  = (exp_estimate - 1) * 100,
    pct_ci_low  = (exp_ci_low - 1) * 100,
    pct_ci_high = (exp_ci_high - 1) * 100,
    
    Significant = p_value < 0.05,
    
    Label = sprintf("%.1f%%", pct_change)
  )

forest_plot <-
  
  ggplot(
    forest_df,
    aes(
      x = pct_change,
      y = Predictor,
      colour = Significant
    )
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  
  geom_errorbar(
    aes(
      xmin = pct_ci_low,
      xmax = pct_ci_high
    ),
    orientation = "y",
    width = 0.18
  ) +
  
  geom_point(size = 3.5) +
  
  geom_text(
    aes(label = Label),
    colour = "black",
    nudge_x = 2,
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
    title = "Effects of watershed characteristics on phosphorus concentrations",
    subtitle = "Points are estimated percent changes; bars show 95% confidence intervals",
    x = "Percent change in median phosphorus concentration",
    y = NULL
  ) +
  
  theme_classic(base_size = 14)

forest_plot
##### OBSERVED VS PREDICTED #####

model_df$predicted <- predict(
  fit,
  type = "response"
)

ggplot(model_df,
       aes(predicted,
           p_value)) +
  
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

# CAFO density
p1 <- ggplot(
  ggpredict(fit, terms = "huc12_upstream_cafo_density"),
  aes(x, predicted)
) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.2
  ) +
  theme_bw() +
  labs(
    x = "CAFO density",
    y = "Predicted phosphorus (mg/L)"
  )

# Landfill density
p2 <- ggplot(
  ggpredict(fit, terms = "huc12_upstream_landfill_density"),
  aes(x, predicted)
) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.2
  ) +
  theme_bw() +
  labs(
    x = "Landfill density",
    y = NULL
  )

# Cropland
p3 <- ggplot(
  ggpredict(fit, terms = "huc12_pct_cropland"),
  aes(x, predicted)
) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.2
  ) +
  theme_bw() +
  labs(
    x = "Cropland (%)",
    y = "Predicted phosphorus (mg/L)"
  )

# Population density
p4 <- ggplot(
  ggpredict(fit, terms = "huc12_pop_density"),
  aes(x, predicted)
) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.2
  ) +
  theme_bw() +
  labs(
    x = "Population density",
    y = NULL
  )

# Combine panels
(p1 | p2) /
  (p3 | p4)

plot(
  model_df$p_value,
  model_df$n_value,
  pch = 16,
  cex = 0.5,
  xlab = "Phosphorus (mg/L)",
  ylab = "Nitrogen (mg/L)"
)

abline(lm(n_value ~ p_value, data = model_df),
       col = "red",
       lwd = 2)

plot(
  df$phosphorusmgl,
  df$no2no3mgl,
  pch = 16,
  cex = 0.5,
  xlab = "Phosphorus (mg/L)",
  ylab = "Nitrogen (mg/L)"
)

abline(lm(n_value ~ p_value, data = model_df),
       col = "red",
       lwd = 2)

