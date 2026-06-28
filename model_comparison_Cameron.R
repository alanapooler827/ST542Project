##### Alana's work but shifted to HUC10 #####

library(lme4)
library(lmerTest)
library(tidyverse)

df <- read_csv("waterdata.csv")

# format date column
df$date <- as.Date(df$date, format = '%m/%d/%y')

# extract month from date as new column
df <- df |>
  mutate(month = month(date))

# turn month into wet/dry season

wet <- c(5,6,7,8,9)
df$szn <- "dry"
df$szn <- ifelse(df$month %in% wet, "wet", "dry")
df$szn <- factor(df$szn, levels = c("dry", "wet"))  # dry is reference

# convert columns to factors where necessary
df <- df |>
  mutate(
    huc8id = as.factor(huc8id),
    huc10id = as.factor(huc10id),
    huc10id = as.factor(huc10id),
    station_id = as.factor(station_id),
    stationcategory = as.factor(stationcategory),
  )

huc10_df <- df |>
  drop_na(huc10_nearest_upstream_cafo_dist)

huc10_df <- huc10_df %>%
  mutate(across(
    c(huc10_upstream_cafo_density, 
      huc10_upstream_cafo_count,
      huc10_weighted_cafo_load, 
      huc10_nearest_upstream_cafo_dist),
    ~ scale(.)[,1],
    .names = "{.col}_scaled"
  ))

max(huc10_df$huc10_upstream)
# scaled cafo predictors to test one at a time
all_cafo_predictors <- c(
  'huc10_upstream_cafo_density_scaled',
  'huc10_upstream_cafo_count_scaled',
  'huc10_weighted_cafo_load_scaled',
  'huc10_nearest_upstream_cafo_dist_scaled'
)

# covariates for backward elimination
covariates <- c(
  'huc10_pop_density',
  'huc10_pct_cropland',
  'huc10_upstream_wastewater_count',
  'huc10_upstream_landfill_count',
  'huc10_upstream_other_npdes_count',
  'stationcategory',
  'szn'
)

nitlog_results <- cafo_model_selection(huc10_df,
                                       response = 'nitlog',
                                       all_cafo_predictors,
                                       covariates,
                                       random_effects = '+ (1 | huc8id / station_id)')

feclog_results <- cafo_model_selection(huc10_df,
                                       response = 'feclog',
                                       all_cafo_predictors,
                                       covariates,
                                       random_effects = '+ (1 | huc8id / station_id)')

phoslog_results <- cafo_model_selection(huc10_df,
                                       response = 'phoslog',
                                       all_cafo_predictors,
                                       covariates,
                                       random_effects = '+ (1 | huc8id / station_id)')
nitlog_results$best_model
feclog_results$best_model
phoslog_results$best_model
