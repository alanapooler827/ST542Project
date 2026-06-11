# ============================================================
# Water Quality EDA Script
# ============================================================
# Covers:
#   1. Data loading & cleaning
#   2. Summary statistics
#   3. Response variable distributions
#   4. Exceedance rates
#   5. Temporal trends
#   6. CAFO predictor relationships
#   7. Spatial/watershed summaries
#   8. Correlation matrix
#   9. Mixed effects model prep checks
# ============================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(GGally)      # ggpairs correlation plots
library(ggridges)    # ridge plots
library(patchwork)   # combine plots
library(scales)      # axis formatting

theme_set(theme_bw(base_size = 12))

# ============================================================
# 1. LOAD & CLEAN
# ============================================================

water <- read.csv("waterdata.csv") |>
  filter(!is.na(station_id) & station_id != "") |>
  mutate(
    exceed_fec  = as.factor(exceed_fec),
    exceed_nit  = as.factor(exceed_nit),
    exceed_phos = as.factor(exceed_phos),
    year        = as.factor(year),
    stationcategory = as.factor(stationcategory),
    streamclass     = as.factor(streamclass),
    subbasin        = as.factor(subbasin)
  )

water$huc10id_f <- as.factor(water$huc10id)
cat("Rows:", nrow(water), "\n")
cat("Unique stations:", n_distinct(water$station_id), "\n")
cat("Unique HUC8s:", n_distinct(water$huc8id), "\n")
cat("Unique HUC10s:", n_distinct(water$huc10id), "\n")
cat("Unique HUC12s:", n_distinct(water$huc12id), "\n")

cat("Years:", paste(sort(unique(water$year)), collapse = ", "), "\n")


# ============================================================
# 2. SUMMARY STATISTICS
# ============================================================

response_vars <- c("fecalcoliformcol100ml", "feclog",
                   "no2no3mgl", "nitlog",
                   "phosphorusmgl", "phoslog")

cat("\n--- Summary statistics: response variables ---\n")
print(summary(water[, response_vars]))

cafo_vars <- c("weighted_cafo_load", "weighted_poultry_cafo_load",
               "huc10_weighted_cafo_load", "huc10_weighted_poultry_cafo_load",
               "huc12_weighted_cafo_load", "huc12_weighted_poultry_cafo_load",
               "upstream_cafo_density", "nearest_upstream_cafo_distance_m")

cat("\n--- Summary statistics: CAFO predictors ---\n")
print(summary(water[, cafo_vars]))

cat("\n--- Missing values ---\n")
miss <- sapply(water[, c(response_vars, cafo_vars)], function(x) sum(is.na(x)))
print(miss[miss > 0])


# ============================================================
# 3. RESPONSE VARIABLE DISTRIBUTIONS
# ============================================================

# Raw vs log-transformed side by side
p_fec_raw <- ggplot(water, aes(x = fecalcoliformcol100ml)) +
  geom_histogram(bins = 50, fill = "#378ADD", color = "white", linewidth = 0.2) +
  scale_x_continuous(labels = comma) +
  labs(title = "Fecal coliform (raw)", x = "col/100mL", y = "Count")

p_fec_log <- ggplot(water, aes(x = feclog)) +
  geom_histogram(bins = 40, fill = "#378ADD", color = "white", linewidth = 0.2) +
  labs(title = "Fecal coliform (log)", x = "log(col/100mL)", y = "Count")

p_nit_raw <- ggplot(water, aes(x = no2no3mgl)) +
  geom_histogram(bins = 50, fill = "#1D9E75", color = "white", linewidth = 0.2) +
  labs(title = "Nitrate (raw)", x = "mg/L", y = "Count")

p_nit_log <- ggplot(water, aes(x = nitlog)) +
  geom_histogram(bins = 40, fill = "#1D9E75", color = "white", linewidth = 0.2) +
  labs(title = "Nitrate (log)", x = "log(mg/L)", y = "Count")

p_phos_raw <- ggplot(water, aes(x = phosphorusmgl)) +
  geom_histogram(bins = 50, fill = "#D85A30", color = "white", linewidth = 0.2) +
  labs(title = "Phosphorus (raw)", x = "mg/L", y = "Count")

p_phos_log <- ggplot(water, aes(x = phoslog)) +
  geom_histogram(bins = 40, fill = "#D85A30", color = "white", linewidth = 0.2) +
  labs(title = "Phosphorus (log)", x = "log(mg/L)", y = "Count")

(p_fec_raw | p_fec_log) /
  (p_nit_raw | p_nit_log) /
  (p_phos_raw | p_phos_log)

ggsave("eda_01_distributions.png", width = 10, height = 10, dpi = 150)


# ============================================================
# 4. EXCEEDANCE RATES
# ============================================================

# Overall rates
exceed_summary <- water |>
  summarise(
    fec_rate  = mean(exceed_fec  == 1, na.rm = TRUE),
    nit_rate  = mean(exceed_nit  == 1, na.rm = TRUE),
    phos_rate = mean(exceed_phos == 1, na.rm = TRUE)
  ) |>
  pivot_longer(everything(), names_to = "pollutant", values_to = "rate") |>
  mutate(pollutant = recode(pollutant,
                            fec_rate  = "Fecal coliform",
                            nit_rate  = "Nitrate",
                            phos_rate = "Phosphorus"
  ))

p_exceed_overall <- ggplot(exceed_summary, aes(x = pollutant, y = rate, fill = pollutant)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("Fecal coliform" = "#378ADD",
                               "Nitrate"         = "#1D9E75",
                               "Phosphorus"      = "#D85A30")) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  labs(title = "Overall exceedance rates", x = NULL, y = "% of samples") +
  theme(legend.position = "none")

# Exceedance by station category
exceed_by_cat <- water |>
  group_by(stationcategory) |>
  summarise(
    fec  = mean(exceed_fec  == 1, na.rm = TRUE),
    nit  = mean(exceed_nit  == 1, na.rm = TRUE),
    phos = mean(exceed_phos == 1, na.rm = TRUE),
    n    = n()
  ) |>
  pivot_longer(c(fec, nit, phos), names_to = "pollutant", values_to = "rate")

p_exceed_cat <- ggplot(exceed_by_cat, aes(x = stationcategory, y = rate, fill = pollutant)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c(fec = "#378ADD", nit = "#1D9E75", phos = "#D85A30"),
                    labels = c("Fecal coliform", "Nitrate", "Phosphorus")) +
  scale_y_continuous(labels = percent) +
  labs(title = "Exceedance by station category", x = NULL, y = "% exceeding", fill = NULL) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

p_exceed_overall | p_exceed_cat
ggsave("eda_02_exceedance.png", width = 12, height = 5, dpi = 150)


# ============================================================
# 5. TEMPORAL TRENDS
# ============================================================

temporal <- water |>
  group_by(year) |>
  summarise(
    feclog_mean  = mean(feclog,  na.rm = TRUE),
    nitlog_mean  = mean(nitlog,  na.rm = TRUE),
    phoslog_mean = mean(phoslog, na.rm = TRUE),
    fec_exceed   = mean(exceed_fec  == 1, na.rm = TRUE),
    nit_exceed   = mean(exceed_nit  == 1, na.rm = TRUE),
    phos_exceed  = mean(exceed_phos == 1, na.rm = TRUE),
    n            = n()
  )

# Mean log levels over time
p_temp_levels <- temporal |>
  pivot_longer(c(feclog_mean, nitlog_mean, phoslog_mean),
               names_to = "pollutant", values_to = "mean_log") |>
  mutate(pollutant = recode(pollutant,
                            feclog_mean  = "Fecal coliform",
                            nitlog_mean  = "Nitrate",
                            phoslog_mean = "Phosphorus"
  )) |>
  ggplot(aes(x = year, y = mean_log, color = pollutant, group = pollutant)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Fecal coliform" = "#378ADD",
                                "Nitrate"         = "#1D9E75",
                                "Phosphorus"      = "#D85A30")) +
  labs(title = "Mean log pollutant levels by year", x = "Year",
       y = "Mean log value", color = NULL)

# Exceedance rates over time
p_temp_exceed <- temporal |>
  pivot_longer(c(fec_exceed, nit_exceed, phos_exceed),
               names_to = "pollutant", values_to = "rate") |>
  mutate(pollutant = recode(pollutant,
                            fec_exceed  = "Fecal coliform",
                            nit_exceed  = "Nitrate",
                            phos_exceed = "Phosphorus"
  )) |>
  ggplot(aes(x = year, y = rate, color = pollutant, group = pollutant)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Fecal coliform" = "#378ADD",
                                "Nitrate"         = "#1D9E75",
                                "Phosphorus"      = "#D85A30")) +
  scale_y_continuous(labels = percent) +
  labs(title = "Exceedance rates by year", x = "Year",
       y = "% exceeding", color = NULL)

p_temp_levels / p_temp_exceed
ggsave("eda_03_temporal.png", width = 10, height = 8, dpi = 150)


# ============================================================
# 6. CAFO PREDICTOR RELATIONSHIPS
# ============================================================

# Weighted CAFO load distributions
p_cafo_dist <- water |>
  select(weighted_cafo_load, huc10_weighted_cafo_load, huc12_weighted_cafo_load) |>
  pivot_longer(everything(), names_to = "scale", values_to = "load") |>
  mutate(scale = recode(scale,
                        weighted_cafo_load       = "Upstream",
                        huc10_weighted_cafo_load = "HUC10",
                        huc12_weighted_cafo_load = "HUC12"
  )) |>
  ggplot(aes(x = load, fill = scale)) +
  geom_histogram(bins = 40, color = "white", linewidth = 0.2, alpha = 0.85) +
  facet_wrap(~scale, scales = "free") +
  scale_fill_manual(values = c(Upstream = "#378ADD", HUC10 = "#7F77DD", HUC12 = "#1D9E75")) +
  labs(title = "Weighted CAFO load distributions by spatial scale",
       x = "Weighted CAFO load", y = "Count") +
  theme(legend.position = "none")

ggsave("eda_04_cafo_distributions.png", p_cafo_dist, width = 12, height = 4, dpi = 150)

# Scatter: weighted CAFO load vs each response
scatter_cafo <- function(cafo_col, cafo_label) {
  p1 <- ggplot(water, aes(x = .data[[cafo_col]], y = feclog)) +
    geom_point(alpha = 0.3, size = 1.2, color = "#378ADD") +
    geom_smooth(method = "lm", color = "black", linewidth = 0.8, se = TRUE) +
    labs(x = cafo_label, y = "Fecal coliform (log)")
  
  p2 <- ggplot(water, aes(x = .data[[cafo_col]], y = nitlog)) +
    geom_point(alpha = 0.3, size = 1.2, color = "#1D9E75") +
    geom_smooth(method = "lm", color = "black", linewidth = 0.8, se = TRUE) +
    labs(x = cafo_label, y = "Nitrate (log)")
  
  p3 <- ggplot(water, aes(x = .data[[cafo_col]], y = phoslog)) +
    geom_point(alpha = 0.3, size = 1.2, color = "#D85A30") +
    geom_smooth(method = "lm", color = "black", linewidth = 0.8, se = TRUE) +
    labs(x = cafo_label, y = "Phosphorus (log)")
  
  p1 | p2 | p3
}

scatter_cafo("weighted_cafo_load", "Upstream weighted CAFO load")
ggsave("eda_05_cafo_scatter_upstream.png", width = 13, height = 4, dpi = 150)

scatter_cafo("huc10_weighted_cafo_load", "HUC10 weighted CAFO load")
ggsave("eda_06_cafo_scatter_huc10.png", width = 13, height = 4, dpi = 150)

scatter_cafo("huc12_weighted_cafo_load", "HUC12 weighted CAFO load")
ggsave("eda_07_cafo_scatter_huc12.png", width = 13, height = 4, dpi = 150)

# Boxplots: exceedance vs weighted CAFO load
p_cafo_exceed <- water |>
  select(exceed_fec, weighted_cafo_load, weighted_poultry_cafo_load) |>
  pivot_longer(c(weighted_cafo_load, weighted_poultry_cafo_load),
               names_to = "cafo_type", values_to = "load") |>
  filter(!is.na(exceed_fec)) |>
  mutate(cafo_type = recode(cafo_type,
                            weighted_cafo_load         = "All CAFOs",
                            weighted_poultry_cafo_load = "Poultry CAFOs"
  )) |>
  ggplot(aes(x = exceed_fec, y = load, fill = exceed_fec)) +
  geom_boxplot(outlier.alpha = 0.2, outlier.size = 0.8) +
  facet_wrap(~cafo_type) +
  scale_fill_manual(values = c("0" = "#B5D4F4", "1" = "#378ADD")) +
  labs(title = "Upstream CAFO load by fecal coliform exceedance",
       x = "Exceeds standard", y = "Weighted CAFO load") +
  theme(legend.position = "none")

ggsave("eda_08a_cafo_exceedance_boxplot.png", p_cafo_exceed, width = 8, height = 5, dpi = 150)

# Exceedance at HUC10

p_cafo_exceed_10 <- water |>
  select(exceed_fec, huc10_weighted_cafo_load, huc10_weighted_poultry_cafo_load) |>
  pivot_longer(c(huc10_weighted_cafo_load, huc10_weighted_poultry_cafo_load),
               names_to = "cafo_type", values_to = "load") |>
  filter(!is.na(exceed_fec)) |>
  mutate(cafo_type = recode(cafo_type,
                            huc10_weighted_cafo_load         = "HUC10 - All CAFOs",
                            huc10_weighted_poultry_cafo_load = "HUC10 - Poultry CAFOs"
  )) |>
  ggplot(aes(x = exceed_fec, y = load, fill = exceed_fec)) +
  geom_boxplot(outlier.alpha = 0.2, outlier.size = 0.8) +
  facet_wrap(~cafo_type) +
  scale_fill_manual(values = c("0" = "#B5D4F4", "1" = "#378ADD")) +
  labs(title = "HUC10 CAFO load by fecal coliform exceedance",
       x = "Exceeds standard", y = "Weighted CAFO load") +
  theme(legend.position = "none")

ggsave("eda_08b_cafo_exceedance_boxplot_10.png", p_cafo_exceed_10, width = 8, height = 5, dpi = 150)


# Exceedance at HUC12

p_cafo_exceed_12 <- water |>
  select(exceed_fec, huc12_weighted_cafo_load, huc12_weighted_poultry_cafo_load) |>
  pivot_longer(c(huc12_weighted_cafo_load, huc12_weighted_poultry_cafo_load),
               names_to = "cafo_type", values_to = "load") |>
  filter(!is.na(exceed_fec)) |>
  mutate(cafo_type = recode(cafo_type,
                            huc12_weighted_cafo_load         = "HUC12 - All CAFOs",
                            huc12_weighted_poultry_cafo_load = "HUC12 - Poultry CAFOs"
  )) |>
  ggplot(aes(x = exceed_fec, y = load, fill = exceed_fec)) +
  geom_boxplot(outlier.alpha = 0.2, outlier.size = 0.8) +
  facet_wrap(~cafo_type) +
  scale_fill_manual(values = c("0" = "#B5D4F4", "1" = "#378ADD")) +
  labs(title = "HUC12 CAFO load by fecal coliform exceedance",
       x = "Exceeds standard", y = "Weighted CAFO load") +
  theme(legend.position = "none")

ggsave("eda_08c_cafo_exceedance_boxplot_12.png", p_cafo_exceed_12, width = 8, height = 5, dpi = 150)

ggsave("eda_08d_cafo_exceedance_boxplot_ALL.png", 
        p_cafo_exceed /
        p_cafo_exceed_10 /
        p_cafo_exceed_12, 
        width = 6, height = 10, dpi = 150)

p_cafo_exceed /
  p_cafo_exceed_10 /
  p_cafo_exceed_12

# ============================================================
# 7. WATERSHED / SPATIAL SUMMARIES
# ============================================================

# Mean response by HUC8
huc8_summary <- water |>
  group_by(huc8id) |>
  summarise(
    feclog_mean  = mean(feclog,  na.rm = TRUE),
    nitlog_mean  = mean(nitlog,  na.rm = TRUE),
    phoslog_mean = mean(phoslog, na.rm = TRUE),
    cafo_load    = mean(weighted_cafo_load, na.rm = TRUE),
    n_obs        = n(),
    n_stations   = n_distinct(station_id)
  )

cat("\n--- HUC8 watershed summary ---\n")
print(huc8_summary)

# Ridge plot: feclog distribution by subbasin
p_ridge <- ggplot(water, aes(x = feclog, y = subbasin, fill = subbasin)) +
  geom_density_ridges(alpha = 0.7, scale = 1.2, color = "white") +
  labs(title = "Fecal coliform distribution by subbasin",
       x = "Fecal coliform (log)", y = NULL) +
  theme(legend.position = "none")

ggsave("eda_09a_ridge_subbasin.png", p_ridge, width = 9, height = 6, dpi = 150)

# Ridge plot 2: feclog distibution by 
p_ridge2 <- ggplot(water, aes(x = feclog, y = huc10id_f, fill = huc10id_f)) +
  geom_density_ridges(alpha = 0.7, scale = 1.2, color = "white") +
  labs(title = "Fecal coliform distribution by HUC10 ID",
       x = "Fecal coliform (log)", y = NULL) +
  theme(legend.position = "none")

ggsave("eda_09b_ridge_huc10.png", p_ridge, width = 9, height = 6, dpi = 150)

# Station-level mean vs CAFO load (shows spatial pattern)
station_summary <- water |>
  group_by(station_id, latitude, longitude) |>
  summarise(
    feclog_mean       = mean(feclog, na.rm = TRUE),
    cafo_load_mean    = mean(weighted_cafo_load, na.rm = TRUE),
    exceed_fec_rate   = mean(exceed_fec == 1, na.rm = TRUE),
    n                 = n(),
    .groups = "drop"
  )

p_station_cafo <- ggplot(station_summary, aes(x = cafo_load_mean, y = feclog_mean,
                                              size = n, color = exceed_fec_rate)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "#B5D4F4", high = "#D85A30", labels = percent) +
  scale_size_continuous(range = c(2, 8)) +
  labs(title = "Station-level mean fecal coliform vs mean CAFO load",
       x = "Mean upstream weighted CAFO load",
       y = "Mean fecal coliform (log)",
       color = "Exceedance\nrate",
       size = "n obs")

ggsave("eda_10_station_cafo.png", p_station_cafo, width = 8, height = 6, dpi = 150)




# ============================================================
# 8. CORRELATION MATRIX
# ============================================================

cor_vars <- water |>
  select(feclog, nitlog, phoslog,
         weighted_cafo_load, weighted_poultry_cafo_load,
         huc10_weighted_cafo_load, huc12_weighted_cafo_load,
         upstream_cafo_density,
         nearest_upstream_cafo_distance_m,
         huc8_pct_cropland, huc8_pop_density) |>
  rename(
    fec              = feclog,
    nit              = nitlog,
    phos             = phoslog,
    cafo_up          = weighted_cafo_load,
    poultry_up       = weighted_poultry_cafo_load,
    cafo_h10         = huc10_weighted_cafo_load,
    cafo_h12         = huc12_weighted_cafo_load,
    cafo_dens        = upstream_cafo_density,
    cafo_dist        = nearest_upstream_cafo_distance_m,
    pct_crop         = huc8_pct_cropland,
    pop_dens         = huc8_pop_density
  )

cat("\n--- Pearson correlations ---\n")
print(round(cor(cor_vars, use = "complete.obs"), 2))

test <- round(cor(cor_vars, use = "complete.obs"), 2)

write.csv(test, "correlation_matrix.csv")

png("eda_11_correlation_matrix.png", width = 900, height = 900, res = 120)
ggpairs(cor_vars,
        upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", alpha = 0.15, size = 0.5)),
        diag  = list(continuous = wrap("densityDiag"))) +
  theme_bw(base_size = 9) +
  labs(title = "Correlation matrix — responses and CAFO predictors")
dev.off()


# ============================================================
# 9. MODEL PREP CHECKS
# ============================================================

# Observations per station (check for sufficient replication)
obs_per_station <- water |> count(station_id, name = "n_obs")
cat("\n--- Observations per station ---\n")
print(summary(obs_per_station$n_obs))

p_obs <- ggplot(obs_per_station, aes(x = n_obs)) +
  geom_histogram(bins = 30, fill = "#7F77DD", color = "white", linewidth = 0.2) +
  labs(title = "Observations per station", x = "n observations", y = "Count")

# Stations per HUC8 (check hierarchy balance)
stations_per_huc8 <- water |>
  group_by(huc8id) |>
  summarise(n_stations = n_distinct(station_id), n_obs = n())

cat("\n--- Stations per HUC8 ---\n")
print(stations_per_huc8)

p_huc8 <- ggplot(stations_per_huc8, aes(x = reorder(huc8id, -n_stations),
                                        y = n_stations)) +
  geom_col(fill = "#1D9E75") +
  geom_text(aes(label = n_stations), vjust = -0.4, size = 3.5) +
  labs(title = "Stations per HUC8 watershed",
       x = "HUC8 ID", y = "Number of stations") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

p_obs | p_huc8
ggsave("eda_12_model_prep.png", width = 12, height = 5, dpi = 150)

cat("\n===== EDA complete. Check your working directory for output PNGs. =====\n")
