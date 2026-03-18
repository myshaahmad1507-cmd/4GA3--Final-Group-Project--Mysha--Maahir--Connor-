
# Import libraries
library(sf)
library(dplyr)
library(ggplot2)
library(spdep)

# Load dataset
ct_ham <- readRDS("data_processed/ct_ham_final.rds")

# Remove NA values
ct_ham <- ct_ham %>%
  filter(!is.na(med_income),
         !is.na(gs_pct),
         !is.na(mean_temp))

# =========================================================
# STEP 10 - Create Exploratory Maps
# =========================================================

# Create a folder for maps (../outputs/maps/)
dir.create("outputs/maps", recursive = TRUE, showWarnings = FALSE)

# Define map variables and parameters
variables <- c("med_income", "gs_pct", "mean_temp")
labels <- c("Income ($)", "%", "Temp (°C)")
classes <- c(3, 6, 9)

# Double for loop to create maps for the variables (variables) with N classes (classes)
for (j in 1:length(variables)) {
  for (i in 1:length(classes)) {
    p <- ggplot(ct_ham) +
      geom_sf(aes(fill = cut_interval(.data[[variables[j]]], classes[i])),
              color = "white",
              size = 0.05) +
      scale_fill_brewer(palette = "Oranges") +
      theme_minimal() +
      labs(fill = labels[j])

    # Save map as .png to output folder ("outputs/maps/")
    ggsave(paste0("outputs/maps/",variables[j],"_",classes[i],"_classes.png"), plot = p, width = 8, height = 6, dpi = 300)
  }
}

# =========================================================
# STEP 11 - Create Scatterplots
# =========================================================

# Create a folder for scatterplots (../outputs/scatterplots/)
dir.create("outputs/scatterplots", recursive = TRUE, showWarnings = FALSE)

# Define neighbours using "Queen" criteria
ct_ham.nb <- poly2nb(pl = ct_ham)
#summary(ct_ham.nb)

# Create spatial weights matrix
ct_ham.w <- nb2listw(ct_ham.nb)

# Calculate spatial moving averages
ct_ham$gs_pct.sma <- lag.listw(ct_ham.w, ct_ham$gs_pct)
ct_ham$med_income.sma <- lag.listw(ct_ham.w, ct_ham$med_income)
ct_ham$mean_temp.sma <- lag.listw(ct_ham.w, ct_ham$mean_temp)

# Scatter plot for Greenspace
s1 <- ggplot(ct_ham, aes(x = gs_pct, y = gs_pct.sma)) +
  geom_point(alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Greenspace (%)",
       y = "Spatial Moving Average")

# Scatter plot for Median Income
s2 <- ggplot(ct_ham, aes(x = med_income, y = med_income.sma)) +
  geom_point(alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Median Household Income ($)",
       y = "Spatial Moving Average")

# Scatter plot for Mean Temperature
s3 <- ggplot(ct_ham, aes(x = mean_temp, y = mean_temp.sma)) +
  geom_point(alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Mean Temperature (°C)",
       y = "Spatial Moving Average")

# Save scatterplots as .png to output folder ("outputs/scatterplots/")
ggsave(filename = paste0("outputs/scatterplots/gs_pct_SMA.png"), plot = s1, width = 8, height = 6, dpi = 300)
ggsave(filename = paste0("outputs/scatterplots/med_income_SMA.png"), plot = s2, width = 8, height = 6, dpi = 300)
ggsave(filename = paste0("outputs/scatterplots/mean_temperature_SMA.png"), plot = s3, width = 8, height = 6, dpi = 300)

# =========================================================
# STEP 12 - Null Landscapes
# =========================================================

# Set rondomization seed
set.seed(123)

# Scramble values for each variable
ct_ham$med_income_random <- sample(ct_ham$med_income)
ct_ham$gs_pct_random <- sample(ct_ham$gs_pct)
ct_ham$mean_temp_random <- sample(ct_ham$mean_temp)

# Compute spatial moving average for each variable
ct_ham$med_income_random.sma <- lag.listw(ct_ham.w, ct_ham$med_income_random)
ct_ham$gs_pct_random.sma <- lag.listw(ct_ham.w, ct_ham$gs_pct_random)
ct_ham$mean_temp_random.sma <- lag.listw(ct_ham.w, ct_ham$mean_temp_random)

# Scatter plot for Greenspace (RANDOM)
s4 <- ggplot(ct_ham, aes(x = gs_pct_random, y = gs_pct_random.sma)) +
  geom_point(alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Random Greenspace (%)",
       y = "Spatial Moving Average")

# Scatter plot for Median Income (RANDOM)
s5 <- ggplot(ct_ham, aes(x = med_income_random, y = med_income_random.sma)) +
  geom_point(alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Random Median Household Income ($)",
       y = "Spatial Moving Average")

# Scatter plot for Mean Temperature (RANDOM)
s6 <- ggplot(ct_ham, aes(x = mean_temp_random, y = mean_temp_random.sma)) +
  geom_point(alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Random Mean Temperature (°C)",
       y = "Spatial Moving Average")

# Save scatterplots as .png to output folder ("outputs/scatterplots/")
ggsave(filename = paste0("outputs/scatterplots/gs_pct_SMA_random.png"), plot = s4, width = 8, height = 6, dpi = 300)
ggsave(filename = paste0("outputs/scatterplots/med_income_SMA_random.png"), plot = s5, width = 8, height = 6, dpi = 300)
ggsave(filename = paste0("outputs/scatterplots/mean_temperature_SMA_random.png"), plot = s6, width = 8, height = 6, dpi = 300)

