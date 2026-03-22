
# Import libraries
library(sf)
library(dplyr)
library(ggplot2)
library(spdep)

# Clear workspace
rm(list = ls())

# Load dataset
ct_ham <- readRDS("data_processed/ct_ham_final.rds")

# Remove NA values
ct_ham <- ct_ham %>%
  filter(!is.na(med_income),
         !is.na(gs_pct),
         !is.na(mean_temp))

# =========================================================
# STEP 13 - Measure spatial autocorrelation
# =========================================================

# Create spatial weights matrix
ct_ham.w <- ct_ham |> poly2nb() |> nb2listw()

# Function that saves a Moran's I plot
moran_plot <- function(VAR, nicename) {

  # do moran test
  var_m <- moran.test(ct_ham[[VAR]], ct_ham.w)

  # start creating moran plots in scatterplots folder
  png(paste0("outputs/scatterplots/", VAR, "_moran.png"), width = 2400, height = 1800, res = 300)

  # generate moran plot
  moran.plot(ct_ham[[VAR]], ct_ham.w, xlab = nicename, ylab = paste0("Spatially Lagged ", nicename))

  # add text describing the moran test
  mtext(paste0("Moran's I Plot for ", nicename),
        side=3, line=2, col="black", cex=1.3)
  mtext(paste0("Statistic: ",
               round(var_m$estimate[1], digits = 5),
               ", p-value: ",
               signif(var_m$p.value, digits = 3)),
        side=3, line=1, col="black", cex=1.0)

  # finish
  dev.off()

}

# Save Moran plots for each variable
moran_plot("gs_pct", "Percent Greenspace")
moran_plot("med_income", "Median Household Income")
moran_plot("mean_temp", "Mean Land Surface Temperature")



# =========================================================
# STEP 14 - Detect local spatial patterns
# =========================================================

# 14.1 Local Moran's I

# Function that creates a map of the local Moran's I - adapted from the function of the same name from the isdas package
localmoran.map <- function(VAR, nicename) {

  # extract relevant variable and calculate local moran and type case
  df_moran <- ct_ham |>
    rename(VAR = as.name(VAR),
           key = "CTUID") |>
    transmute(key,
              VAR,
              Z = (VAR - mean(VAR)) / var(VAR),
              SMA = lag.listw(ct_ham.w, Z),
              Type = case_when(Z < 0 & SMA < 0 ~ "LL",
                               Z > 0 & SMA > 0 ~ "HH",
                               TRUE ~ "HL/LH"))

  # calculate p-values
  local_I <- localmoran(df_moran$VAR, ct_ham.w)
  colnames(local_I) <- c("Ii", "E.Ii", "Var.Ii", "Z.Ii", "p.val")
  df_moran <- left_join(df_moran,
                        data.frame(key = df_moran$key,
                                   local_I),
                        by = "key")

  # split type case by significant/not significant
  df_moran <- mutate(df_moran,
                     Type = paste(Type, ifelse(p.val < 0.05, "significant", "not significant")))

  # display colors for each type case
  type_colors <- c(
    "HH significant"       = "red",
    "HH not significant"   = "indianred",
    "LL significant"       = "dodgerblue",
    "LL not significant"   = "lightskyblue2",
    "HL/LH significant"    = "yellow",
    "HL/LH not significant"= "lemonchiffon1"
  )

  # dummy rows for options that may not exist in the data - makes sure they're still on the legend
  dummy <- data.frame(Type = names(type_colors))
  df_moran <- bind_rows(df_moran, dummy)

  # turn Type into a factor vector so that it's ordered nicely in the legend, instead of alphabetically
  df_moran$Type <- factor(df_moran$Type, levels = names(type_colors))

  # create map
  p <- ggplot(df_moran) +
    geom_sf(aes(fill = Type),
            size = 0.05) +
    theme_minimal() +
    scale_fill_manual(values = type_colors) +
    labs(title = paste0("Local Moran's I Plot for ", nicename),
         fill = "Type")

  # save to maps folder
  ggsave(paste0("outputs/maps/", VAR ,"_local_moran.png"), plot = p, width = 8, height = 6, dpi = 300)

}

# Save local Moran maps for each variable
localmoran.map("gs_pct", "Percent Greenspace")
localmoran.map("med_income", "Median Household Income")
localmoran.map("mean_temp", "Mean Land Surface Temperature")


# 14.2 Local G

# Function that creates a map of the Gi* - adapted from the function of the same name from the isdas package
gistar.map <- function(VAR, nicename, bonferroni = FALSE) {

  # extract relevant variable and calculate gistar and type case
  df_gistar <- ct_ham
  df.lg <- localG(df_gistar[[VAR]], ct_ham.w)
  df.lg <- as.numeric(df.lg)
  df.lg <- data.frame(Gstar = df.lg, p.val = 2 * pnorm(abs(df.lg), lower.tail = FALSE))

  # bonferroni correction makes the required alpha smaller
  alpha = ifelse(bonferroni, 0.05/nrow(df_gistar), 0.05)

  df.lg <- mutate(df.lg,
                  Type = case_when(Gstar < 0 & p.val <= alpha ~ "Low Concentration",
                                   Gstar > 0 & p.val <= alpha ~ "High Concentration",
                                   TRUE ~ "Not signicant"))

  df_gistar <- left_join(df_gistar,
                         data.frame(CTUID = df_gistar$CTUID, df.lg))

  # display colors for each type case
  type_colors <- c(
    "High Concentration"   = "red",
    "Low Concentration"    = "dodgerblue",
    "Not signicant"        = "lightgray"
  )

  # dummy rows for options that may not exist in the data - makes sure they're still on the legend
  dummy <- data.frame(Type = names(type_colors))
  df_gistar <- bind_rows(df_gistar, dummy)

  # turn Type into a factor vector so that it's ordered nicely in the legend, instead of alphabetically
  df_gistar$Type <- factor(df_gistar$Type, levels = names(type_colors))

  # change plot name for bonferroni
  if (bonferroni) {
    nicename = paste0(nicename, " (Bonferroni Correction)")
    VAR = paste0(VAR, "_bonferroni")
  }

  # create map
  p <- ggplot(df_gistar) +
    geom_sf(aes(fill = Type),
            size = 0.05) +
    theme_minimal() +
    scale_fill_manual(values = type_colors) +
    labs(title = paste0("Gi* Plot for ", nicename),
         fill = "Type")

  # save to maps folder
  ggsave(paste0("outputs/maps/", VAR ,"_gistar.png"), plot = p, width = 8, height = 6, dpi = 300)

}

# Save Gi* maps for each variable
gistar.map("gs_pct", "Percent Greenspace")
gistar.map("med_income", "Median Household Income")
gistar.map("mean_temp", "Mean Land Surface Temperature")


# 14.3 Local G with Bonferroni Correction

# Save Gi* maps for each variable with Bonferroni correction
gistar.map("gs_pct", "Percent Greenspace", TRUE)
gistar.map("med_income", "Median Household Income", TRUE)
gistar.map("mean_temp", "Mean Land Surface Temperature", TRUE)



# =========================================================
# STEP 15 - Regression analysis
# =========================================================

# 15.1 Does greenspace affect urban heat?

# We want to know if greenspace percent and mean temperature are related.
# First, we test a linear model.

lm_model1 <- lm(formula = mean_temp ~ gs_pct,
                data = ct_ham)
summary(lm_model1)

moran.test(lm_model1$residuals, ct_ham.w)

# The residuals are not random, so this model is no good.
# Next, we'll try spatial regression models. We learned about these in Ch 29.
# They can correct for missing model variables or incorrect functional forms.
# These models are commonly used in urban heat island temperature studies.
# This is because they can compensate for the fact that temperature in one
# location will influence temperature of nearby areas, or for the missing
# model variables that represent local environmental factors (distance to
# nearby water bodies, wind speed, etc).

# Examples of studies that use these models:
# https://doi.org/10.1016/j.scitotenv.2018.01.165
# https://doi.org/10.1016/j.landurbplan.2014.01.016
# https://doi.org/10.1016/j.scs.2020.102286
# https://doi.org/10.1016/j.envpol.2025.126025
# https://doi.org/10.1016/j.uclim.2022.101213
# https://doi.org/10.1038/s41598-025-07980-w
# https://doi.org/10.1016/j.landurbplan.2013.11.014

# Lag model

lag_model1 <- lagsarlm(formula = mean_temp ~ gs_pct,
                       data = ct_ham,
                       listw = ct_ham.w)
summary(lag_model1)

moran.test(lag_model1$residuals, ct_ham.w)

# Error model

err_model1 <- errorsarlm(formula = mean_temp ~ gs_pct,
                         data = ct_ham,
                         listw = ct_ham.w)
summary(err_model1)

moran.test(err_model1$residuals, ct_ham.w)

# Mixed/General model

mixed_model1 <- lagsarlm(formula = mean_temp ~ gs_pct,
                         data = ct_ham,
                         listw = ct_ham.w,
                         type = "mixed")
summary(mixed_model1)

moran.test(mixed_model1$residuals, ct_ham.w)

# The residuals for the moran test are random for all three of these models.
# That's good! But which one do we use?

# Compare the models using Akaike Information Criterion (AIC).
# The lowest value is the best model.
# The first source above also used this.
AIC(lm_model1, lag_model1, err_model1, mixed_model1)
# The result is that the mixed model is the best.
# This agrees with the sources above.

# Plot the predicted vs observed temperature for the mixed model
ct_ham$mean_temp_pred <- as.numeric(predict(mixed_model1))

mixed_model1.p <- ggplot(ct_ham, aes(x = mean_temp_pred, y = mean_temp)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "blue") +
  labs(title = "Mixed Spatial Regression Model - Model Performance",
       subtitle = "Model of Greenspace Percent vs Mean Land Surface Temperature",
       x = "Predicted Mean Land Surface Temperature",
       y = "Observed Mean Land Surface Temperature")

# save to scatterplots folder
ggsave(paste0("outputs/scatterplots/gs_pct_vs_mean_temp_mixed.png"), plot = mixed_model1.p, width = 8, height = 6, dpi = 300)

# Conclusion:
# We can predict mean temperature using greenspace percent.
# There is clearly a relationship between urban heat and greenspace.


# 15.2 Do lower income neighbourhoods have less greenspace?

# We want to know if median income and greenspace percent are related.
# First, we test a linear model.

lm_model2 <- lm(formula = gs_pct ~ med_income,
                data = ct_ham)
summary(lm_model2)

moran.test(lm_model2$residuals, ct_ham.w)

# The residuals are not random, so this model is no good.
# Next, we'll try spatial regression models.

# Lag model

lag_model2 <- lagsarlm(formula = gs_pct ~ med_income,
                       data = ct_ham,
                       listw = ct_ham.w)
summary(lag_model2)

moran.test(lag_model2$residuals, ct_ham.w)

# Error model

err_model2 <- errorsarlm(formula = gs_pct ~ med_income,
                         data = ct_ham,
                         listw = ct_ham.w)
summary(err_model2)

moran.test(err_model2$residuals, ct_ham.w)

# Mixed/General model

mixed_model2 <- lagsarlm(formula = gs_pct ~ med_income,
                         data = ct_ham,
                         listw = ct_ham.w,
                         type = "mixed")
summary(mixed_model2)

moran.test(mixed_model2$residuals, ct_ham.w)

# The residuals for the moran test are random for all three of these models.
# That's good! But which one do we use?

# Compare the models using Akaike Information Criterion (AIC).
# The lowest value is the best model.
AIC(lm_model2, lag_model2, err_model2, mixed_model2)
# The result is that the mixed model is the best.

# Plot the predicted vs observed greenspace percent for the mixed model
ct_ham$gs_pct_pred <- as.numeric(predict(mixed_model2))

mixed_model2.p <- ggplot(ct_ham, aes(x = gs_pct_pred, y = gs_pct)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "blue") +
  labs(title = "Mixed Spatial Regression Model - Model Performance",
       subtitle = "Model of Median Income vs Greenspace Percent",
       x = "Predicted Greenspace Percent",
       y = "Observed Greenspace Percent")

# save to scatterplots folder
ggsave(paste0("outputs/scatterplots/med_income_vs_gs_pct_mixed.png"), plot = mixed_model2.p, width = 8, height = 6, dpi = 300)

# Conclusion:
# We can predict greenspace percent using median income.
# There is clearly a relationship between greenspace and median income.
