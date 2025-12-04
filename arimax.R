# import libraries
# install.packages("forecast")
library(forecast)
# install.packages("ggplot2")
library(ggplot2)
# install.packages("patchwork")
library(patchwork)
# install.packages("MASS")
library(MASS)

# load in sp500 returns + macro surprises
spx_data <- read.csv("data/sp500_hist_returns.csv")
spx_data <- na.omit(spx_data)
macro_data <- read.csv("data/macro_indicators_data.csv")

# reformat date columns
spx_data$date <- as.Date(spx_data$date, format = "%Y-%m-%d")
macro_data$date <- as.Date(macro_data$date, format = "%Y-%m-%d")

# merge datasets
data_full <- merge(spx_data, macro_data, by = "date", all.x = TRUE)

data <- data_full[c("date",
                    "ln_return",
                    "NFP_surprise",
                    "core_PCE_surprise",
                    "core_retail_sales_surprise")]


corr_df <- data[c("ln_return", 
                  "NFP_surprise", 
                  "core_PCE_surprise", 
                  "core_retail_sales_surprise")]

# get correlation matrix while handling NAs
corr_mat <- cor(corr_df, use = "pairwise.complete.obs")

print(corr_mat)

macro_cols <- c("NFP_surprise", "core_PCE_surprise", "core_retail_sales_surprise")

# define X matrix of macro surprises
X <- data[macro_cols]

# extract means and sds for interpretation + scenario scaling
macro_means <- colMeans(X, na.rm = TRUE)
macro_sds <- apply(X, 2, sd, na.rm = TRUE)

print("Macro indicator means:")
print(macro_means)
print("Macro indicator sds:")
print(macro_sds)

# scale all X
scaled_X <- scale(X)

# replace NA surprises with 0 (no surprise on that day)
scaled_X[is.na(scaled_X)] <- 0

# define xreg
xreg <- as.matrix(scaled_X)

# define target variable
y <- data$ln_return

# fit t-distribution using MLE
df_fit <- fitdistr(y, "t")
deg_freedom <- df_fit$estimate["df"]

# train-test split validation
cut_date <- as.Date("2022-12-31")

train_idx <- data$date <= cut_date
test_idx  <- data$date >  cut_date

y_train <- y[train_idx]
y_test <- y[test_idx]
xreg_train <- xreg[train_idx, , drop = FALSE]
xreg_test <- xreg[test_idx, , drop = FALSE]

cat("Train size:", length(y_train), " Test size:", length(y_test), "\n")

# fit ARIMAX on training data only
fit_train <- auto.arima(y_train, xreg = xreg_train, seasonal = FALSE)

cat("Training model summary:\n")
print(summary(fit_train))

# out-of-sample forecast on test set
fc_test <- forecast(fit_train, xreg = xreg_test, h = nrow(xreg_test))

y_hat <- fc_test$mean

# basic error metrics
rmse <- sqrt(mean((y_hat - y_test)^2, na.rm = TRUE))
mae <- mean(abs(y_hat - y_test), na.rm = TRUE)
dir_acc <- mean(sign(y_hat) == sign(y_test), na.rm = TRUE)

cat("Out-of-sample performance;\n")
cat("RMSE:", rmse, "\n")
cat("MAE: ", mae,  "\n")
cat("Directional accuracy:", round(100 * dir_acc, 2), "%\n")

# refit on full data after satisfactory training performance
train_specs <- arimaorder(fit_train)
print(train_specs)
fit <- Arima(y, order = train_specs, xreg = xreg, seasonal = FALSE)

cat("\nFull-sample model:\n")
print(summary(fit))

# extract and print coefficients
cat("\nCoefficient Table:\n")
coef_table <- summary(fit)$coef
print(coef_table)

coef_nfp <- coef_table["NFP_surprise"]
coef_pce <- coef_table["core_PCE_surprise"]
coef_retail <- coef_table["core_retail_sales_surprise"]

bps_per_sd <- 100 * 100 * c(NFP = coef_nfp, PCE = coef_pce, Retail = coef_retail) 
print(bps_per_sd)

# check residuals
print(checkresiduals(fit))

scenario_analysis <- function(nfp_surprise = 0,
                              core_PCE_surprise = 0,
                              core_retail_sales_surprise = 0,
                              scaled = FALSE) {

  # collect scenario inputs
  scenario <- c(
    NFP_surprise = nfp_surprise,
    core_PCE_surprise = core_PCE_surprise,
    core_retail_sales_surprise = core_retail_sales_surprise
  )
  
  # if not scaled, convert raw values to standardized form
  if (!scaled) {
    scenario <- (scenario - macro_means) / macro_sds
  }
  
  # convert matrix for forecasting
  future_xreg <- matrix(scenario, nrow = 1)
  colnames(future_xreg) <- colnames(xreg)
  
  # forecast 1-day ahead with scenario
  fc <- forecast(fit, xreg = future_xreg, h = 1, level = 95)
  mu <- fc$mean[1]
  
  # back out standard deviation from 95% interval
  z  <- qnorm(0.975)
  lo95 <- fc$lower[1, "95%"]
  hi95 <- fc$upper[1, "95%"]
  
  sd_from_hi <- (hi95 - mu) / z
  sd_from_lo <- (mu - lo95) / z
  sd_hat <- mean(c(sd_from_hi, sd_from_lo))
  
  # return predictive mean and sd
  return(list(mu = mu, sd = sd_hat))
}

plot_scenario <- function(nfp_surprise = 0,
                          core_PCE_surprise = 0,
                          core_retail_sales_surprise = 0,
                          scaled = FALSE) {
  
  # call scenario analysis function
  scenario_res <- scenario_analysis(
    nfp_surprise = nfp_surprise,
    core_PCE_surprise = core_PCE_surprise,
    core_retail_sales_surprise = core_retail_sales_surprise,
    scaled = scaled
  )
  
  # extract degrees of freedom from earlier fit
  deg_freedom <- df_fit$estimate["df"]

  mu <- scenario_res$mu
  sd <- scenario_res$sd

  # adjust sd for t-distribution
  sd_t <- sd * (qnorm(0.975) / qt(0.975, deg_freedom))
  
  # sample from t-student predictive distribution
  set.seed(42)
  N <- 10000
  sim_ret <- mu + sd_t * rt(N, df = deg_freedom)
  
  # Plot histogram + density of simulated returns (in %)
  sim_df <- data.frame(ret = 100 * sim_ret)
  
  p <- ggplot(sim_df, aes(x = ret)) +
    geom_histogram(aes(y = ..density..), bins = 50) +
    geom_density(linewidth = 1) +
    geom_vline(xintercept = mean(sim_df$ret),
               linetype = "dashed", linewidth = 1) +
    labs(
      title = sprintf(
        "Predictive Distribution of Daily SPX Returns\nunder Macro Surprise Scenario\nNFP: %s, Core PCE: %s, Core Retail Sales: %s\n(scaled = %s)",
        nfp_surprise, core_PCE_surprise, core_retail_sales_surprise, scaled
      ),
      x = "Simulated daily return (%)",
      y = "Density"
    )
  return(p)
}

# define base plots
base_plot_scaled <- plot_scenario(scaled = TRUE)
base_plot_unscaled <- plot_scenario(scaled = FALSE)

# define scenario plots at +/- 2 sd surprises
nfp_2sd_plot <- plot_scenario(nfp_surprise = 2, scaled = TRUE)
nfp_neg2sd_plot <- plot_scenario(nfp_surprise = -2, scaled = TRUE)
PCE_2sd_plot <- plot_scenario(core_PCE_surprise = 2, scaled = TRUE)
PCE_neg2sd_plot <- plot_scenario(core_PCE_surprise = -2, scaled = TRUE)
retail_2sd_plot <- plot_scenario(core_retail_sales_surprise = 2, scaled = TRUE)
retail_neg2sd_plot <- plot_scenario(core_retail_sales_surprise = -2, scaled = TRUE)
# print scenario plots
print(base_plot_unscaled / nfp_2sd_plot / nfp_neg2sd_plot)
print(scenario_analysis(nfp_surprise = 2, scaled = TRUE)$mu)
print(scenario_analysis(nfp_surprise = -2, scaled = TRUE)$mu)
print(scenario_analysis(scaled = TRUE)$mu)
print(scenario_analysis(scaled = FALSE)$mu)
print(base_plot_scaled / PCE_2sd_plot / PCE_neg2sd_plot)
print(scenario_analysis(core_PCE_surprise = 2, scaled = TRUE)$mu)
print(scenario_analysis(core_PCE_surprise = -2, scaled = TRUE)$mu)
print(base_plot_scaled / retail_2sd_plot / retail_neg2sd_plot)
print(scenario_analysis(core_retail_sales_surprise = 2, scaled = TRUE)$mu)
print(scenario_analysis(core_retail_sales_surprise = -2, scaled = TRUE)$mu)

# define range of z-scores for reaction curves
z_vals <- seq(-3, 3, by = 0.1)

# helper to build one reaction curve
build_reaction_df <- function(which_var = c("NFP", "PCE", "Retail")) {
  which_var <- match.arg(which_var)
  mu_vals <- numeric(length(z_vals))

  for (i in seq_along(z_vals)) {
    res <- scenario_analysis(
      nfp_surprise = if (which_var == "NFP") z_vals[i] else 0,
      core_PCE_surprise = if (which_var == "PCE") z_vals[i] else 0,
      core_retail_sales_surprise = if (which_var == "Retail") z_vals[i] else 0,
      scaled = TRUE
    )
    mu_vals[i] <- 100 * res$mu   # convert to %
  }

  df <- data.frame(z = z_vals, mu_pct = mu_vals)

  # baseline at neutral surprise (z = 0)
  mu0 <- df$mu_pct[df$z == 0]

  # change vs baseline, in basis points
  df$delta_bp <- (df$mu_pct - mu0) * 100

  return(df)
}

# build each curve
df_nfp <- build_reaction_df("NFP")
df_PCE <- build_reaction_df("PCE")
df_retail <- build_reaction_df("Retail")

# plots of change in expected return (bps) instead of raw %
nfp_reaction_plot <- ggplot(df_nfp, aes(x = z, y = delta_bp)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 1) +
  coord_cartesian(ylim = c(-60, 60)) +
  labs(
    title = "SPX Reaction to NFP Surprise",
    x = "NFP Surprise (z-score)",
    y = "Change in expected daily return (bps)"
  )

PCE_reaction_plot <- ggplot(df_PCE, aes(x = z, y = delta_bp)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 1) +
  coord_cartesian(ylim = c(-60, 60)) +
  labs(
    title = "SPX Reaction to Core PCE Surprise",
    x = "Core PCE Surprise (z-score)",
    y = "Change in expected daily return (bps)"
  )

retail_reaction_plot <- ggplot(df_retail, aes(x = z, y = delta_bp)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 1) +
  coord_cartesian(ylim = c(-60, 60)) +
  labs(
    title = "SPX Reaction to Retail Sales Surprise",
    x = "Retail Sales Surprise (z-score)",
    y = "Change in expected daily return (bps)"
  )

# show all plots together
print(nfp_reaction_plot / PCE_reaction_plot / retail_reaction_plot)
