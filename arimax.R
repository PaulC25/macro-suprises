# import libraries
# install.packages("forecast")
library(forecast)
library(ggplot2)
# install.packages("patchwork")
library(patchwork)


# load in sp500 returns + macro surprises
spx_data <- read.csv("data/sp500_hist_returns.csv")
macro_data <- read.csv("data/macro_indicators_data.csv")

# reformat date columns
spx_data$date <- as.Date(spx_data$date, format = "%Y-%m-%d")
macro_data$date <- as.Date(macro_data$date, format = "%Y-%m-%d")

# merge datasets
data_full <- merge(spx_data, macro_data, by = "date", all.x = TRUE)

data <- data_full[c("date",
                    "ln_return",
                    "NFP_surprise",
                    "core_CPE_surprise",
                    "core_retail_sales_surprise")]

macro_cols <- c("NFP_surprise", "core_CPE_surprise", "core_retail_sales_surprise")

# define X matrix of macro surprises
X <- data[macro_cols]

# extacrt means and sds for interpretation + scenario scaling
macro_means <- colMeans(X, na.rm = TRUE)
macro_sds <- apply(X, 2, sd, na.rm = TRUE)

# replace NA surprises with 0 (no surprise on that day)
X[is.na(X)] <- 0

# scale all X
scaled_X <- scale(X)
xreg <- as.matrix(scaled_X)

# define target variable
y <- data$ln_return

# train-test split validation
cut_date <- as.Date("2022-12-31")

train_idx <- data$date <= cut_date
test_idx  <- data$date >  cut_date

y_train <- y[train_idx]
y_test <- y[test_idx]
xreg_train <- xreg[train_idx, , drop = FALSE]
xreg_test <- xreg[test_idx,  , drop = FALSE]

cat("Train size:", length(y_train), " Test size:", length(y_test), "\n")

# fit ARIMAX on training data only
fit_train <- auto.arima(y_train, xreg = xreg_train, seasonal = FALSE)

cat("Training model summary:\n")
print(summary(fit_train))

# out-of-sample forecast on test set
fc_test <- forecast(fit_train, xreg = xreg_test, h = nrow(xreg_test))

# y_hat <- as.numeric(fc_test$mean)
y_hat <- fc_test$mean

# basic error metrics
rmse <- sqrt(mean((y_hat - y_test)^2, na.rm = TRUE))
mae  <- mean(abs(y_hat - y_test), na.rm = TRUE)
dir_acc <- mean(sign(y_hat) == sign(y_test), na.rm = TRUE)

cat("Out-of-sample performance;\n")
cat("RMSE:", rmse, "\n")
cat("MAE: ", mae,  "\n")
cat("Directional accuracy:", round(100 * dir_acc, 2), "%\n")

# refit on full data after satisfactory validation performance
fit <- auto.arima(y, xreg = xreg, seasonal = FALSE)

cat("\nFull-sample model (used for scenario analysis):\n")
print(summary(fit))
print(checkresiduals(fit))

scenario_analysis <- function(nfp_surprise = 0,
                              core_CPE_surprise = 0,
                              core_retail_sales_surprise = 0,
                              scaled = FALSE) {

  # collect scenario inputs
  scenario <- c(
    NFP_surprise               = nfp_surprise,
    core_CPE_surprise          = core_CPE_surprise,
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
  z    <- qnorm(0.975)
  lo95 <- fc$lower[1, "95%"]
  hi95 <- fc$upper[1, "95%"]
  
  sd_from_hi <- (hi95 - mu) / z
  sd_from_lo <- (mu - lo95) / z
  sd_hat     <- mean(c(sd_from_hi, sd_from_lo))
  
  # return predictive mean and sd
  return(list(mu = mu, sd = sd_hat))
}


# print("function test:")
# print(scenario_analysis(nfp_surprise = 100000,
#                                core_CPE_surprise = 0,
#                                core_retail_sales_surprise = 0))

plot_scenario <- function(nfp_surprise = 0,
                          core_CPE_surprise = 0,
                          core_retail_sales_surprise = 0,
                          scaled = FALSE) {
  
  # Call your updated scenario_analysis()
  scenario_res <- scenario_analysis(
    nfp_surprise = nfp_surprise,
    core_CPE_surprise = core_CPE_surprise,
    core_retail_sales_surprise = core_retail_sales_surprise,
    scaled = scaled
  )
  
  mu <- scenario_res$mu
  sd <- scenario_res$sd
  
  # Monte Carlo sample from predictive distribution
  set.seed(42)
  N <- 10000
  sim_ret <- rnorm(N, mean = mu, sd = sd)
  
  # Plot histogram + density of simulated returns (in %)
  sim_df <- data.frame(ret = 100 * sim_ret)
  
  p <- ggplot(sim_df, aes(x = ret)) +
    geom_histogram(aes(y = ..density..), bins = 50) +
    geom_density(linewidth = 1) +
    geom_vline(xintercept = mean(sim_df$ret),
               linetype = "dashed", linewidth = 1) +
    labs(
      title = sprintf(
        "Predictive Distribution of Daily SPX Returns\nunder Macro Surprise Scenario\nNFP: %s, Core CPE: %s, Core Retail Sales: %s\n(scaled = %s)",
        nfp_surprise, core_CPE_surprise, core_retail_sales_surprise, scaled
      ),
      x = "Simulated daily return (%)",
      y = "Density"
    )
  return(p)
}

# Sequence of surprises in SD units (z-scores)
z_vals <- seq(-3, 3, by = 0.1)

# Store expected returns here
mu_vals <- numeric(length(z_vals))

# loop through z values, holding other indicators at 0
for (i in seq_along(z_vals)) {
  res <- scenario_analysis(
    nfp_surprise = z_vals[i],
    core_CPE_surprise = 0,
    core_retail_sales_surprise = 0,
    scaled = TRUE 
  )
  mu_vals[i] <- 100 * res$mu  # convert to %
}

# plot nfp reaction curve
df_nfp <- data.frame(z = z_vals, mu_pct = mu_vals)

nfp_reaction_plot <- ggplot(df_nfp, aes(x = z, y = mu_pct)) + geom_line(linewidth = 1) + coord_cartesian(ylim = c(0, 0.25)) +
  labs(
    title = "SPX Expected Return vs. NFP Surprise (in SD units)",
    x = "NFP Surprise (z‐score)",
    y = "Expected daily return (%)"
  )


print(df_nfp)

for (i in seq_along(z_vals)) {
  res <- scenario_analysis(
    nfp_surprise = 0,
    core_CPE_surprise = z_vals[i],
    core_retail_sales_surprise = 0,
    scaled = TRUE
  )
  mu_vals[i] <- 100 * res$mu
}

# plot cpe reaction curve
df_cpe <- data.frame(z = z_vals, mu_pct = mu_vals)

cpe_reaction_plot <- ggplot(df_cpe, aes(x = z, y = mu_pct)) + geom_line(linewidth = 1) + coord_cartesian(ylim = c(0, 0.25)) +
  labs(
    title = "SPX Expected Return vs. Core CPE Surprise (in SD units)",
    x = "Core CPI Surprise (z‐score)",
    y = "Expected daily return (%)"
  )


print(df_cpe)


for (i in seq_along(z_vals)) {
  res <- scenario_analysis(
    nfp_surprise = 0,
    core_CPE_surprise = 0,
    core_retail_sales_surprise = z_vals[i],
    scaled = TRUE
  )
  mu_vals[i] <- 100 * res$mu
}

# plot retail sales reaction curve
df_retail <- data.frame(z = z_vals, mu_pct = mu_vals)

retail_reaction_plot <- ggplot(df_retail, aes(x = z, y = mu_pct)) + geom_line(linewidth = 1) + coord_cartesian(ylim = c(0, 0.25)) +
  labs(
    title = "SPX Expected Return vs. Retail Sales Surprise (in SD units)",
    x = "Retail Surprise (z‐score)",
    y = "Expected daily return (%)"
  )

print(df_retail)

print(nfp_reaction_plot / cpe_reaction_plot / retail_reaction_plot)