## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message = FALSE, warning = FALSE----------------------------------
library(ReliaGrowR)
library(WeibullR)

## ----fleet-data---------------------------------------------------------------
set.seed(8)
n_units <- 100
test_end <- 1000
beta_mech <- 4 # wear-out failure mechanism (shape), constant across builds
eta_matured <- 1456 # characteristic life of the final, matured design
fleet_life <- eta_matured * (-log(runif(n_units)))^(1 / beta_mech)

failures <- sort(fleet_life[fleet_life <= test_end])
suspensions <- rep(test_end, sum(fleet_life > test_end))

n_test_failures <- length(failures)
n_surviving <- length(suspensions)

## ----weibull-fit--------------------------------------------------------------
obj <- wblr(failures, suspensions,
  col = "steelblue", label = "Fielded Fleet", is.plot.cb = FALSE
)
obj <- wblr.fit(obj, method.fit = "mle")
sim_beta <- obj$fit[[1]]$beta
eta_orig <- obj$fit[[1]]$eta
plot(obj,
  main = "Weibull Fit to Fielded Fleet (No-Growth Baseline)",
  xlab = "Operating Time", ylab = "Unreliability"
)

## ----dev-test-----------------------------------------------------------------
set.seed(3)
N_dev <- 30 # prototypes built and tested in sequence
g_exp <- 0.8 # strength of design-life improvement with build index

# Design life improves with build index; the last build is the matured design.
eta_build <- eta_matured * (seq_len(N_dev) / N_dev)^g_exp
dev_life <- eta_build * (-log(runif(N_dev)))^(1 / beta_mech)
dev_cum <- cumsum(dev_life)

## ----dev-fit, results = "asis"------------------------------------------------
dev_fit <- rga(
  data.frame(times = dev_cum, failures = rep(1, N_dev)),
  times_type = "cumulative_failure_times"
)
growth_beta <- as.numeric(dev_fit$betas)
dev_r2 <- summary(dev_fit$model)$r.squared

plot(dev_fit,
  main = "Developmental Reliability Growth Test",
  xlab = "Cumulative Test Time", ylab = "Cumulative Failures"
)

## ----forecast-calc------------------------------------------------------------
# Forecast window: a modest forward increment (~40% of the observation period) applied to
# each surviving unit, chosen rather than solved for.
forecast_window <- 400
test_end_cum_time <- sum(c(failures, suspensions))

# No-growth count the fitted fleet Weibull expects over the window (forward hazard).
expected_failures <- function(w, eta, r, beta) {
  sum(1 - exp(-(((r + w) / eta)^beta - (r / eta)^beta)))
}
# Scale that reproduces a target expected count over the window.
solve_eta <- function(target, r, beta, window) {
  uniroot(function(eta) expected_failures(window, eta, r, beta) - target,
    interval = c(1e-3, 1e9)
  )$root
}
# Crow-AMSAA growth count-ratio: growth lowers the expected count; equals 1 at beta_growth = 1.
count_ratio <- function(bg, rho) ((1 + rho)^bg - 1) / rho

f0 <- expected_failures(forecast_window, eta_orig, rep(test_end, n_surviving), sim_beta)
rho <- (n_surviving * forecast_window) / test_end_cum_time
n_expected <- f0 * count_ratio(growth_beta, rho)
n_forecast <- round(n_expected)
n_target <- n_test_failures + n_forecast

## ----forecast-plot, warning = FALSE-------------------------------------------
dev_target <- N_dev + n_forecast
t_f_dev <- (dev_target / as.numeric(dev_fit$lambdas))^(1 / growth_beta)
dev_forecast_times <- seq(max(dev_cum) * 1.001, t_f_dev, length.out = 50)

fc <- predict_rga(dev_fit, times = dev_forecast_times)
plot(fc,
  main = "Reliability Growth Forecast",
  xlab = "Cumulative Test Time", ylab = "Cumulative Failures"
)

## -----------------------------------------------------------------------------
set.seed(123)
sim_eta <- solve_eta(n_expected, rep(test_end, n_surviving), sim_beta, forecast_window)
sim_result <- sim_failures(
  n = n_forecast,
  runtimes = rep(1000, n_surviving),
  window = forecast_window,
  beta = sim_beta,
  eta = sim_eta
)

## -----------------------------------------------------------------------------
sim_fail_times <- sim_result$runtime[sim_result$type == "Failure"]
sim_susp_times <- sim_result$runtime[sim_result$type == "Suspension"]
combined_failures <- c(failures, sim_fail_times)
combined_suspensions <- sim_susp_times

## -----------------------------------------------------------------------------
# Headline comparison is like-for-like: a 2-parameter Weibull fitted by MLE to both the
# fielded fleet (no-growth baseline) and the growth-augmented data, so the two fits differ
# only in their data, not in parameterization or estimation method.
obj_nogrowth <- wblr(failures, suspensions,
  col = "grey40", label = "Without Growth", is.plot.cb = FALSE
)
obj_nogrowth <- wblr.fit(obj_nogrowth, method.fit = "mle")

obj_growth <- wblr(combined_failures, combined_suspensions,
  col = "tomato", label = "With Growth (2P)", is.plot.cb = FALSE
)
obj_growth <- wblr.fit(obj_growth, method.fit = "mle")

plot.wblr(list(obj_nogrowth, obj_growth),
  main = "Weibull Comparison: With vs. Without Reliability Growth",
  is.plot.legend = TRUE
)

## -----------------------------------------------------------------------------
growth_wb_beta <- obj_growth$fit[[1]]$beta
growth_wb_eta <- obj_growth$fit[[1]]$eta
nogrowth_wb_beta <- obj_nogrowth$fit[[1]]$beta
nogrowth_wb_eta <- obj_nogrowth$fit[[1]]$eta

# B10 life (10th-percentile life): t such that F(t) = 0.10.
b10 <- function(beta, eta, t0 = 0) t0 + eta * (-log(0.90))^(1 / beta)
nogrowth_b10 <- b10(nogrowth_wb_beta, nogrowth_wb_eta)
growth_b10 <- b10(growth_wb_beta, growth_wb_eta)

params <- data.frame(
  Scenario = c("Without Growth (2P)", "With Growth (2P)"),
  Beta = round(c(nogrowth_wb_beta, growth_wb_beta), 3),
  Eta = round(c(nogrowth_wb_eta, growth_wb_eta), 1),
  B10 = round(c(nogrowth_b10, growth_b10), 1)
)

knitr::kable(params,
  caption = "Weibull parameters and B10 life: with vs. without reliability growth (two-parameter MLE)",
  col.names = c("Scenario", "β", "η", "B10"),
  row.names = FALSE
)

## ----mc-helpers, include = FALSE----------------------------------------------
# Given a target expected count over the window, calibrate the scale, draw round(target)
# failures, combine with the observed failures, and refit a two-parameter Weibull.
sim_and_fit <- function(test_failures, target_exp, runtimes, window, beta) {
  eta_c <- solve_eta(target_exp, runtimes, beta, window)
  sim_i <- sim_failures(
    n = round(target_exp), runtimes = runtimes, window = window,
    beta = beta, eta = eta_c
  )
  sim_f <- sim_i$runtime[sim_i$type == "Failure"]
  sim_s <- sim_i$runtime[sim_i$type == "Suspension"]
  obj_i <- wblr(c(test_failures, sim_f), sim_s, is.plot.cb = FALSE)
  obj_i <- wblr.fit(obj_i, method.fit = "mle")
  data.frame(beta = obj_i$fit[[1]]$beta, eta = obj_i$fit[[1]]$eta)
}

## ----mc-loop, message = FALSE, warning = FALSE, results = "hide"--------------
set.seed(99)
n_mc <- 500
mc_results <- vector("list", n_mc)

# Sampling distribution of the developmental growth parameter.
betag_hat <- as.numeric(dev_fit$betas)
betag_se <- as.numeric(dev_fit$betas_se)

for (i in seq_len(n_mc)) {
  tryCatch(
    {
      # Draw a growth parameter and recompute the growth-adjusted count over the window.
      beta_i <- rnorm(1, betag_hat, betag_se)
      if (beta_i <= 0) stop("non-positive growth parameter draw")
      target_i <- f0 * count_ratio(beta_i, rho)
      if (!is.finite(target_i) || target_i <= 0) stop("non-positive count draw")

      mc_results[[i]] <- sim_and_fit(
        failures, target_i, rep(1000, n_surviving), forecast_window, sim_beta
      )
    },
    error = function(e) NULL
  )
}

mc_df <- do.call(rbind, Filter(Negate(is.null), mc_results))

## ----mc-histograms, fig.height = 4--------------------------------------------
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

hist(mc_df$beta,
  breaks = "Sturges", col = "steelblue", border = "white",
  main = "MC Distribution of β",
  xlab = expression(hat(beta)), ylab = "Count", freq = TRUE
)
abline(v = median(mc_df$beta), lty = 2, lwd = 2)
abline(v = nogrowth_wb_beta, lty = 3, lwd = 2, col = "grey40")

hist(mc_df$eta,
  breaks = "Sturges", col = "tomato", border = "white",
  main = "MC Distribution of η",
  xlab = expression(hat(eta)), ylab = "Count", freq = TRUE
)
abline(v = median(mc_df$eta), lty = 2, lwd = 2)
abline(v = nogrowth_wb_eta, lty = 3, lwd = 2, col = "grey40")

par(mfrow = c(1, 1))

## ----mc-table-----------------------------------------------------------------
mc_df$b10 <- b10(mc_df$beta, mc_df$eta)

mc_summary <- data.frame(
  Parameter = c("β", "η", "B10"),
  NoGrowth = round(c(nogrowth_wb_beta, nogrowth_wb_eta, nogrowth_b10), 3),
  Median = round(c(median(mc_df$beta), median(mc_df$eta), median(mc_df$b10)), 3),
  Mean = round(c(mean(mc_df$beta), mean(mc_df$eta), mean(mc_df$b10)), 3),
  CI_lo = round(c(
    quantile(mc_df$beta, 0.025),
    quantile(mc_df$eta, 0.025),
    quantile(mc_df$b10, 0.025)
  ), 3),
  CI_hi = round(c(
    quantile(mc_df$beta, 0.975),
    quantile(mc_df$eta, 0.975),
    quantile(mc_df$b10, 0.975)
  ), 3)
)

knitr::kable(mc_summary,
  caption = paste0(
    "Monte Carlo summary of Weibull parameters and B10 life (",
    nrow(mc_df), " valid iterations)"
  ),
  col.names = c(
    "Parameter", "No Growth", "Median", "Mean", "2.5%", "97.5%"
  ),
  row.names = FALSE
)

## ----sens-fail-setup----------------------------------------------------------
assessment_times <- c(850, 950, 1050, 1150)
fleet_obs <- vapply(assessment_times, function(te) sum(fleet_life <= te), integer(1))
fleet_labels <- paste0(fleet_obs, " failures")

## ----sens-fail-mc, message = FALSE, warning = FALSE, results = "hide"---------
set.seed(42)
n_mc_sens <- 200

sens_fail_list <- lapply(seq_along(assessment_times), function(k) {
  te <- assessment_times[k]
  fail_te <- sort(fleet_life[fleet_life <= te])
  susp_te <- rep(te, sum(fleet_life > te))
  n_obs <- length(fail_te)
  n_surv <- length(susp_te)

  obj_te <- wblr.fit(wblr(fail_te, susp_te, is.plot.cb = FALSE), method.fit = "mle")
  beta_te <- obj_te$fit[[1]]$beta
  eta_te <- obj_te$fit[[1]]$eta

  if (n_surv <= 0) {
    return(NULL)
  }
  # No-growth count over the fixed window from this assessment's own fitted baseline,
  # lowered by the growth count-ratio at this assessment's cumulative exposure.
  f0_te <- expected_failures(forecast_window, eta_te, rep(te, n_surv), beta_te)
  rho_te <- (n_surv * forecast_window) / sum(c(fail_te, susp_te))
  target_te <- f0_te * count_ratio(growth_beta, rho_te)
  if (!is.finite(target_te) || target_te <= 0) {
    return(NULL)
  }

  rows <- lapply(seq_len(n_mc_sens), function(i) {
    tryCatch(
      {
        r <- sim_and_fit(fail_te, target_te, rep(te, n_surv), forecast_window, beta_te)
        cbind(scenario = fleet_labels[k], r)
      },
      error = function(e) NULL
    )
  })
  do.call(rbind, Filter(Negate(is.null), rows))
})

sens_fail_df <- do.call(rbind, Filter(Negate(is.null), sens_fail_list))
sens_fail_df$scenario <- factor(sens_fail_df$scenario, levels = fleet_labels)

## ----sens-fail-boxplot, fig.height = 4----------------------------------------
par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))

fail_cols <- c("#2E86C1", "#27AE60", "#F39C12", "#C0392B")
boxplot(beta ~ scenario,
  data = sens_fail_df,
  col = fail_cols, outline = FALSE,
  main = "Fitted β by Failures Observed",
  xlab = "Failures observed at assessment", ylab = expression(hat(beta)),
  las = 2, cex.axis = 0.8
)

boxplot(eta ~ scenario,
  data = sens_fail_df,
  col = fail_cols, outline = FALSE,
  main = "Fitted η by Failures Observed",
  xlab = "Failures observed at assessment", ylab = expression(hat(eta)),
  las = 2, cex.axis = 0.8
)

par(mfrow = c(1, 1))

## ----growth-setup-------------------------------------------------------------
growth_scenarios <- c(0.4, 0.6, 0.8, 1.0)
growth_labels <- c(
  "βg = 0.4 (strong)",
  "βg = 0.6 (moderate)",
  "βg = 0.8 (mild)",
  "βg = 1.0 (none)"
)

## ----growth-mc, message = FALSE, warning = FALSE, results = "hide"------------
set.seed(77)

sens_growth_list <- lapply(seq_along(growth_scenarios), function(k) {
  gb <- growth_scenarios[k]

  # Lower the fixed-window no-growth count by the growth ratio for this scenario.
  target_k <- f0 * count_ratio(gb, rho)

  if (target_k <= 0) {
    return(NULL)
  }

  rows <- lapply(seq_len(n_mc_sens), function(i) {
    tryCatch(
      {
        r <- sim_and_fit(
          failures, target_k, rep(1000, n_surviving), forecast_window, sim_beta
        )
        cbind(scenario = growth_labels[k], r)
      },
      error = function(e) NULL
    )
  })
  do.call(rbind, Filter(Negate(is.null), rows))
})

sens_growth_df <- do.call(rbind, Filter(Negate(is.null), sens_growth_list))
sens_growth_df$scenario <- factor(sens_growth_df$scenario, levels = growth_labels)

## ----growth-boxplot, fig.height = 4-------------------------------------------
par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))

growth_cols <- c("#2E86C1", "#27AE60", "#F39C12", "#C0392B")
boxplot(beta ~ scenario,
  data = sens_growth_df,
  col = growth_cols, outline = FALSE,
  main = "Fitted β by Growth Strength",
  xlab = "Growth scenario", ylab = expression(hat(beta)),
  las = 2, cex.axis = 0.75
)

boxplot(eta ~ scenario,
  data = sens_growth_df,
  col = growth_cols, outline = FALSE,
  main = "Fitted η by Growth Strength",
  xlab = "Growth scenario", ylab = expression(hat(eta)),
  las = 2, cex.axis = 0.75
)

par(mfrow = c(1, 1))

