## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message = FALSE, warning = FALSE----------------------------------
# Prefer the local source checkout when knitting from the package repo.
pkg_root <- NULL
if (requireNamespace("pkgload", quietly = TRUE)) {
  candidate_roots <- c(".", "..")
  for (path in candidate_roots) {
    desc_path <- file.path(path, "DESCRIPTION")
    if (!file.exists(desc_path)) {
      next
    }

    pkg_name <- tryCatch(
      as.character(read.dcf(desc_path, fields = "Package")[1, 1]),
      error = function(e) ""
    )
    if (identical(pkg_name, "ReliaGrowR")) {
      pkg_root <- normalizePath(path, winslash = "/", mustWork = TRUE)
      break
    }
  }
}

if (!is.null(pkg_root)) {
  pkgload::load_all(
    pkg_root,
    export_all = FALSE,
    helpers = FALSE,
    attach_testthat = FALSE,
    quiet = TRUE
  )
} else {
  library(ReliaGrowR)
}

library(WeibullR)

## -----------------------------------------------------------------------------
failures <- c(
  500, 600, 650, 700, 730, 760, 790, 810, 830, 845,
  860, 875, 890, 905, 920, 935, 945, 960, 975, 990
)
suspensions <- rep(1000, 80)

n_test_failures <- length(failures)
n_surviving <- length(suspensions)

## -----------------------------------------------------------------------------
obj <- wblr(failures, suspensions,
  col = "steelblue", label = "Test Data", is.plot.cb = FALSE
)
obj <- wblr.fit(obj)
sim_beta <- obj$fit[[1]]$beta
eta_orig <- obj$fit[[1]]$eta
plot(obj,
  main = "Weibull Fit to Test Data",
  xlab = "Operating Time", ylab = "Unreliability"
)

## -----------------------------------------------------------------------------
rga_data <- weibull_to_rga(failures, suspensions)
rga_input <- data.frame(
  times = rga_data$CumulativeTime,
  failures = rga_data$Failures
)

## -----------------------------------------------------------------------------
test_end_cum_time <- max(rga_input$times)
fit <- rga(rga_input, times_type = "cumulative_failure_times")
growth_beta <- fit$betas
growth_lambda <- fit$lambdas

## ----results = "asis"---------------------------------------------------------
plot(fit,
  main = "Reliability Growth Analysis",
  xlab = "Cumulative Time", ylab = "Cumulative Failures"
)

## -----------------------------------------------------------------------------
n_forecast <- 60
n_target <- n_test_failures + n_forecast

t_forecast <- (n_target / growth_lambda)^(1 / growth_beta)
delta_cum_time <- t_forecast - test_end_cum_time
effective_fleet <- n_surviving - n_forecast / 2
per_unit_time <- delta_cum_time / effective_fleet

## ----warnings = FALSE---------------------------------------------------------
forecast_times <- seq(test_end_cum_time, t_forecast, length.out = 50)
fc <- predict_rga(fit, forecast_times)
plot(fc,
  main = "Reliability Growth Forecast",
  xlab = "Cumulative Time", ylab = "Cumulative Failures"
)

## -----------------------------------------------------------------------------
set.seed(123)
sim_result <- sim_failures(
  n = n_forecast,
  runtimes = rep(1000, n_surviving),
  window = per_unit_time,
  beta = sim_beta
)
sim_eta <- attr(sim_result, "weibull_eta")

## -----------------------------------------------------------------------------
sim_fail_times <- sim_result$runtime[sim_result$type == "Failure"]
sim_susp_times <- sim_result$runtime[sim_result$type == "Suspension"]
combined_failures <- c(failures, sim_fail_times)
combined_suspensions <- sim_susp_times

## -----------------------------------------------------------------------------
obj_growth <- wblr(combined_failures, combined_suspensions,
  col = "tomato", label = "With Growth", is.plot.cb = FALSE
)
obj_growth <- wblr.fit(obj_growth)

obj_nogrowth <- wblr(failures, suspensions,
  col = "grey40", label = "Without Growth", is.plot.cb = FALSE
)
obj_nogrowth <- wblr.fit(obj_nogrowth)

plot.wblr(list(obj_nogrowth, obj_growth),
  main = "Weibull Comparison: With vs. Without Reliability Growth",
  is.plot.legend = TRUE
)

## -----------------------------------------------------------------------------
growth_wb_beta <- obj_growth$fit[[1]]$beta
growth_wb_eta <- obj_growth$fit[[1]]$eta
nogrowth_wb_beta <- obj_nogrowth$fit[[1]]$beta
nogrowth_wb_eta <- obj_nogrowth$fit[[1]]$eta

params <- data.frame(
  Scenario = c("Without Growth", "With Growth"),
  Beta = round(c(nogrowth_wb_beta, growth_wb_beta), 3),
  Eta = round(c(nogrowth_wb_eta, growth_wb_eta), 1)
)

knitr::kable(params,
  caption = "Weibull parameters: with vs. without reliability growth",
  col.names = c("Scenario", "\u03b2", "\u03b7"),
  row.names = FALSE
)

## ----mc-helpers, include = FALSE----------------------------------------------
sim_and_fit <- function(test_failures, n_sim, runtimes, window, beta) {
  sim_i <- sim_failures(
    n = n_sim, runtimes = runtimes, window = window, beta = beta
  )
  sim_f <- sim_i$runtime[sim_i$type == "Failure"]
  sim_s <- sim_i$runtime[sim_i$type == "Suspension"]
  obj_i <- wblr(c(test_failures, sim_f), sim_s, is.plot.cb = FALSE)
  obj_i <- wblr.fit(obj_i)
  data.frame(beta = obj_i$fit[[1]]$beta, eta = obj_i$fit[[1]]$eta)
}

## ----mc-loop, message = FALSE, warning = FALSE, results = "hide"--------------
set.seed(99)
n_mc <- 500
mc_results <- vector("list", n_mc)

for (i in seq_len(n_mc)) {
  tryCatch(
    {
      mc_results[[i]] <- sim_and_fit(
        failures, n_forecast, rep(1000, n_surviving), per_unit_time, sim_beta
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
  main = "MC Distribution of \u03b2",
  xlab = expression(hat(beta)), ylab = "Count", freq = TRUE
)
abline(v = median(mc_df$beta), lty = 2, lwd = 2)
abline(v = nogrowth_wb_beta, lty = 3, lwd = 2, col = "grey40")

hist(mc_df$eta,
  breaks = "Sturges", col = "tomato", border = "white",
  main = "MC Distribution of \u03b7",
  xlab = expression(hat(eta)), ylab = "Count", freq = TRUE
)
abline(v = median(mc_df$eta), lty = 2, lwd = 2)
abline(v = nogrowth_wb_eta, lty = 3, lwd = 2, col = "grey40")

par(mfrow = c(1, 1))

## ----mc-table-----------------------------------------------------------------
mc_summary <- data.frame(
  Parameter = c("\u03b2", "\u03b7"),
  NoGrowth = round(c(nogrowth_wb_beta, nogrowth_wb_eta), 3),
  Median = round(c(median(mc_df$beta), median(mc_df$eta)), 3),
  Mean = round(c(mean(mc_df$beta), mean(mc_df$eta)), 3),
  CI_lo = round(c(
    quantile(mc_df$beta, 0.025),
    quantile(mc_df$eta, 0.025)
  ), 3),
  CI_hi = round(c(
    quantile(mc_df$beta, 0.975),
    quantile(mc_df$eta, 0.975)
  ), 3)
)

knitr::kable(mc_summary,
  caption = paste0(
    "Monte Carlo summary of Weibull parameters (",
    nrow(mc_df), " valid iterations)"
  ),
  col.names = c(
    "Parameter", "No Growth", "Median", "Mean", "2.5%", "97.5%"
  ),
  row.names = FALSE
)

## ----sens-fail-setup----------------------------------------------------------
# Build failure scenarios by subsampling or extending the base failures
fail_10 <- failures[seq(1, 20, by = 2)]
fail_20 <- failures
fail_30 <- sort(c(failures, seq(505, 995, length.out = 10)))
fail_40 <- sort(c(failures, seq(510, 998, length.out = 20)))

test_fail_scenarios <- list(
  "10 failures" = list(f = fail_10, s = rep(1000, 90)),
  "20 failures" = list(f = fail_20, s = rep(1000, 80)),
  "30 failures" = list(f = fail_30, s = rep(1000, 70)),
  "40 failures" = list(f = fail_40, s = rep(1000, 60))
)

## ----sens-fail-mc, message = FALSE, warning = FALSE, results = "hide"---------
set.seed(42)
n_mc_sens <- 200

sens_fail_list <- lapply(names(test_fail_scenarios), function(sc_name) {
  sc <- test_fail_scenarios[[sc_name]]
  n_f <- length(sc$f)
  n_s <- length(sc$s)
  n_fc <- round(n_s * 0.75) # forecast 75% of survivors

  tryCatch(
    {
      # Fit Weibull
      obj_sc <- wblr(sc$f, sc$s, is.plot.cb = FALSE)
      obj_sc <- wblr.fit(obj_sc)
      beta_sc <- obj_sc$fit[[1]]$beta

      # Fit RGA
      rga_sc <- weibull_to_rga(sc$f, sc$s)
      rga_in <- data.frame(
        times = rga_sc$CumulativeTime,
        failures = rga_sc$Failures
      )
      fit_sc <- rga(rga_in, times_type = "cumulative_failure_times")

      # Forecast
      t_end_sc <- max(rga_in$times)
      n_tgt <- n_f + n_fc
      t_fc <- (n_tgt / fit_sc$lambdas)^(1 / fit_sc$betas)
      delta_sc <- t_fc - t_end_sc
      effective_sc <- n_s - n_fc / 2
      put_sc <- delta_sc / effective_sc

      if (put_sc <= 0) {
        return(NULL)
      }

      # MC loop
      rows <- lapply(seq_len(n_mc_sens), function(i) {
        tryCatch(
          {
            r <- sim_and_fit(sc$f, n_fc, rep(1000, n_s), put_sc, beta_sc)
            cbind(scenario = sc_name, r)
          },
          error = function(e) NULL
        )
      })
      do.call(rbind, Filter(Negate(is.null), rows))
    },
    error = function(e) NULL
  )
})

sens_fail_df <- do.call(rbind, Filter(Negate(is.null), sens_fail_list))
sens_fail_df$scenario <- factor(
  sens_fail_df$scenario,
  levels = names(test_fail_scenarios)
)

## ----sens-fail-boxplot, fig.height = 4----------------------------------------
par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))

fail_cols <- c("#2E86C1", "#27AE60", "#F39C12", "#C0392B")
boxplot(beta ~ scenario,
  data = sens_fail_df,
  col = fail_cols, outline = FALSE,
  main = "Fitted \u03b2 by Test Failures",
  xlab = "Scenario", ylab = expression(hat(beta)),
  las = 2, cex.axis = 0.8
)

boxplot(eta ~ scenario,
  data = sens_fail_df,
  col = fail_cols, outline = FALSE,
  main = "Fitted \u03b7 by Test Failures",
  xlab = "Scenario", ylab = expression(hat(eta)),
  las = 2, cex.axis = 0.8
)

par(mfrow = c(1, 1))

## ----sens-helpers, include = FALSE--------------------------------------------
summarize_mc <- function(df) {
  do.call(rbind, lapply(levels(df$scenario), function(lbl) {
    sub <- df[df$scenario == lbl, ]
    if (nrow(sub) == 0) return(NULL)
    data.frame(
      Scenario = lbl, Valid = nrow(sub),
      Beta_med = round(median(sub$beta), 3),
      Eta_med = round(median(sub$eta), 1),
      Eta_lo = round(quantile(sub$eta, 0.025), 1),
      Eta_hi = round(quantile(sub$eta, 0.975), 1)
    )
  }))
}

## ----sens-fail-table----------------------------------------------------------
sens_fail_tbl <- summarize_mc(sens_fail_df)

knitr::kable(sens_fail_tbl,
  caption = paste0(
    "Growth-adjusted Weibull by test-failure scenario (",
    n_mc_sens, " iterations each)"
  ),
  col.names = c(
    "Scenario", "Valid", "Med. \u03b2",
    "Med. \u03b7", "2.5% \u03b7", "97.5% \u03b7"
  ),
  row.names = FALSE
)

## ----growth-setup-------------------------------------------------------------
growth_scenarios <- c(0.4, 0.6, 0.8, 1.0)
growth_labels <- c(
  "\u03b2g = 0.4 (strong)",
  "\u03b2g = 0.6 (moderate)",
  "\u03b2g = 0.8 (mild)",
  "\u03b2g = 1.0 (none)"
)

## ----growth-mc, message = FALSE, warning = FALSE, results = "hide"------------
set.seed(77)

sens_growth_list <- lapply(seq_along(growth_scenarios), function(k) {
  gb <- growth_scenarios[k]

  # Re-anchor lambda so model matches observed data, then forecast
  lambda_k <- n_test_failures / test_end_cum_time^gb
  t_fc_k <- (n_target / lambda_k)^(1 / gb)
  delta_k <- t_fc_k - test_end_cum_time
  put_k <- delta_k / effective_fleet

  if (put_k <= 0) {
    return(NULL)
  }

  rows <- lapply(seq_len(n_mc_sens), function(i) {
    tryCatch(
      {
        r <- sim_and_fit(
          failures, n_forecast, rep(1000, n_surviving), put_k, sim_beta
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
  main = "Fitted \u03b2 by Growth Strength",
  xlab = "Growth scenario", ylab = expression(hat(beta)),
  las = 2, cex.axis = 0.75
)

boxplot(eta ~ scenario,
  data = sens_growth_df,
  col = growth_cols, outline = FALSE,
  main = "Fitted \u03b7 by Growth Strength",
  xlab = "Growth scenario", ylab = expression(hat(eta)),
  las = 2, cex.axis = 0.75
)

par(mfrow = c(1, 1))

## ----growth-table-------------------------------------------------------------
sens_growth_tbl <- summarize_mc(sens_growth_df)

knitr::kable(sens_growth_tbl,
  caption = paste0(
    "Growth-adjusted Weibull by growth parameter scenario (",
    n_mc_sens, " iterations each)"
  ),
  col.names = c(
    "Scenario", "Valid", "Med. \u03b2",
    "Med. \u03b7", "2.5% \u03b7", "97.5% \u03b7"
  ),
  row.names = FALSE
)

