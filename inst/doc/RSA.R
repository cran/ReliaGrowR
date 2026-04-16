## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(ReliaGrowR)

id   <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3)
time <- c(100, 350, 500, 80, 300, 600, 150, 250, 400, 700)

## -----------------------------------------------------------------------------
result <- mcf(id, time)
plot(result, main = "Mean Cumulative Function",
     xlab = "Time", ylab = "MCF")

## -----------------------------------------------------------------------------
df <- data.frame(id = id, time = time)
result2 <- mcf(data = df)

## -----------------------------------------------------------------------------
id    <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
time  <- c(100, 350, 500, 80, 300, 400, 150, 250, 700)
event <- c(  1,   1,   0,  1,   1,   0,   1,   1,   1)

result <- mcf(id, time, event)

## -----------------------------------------------------------------------------
id   <- c(1, 1, 2, 2, 3, 3)
time <- c(100, 300, 150, 400, 200, 350)

# Without end_time: system observation ends at last event
mcf_basic <- mcf(id, time)

# With end_time: all systems observed until time 800
mcf_adj <- mcf(id, time, end_time = c("1" = 800, "2" = 800, "3" = 800))

## -----------------------------------------------------------------------------
par(mfrow = c(1, 2))
plot(mcf_basic, main = "MCF (inferred exposure)",
     xlab = "Time", ylab = "MCF")
plot(mcf_adj, main = "MCF (explicit exposure)",
     xlab = "Time", ylab = "MCF")

## -----------------------------------------------------------------------------
id    <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
time  <- c(100, 350, 500, 80, 300, 400, 150, 250, 700)
event <- c(  1,   1,   0,  1,   1,   0,   1,   1,   1)

exp_result <- exposure(id, time, event)
mcf_result <- mcf(id, time, event, end_time = exp_result$end_times)

## -----------------------------------------------------------------------------
id   <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3)
time <- c(100, 350, 500, 80, 300, 600, 150, 250, 400, 700)
result <- exposure(id, time)

## ----fig.height=8-------------------------------------------------------------
plot(result)

## -----------------------------------------------------------------------------
plot(result, which = "exposure")

## -----------------------------------------------------------------------------
plot(result, which = "event_rate")

## -----------------------------------------------------------------------------
id    <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
time  <- c(100, 350, 500, 80, 300, 400, 150, 250, 700)
event <- c(  1,   1,   0,  1,   1,   0,   1,   1,   1)

result <- exposure(id, time, event)

## -----------------------------------------------------------------------------
time  <- c(200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000)
event <- c(  3,   5,   4,   7,    6,    8,    5,    9,    7,   10)

## -----------------------------------------------------------------------------
fit_mle <- nhpp(time, event, method = "MLE")
plot(fit_mle, main = "Power Law NHPP (MLE)",
     xlab = "Cumulative Time", ylab = "Cumulative Events")

## -----------------------------------------------------------------------------
fit_ls <- nhpp(time, event, method = "LS")

## -----------------------------------------------------------------------------
result_ll <- nhpp(time, event, model_type = "Log-Linear")
plot(result_ll, main = "Log-Linear NHPP",
     xlab = "Cumulative Time", ylab = "Cumulative Events")

## -----------------------------------------------------------------------------
time  <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
           1100, 1200, 1300, 1400, 1500)
event <- c(  1,   1,   2,   4,   4,   1,   1,   2,   1,    4,
             1,   1,   3,   3,   4)

result_pw <- nhpp(time, event, breaks = c(500), method = "LS")
plot(result_pw, main = "Piecewise Power Law NHPP",
     xlab = "Cumulative Time", ylab = "Cumulative Events")

## -----------------------------------------------------------------------------
time  <- c(200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000)
event <- c(  3,   5,   4,   7,    6,    8,    5,    9,    7,   10)

fit <- nhpp(time, event, method = "MLE")
fc <- predict_nhpp(fit, time = c(2001, 2500, 3000, 4000, 5000))

## -----------------------------------------------------------------------------
plot(fc, main = "NHPP Forecast",
     xlab = "Cumulative Time", ylab = "Cumulative Events")

## -----------------------------------------------------------------------------
fit_ll <- nhpp(time, event, model_type = "Log-Linear")
fc_ll <- predict_nhpp(fit_ll, time = c(2500, 3000))

