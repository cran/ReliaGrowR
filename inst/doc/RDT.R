## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----setup--------------------------------------------------------------------
library(ReliaGrowR)

## ----solve-time---------------------------------------------------------------
plan <- rdt(
  target       = 0.90,   # 90% reliability
  mission_time = 500,    # hours
  conf_level   = 0.90,   # 90% confidence
  n            = 20      # 20 test units
)
print(plan)

## ----solve-n------------------------------------------------------------------
plan2 <- rdt(
  target       = 0.90,
  mission_time = 500,
  conf_level   = 0.90,
  test_time    = 800
)
print(plan2)

## ----weibull-shape------------------------------------------------------------
betas <- c(0.8, 1.0, 1.5)
for (b in betas) {
  p <- rdt(target = 0.90, mission_time = 500, conf_level = 0.90, n = 20, beta = b)
  cat(sprintf("beta = %.1f  ->  required test time = %.1f hours\n",
              b, p$Required_Test_Time))
}

## ----allow-failures-----------------------------------------------------------
for (f in 0:3) {
  p <- rdt(target = 0.90, mission_time = 500, conf_level = 0.90, n = 20, f = f)
  cat(sprintf("f = %d  ->  required test time = %.1f hours\n",
              f, p$Required_Test_Time))
}

