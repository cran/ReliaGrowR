% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/duane_plot.R
\name{duane_plot}
\alias{duane_plot}
\title{Duane Analysis.}
\usage{
duane_plot(
  times,
  failures,
  plot = TRUE,
  log = TRUE,
  point_col = "black",
  line_col = "black",
  xlab = "Cumulative Time",
  ylab = "Cumulative MTBF",
  main = "Duane Plot with Cumulative MTBF"
)
}
\arguments{
\item{times}{A vector of cumulative times at which failures occurred.}

\item{failures}{A vector of the number of failures at each corresponding time in times.}

\item{plot}{Show Duane plot (TRUE) or hide plot (FALSE).}

\item{log}{Logical indicating whether to use logarithmic scale for the plot (default: TRUE).}

\item{point_col}{Color for the data points (default: "black").}

\item{line_col}{Color for the fitted line (default: "black").}

\item{xlab}{Label for the x-axis (default: "Cumulative Time").}

\item{ylab}{Label for the y-axis (default: "Cumulative MTBF").}

\item{main}{Title for the plot (default: "Duane Plot with Cumulative MTBF").}
}
\value{
A list containing the fitted model, AIC, and BIC.
}
\description{
This function performs a Duane analysis (1962) \url{doi:10.1109/TA.1964.4319640} on the provided failure data.
}
\examples{
times <- c(100, 200, 300, 400, 500)
failures <- c(1, 2, 1, 3, 2)
fit <- duane_plot(times, failures)
print(fit)
}
