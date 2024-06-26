% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trendplot.R
\name{trendplot}
\alias{trendplot}
\title{Time Series Visualizations for Water Quality Data}
\usage{
trendplot(x, y, type = "ts")
}
\arguments{
\item{x}{Vector of the dates corresponding to \code{y}.}

\item{y}{Vector of the response variable.}

\item{type}{One of:
\itemize{
\item \code{'ts'} = (\code{DEFAULT}) plots the response variable \code{(y)} with time \code{(x)}.
\item \code{'kw'} = Kruskal-Wallis test visualized through box plots across months.
\item \code{'vbh'} = van Belle & Hughes test visualized through a line graph grouped by month.
\item \code{'levene'} = Levene's test visualized through box plots across years.
}}
}
\description{
This is used to visualize the results from the \code{`trendtest`} function.
It is a wrapper for various functions in \pkg{tidyverse} for data management and visualizations. Please cite the that package.
}
\examples{
data <- data.frame(DATE = c("1988-06-27", "1988-08-25", "1990-05-03", "1990-05-31", "1990-06-26", "1990-08-21", "1992-07-16", "1992-08-10", "1994-05-25", "1994-08-03"),
RESPONSE = c(22, 49, 75, 41, 36, 45, 2, 20, 31, 40))
trendplot(data$DATE, data$RESPONSE)

}
