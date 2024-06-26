% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trendplot.dbrda.R
\name{trendplot.dbrda}
\alias{trendplot.dbrda}
\title{Time Series Visualizations for Community Data}
\usage{
trendplot.dbrda(object, type = "ts")
}
\arguments{
\item{object}{A \code{trendtest.dbrda} object.}

\item{type}{One of:
\itemize{
\item \code{'ts'} = Default db-RDA plot. Fitted explanatory variables are shown as vectors while species are labeled points.
}}
}
\description{
This is used to visualize the results from the \code{`trendtest.dbrda`} function.
It is a wrapper for various functions in \pkg{tidyverse} for data management and visualizations. Please cite the those packages.
For saving the plot, I recommend using \code{ggsave("path/filename.pdf", width = 9, height = 6, units = "in", dpi = 300)}.
}
\examples{
data <- data.frame(DATE = c("1988-06-27", "1988-08-25", "1990-05-03", "1990-05-31", "1990-06-26", "1990-08-21", "1992-07-16", "1992-08-10", "1994-05-25", "1994-08-03"),
Sp1 = c(22, 49, 75, 41, 36, 45, 2, 20, 31, 40),
Sp2 = c(70, 58, 96, 28, 78, 79, 53, 66, 84, 80))
trendplot.dbrda(trendtest.dbrda(data$DATE, data[-1]), type = "ts")

}
