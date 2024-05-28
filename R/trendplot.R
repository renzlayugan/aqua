#' Time Series Visualizations for Water Quality Data
#'
#' @description
#' This is used to visualize the results from the \code{`trendtest`} function.
#' It is a wrapper for various functions in \pkg{tidyverse} for data management and visualizations. Please cite the that package.
#'
#' @import tidyverse
#'
#' @param x Vector of the dates corresponding to \code{y}.
#' @param y Vector of the response variable.
#' @param type One of:
#' \itemize{
#' \item \code{'ts'} = (\code{DEFAULT}) plots the response variable \code{(y)} with time \code{(x)}.
#' \item \code{'kw'} = Kruskal-Wallis test visualized through box plots across months.
#' \item \code{'vbh'} = van Belle & Hughes test visualized through a line graph grouped by month.
#' \item \code{'levene'} = Levene's test visualized through box plots across years.
#' }
#'
#' @examples
#' data <- data.frame(DATE = c("1988-06-27", "1988-08-25", "1990-05-03", "1990-05-31", "1990-06-26", "1990-08-21", "1992-07-16", "1992-08-10", "1994-05-25", "1994-08-03"),
#' RESPONSE = c(22, 49, 75, 41, 36, 45, 2, 20, 31, 40))
#' trendplot(data$DATE, data$RESPONSE)
#'
#' @export

trendplot <- function(x, y, type = "ts") {
  # required packages
  reqpkgs <- c("tidyverse")
  for (pkg in reqpkgs) {
    # install if required
    if (!pkg %in% installed.packages()) {install.packages(pkg, dependencies = T)}

    # load when ready
    if (pkg %in% installed.packages()) {suppressMessages(library(pkg, character.only = T))}
  }

  # load data and add MONTH and YEAR
  data1 <- data.frame(RESPONSE = y,  DATE = as_date(x))
  data1$MONTH = month(data1$DATE)
  data1$YEAR = year(data1$DATE)

  # plot using ggplot2
  if (type == "ts") {
    plot1 <- ggplot(data = data1, aes(DATE, RESPONSE)) +
      geom_line() +
      geom_point() +
      theme(panel.background = element_rect(fill = F),
            panel.border = element_rect(fill = F))
    print(plot1)
  } else if (type == "kw") {
    plot1 <- ggplot(data = data1, aes(as.factor(MONTH), RESPONSE)) +
      geom_boxplot() +
      theme(panel.background = element_rect(fill = F),
            panel.border = element_rect(fill = F))
    print(plot1)
  } else if (type == "vbh") {
    plot1 <- ggplot(data = data1, aes(YEAR, RESPONSE, colour = as.factor(MONTH))) +
      geom_line() +
      geom_point() +
      theme(panel.background = element_rect(fill = F),
            panel.border = element_rect(fill = F),
            legend.position = "top") +
      guides(colour = guide_legend(nrow = 1))
    print(plot1)
  } else if (type == "levene") {
    plot1 <- ggplot(data = data1, aes(YEAR, RESPONSE, group = as.factor(YEAR))) +
      geom_boxplot() +
      theme(panel.background = element_rect(fill = F),
            panel.border = element_rect(fill = F))
    print(plot1)
  }
}
