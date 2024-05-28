#' Time Series Analyses for Water Quality Data
#'
#' @description
#' The set of analyses in this are adapted from Alberta Lake Management Society (ALMS) 2018: Guide to trend analysis on Alberta lakes. \cr
#' (\url{https://alms.ca/wp-content/uploads/2018/04/ALMS-Trend-Analysis-Methods-04-2018.pdf}) \cr\cr
#'
#' This is primarily a wrapper for the functions in the following packages, including various functions in \pkg{tidyverse} for data management. Please cite the listed packages.
#' \itemize{
#' \item \pkg{FSA}: \code{`dunnTest`}
#' \item \pkg{EnvStats}: \code{`kendallSeasonalTrendTest`} and \code{`kendallTrendTest`}
#' \item \pkg{car}: \code{`leveneTest`}
#' }
#'
#' @import tidyverse
#' @import FSA
#' @import EnvStats
#' @import car
#'
#' @param x Vector of the dates corresponding to \code{y}.
#' @param y Vector of the response variable.
#' @param summary If \code{TRUE (DEFAULT)}, it returns the summary output for the results.
#' If \code{FALSE}, it returns the full native results in a \code{list} object
#' @param median If \code{TRUE}, the Mann-Kendall test is performed on annual median.
#' If \code{FALSE (DEFAULT)}, it is performed on all values while correcting for ties in the estimation.
#'
#' @examples
#' data <- data.frame(DATE = c("1988-06-27", "1988-08-25", "1990-05-03", "1990-05-31", "1990-06-26", "1990-08-21", "1992-07-16", "1992-08-10", "1994-05-25", "1994-08-03"),
#' RESPONSE = c(22, 49, 75, 41, 36, 45, 2, 20, 31, 40))
#' trendtest(data$DATE, data$RESPONSE)
#'
#' @export

trendtest <- function(x, y, summary = TRUE, median = FALSE) {

  # required packages
  reqpkgs <- c("tidyverse", "FSA", "EnvStats", "car")
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

  # Kruskal-Wallis — a test for seasonality or differences among months
  kw <- kruskal.test(RESPONSE ~ as.factor(MONTH), data = data1)

  # Dunn's multiple pairwise test (post-hoc for Kruskal-Wallis) - which months are different
  if(kw$p.value < 0.05){dunn <- FSA::dunnTest(RESPONSE ~ as.factor(MONTH), data = data1, method = "bonferroni")}

  # Seasonal Kendall - test for monotonic trends, while considering seasonal variation
  sk <- EnvStats::kendallSeasonalTrendTest(RESPONSE ~ as.factor(MONTH) + YEAR, data = data1)

  # Mann-Kendall - test for monotonic trend
  if (median == TRUE) {
    data1.median <- data1 %>% group_by(YEAR) %>% summarize(RESPONSE = median(RESPONSE))
    mk <- EnvStats::kendallTrendTest(RESPONSE ~ YEAR, data = data1.median)
  } else {mk <- EnvStats::kendallTrendTest(RESPONSE ~ YEAR, data = data1)}

  # Levene's test for heterogeneity of variance (test for seasonal variability)
  levene.y <- car::leveneTest(RESPONSE ~ as.factor(YEAR), data = data1)

  # output
  if (summary == T) {
    cat("================================================================== \n")
    cat("Time series analysis for water quality \n")
    cat("(adapted from ALMS 2018: Guide to trend analysis on Alberta lakes) \n\n")

    cat("---------------------------------------------------- \n")
    cat("Kruskal-Wallis rank sum test (test for seasonality): \n\n")
    kwres <- matrix(nrow = 1, ncol = 3)
    rownames(kwres) <- c(" ")
    colnames(kwres) <- c("Chi-squared", "df", "p-value")
    kwres[1, "Chi-squared"] <- kw$statistic
    kwres[1, "df"] <- kw$parameter
    kwres[1, "p-value"] <- kw$p.value
    print(kwres)
    cat("\n")

    if (kw$p.value < 0.05) {
      cat("Dunn's multiple pairwise test (Bonferroni-adjusted, α < 0.10): \n\n")
      print(dunn$res[dunn$res$P.adj < 0.10,])
      cat("\n")
    }

    cat("---------------------------------------------------- \n")
    cat("van Belle Hughes homogeneity of trends test: \n\n")
    vbhres <- matrix(nrow = 1, ncol = 3)
    rownames(vbhres) <- c(" ")
    colnames(vbhres) <- c("Chi-squared", "df", "p-value")
    vbhres[1, "Chi-squared"] <- sk$statistic[1]
    vbhres[1, "df"] <- sk$parameters
    vbhres[1, "p-value"] <- sk$p.value[1]
    print(vbhres)
    cat("\n")

    cat("---------------------------------------------------- \n")
    cat("Seasonal Kendall test: \n\n")
    cat(" Estimates: \n")
    sk1res <- matrix(nrow = 1, ncol = 4)
    rownames(sk1res) <- c(" ")
    colnames(sk1res) <- c("intercept", "slope", "lower", "upper")
    sk1res[1, "intercept"] <- sk$estimate["intercept"]
    sk1res[1, "slope"] <- sk$estimate["slope"]
    sk1res[1, "lower"] <- sk$interval$limits[1]
    sk1res[1, "upper"] <- sk$interval$limits[2]
    print(sk1res)
    cat("\n")

    cat(" Hypothesis testing: \n")
    sk2res <- matrix(nrow = 1, ncol = 4)
    rownames(sk2res) <- c(" ")
    colnames(sk2res) <- c("tau", "z", "n", "p-value")
    sk2res[1, "tau"] <- sk$estimate["tau"]
    sk2res[1, "z"] <- sk$statistic[2]
    sk2res[1, "n"] <- sk$sample.size["Total"]
    sk2res[1, "p-value"] <- sk$p.value[2]
    print(sk2res)
    cat("\n")

    cat("---------------------------------------------------- \n")
    if (median == TRUE) {cat("Mann-Kendall test (annual median): \n\n")}
    else {cat("Mann-Kendall test: \n\n")}
    cat(" Estimates: \n")
    mk1res <- matrix(nrow = 1, ncol = 4)
    rownames(mk1res) <- c(" ")
    colnames(mk1res) <- c("intercept", "slope", "lower", "upper")
    mk1res[1, "intercept"] <- mk$estimate["intercept"]
    mk1res[1, "slope"] <- mk$estimate["slope"]
    mk1res[1, "lower"] <- mk$interval$limits[1]
    mk1res[1, "upper"] <- mk$interval$limits[2]
    print(mk1res)
    cat("\n")

    cat(" Hypothesis testing: \n")
    mk2res <- matrix(nrow = 1, ncol = 4)
    rownames(mk2res) <- c(" ")
    colnames(mk2res) <- c("tau", "z", "n", "p-value")
    mk2res[1, "tau"] <- mk$estimate["tau"]
    mk2res[1, "z"] <- mk$statistic
    mk2res[1, "n"] <- mk$sample.size
    mk2res[1, "p-value"] <- mk$p.value
    print(mk2res)
    cat("\n")

    cat("---------------------------------------------------- \n")
    cat("Levene's test for heterogeneity of variance (by YEAR): \n\n")
    leveneres <- matrix(nrow = 1, ncol = 4)
    rownames(leveneres) <- c(" ")
    colnames(leveneres) <- c("F", "df1", "df2", "p-value")
    leveneres[1, "F"] <- levene.y$`F value`[1]
    leveneres[1, "df1"] <- levene.y$Df[1]
    leveneres[1, "df2"] <- levene.y$Df[2]
    leveneres[1, "p-value"] <- levene.y$`Pr(>F)`[1]
    print(leveneres)
    cat("\n")

    cat("================================================================== \n")

  } else {return(list(KW = kw, DUNN = dunn, SK = sk, MK = mk, LEVENE.Y = levene.y))}
}
