#' Pairwise Synchrony Analysis
#'
#' @description
#' This function tests for lake synchrony or can be used for multiple variables using Kendall tau-b's coefficient. \cr\cr
#'
#' This is primarily a wrapper for the \code{`Kendall`} function in \pkg{Kendall} and various functions in \pkg{tidyverse} for data management. Please cite the listed packages.
#'
#' @import tidyverse
#' @import Kendall
#'
#' @param data A wide-form \code{data.frame} object where columns are variables and rows are samples.
#'
#' @examples
#' data <- data.frame(DATE = c("1988-06-27", "1988-08-25", "1990-05-03", "1990-05-31", "1990-06-26", "1990-08-21", "1992-07-16", "1992-08-10", "1994-05-25", "1994-08-03"),
#' var1 = c(22, 49, 75, 41, 36, 45, 2, 20, 31, 40),
#' var2 = c(70, 58, 96, 28, 78, 79, 53, 66, 84, 80))
#' pairwise.kendall(data[-1])
#'
#' @export

pairwise.kendall <- function(data) {

  # required packages
  reqpkgs <- c("tidyverse", "Kendall")
  for (pkg in reqpkgs) {
    # install if required
    if (!pkg %in% installed.packages()) {install.packages(pkg, dependencies = T)}

    # load when ready
    if (pkg %in% installed.packages()) {suppressMessages(library(pkg, character.only = T))}
  }

  # load data
  data1 <- data

  ## get new data comparing var a to var b
  # get a matrix of variable names
  data1.mtrx <- matrix(nrow = length(names(data1)), ncol = length(names(data1)))
  rownames(data1.mtrx) <- names(data1)
  colnames(data1.mtrx) <- names(data1)
  # only get data for the upper triangle = non-repeating pairwise combinations
  data1.mtrx[upper.tri(data1.mtrx)] <- 1
  # revert matrix to data frame - make it easy for analysis
  data1.df <- as.data.frame(data1.mtrx)
  # add var a from rownames to a column
  data1.df$VAR1 <- rownames(data1.df)
  # convert data frame into long form (non-repeating pairwise combinations)
  data1.df <- gather(data1.df, key = "VAR2", value = "YES", -VAR1, na.rm = T)

  ## pairwise Kendall tau-b correlations
  data1.results <- data1.df %>% select(-YES)
  data1.results$n <- NA
  data1.results$tau <- NA
  data1.results$S <- NA
  data1.results$D <- NA
  data1.results$varS <- NA
  data1.results$p <- NA

  for (i in rownames(data1.results)) {
    x <- data1.results$VAR1[rownames(data1.results) == i]
    y <- data1.results$VAR2[rownames(data1.results) == i]

    n <- sum(!is.na(data1[data1.results$VAR1[rownames(data1.results) == i]]) & !is.na(data1[data1.results$VAR2[rownames(data1.results) == i]]))
    out1 <- Kendall(as.matrix(data1[,x]), as.matrix(data1[,y]))

    data1.results["n"][rownames(data1.results) == i,] <- n
    data1.results["tau"][rownames(data1.results) == i,] <- out1$tau
    data1.results["S"][rownames(data1.results) == i,] <- out1$S
    data1.results["D"][rownames(data1.results) == i,] <- out1$D
    data1.results["varS"][rownames(data1.results) == i,] <- out1$varS
    data1.results["p"][rownames(data1.results) == i,] <- out1$sl
  }

  # calculate bonferroni-corrected p-values: p.adj = p * [k(k-1) / 2], where k = the number of lakes
  k <- length(names(data1))
  data1.results$p.adj <- data1.results$p * ((k * (k - 1)) / 2)
  data1.results$p.adj[data1.results$p.adj > 1] <- 1

  # remove values for pairs with no observations
  data1.results$tau[data1.results$n == 0] <- NA
  data1.results$S[data1.results$n == 0] <- NA
  data1.results$D[data1.results$n == 0] <- NA
  data1.results$varS[data1.results$n == 0] <- NA
  data1.results$p[data1.results$n == 0] <- NA
  data1.results$p.adj[data1.results$n == 0] <- NA

  # output
  return(data1.results)
}
