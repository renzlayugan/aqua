#' Data Preparation for `trendplot.gam`
#'
#' @description
#' This is used to prepare the data for the \code{`trendplot.gam`} function.
#' It is a wrapper for various functions in \pkg{tidyverse}, \pkg{mgcv}, and \pkg{gratia}. Please cite those packages.
#'
#' @import tidyverse
#' @import mgcv
#' @import gratia
#'
#' @param object A `mgcv::gam` or `mgcv::gamm` object.
#' @param data Data frame source of the \code{object} argument.
#' @param fitted A logical argument. \code{TRUE} dictates that derivatives are calculated from equally separated, weekly fitted values \code{(DEFAULT)}. \code{FALSE} dictates that derivatives are calculated from predictor values in \code{data}.
#' @param ... Additional arguments passed onto `gratia::derivatives`.
#' However, you cannot change the object, data, select, partial_match, type, and interval arguments as this function have preset options.
#'
#' @examples
#' data1 <- data.frame(
#' DATE = c("1988-06-27", "1988-08-25", "1990-05-03", "1990-05-31", "1990-06-26", "1988-06-27", "1988-08-25", "1990-05-03", "1990-05-31", "1990-06-26"),
#' NUMDATE = c(6752, 6811, 7427, 7455, 7481, 6752, 6811, 7427, 7455, 7481),
#' JULIAN = c(179, 238, 123, 151, 177, 179, 238, 123, 151, 177),
#' GROUP = c('A', 'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'B'),
#' RESPONSE = c(22, 49, 75, 41, 36, 45, 2, 20, 31, 40)
#' )
#'
#' mod1 <- mgcv::gam(RESPONSE ~ GROUP + s(NUMDATE, by = GROUP, k = kyears, bs = 'tp') +
#' s(JULIAN, by = GROUP, k = kmonths, bs = 'cc'),
#' data = data1, method = 'REML', family = 'nb',
#' knots = list(JULIAN = c(1,365)))
#'
#' dataprep.derivatives(mod1, data1, fitted = TRUE)
#'
#' @export

dataprep.derivatives <- function(object, data, fitted = TRUE, ...) {

  # load required packages
  reqpkgs <- c("tidyverse", "mgcv", "gratia")
  for (pkg in reqpkgs) {
    # install if required
    if (!pkg %in% installed.packages()) {install.packages(pkg, dependencies = T)}
    # load when ready
    if (pkg %in% installed.packages()) {suppressMessages(library(pkg, character.only = T))}
  }

  # return error if object is not supplied with a gam / gamm object
  if(!any(class(object) %in% c('gam', 'gamm'))) {stop('Ensure the object argument is supplied with a `mgcv::gam` or `mgcv::gamm` object.')}

  # determine if we based it off fitted values
  if(fitted == TRUE) {
    # get number of years and months, respectively
    kyears <- length(unique(data$YEAR))
    kmonths <- length(unique(data$MONTH))

    # get the derivatives from the model
    numdate <- gratia::derivatives(object = object, data = NULL, select = 's(NUMDATE)', partial_match = T,
                                   type = 'central', interval = 'confidence', n = kyears * 52, ...)
    julian <- gratia::derivatives(object = object, select = 's(JULIAN)', partial_match = T,
                                  type = 'central', data = NULL, n = kmonths * 4, interval = 'confidence', ...)
  } else {
    data0 <- data

    # get the derivatives from the model
    numdate <- gratia::derivatives(object = object, data = data0, select = 's(NUMDATE)', partial_match = T,
                                   type = 'central', interval = 'confidence', ...)
    julian <- gratia::derivatives(object = object, select = 's(JULIAN)', partial_match = T,
                                  type = 'central', data = data0, interval = 'confidence', ...)

    # add observed response values
    numdate$RESPONSE <- data$RESPONSE
  }

  # exclude extrapolated derivatives from 's(NUMDATE')'
  numdate <- numdate[month(as_date(numdate$NUMDATE)) >= min(data$MONTH) & month(as_date(numdate$NUMDATE)) <= max(data$MONTH),]

  # get numerical indicators of significant change
  numdate <- gratia::add_sizer(numdate, type = 'change')
  julian <- gratia::add_sizer(julian, type = 'change')


  # get factorial indicators of significant change
  numdate$.trend <- ifelse(is.na(numdate$.change), 'no change',
                           ifelse(numdate$.change > 0, 'increasing', 'decreasing'))
  numdate$.trend <- as.factor(numdate$.trend)

  julian$.trend <- ifelse(is.na(julian$.change), 'no change',
                           ifelse(julian$.change > 0, 'increasing', 'decreasing'))
  julian$.trend <- as.factor(julian$.trend)

  # replace NAs as zeros for the numerical indicator of significant change.
  numdate$.change[is.na(numdate$.change)] <- 0
  julian$.change[is.na(julian$.change)] <- 0

  # add a date column
  numdate$DATE <- as_date(numdate$NUMDATE)

  class(numdate)[5] <- 'aqua.derivatives'
  class(julian)[5] <- 'aqua.derivatives'

  output <- list(numdate = numdate, julian = julian)
  class(output) <- 'aqua.derivatives'
  return(output)
}
