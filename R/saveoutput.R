#' Save Output
#'
#' @description
#' This is used to save the output from a `trendtest()` object.
#'
#' @param object A `trendtest()` object.
#' @param filename Your path and file name in a string form. See example below.
#' @param description A short description for the analyses.
#'
#' @examples
#' data <- data.frame(DATE = c("1988-06-27", "1988-08-25", "1990-05-03", "1990-05-31", "1990-06-26", "1990-08-21", "1992-07-16", "1992-08-10", "1994-05-25", "1994-08-03"),
#' RESPONSE = c(22, 49, 75, 41, 36, 45, 2, 20, 31, 40))
#' saveoutput(trendtest(data$DATE, data$RESPONSE), file = "path/filename.txt", description = "your description")
#'
#' @export

saveoutput <- function(object, filename = NULL, description = NULL) {

  # create the output file
  if (!is.null(filename)){filename = filename}
  else {filename = paste("output", Sys.time(), ".txt", sep = "")}
  sink(file = filename)

  # create, if any, the description for the set of analyses
  if (!is.null(description)) {cat(description, "\n\n")}
  else {cat("ADD DOCUMENT TITLE/DESCRIPTION HERE \n\n")}

  # the output object
  object

  # end and save the created output file
  sink()
}
