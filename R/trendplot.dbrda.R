#' Time Series Visualizations for Community Data
#'
#' @description
#' This is used to visualize the results from the \code{`trendtest.dbrda`} function.
#' It is a wrapper for various functions in \pkg{tidyverse} for data management and visualizations. Please cite the those packages.
#' For saving the plot, I recommend using \code{ggsave("path/filename.pdf", width = 9, height = 6, units = "in", dpi = 300)}.
#'
#' @import tidyverse
#' @import ggrepel
#'
#' @param object A \code{trendtest.dbrda} object.
#' @param type One of:
#' \itemize{
#' \item \code{'ts'} = Default db-RDA plot. Fitted explanatory variables are shown as vectors while species are labeled points.
#' }
#'
#' @examples
#' data <- data.frame(DATE = c("1988-06-27", "1988-08-25", "1990-05-03", "1990-05-31", "1990-06-26", "1990-08-21", "1992-07-16", "1992-08-10", "1994-05-25", "1994-08-03"),
#' Sp1 = c(22, 49, 75, 41, 36, 45, 2, 20, 31, 40),
#' Sp2 = c(70, 58, 96, 28, 78, 79, 53, 66, 84, 80))
#' trendplot.dbrda(trendtest.dbrda(data$DATE, data[-1]), type = "ts")
#'
#' @export

trendplot.dbrda <- function(object, type = "ts") {

  # required packages
  reqpkgs <- c("tidyverse", "ggrepel")
  for (pkg in reqpkgs) {
    # install if required
    if (!pkg %in% installed.packages()) {install.packages(pkg, dependencies = T)}

    # load when installed
    if (pkg %in% installed.packages()) {suppressMessages(library(pkg, character.only = T))}
  }

  # load data
  results <- object

  # plots
  if (type == "ts") {
    # set up the plot
    plot1 <- ggplot() +
      # add explanatory vectors
      geom_segment(data = results$summary$biplot, aes(x = 0, xend = dbRDA1, y = 0, yend = dbRDA2),
                   arrow = arrow(length = unit(0.25, "cm")), linewidth = 1, colour = "black") +
      ggrepel::geom_text_repel(data = results$summary$biplot, aes(x = dbRDA1, y = dbRDA2, label = rownames(results$summary$biplot))) +
      # add species as points
      geom_point(data = results$summary$species, aes(dbRDA1, dbRDA2)) +
      ggrepel::geom_text_repel(data = results$summary$species, aes(dbRDA1, dbRDA2, label = rownames(results$summary$species), fontface = "italic"), max.overlaps = 20) +
      # add theme
      labs(x = paste("dbRDA1 (", round(results$summary$cont$importance[2, "dbRDA1"] * 100, digits = 1), "%)", sep = ""),
           y = paste("dbRDA2 (", round(results$summary$cont$importance[2, "dbRDA2"] * 100, digits = 1), "%)", sep = "")) +
      theme(panel.background = element_rect(fill = F),
            panel.border = element_rect(fill = F),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14))

    return(plot1)
  }
}
