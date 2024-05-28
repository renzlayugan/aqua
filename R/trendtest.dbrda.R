#' Time Series Analyses for Community Data
#'
#' @description
#' This function tests for the interannual and seasonal shifts in a given community composition.
#' It is based on Bray-Curtis dissimilarity distance-based redundancy analysis. \cr\cr
#'
#' This is primarily a wrapper for the following functions in \pkg{vegan} and various functions in \pkg{tidyverse} for data management. Please cite the listed packages.
#' \itemize{
#' \item \pkg{vegan}:
#' \code{`decostand`} for Hellinger standardization,
#' \code{`vegdist`} for calculating Bray-Curtis dissimilarity matrix,
#' \code{`adonis2`} for perMANOVA,
#' \code{`dbrda`} for distanced-based RDA,
#' \code{`sppscores`} for generating species scores,
#' \code{`score`} for extracting scores,
#' \code{`ordiR2step`} for forward selection,
#' \code{`varpart`} for variance partitioning.
#' }
#'
#' @import tidyverse
#' @import vegan
#'
#' @param x A vector of dates or a data frame of the explanatory variables.
#' If \code{x} is a vector of dates (determined by `is.Date` function of \pkg{lubridate}),
#' then this function automatically generates a data frame of the DATE, YEAR, MONTH, and YDAY (Julian day).
#' @param y Data frame of response variables, where rows represent sites and columns represent species.
#' @param summary \code{TRUE} returns a summarized version of the db-RDA results. \code{FALSE} returns the full uncondensed results in a \code{list} object.
#' @param nperm Number of permutations. \code{999} is default.
#' @param seed \code{FALSE} is default. Otherwise, it can be any \code{numeric} to set up the seed using `set.seed(seed)`, which ensure the results from permutation steps are replicable.
#'
#' @examples
#' data <- data.frame(DATE = c("1988-06-27", "1988-08-25", "1990-05-03", "1990-05-31", "1990-06-26", "1990-08-21", "1992-07-16", "1992-08-10", "1994-05-25", "1994-08-03"),
#' Sp1 = c(22, 49, 75, 41, 36, 45, 2, 20, 31, 40),
#' Sp2 = c(70, 58, 96, 28, 78, 79, 53, 66, 84, 80))
#' trendtest.dbrda(data$DATE, data[-1])
#'
#' @export

trendtest.dbrda <- function(x, y, summary = TRUE, nperm = 999, seed = FALSE) {
  # required packages
  reqpkgs <- c("tidyverse", "vegan")
  for (pkg in reqpkgs) {
    # install if required
    if (!pkg %in% installed.packages()) {install.packages(pkg, dependencies = T)}

    # load when installed
    if (pkg %in% installed.packages()) {suppressMessages(library(pkg, character.only = T))}
  }

  # set seed - to ensure permutations are consistent / replicable
  if (seed != F) {set.seed(seed)}

  # load data
  comm <- y
  comm.hell <- decostand(comm, method = "hellinger")
  if (is.Date(x)) {
    expl <- data.frame(DATE = x, YEAR = year(x), MONTH = month(x), YDAY = yday(x))
  } else {expl <- x}

  ## perform dbRDA with a Bray-Curtis dissimilarity matrix
  mod1 <- dbrda(comm.hell ~ YEAR + YDAY, expl, distance = "bray")
  # add species scores to dbRDA - superimposed
  sppscores(mod1) <- comm.hell
  # extract all scores from the dbRDA
  scores <- scores(mod1, display = "all")

  ## conduct type III Sum Sq PERMANOVA
  #* identical results to `anova.cca`
  permanova <- adonis2(comm.hell ~ YEAR + YDAY, expl, permutations = nperm, method = "bray", by = "margin")

  ## perform forward selection of dbRDA
  # set-up intercept-only model
  mod0 <- dbrda(comm.hell ~ 1, expl, distance = "bray")
  # get full model's R-squared
  adjR2 <- RsquareAdj(mod1)$adj.r.squared
  # run forward selection
  fwdsel <- ordiR2step(mod0, scope = formula(mod1), R2scope = adjR2, direction = 'forward', permutations = nperm, trace = F)
  # adjust p-values according Bonferroni's method
  fwdsel2 <- fwdsel
  fwdsel2$anova$`Pr(>F)` <- p.adjust (fwdsel2$anova$`Pr(>F)`, method = 'bonferroni', n = ncol (expl))

  ## perform variance partioning analysis
  #* Individual explanatory power calculations (R2)
  #* Variance explained by A = Combined R2 of A&B — Individual R2 of B
  #* Variance explained by B = Combined R2 of A&B — Individual R2 of A
  #* Variance explained by A&B = 1 — Residual R2 - Variance explained by A - Variance explained by B
  varpartition <- varpart(vegdist(comm.hell, method = "bray"), ~YEAR, ~YDAY, data = expl)

  # add TITAN2:titan REMEMBER

  if (summary == T) {
    cat("================================================================== \n")
    cat("Time series analysis for a multivariate community \n")
    cat("Testing for interannual and seasonal trends of a community composition \n\n")

    cat("---------------------------------------------------- \n")
    cat("db-RDA Summary \n\n")
    cat("Variance Partitioned: \n")
    inertia <- matrix(nrow = 3, ncol = 2)
    rownames(inertia) <- c("Total", "Constrained", "Unconstrained")
    colnames(inertia) <- c("Inertia", "Proportion")
    inertia["Total", "Inertia"] <- summary(mod1)$tot.chi
    inertia["Constrained", "Inertia"] <- summary(mod1)$constr.chi
    inertia["Unconstrained", "Inertia"] <- summary(mod1)$unconst.chi
    inertia["Total", "Proportion"] <- 1
    inertia["Constrained", "Proportion"] <- summary(mod1)$constr.chi / summary(mod1)$tot.chi
    inertia["Unconstrained", "Proportion"] <- summary(mod1)$unconst.chi / summary(mod1)$tot.chi
    print(inertia)
    cat("\n")

    cat("Importance of components: \n")
    compimportance <- matrix(nrow = 3, ncol = 3)
    rownames(compimportance) <- c("Eigenvalue", "Propotion (total)", "Propotion (constr.)")
    colnames(compimportance) <- c("dbRDA1", "dbRDA2", "Sum")
    compimportance["Eigenvalue", "dbRDA1"] <- summary(mod1)$cont$importance["Eigenvalue", "dbRDA1"]
    compimportance["Eigenvalue", "dbRDA2"] <- summary(mod1)$cont$importance["Eigenvalue", "dbRDA2"]
    compimportance["Eigenvalue", "Sum"] <- summary(mod1)$cont$importance["Eigenvalue", "dbRDA1"] + summary(mod1)$cont$importance["Eigenvalue", "dbRDA2"]
    compimportance["Propotion (total)", "dbRDA1"] <- summary(mod1)$cont$importance["Proportion Explained", "dbRDA1"]
    compimportance["Propotion (total)", "dbRDA2"] <- summary(mod1)$cont$importance["Proportion Explained", "dbRDA2"]
    compimportance["Propotion (total)", "Sum"] <- summary(mod1)$cont$importance["Cumulative Proportion", "dbRDA2"]
    compimportance["Propotion (constr.)", "dbRDA1"] <- summary(mod1)$concont$importance["Proportion Explained", "dbRDA1"]
    compimportance["Propotion (constr.)", "dbRDA2"] <- summary(mod1)$concont$importance["Proportion Explained", "dbRDA2"]
    compimportance["Propotion (constr.)", "Sum"] <- summary(mod1)$concont$importance["Cumulative Proportion", "dbRDA2"]
    print(compimportance)
    cat("\n")

    cat("Constraint loadings: \n")
    print(summary(mod1)$biplot[,c(1:2)])
    cat("\n")

    cat("---------------------------------------------------- \n")
    cat("perMANOVA Table (Type III SS,", nperm, "permutations):", "\n")
    print(as.matrix(permanova))
    cat("\n")

    cat("---------------------------------------------------- \n")
    cat("Forward selection of constrained variables:", "\n")
    print(as.matrix(fwdsel2$anova))
    cat("\n")

    cat("---------------------------------------------------- \n")
    cat("Variance partitioning analysis:", "\n")
    varpart.summary <- varpartition$part$indfract[,-c(2,4)]
    varpart.summary[1,"Comment"] <- "unique effect of YEAR"
    varpart.summary[2,"Comment"] <- "unique effect of YDAY"
    varpart.summary[3,"Comment"] <- "combined effect"
    varpart.summary[4,"Comment"] <- "unaccounted effects"
    print(varpart.summary)
    cat("\n")

    cat("================================================================== \n")

  } else {
    results1 <- list(summary = summary(mod1), permanova = permanova, forward.selection = fwdsel2, variance.partitioning = varpartition)
    return(results1)
  }
}
