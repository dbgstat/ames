#' ames: statistical solutions for Ames test data
#'
#' A complete functional tool for Salmonella/microsome (Ames test) data analysis:
#' model fitting and selection; residual and diagnostics analyses; plotting and envelope functions.
#'
#' \itemize{
#'   \item \code{\link{fit.ames}}: fits several different possible modes to Ames test data
#'   \item \code{\link{aic.fitames}}: given multiple models, helps to select the best fit
#'   \item \code{\link{dtest.ames}}: performs overdispersion tests for fitted Poisson models or raw data
#'   \item \code{\link{diagnostics.fitames}}: provides diagnostic measures such as leverage and influence for fitted models
#'   \item \code{\link{plot.fitames}}: offers different type of plots for fitted values and residual analyses
#'   \item \code{\link{mmd.fitames}}: uses one or multiple models to estimate the MMD (maximum mutagenic dose)
#'   \item \code{\link{rames}}: usefull for building confidence intervals, generates random numbers and/or simulates an assay
#' }
#'
#' @author Davi Butturi-Gomes
#'
#' Silvio S. Zocchi
#'
#'
#'
#' @docType package
#' @keywords packages ames
#' @name ames
NULL
