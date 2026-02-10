#' Parallel Rarefaction Curve
#'
#' @description Wrapper around calculate_alpha_diversity to generate rarefaction curve data with multiple rarefactions of each sample using parallel processing
#' @param pobject A \code{phyloseq} object.
#' @param ntables Number or rarefactions to perform. Default 10
#' @param step Resolution of sequencing depth (interval between rarefactions). Default 250]
#' @param maxdepth Maximum sequencing depth to include Default 100000 - suitable for Nanopore sequencing
#' @param metrics Which alpha diversity metrics to calculate. Default c("Observed","Chao1","FaithPD","Shannon")
#' @param seedstart Random seed to use. Default 500
#' @param verbose Logical, indicating if information should be reported during calculation. Default FALSE
#' @param nthreads Integer. Number of parallel threads. 0 = auto detect maxcores available (all cores - 1), 1 = no parallelization, >1 = specific number. Default 0
#'
#' @return dataframe with average alpha diversity per sample at each sequencing depth (from 1 to maxdepth)
#'
#' @export
#'
#' @importFrom future plan multisession availableCores
#' @importFrom future.apply future_lapply
#'
#' @examples
#' \dontrun{
#' # Auto detect max cores available (uses all -1)
#' rdat <- Rcurve_data_paral(phy)
#'
#' # Use a specific number of cores
#' rdat <- Rcurve_data_paral(phy, nthreads = 4)
#'
#' Do not use parallelisation - could also just run Rcurve_data()
#' rdat <- Rcurve_data_paral(phy, nthreads = 1)
#'  }



Rcurve_data_paral <- function(pobject, ntables = 10, step = 250, maxdepth = 100000,
                              metrics = c("Observed", "Chao1", "FaithPD", "Shannon"),
                              seedstart = 500, verbose = FALSE, nthreads = 0) {


  require("vegan")
  loadNamespace("phyloseq")

  # Input check

  if (!inherits(pobject, "phyloseq")){
    stop("pobject must be a phyloseq object")
  }

  if (nthreads < 0){
    stop()
  }

  # Determine number of cores
  if(nthreads == 0) {
    ncores <- max(1, future::availableCores() - 1)
  } else if(nthreads == 1) {
    ncores <- 1
  } else {
    ncores <- nthreads
  }

  # Report on number of cores being used
  if (verbose){
    message("Using ", ncores, " cores for parallel processing")
  }

  # Set up parallel plan if using more than 1 core
  if (ncores > 1){
    oplan <- future::plan(future::multisession, workers  = ncores)
    on.exit(future::plan(oplan), add = TRUE) # reset plan on exit
  }

  # Create sequence of depths
  step.seq <- seq(from = 10, to = maxdepth, by = step)

  if (verbose){
    message("Using ", length(step.seq), " depth values from 10 to ", maxdepth)
  }

  # Run rarefaction (parallel or sequential based on ncores)
  if (ncores > 1) {
    rare_tab <- future.apply::future_lapply(
      step.seq,
      function(k) {
        library(phyloseq)
        library(vegan)
        calculate_alpha_diversity(
          pobject = pobject,
          ntables = ntables,
          depth = k,
          INDECES = metrics,
          seedstart = seedstart,
          verbose = FALSE  # Suppress verbose in parallel to avoid clutter
        )
      },
      future.seed = TRUE
    )
  } else {
    # Sequential processing
    rare_tab <- lapply(
      step.seq,
      function(k) {
        calculate_alpha_diversity(
          pobject = pobject,
          ntables = ntables,
          depth = k,
          INDECES = metrics,
          seedstart = seedstart,
          verbose = verbose)
      }
    )
  }

  # Combine results
  rare_tab <- do.call(rbind, rare_tab)

  if (verbose) {
    message("Rarefaction curve data generated successfully")
  }
  return(rare_tab)
}
