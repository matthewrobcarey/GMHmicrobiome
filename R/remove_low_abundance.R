#' Remove low abundance taxa
#'
#' @description
#' Wrapper around phyloseq filter_taxa to allow users to remove (prune) low abundance taxa
#'
#' @param pobject Phyloseq object to remove low abundance taxa
#' @param relative remove taxa based on relative abundance (TRUE) or raw reads (FALSE). Default = FALSE remove based on raw reads.
#' @param read_count value to filter taxa by before removal when using reads. Default = 100
#' @param rel_abun value to filter taxa by before removal when using relative abundance. Default = 0.0001 (0.01 percent).
#' @param min_prev prevalence (the proportion of samples) that must contain the taxa at the minimum level to be kept. Default = 0.01 (1 percent).
#' @param rel_fun equation used to work out relative abundance, proportion used rather than percent. Default = function(x) x/sum(x).
#' @param verbose logical, prints number of taxa removed. Default = TRUE
#'
#' @return phyloseq object with taxa removed based on rules
#'
#' @export
#'
#' @import phyloseq
#'
#' @examples
#' \dontrun{
#' # Remove taxa where 1% of the samples does not have more than 100 reads
#' remove_low_abundance(phy)
#'
#' # Remove taxa where 10% of the samples do not contain 0.25 percent relative abundance
#' remove_low_abundance(phy, relative = TRUE, rel_abun = 0.0025, min_prev = 0.1)
#'
#' # Remove taxa where 50% of the samples do not contain 0.001 percent relative abundance
#' remove_low_abundance(phy, relative = TRUE, rel_abun = 0.0001, min_prev = 0.5)
#' }
#'


remove_low_abundance <- function(pobject, relative  = FALSE, read_count = 100, rel_abun = 0.0001, min_prev = 0.01,
                                 rel_fun = function(x) x/sum(x), verbose = TRUE){

  loadNamespace("phyloseq")

  # Check input

  if (!inherits(pobject, "phyloseq")){
    stop("pobject must be a phyloseq object")
  }

  ntaxa_before <- phyloseq::ntaxa(pobject)


  if (relative){
    phy_prop <- phyloseq::transform_sample_counts(physeq = pobject, fun = rel_fun)
    phy_filtered <- phyloseq::filter_taxa(phy_prop, function(x){mean(x > rel_abun) >= min_prev}, prune = TRUE)

  } else {
    phy_filtered <- phyloseq::filter_taxa(pobject, function(x){mean(x > read_count) >= min_prev}, prune = TRUE)
  }

  ntaxa_after <- phyloseq::ntaxa(phy_filtered)

  if (verbose){
    message("Removed ", ntaxa_before - ntaxa_after, " taxa ( ",
            round((ntaxa_before - ntaxa_after)/ntaxa_before * 100, 1),
            "%). ", ntaxa_after, " taxa were retained.")
  }

  return(phy_filtered)

}









