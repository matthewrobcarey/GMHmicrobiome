#' Collapse closely related taxa
#'
#' @description
#' A wrapper function that uses tip_glom to group closely related species and/or genus.
#'
#' @param pobject Phyloseq object to cluster closely realted species.
#' @param height Numeric scalar of the height where the tree should be cut. Default = 0.05
#' @param h_method Agglomeration method used, taken from hclust (from stats package). Options "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid". Default = "average"
#' @param rename Logical, rename the tips the name of the most abundant taxa in the cluster (TRUE), or leave the tip names as the most common taxonomic level (FALSE). Default = TRUE
#' @param export_csv Logical, export the csv of species clusters. Default = TRUE
#' @param export_path Character, path to export cluster.csv file. Default = "output/cluster.csv"
#' @param verbose logical, reports change in the number of taxa as a results of clustering. Default = TRUE
#'
#'
#' @return A phyloseq object with closely related species clustered
#'
#' @export
#'
#' @import phyloseq
#' @importFrom tibble enframe
#' @importFrom tidyr unnest_longer
#' @importFrom dplyr mutate
#' @importFrom utils write.csv
#'
#'
#' @examples
#' \dontrun{
#' # Cluster species with a height of 0.2 and don't rename
#' collapse_taxa(phy, height = 0.2, rename = FALSE)
#'
#' # Cluster species with a height of 0.01, rename taxa and use the ward.D2 method
#' collapse_taxa(phy, height = 0.01, method = "ward.D2" )
#' }
#'


collapse_taxa <- function(pobject, height = 0.05, h_method = "average", rename = TRUE, export_csv = TRUE, export_path = "output/clusters.csv", verbose = TRUE){


  # Check inputs
  if (!inherits(pobject, "phyloseq")){
    stop("pobject must be a phyloseq object")
  }

  if (is.null(phy_tree(pobject, FALSE))){
    stop("pobject must contain phylogenetic tree")
  }

  # Count taxa before
  ntaxa_before <- phyloseq::ntaxa(pobject)

  if(!rename){
    glom_phy <- phyloseq::tip_glom(pobject, h = height, hcfun = hclust, method = h_method)
  } else {

    # Extract original taxonomy for renaming
    taxa_org <- as.data.frame(phyloseq::tax_table(pobject))
    # Carry out tip_glom
    glom_phy <- phyloseq::tip_glom(pobject, h = height, hcfun = hclust, method = h_method)

    # Filter original taxonomy to match clusters
    taxa_updated <- taxa_org[phyloseq::taxa_names(glom_phy), ]

    # Change NA species to genus + "spp"
    taxa_updated$Species <- ifelse(
      is.na(taxa_updated$Species) | taxa_updated$Species == "",
      paste0(taxa_updated$Genus, "_spp"),
      taxa_updated$Species
    )

    # Update glom_phy object
    phyloseq::tax_table(glom_phy) <- as.matrix(taxa_updated)
  }

  ntaxa_after <- phyloseq::ntaxa(glom_phy)

  if (verbose){
    if ( ntaxa_after == ntaxa_before){
      message("No taxa were collapsed at height ", height, ". Consider increasing height parameter.")
    }else{
      message("Collapsed ", ntaxa_before - ntaxa_after, " taxa (",
              round((ntaxa_before - ntaxa_after) / ntaxa_before * 100, 1),
              "%). Retained ", ntaxa_after, " taxa.")
    }


  }

  if (export_csv){

    # Create output directory if it doesn't exist
    export_dir <- dirname(export_path)
    if (!dir.exists(export_dir)){
      dir.create(export_dir, recursive = TRUE)
    }

    # Compute distances a tips using h_clust
    # Note: This re-does clustering for export - slight redundancy with tip_glom()
    # but was unable to access clusters from tip_glom()
    tree <- phyloseq::phy_tree(pobject)
    D <- stats::cophenetic(tree)
    # Cluster tips and cut at h_value
    hc  <- stats::hclust(as.dist(D), method = h_method)
    grp <- stats::cutree(hc, h = height)
    # Turn group assignments into clusters of tip labels
    clusters <- split(names(grp), grp)

    # Get representative names
    rep_names <- phyloseq::taxa_names(glom_phy)
    # Extract original taxonomy for renaming
    taxa_org <- as.data.frame(phyloseq::tax_table(pobject))


    # Build a named character vector: names = OTU IDs, values = Species
    species_lookup <- stats::setNames(
      as.vector(taxa_org[, "Species"]),
      taxa_names(pobject)
    )

    # Change NA species to genus plus "_spp"
    species_lookup <- ifelse(
      is.na(species_lookup) | species_lookup == "",
      paste0(taxa_org[names(species_lookup), "Genus"], "_spp"),
      species_lookup
    )

    # Create dataframe of species combinations
    df_clusters <- tibble::enframe(clusters, name = "cluster_id", value = "otus") %>%
      tidyr::unnest_longer(otus, values_to = "OTU") %>%
      dplyr::mutate(Species = species_lookup[OTU],
                    is_reprentative = OTU %in% rep_names,
                    cluster_size = lengths(clusters)[cluster_id]
                    )


    utils::write.csv(df_clusters, file = export_path, row.names = FALSE)

    if (verbose){
      message("Cluster mapping saved to ", export_path)
    }
  }

  return(glom_phy)

}
