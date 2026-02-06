#' Initiate GMH project
#'
#' @description Creates the folders and/or import the base Rmarkdown templates needed to run basic microbiome analysis using the \code{GMHmicrobiome} pipeline.
#' @param folders Logical. indicator if folder should be created. Default = TRUE
#' @param files Logical. indicator if Rmarkdown templates should be created. Default = TRUE
#' @param overwrite Logical. indicator if files should be overwritten if present. Default = FALSE
#' @param platform Character. Sequencing platform used: "nanopore" or "iontorrent". Default = "nanopore"
#' @return Create system files and folders
#' @export
#'
#' @examples
#' \dontrun{
#' # Create project with nanopore templates (default)
#' initiate_GMH_project()
#'
#' # Create project with iontorrent templates
#'  initiate_GMH_project(platform = "iontorrent")
#'
#' # Only create the folders, no template files
#' initiate_GMH_project(files = FALSE)
#'
#' Overwrite existing template files - *WARNING this will delete any existing scripts* -
#' initiate_GMH_project(overwrite = TRUE)
#' }

initiate_GMH_project <- function(folders = TRUE, files = TRUE, overwrite = FALSE, platform = "nanopore"){

  platform <- match.arg(platform, choices = c("nanopore", "iontorrent"))

  # Define template mappings based on platform selection
  import_template <- paste0("gmh_import_", platform)
  import_file <- paste0("GMH_1_import_", platform, ".Rmd")

  # Create used folders if missing
  if (isTRUE(folders)) {
    if (!file.exists("R_objects")) dir.create(file.path(getwd(), "R_objects"))
    if (!file.exists("plots")) dir.create(file.path(getwd(), "plots"))
    if (!file.exists("tables")) dir.create(file.path(getwd(), "tables"))
    if (!file.exists("output")) dir.create(file.path(getwd(), "output"))
    if (!file.exists("input")) dir.create(file.path(getwd(), "input"))
  }

  # Create the initial files if not found
  if (isTRUE(files) & isFALSE(overwrite)) {
    if (!file.exists(import_file)) file.copy(system.file("rmarkdown", "templates", import_template, "skeleton", "skeleton.Rmd", package = "GMHmicrobiome"), import_file)
    if (!file.exists("GMH_2_metadata_clean.Rmd")) file.copy(system.file("rmarkdown", "templates", "gmh_metadata_import", "skeleton", "skeleton.Rmd", package = "GMHmicrobiome"), "GMH_2_metadata_clean.Rmd")
    if (!file.exists("GMH_3_description.Rmd")) file.copy(system.file("rmarkdown", "templates", "gmh_description", "skeleton", "skeleton.Rmd", package = "GMHmicrobiome"), "GMH_3_description.Rmd")
    if (!file.exists("GMH_4_test_variables.Rmd")) file.copy(system.file("rmarkdown", "templates", "gmh_test_variables", "skeleton", "skeleton.Rmd", package = "GMHmicrobiome"), "GMH_4_test_variables.Rmd")
    if (!file.exists("GMH_4x_test_variables_code.Rmd")) file.copy(system.file("rmarkdown", "templates", "gmh_test_variables_code", "skeleton", "skeleton.Rmd", package = "GMHmicrobiome"), "GMH_4x_test_variables_code.Rmd")
    if (!file.exists("GMH_5_betadiversity.Rmd")) file.copy(system.file("rmarkdown", "templates", "gmh_beta_diversity", "skeleton", "skeleton.Rmd", package = "GMHmicrobiome"), "GMH_5_betadiversity.Rmd")
    if (!file.exists("GMH_6_differential_abundance.Rmd")) file.copy(system.file("rmarkdown", "templates", "gmh_differential_abundance", "skeleton", "skeleton.Rmd", package = "GMHmicrobiome"), "GMH_6_differential_abundance.Rmd")
    if (!file.exists("GMH_7_multiomics.Rmd")) file.copy(system.file("rmarkdown", "templates", "gmh_multiomics", "skeleton", "skeleton.Rmd", package = "GMHmicrobiome"), "GMH_7_multiomics.Rmd")
    if (!file.exists("GMH_8_machine_learning.Rmd")) file.copy(system.file("rmarkdown", "templates", "gmh_machine_learning", "skeleton", "skeleton.Rmd", package = "GMHmicrobiome"), "GMH_8_machine_learning.Rmd")
  }

  # Create and/or overwrite the initial files
  if (isTRUE(files) & isTRUE(overwrite)) {
    file.copy(system.file("rmarkdown", "templates", import_template, "skeleton", "skeleton.Rmd", package = "GMHmicrobiome"), import_file, overwrite = TRUE)
    file.copy(system.file("rmarkdown", "templates", "gmh_metadata_import", "skeleton", "skeleton.Rmd", package = "GMHmicrobiome"), "GMH_2_metadata_clean.Rmd", overwrite = TRUE)
    file.copy(system.file("rmarkdown", "templates", "gmh_description", "skeleton", "skeleton.Rmd", package = "GMHmicrobiome"), "GMH_3_description.Rmd", overwrite = TRUE)
    file.copy(system.file("rmarkdown", "templates", "gmh_test_variables", "skeleton", "skeleton.Rmd", package = "GMHmicrobiome"), "GMH_4_test_variables.Rmd", overwrite = TRUE)
    file.copy(system.file("rmarkdown", "templates", "gmh_test_variables_code", "skeleton", "skeleton.Rmd", package = "GMHmicrobiome"), "GMH_4x_test_variables_code.Rmd", overwrite = TRUE)
    file.copy(system.file("rmarkdown", "templates", "gmh_beta_diversity", "skeleton", "skeleton.Rmd", package = "GMHmicrobiome"), "GMH_5_betadiversity.Rmd", overwrite = TRUE)
    file.copy(system.file("rmarkdown", "templates", "gmh_differential_abundance", "skeleton", "skeleton.Rmd", package = "GMHmicrobiome"), "GMH_6_differential_abundance.Rmd", overwrite = TRUE)
    file.copy(system.file("rmarkdown", "templates", "gmh_multiomics", "skeleton", "skeleton.Rmd", package = "GMHmicrobiome"), "GMH_7_multiomics.Rmd", overwrite = TRUE)
    file.copy(system.file("rmarkdown", "templates", "gmh_machine_learning", "skeleton", "skeleton.Rmd", package = "GMHmicrobiome"), "GMH_8_machine_learning.Rmd", overwrite = TRUE)
  }
}

