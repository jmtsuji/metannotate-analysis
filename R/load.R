# load.R
# Load and map naming info for metannotate files
# Copyright Jackson M. Tsuji, 2020 (Neufeld Lab)

#' Loads MetAnnotate data table
#'
#' @aliases load_metannotate_data load read
#' @param metannotate_table_filename filepath to the metannotate table (TSV format; can be compressed)
#' @return tibble of metannotate data
#' @export
read_metannotate_data <- function(metannotate_table_filename) {
  
  # Load the data
  metannotate_data <- read.table(metannotate_table_filename, sep = "\t", header = TRUE, 
                                comment.char = "", stringsAsFactors = FALSE) %>%
    tibble::as_tibble()

  # Check for the key required cols
  key_columns <- c("Dataset", "HMM.Family", "ORF", "HMM.E.val", "Closest.Homolog.Species", "Closest.Homolog.Genus",
                   "Closest.Homolog.Family", "Closest.Homolog.Order", "Closest.Homolog.Class",
                   "Closest.Homolog.Phylum", "Closest.Homolog.Superkingdom")
  if (all(key_columns %in% colnames(metannotate_data)) == FALSE) {
    missing_cols <- glue::glue_collapse(key_columns[!key_columns %in% colnames(metannotate_data)], sep = ", ")
    stop("The provided metannotate data table is missing key columns: '", missing_cols,
         "'. Please check carefully. It's possible that they are mis-named.")
  }

  return(metannotate_data)
}

#' Creates templates for user-specified HMM.Family and Dataset naming data
#'
#' @param metannotate_data Tibble of metannotate data
#' @param write_tables Logical; write the tables to file?
#' @param hmm_info_filename File to optionally write HMM information to
#' @param dataset_info_filename File to optionally write dataset information to
#' @return List of two tibbles: hmm_info and dataset_info; user can fill these out
#' @keywords internal
create_setup_templates <- function(metannotate_data, write_tables = FALSE,
                                   hmm_info_filename = "hmm_info_template.tsv",
                                   dataset_info_filename = "dataset_info_template.tsv") {
  # Make template for HMM naming
  # Columns will be: c("HMM.Family", "raw_name", "HMM_length", "notes")
  hmm_info <- tibble::tibble(raw_name = unique(metannotate_data$HMM.Family),
                             HMM.Family = "", HMM_length = "", notes = "")
  
  if (write_tables == TRUE) {
    write.table(hmm_info, file = hmm_info_filename, sep = "\t",
                row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
    
  # Make template for sample naming
  # Columns will be: c("raw_name", "Dataset")
  dataset_info <- tibble::tibble(raw_name = unique(metannotate_data$Dataset),
                                 Dataset = "")

  if (write_tables == TRUE) {
    write.table(dataset_info, file = dataset_info_filename, sep = "\t",
                row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
  
  output_list <- list(hmm_info, dataset_info)
  names(output_list) <- c("hmm_info", "dataset_info")
  return(output_list)
}

#' Map user information onto metannotate data
#' 
#' @description Applies user-provided HMM.Family/Dataset naming/length information onto the metannotate data.
#' A few interesting quirks about this function:
#' - The order of the elements in the provided tables is meaningful. It determines the order of the factor for HMM
#'   Dataset names. This dictates the order during plotting
#' - You can give two elements the same plotting name (e.g., R1 and R2 reads). They will be merged together downstream.
#' - Only elements in these tables will be used in the script. Any elements you omit from the info tables (e.g., R2 reads,
#'   if you don't want to plot these) will be deleted from the metannotate table.
#' @param metannotate_data Tibble of metannotate data
#' @param hmm_naming_info_filename Filepath to HMM naming template exported from create_setup_templates()
#' @param dataset_naming_info_filename Filepath to Dataset naming template exported from create_setup_templates()
#' @return Tibble of metannotate data (renamed with factors)
#' @export
map_naming_information <- function(metannotate_data, hmm_naming_info_filename,
                                   dataset_naming_info_filename) {
  
  # Read the HMM info
  hmm_info <- read.table(hmm_naming_info_filename, sep = "\t", header = TRUE,
                         stringsAsFactors = FALSE, comment.char = "") %>%
    tibble::as_tibble()

  key_columns <- c("HMM.Family", "raw_name", "HMM_length")
  if (all(key_columns %in% colnames(hmm_info)) == FALSE) {
    missing_cols <- glue::glue_collapse(key_columns[!key_columns %in% colnames(hmm_info)], sep = ", ")
    stop("The provided HMM naming file is missing key columns: '", missing_cols,
         "'. Please check carefully. It's possible that they are mis-named.")
  }

  # Read the dataset info
  dataset_info <- read.table(dataset_naming_info_filename, sep = "\t", header = TRUE, 
                             stringsAsFactors = FALSE, comment.char = "") %>%
    tibble::as_tibble()

  key_columns <- c("raw_name", "Dataset")
  if (all(key_columns %in% colnames(dataset_info)) == FALSE) {
    missing_cols <- glue::glue_collapse(key_columns[!key_columns %in% colnames(dataset_info)], sep = ", ")
    stop("The provided dataset naming file is missing key columns: '", missing_cols,
         "'. Please check carefully. It's possible that they are mis-named.")

  }

  # Remove elements of the metannotate table not in the naming info tables
  metannotate_data <- dplyr::filter(metannotate_data, (HMM.Family %in% hmm_info$raw_name) &
                                      (Dataset %in% dataset_info$raw_name))
  
  # Order elements based on row names and rename
  # N.B., levels = unique(...) was used for applying the factor in case the user wanted to combine two types of data
  #       (e.g., R1 and R2 reads) to have the same name.
  metannotate_data$HMM.Family <- plyr::mapvalues(x = metannotate_data$HMM.Family,
                                                 from = hmm_info$raw_name, to = hmm_info$HMM.Family)
  metannotate_data$HMM.Family <- factor(metannotate_data$HMM.Family, levels = unique(hmm_info$HMM.Family), ordered = TRUE)
  
  metannotate_data$Dataset <- plyr::mapvalues(x = metannotate_data$Dataset,
                                              from = dataset_info$raw_name, to = dataset_info$Dataset)
  metannotate_data$Dataset <- factor(metannotate_data$Dataset, levels = unique(dataset_info$Dataset), ordered = TRUE)
  
  # Lastly, add the HMM lengths
  hmm_info <- dplyr::select(hmm_info, HMM.Family, HMM_length)
  hmm_info$HMM.Family <- factor(hmm_info$HMM.Family, levels = unique(hmm_info$HMM.Family), ordered = TRUE)
  metannotate_data <- dplyr::left_join(metannotate_data, hmm_info, by = "HMM.Family")
  
  return(metannotate_data)
}
