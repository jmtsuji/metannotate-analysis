# filter.R
# MetAnnotate barplot generator
# Copyright Jackson M. Tsuji, 2020 (Neufeld Lab)

#' Counts all HMM hits for the gene of interest
#' 
#' @param metannotate_data Tibble of metannotate data
#' @param gene Character vector (length 1) of a gene to summarize
#' @param collapsed Logical (length 1) - has the table already been collapsed to a taxonomy rank?
#' If TRUE, then the function will sum the "hits" column instead of a simple count.
#' @return Tibble of summarized hit counts
#' @keywords internal
summarize_total_reads <- function(metannotate_data, gene = "rpoB", collapsed = FALSE) {
  # Filter down to gene of interest
  metannotate_summ <- dplyr::filter(metannotate_data, HMM.Family %in% gene)
  
  # Make summary table
  metannotate_summ <- dplyr::group_by(metannotate_summ, Dataset)
  
  if (collapsed == FALSE) {
    
    flog.debug("Assuming table has not been collapsed by taxonomy. Summing rows.")
    metannotate_summ <- dplyr::summarise(metannotate_summ, hits = n())
    
  } else if (collapsed == TRUE) {
    
    flog.debug("Assuming table has been collapsed by taxonomy. Summing 'hits' column.")
    metannotate_summ <- dplyr::summarise(metannotate_summ, hits = sum(hits))
    
  }
  
  return(metannotate_summ)
}

#' Counts all HMM hits for all genes in the metannotate dataset
#' 
#' @param metannotate_data Tibble of metannotate data
#' @param format Character vector (length 1) of either 'wide' or 'long' (tidyr terminology)
#' for the style of output table
#' @param collapsed Logical (length 1) - has the table already been collapsed to a taxonomy rank?
#' If TRUE, then the function will sum the "hits" column instead of a simple count.
#' @return Tibble of summarized hit counts, wide format
#' @keywords internal
summarize_total_reads_all_genes <- function(metannotate_data, format = "wide", collapsed = FALSE) {
  
  # Generate summary for all genes
  gene_summaries <- lapply(unique(metannotate_data$HMM.Family), function(x) {
    summarize_total_reads(metannotate_data, gene = x, collapsed = collapsed) %>%
      tibble::add_column(gene = x)
  }) %>%
    dplyr::bind_rows()

  # If genes have been assigned a factor, then re-assign
  if (is.factor(metannotate_data$HMM.Family) == TRUE) {
    flog.debug("Re-assigning factors to output table")
    gene_summaries$gene <- factor(gene_summaries$gene, levels = levels(metannotate_data$HMM.Family),
                                     ordered = TRUE)
  }
  
  # Make wide format if desired
  if (tolower(format) == "wide") {
    
    flog.debug("Returning wide-format data frame")
    gene_summaries <- gene_summaries %>%
      tidyr::pivot_wider(names_from = gene, values_from = hits)
    
  } else if (tolower(format) == "long") {

    flog.debug("Returning long-format data frame")

  } else {

    stop(paste0("'format' must be either 'wide' or 'long', ",
                "but you specified '", format, "'."))

  }

  return(gene_summaries)
}

#' Filter metannotate data by e-value
#'
#' @aliases filter
#' @description Filters the hits to the desired e-value cutoff and reports stats
#' @param metannotate_data Tibble of metannotate data
#' @param evalue Numeric vector of the evalue cutoff
#' @return List of two:
#' - Tibble of metannotate data, e-value filtered
#' - List of three read count tables. First and second are from before and after e-value filtration.
#'   See summarize_total_reads_all_genes(). The third is the % change between the above two tables.
#' @export
filter_by_evalue <- function(metannotate_data, evalue = 1e-10) {
  # Check initial stats for the HMMs
  read_counts <- list()
  read_counts$original_data <- summarize_total_reads_all_genes(metannotate_data)
  
  metannotate_data <- dplyr::filter(metannotate_data, HMM.E.val <= evalue)
  
  # Then see how it looks after filtering
  read_counts$evalue_filtered_data <- summarize_total_reads_all_genes(metannotate_data)
  
  # Calculate % change
  pseudo_count <- 1e-10 # To prevent divide-by-zero errors
  read_counts$percent_change <- ((dplyr::select(read_counts$evalue_filtered_data, -Dataset) - 
                                    dplyr::select(read_counts$original_data, -Dataset)) / 
                                   (dplyr::select(read_counts$original_data, -Dataset) + pseudo_count) * 100) %>%
    round(digits = 1) %>%
    tibble::add_column(Dataset = dplyr::pull(read_counts$original_data, Dataset), .before = 1)
  
  output_list <- list(metannotate_data, read_counts)
  names(output_list) <- c("metannotate_data", "read_counts")
  
  return(output_list)
}
