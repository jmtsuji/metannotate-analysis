# metannotate_barplots.R
# MetAnnotate barplot generator
# Copyright Jackson M. Tsuji, 2020 (Neufeld Lab)

##### Load libraries ######
library(futile.logger)
library(tools)
library(glue)
library(tibble)
library(tidyselect)
library(plyr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(scales)
library(ggplot2)
#######

# GLOBAL VARIABLES
# Table relating canonical taxonomy values to the colnames of the metannotate table
TAXONOMY_NAMING <- tibble::tibble(taxonomy = c("domain", "phylum", "class", "order", "family", "genus", "species"),
                                  metannotate_colnames = c("Closest.Homolog.Superkingdom", "Closest.Homolog.Phylum",
                                                           "Closest.Homolog.Class", "Closest.Homolog.Order",
                                                           "Closest.Homolog.Family", "Closest.Homolog.Genus",
                                                           "Closest.Homolog.Species"))

#' Loads MetAnnotate data table
#' 
#' @param metannotate_table_filename filepath to the metannotate table (TSV format; can be compressed)
#' @return tibble of metannotate data
#' @export
read_metannotate_data <- function(metannotate_table_filename) {
  
  # Load the data
  metannotate_data <- read.table(metannotate_table_filename, sep = "\t", header = TRUE, 
                                comment.char = "", quote = "", stringsAsFactors = FALSE) %>%
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


#' Counts all HMM hits for the gene of interest
#' 
#' @param metannotate_data Tibble of metannotate data
#' @param gene Character vector (length 1) of a gene to summarize
#' @param collapsed Logical (length 1) - has the table already been collapsed to a taxonomy rank?
#' If TRUE, then the function will sum the "hits" column instead of a simple count.
#' @return Tibble of summarized hit counts
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

#' Collapses metannotate table to the given taxonomic rank
#' 
#' @param metannotate_data Tibble of metannotate data
#' @param taxon Character vector (length 1) giving the taxon name to collapse to
#' Can be: domain, phylum, class, order, family, genus, species (case insensitive)
#' @return Tibble of metannotate data, collapsed with a 'hits' column
#' @export
collapse_metannotate_table_by_taxon <- function(metannotate_data, taxon = "Family") {
  # # N.B., Expected column names of a metannotate table:
  # [1] "Dataset"                      "HMM.Family"                   "ORF"                         
  # [4] "HMM.E.val"                    "Aligned.Length"               "Closest.Homolog"             
  # [7] "X.id.of.Closest.Homolog"      "Closest.Homolog.Species"      "Closest.Homolog.Genus"       
  # [10] "Closest.Homolog.Family"       "Closest.Homolog.Order"        "Closest.Homolog.Class"       
  # [13] "Closest.Homolog.Phylum"       "Closest.Homolog.Superkingdom" "HMM_length"   

  # Check if user-supplied taxonomy is correct
  if ((tolower(taxon) %in% dplyr::pull(TAXONOMY_NAMING, taxonomy)) == FALSE) {
    stop(paste0("Taxon must be a standard seven-rank taxonomy entry; you provided '",
                taxon, "'."))
  }

  # Determine the taxonomy columns to keep while summarizing (should be everything above the chosen rank)
  taxon_number <- match(x = tolower(taxon), table = TAXONOMY_NAMING$taxonomy)
  metannotate_taxa_to_keep <- TAXONOMY_NAMING$metannotate_colnames[1:taxon_number]
  colnames_to_keep <- append(c("Dataset", "HMM.Family", "HMM_length"),
                             metannotate_taxa_to_keep)

  # Summarize the table, keeping all taxonomic ranks above the desired cutoff
  metannotate_summ <- dplyr::group_by_at(metannotate_data, colnames_to_keep) %>%
    dplyr::summarise(hits = n())
  
  return(metannotate_summ)
}

#' Normalize metannotate data
#'
#' @description Double-normlizes metannotate data by HMM length and a single-copy gene marker
#' @param metannotate_data_collapsed Tibble of taxonomy-collapsed metannotate data - see collapse_metannotate_table_by_taxon()
#' @param normalizing_HMM Character vector (length 1) giving the plotting name of an HMM.Family in your table
#' that you want to normalize all genes within each sample to. This should be a reliable single-copy taxonomic marker gene, 
#' e.g., rpoB, dnaK, and so on.
#' @return List of two:
#' - Tibble of metannotate data; 'hits' has now been changed to 'percent_abundance' and is double-normalized
#' - Tibble summarizing total normalized counts for genes
#' @export
normalize_collapsed_metannotate_data <- function(metannotate_data_collapsed, normalizing_HMM = "rpoB") {
  
  ### First normalize by HMM length (for between-HMM comparison)
  metannotate_data_collapsed$hits <- metannotate_data_collapsed$hits / 
    metannotate_data_collapsed$HMM_length
  
  ### Then normalize relative to the normalizing HMM (allows for comparison between samples)
  # Start by getting the total HMM-length normalized counts for the normalizing HMM
  total_hits <- summarize_total_reads_all_genes(metannotate_data_collapsed, format = "wide", collapsed = TRUE)
  
  # Join the normalizing HMM column onto the main table
  total_hits_join <- dplyr::select(total_hits, Dataset, all_of(normalizing_HMM)) %>%
    dplyr::rename(normalizing_HMM_total_hits = normalizing_HMM)
  metannotate_data_collapsed <- dplyr::left_join(metannotate_data_collapsed, total_hits_join, by = "Dataset")
  
  # Normalize metannotate table
  metannotate_data_collapsed$hits <- metannotate_data_collapsed$hits / 
    metannotate_data_collapsed$normalizing_HMM_total_hits * 100
  
  ### Clean up the output
  metannotate_data_collapsed <- dplyr::ungroup(metannotate_data_collapsed) %>%
    dplyr::select(-HMM_length, -normalizing_HMM_total_hits) %>%
    dplyr::rename(percent_abundance = hits)
  
  ### Collect overall summary stats for the normalization
  total_hits_normalized <- (dplyr::select(total_hits, -Dataset) /
    dplyr::pull(total_hits, normalizing_HMM) * 100) %>%
    tibble::add_column(Dataset = dplyr::pull(total_hits, Dataset), .before = 1)
  
  output_list <- list(metannotate_data_collapsed, total_hits_normalized)
  names(output_list) <- c("metannotate_data_normalized", "total_normalized_hits")

  return(output_list)
}

#' Subset high-abundance taxa from metannotate data
#'
#' @description Subsets normalized metannotate data to some desired threshold of most abundance organisms
#' @param metannotate_data_normalized Tibble of normalized metannotate data - see normalize_collapsed_metannotate_data()
#' @param top_x Numeric vector (length 1) giving the subsetting amount you desire.
#' If top_x >=1, the script will return the "top_x most abundant taxa" for each Dataset/HMM.Family
#' If top_x <1, the script will return "all taxa of (top_x * 100%) abundance or greater for each Dataset/HMM.Family - but see below.
#' @param percent_mode If top_x <1, there are two different methods for keeping the most abundant organisms:
#' - "within_sample" -- the normalized % abundance relative to rpoB is used
#' - "within_HMM" -- the percent abundance of that taxon within the specific HMM gene hits is used.
#' You won't notice much of a different between these modes unless one of your HMMs has very few hits and you want to
#' show some of the taxa that were hit. This would be a good time to use 'within_HMM'.
#' @return List of two:
#' - Tibble of metannotate data; 'hits' has now been changed to 'percent_abundance' and is double-normalized
#' - Tibble summarizing total normalized counts for genes
#' @export
subset_normalized_metannotate_data <- function(metannotate_data_normalized, top_x, percent_mode = "within_sample") {
  # Parse the top_x value to determine what the user wants
  if (top_x >= 1) {
    flog.info(paste0("Subsetting to the top ", top_x, " most abundant taxa per sample"))
    
    metannotate_data_summ <- dplyr::group_by(metannotate_data_normalized, Dataset, HMM.Family) %>%
      dplyr::top_n(n = top_x, wt = percent_abundance)
    
  } else if (top_x < 1 && top_x > 0) {
    if (percent_mode == "within_sample") {
      flog.info(paste0("Subsetting all taxa of at least ", top_x * 100, "% normalized relative abundance"))
      
      metannotate_data_summ <- dplyr::filter(metannotate_data_normalized, percent_abundance >= (top_x * 100))
      
    } else if (percent_mode == "within_HMM") {
      flog.info(paste0("Subsetting all taxa of at least ", top_x * 100, "% relative abundance per HMM"))
      
      # Re-calculate the total hits for each HMM
      # TODO - generalize the summary function so I don't have to rename so many columns here
      total_hits <- dplyr::rename(metannotate_data_normalized, hits = percent_abundance) %>%
        summarize_total_reads_all_genes(format = "long", collapsed = TRUE) %>%
        dplyr::rename(HMM.Family = gene, total_percent_abundance_within_HMM = hits)
      
      # Calculate the relative % contribution that the HMM makes to its own total counts
      metannotate_data_normalized <- dplyr::left_join(metannotate_data_normalized, total_hits, 
                                                      by = c("Dataset", "HMM.Family"))
      metannotate_data_normalized$percent_abudance_within_HMM <- metannotate_data_normalized$percent_abundance /
        metannotate_data_normalized$total_percent_abundance_within_HMM * 100
      
      # Filter
      metannotate_data_summ <- dplyr::filter(metannotate_data_normalized, percent_abudance_within_HMM >= (top_x * 100))
      
      # Clean up tibble
      metannotate_data_summ <- dplyr::select(metannotate_data_summ, -total_percent_abundance_within_HMM, 
                                             -percent_abudance_within_HMM)
    }
  }
  
  return(metannotate_data_summ)
}

#' Chooses a nice colour scale for a discrete number of entries
#' 
#' @description Will try to choose a colour-blind friendly colour set if possible, but has to move to
#' less and less ideal colour schemes as the number_of_entries increases
#' @param number_of_entries input numeric vector (integer) with the number of entries that need colours
#' @return a character vector of colours of the same length as the number provided for the number_of_entries
#' @export
choose_discrete_colour_scale <- function(number_of_entries) {
  
  flog.debug(paste0("Generating ", number_of_entries, " colours"))
  
  # Choose the best colour scale based on the number of entries to be plotted
  if ( number_of_entries == 1 ) {
    colour_palette <- "#000000"
  } else if ( number_of_entries == 2 ) {
    colour_palette <- RColorBrewer::brewer.pal(n = 3, name = "Dark2")[1:2]
  } else if ( number_of_entries <= 8 ) {
    colour_palette <- RColorBrewer::brewer.pal(n = number_of_entries, name = "Dark2")
  } else if ( number_of_entries <= 12 ) {
    colour_palette <- RColorBrewer::brewer.pal(n = number_of_entries, name = "Set3")
  } else if ( number_of_entries > 12 ) {
    colour_palette <- scales::hue_pal(h = c(20,290))(number_of_entries)
  } else {
    stop("Something is wrong with the number_of_entries ('", number_of_entries, "'). Is it non-numeric?")
  }
  
  return(colour_palette)
}

#' Generates colour for a metannotate barplot
#' 
#' @param metannotate_data_normalized Tibble of normalized metannotate data - see normalize_collapsed_metannotate_data()
#' @return a tibble of unique taxa with HTML colour codes in the 'colour' column
generate_plotting_colours <- function(metannotate_data_normalized) {
  
  plotting_colour_data <- dplyr::ungroup(metannotate_data_normalized) %>%
    dplyr::select(-percent_abundance, -Dataset, -HMM.Family) %>%
    unique() %>%
    dplyr::arrange_all() # Sort by domain, phylum, and so on (in that order)

  # Pull the lowest taxon rank by detecting the taxonomy that the data has been collapsed to
  plotting_taxon_colname <- TAXONOMY_NAMING$metannotate_colnames[
    TAXONOMY_NAMING$metannotate_colnames %in% colnames(plotting_colour_data)] %>%
    tail(n = 1)
  flog.debug(paste0("Plotting column auto-detected as '", plotting_taxon_colname, "'."))

  unique_taxa <- dplyr::pull(plotting_colour_data, plotting_taxon_colname)

  flog.info(paste0("Generating automatic colour scheme for ", length(unique_taxa), " unique taxa"))
  
  plotting_colour_data$colour <- choose_discrete_colour_scale(length(unique_taxa))
  
  return(plotting_colour_data)
}

#' Generates colour for a metannotate barplot or loads user-defined colours
#' 
#' @description A wrapper for generate_plotting_colours() to handle bigger-picture decision making
#' @param metannotate_data_normalized Tibble of normalized metannotate data - see normalize_collapsed_metannotate_data()
#' @param colouring_template_filename Filename of the colouring template you want to load
#' If the file does not exist, then this function will write a template to that file
#' If 'NA' is entered, then the function will auto-generate colours and continue on
#' @return a tibble of unique taxa with HTML colour codes in the 'colour' column
process_plotting_colours <- function(metannotate_data_normalized, colouring_template_filename = NA) {

  plotting_colour_data <- generate_plotting_colours(metannotate_data_normalized)

  if (is.na(colouring_template_filename) == TRUE) {
    
    flog.debug("Generated plotting colour template without writing to file")
    
  } else if (file.exists(colouring_template_filename) == FALSE) {
    
    flog.info(paste0("Saving plot colour template to '", colouring_template_filename, "'"))
    write.table(plotting_colour_data, colouring_template_filename, sep = "\t", row.names = FALSE,
                col.names = TRUE, quote = FALSE)
    
  } else if (file.exists(colouring_template_filename) == TRUE) {
    flog.info(paste0("Loading plot colour template from '", colouring_template_filename, "'"))

    user_plotting_colour_data <- read.table(colouring_template_filename, sep = "\t", header = TRUE,
                                   comment.char = "", stringsAsFactors = FALSE) %>%
      tibble::as_tibble()

    # Check it matches the auto-generated table
    template_colnum <- length(colnames(plotting_colour_data))
    main_colnames <- colnames(user_plotting_colour_data)[1:template_colnum]

    if (identical(main_colnames, colnames(plotting_colour_data)) == FALSE) {
      stop("Your custom plotting colour table does not have the same first columns as the template. ",
           "Please check carefully.")
    }

    # TODO - this might be a bit aggressive. I could just check the column that will be plotted, for example.
    if (identical(dplyr::arrange_all(plotting_colour_data[,1:(template_colnum-1)]),
                  dplyr::arrange_all(user_plotting_colour_data[,1:(template_colnum-1)])) == FALSE) {
      stop("Your custom plotting colour table does not contain the exact same taxon entries as the template. ",
           "Please check carefully")
    }

    plotting_colour_data <- user_plotting_colour_data

  }
  
  return(plotting_colour_data)
}

#' Generate a ggplot of the MetAnnotate data
#'
#' @description Generates the ggplot of MetAnnotate data
#' @param metannotate_data_normalized Tibble of normalized metannotate data - see normalize_collapsed_metannotate_data()
#' @param hit_totals Tibble of total normalized hits per HMM/sample
#' This is generated in normalize_collapsed_metannotate_data()
#' @param plotting_colour_data Tibble of plot colours generated by process_plotting_colours()
#' @param plotting_taxon The taxonomy rank to summarize bars/bubbles to; MUST MATCH the collapse taxon of the MetAnnotate table
#' @param normalizing_HMM the name of the normalizing HMM (e.g., 'rpoB')
#' @param plot_type Either 'bar' or 'bubble'
#' @param space ggplot setting; 'fixed' or 'free'
#' @param bubble_size_range numeric vector of length two with the small and large bubble sizes
#' @param alpha ggplot value; transparency (for bubble plots)
#' @param bubble_labels logical; show percent labels on the bubbles for bubble plots?
#' @return A ggplot of MetAnnotate data
#' @export
metannotate_ggplot <- function(metannotate_data_normalized, hit_totals, plotting_colour_data,
                               plotting_taxon, normalizing_HMM, plot_type = "bar",
                               space = "free", bubble_size_range = c(1,20), alpha = 0.8,
                               bubble_labels = TRUE) {

  # Check taxon rank
  if ((tolower(plotting_taxon) %in% dplyr::pull(TAXONOMY_NAMING, taxonomy)) == FALSE) {
    stop(paste0("Taxon must be a standard seven-rank taxonomy entry; you provided '", taxon, "'."))
  }
  plotting_taxon_colname <- TAXONOMY_NAMING$metannotate_colnames[match(tolower(plotting_taxon),
                                                                       TAXONOMY_NAMING$taxonomy)]
  plotting_taxon_label <- paste0(substring(plotting_taxon, 1, 1) %>% toupper(),
                                 substring(plotting_taxon, 2, nchar(plotting_taxon)))

  # TODO - can this be moved down to the 'bubble' area without breaking anything?
  if (plot_type == "bubble") {
    metannotate_data_normalized$label <- round(metannotate_data_normalized$percent_abundance, digits = 1)
  }
  
  # Make the base plot
  metannotate_plot <- ggplot(metannotate_data_normalized) +
    theme_bw() +
    # TODO - expose some theme elements to user or make fonts scalable automatically
    theme(axis.title = element_text(size = 12),
          strip.text = element_text(size = 11, face = "italic"),
          strip.background = element_rect(fill = "#e6e6e6"),
          axis.text = element_text(size = 10, colour = "black"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.ticks = element_line(size = 0.5, colour = "black"),
          axis.line = element_line(size = 0.5),
          legend.text = element_text(size = 10, face = "italic"),
          legend.title = element_text(size = 8),
          legend.key = element_rect(colour = "transparent"),
          legend.key.size = unit(6, "mm"),
          legend.spacing = unit(1, "mm"),
          legend.box.just = "left") +
    xlab("Sample")
  
  if (space == "fixed") {
    metannotate_plot <- metannotate_plot +
      facet_grid(HMM.Family ~ ., scales = "free", space = "fixed")
  } else if (space == "free") {
    metannotate_plot <- metannotate_plot +
      facet_grid(HMM.Family ~ ., scales = "free", space = "free")
  } else {
    stop(paste0("'space' must be either 'free' or 'fixed'. You provided '", space, "'."))
  }
  
  if (plot_type == "bar") {
    flog.debug("Generating barplot")

    metannotate_plot <- metannotate_plot +
      geom_bar(data = hit_totals, aes(x = Dataset, weight = percent_abundance), fill = "#808080") +
      geom_bar(aes_string(x = "Dataset", weight = "percent_abundance",
                          fill = plotting_taxon_colname)) +
      scale_fill_manual(values = plotting_colour_data$colour) +
      theme(panel.grid = element_blank(),
            panel.border = element_rect(colour = "transparent"),
            panel.spacing.y = unit(3, "mm")) +
      guides(fill = guide_legend(ncol = 1, title = element_blank())) +
      ylab(paste0("Gene hits relative to ", normalizing_HMM, " (%; normalized)"))
    
  } else if (plot_type == "bubble") {
    flog.debug("Generating bubble plot")

    # TODO - allow user to toggle which taxon rank to use
    legend_taxon_colname <- "Closest.Homolog.Phylum"
    legend_taxon <- "Phylum"
    fill_colours <- dplyr::pull(metannotate_data_normalized, legend_taxon_colname) %>%
        unique() %>%
        length() %>%
        choose_discrete_colour_scale()

    metannotate_plot <- metannotate_plot +
      geom_point(aes_string(x = "Dataset", y = plotting_taxon_colname,
                            size = "percent_abundance", fill = legend_taxon_colname),
                 shape = 21, alpha = alpha) +
      scale_size_continuous(range = bubble_size_range) +
      scale_fill_manual(values = fill_colours) +
      theme(axis.text.y = element_text(size = 5, face = "italic")) +
      guides(fill = guide_legend(title = legend_taxon)) +
      ylab(paste0(plotting_taxon_label, " of closest homologue"))
    
    if (bubble_labels == TRUE) {
      metannotate_plot <- metannotate_plot +
        geom_text(aes_string(x = "Dataset", y = plotting_taxon_colname, label = "label"),
                  size = 2) +
        guides(size = FALSE)
    } else if (bubble_labels == FALSE) {
      metannotate_plot <- metannotate_plot +
        guides(size = guide_legend(title = paste0("Gene hits relative to \n", normalizing_HMM, " (%; normalized)"),
                                   override.aes = list(fill = "#4d4d4d")))
    }
    
  } else {

    stop(paste0("'plot_type' must be either 'bar' or 'bubble'; ",
                "you provided '", plot_type, "'."))

  }
  
  return(metannotate_plot)
}

#' Wrapper for convenient generation of MetAnnotate ggplot data
#'
#' @description Wrapper to generate a ggplot of MetAnnotate data with subset, colours, labels, and so on
#' @param metannotate_data_normalized_list List output of normalize_collapsed_metannotate_data()
#' @param colouring_template_filename Filename of the colouring template you want to load
#' If the file does not exist, then this function will write a template to that file
#' If 'NA' is entered, then the function will auto-generate colours and continue on
#' @param top_x Numeric vector (length 1) giving the subsetting amount you desire.
#' If top_x >=1, the script will return the "top_x most abundant taxa" for each Dataset/HMM.Family
#' If top_x <1, the script will return "all taxa of (top_x * 100%) abundance or greater for each Dataset/HMM.Family - but see below.
#' @param percent_mode If top_x <1, there are two different methods for keeping the most abundant organisms:
#' - "within_sample" -- the normalized % abundance relative to rpoB is used
#' - "within_HMM" -- the percent abundance of that taxon within the specific HMM gene hits is used.
#' You won't notice much of a different between these modes unless one of your HMMs has very few hits and you want to
#' show some of the taxa that were hit. This would be a good time to use 'within_HMM'.
#' @param normalizing_HMM Name of the normalizing HMM (e.g., 'rpoB')]; specify 'auto' to attempt auto-detection
#' @param plot_normalizing_HMM Retain the normalizing_HMM in the final ggplot?
#' @param ... Other fine-tuned plotting options controlled by metannotate_ggplot()
#' @param dump_raw_data Return the normalized and subsetted table in lieu of a ggplot
#' @return A ggplot of MetAnnotate data (or raw data; see above)
#' @export
metannotate_plotter <- function(metannotate_data_normalized_list, colouring_template_filename = NA,
                                top_x = NA, percent_mode = "within_sample", normalizing_HMM = "auto",
                                plot_normalizing_HMM = TRUE, dump_raw_data = FALSE, ...) {
  # # Example column names of the plotting table, if collapsed to family
  # [1] "Dataset"                      "HMM.Family"                   "Closest.Homolog.Superkingdom"
  # [4] "Closest.Homolog.Phylum"       "Closest.Homolog.Class"        "Closest.Homolog.Order"       
  # [7] "Closest.Homolog.Family"       "percent_abundance" 

  # Extract list components
  metannotate_data <- metannotate_data_normalized_list$metannotate_data_normalized
  hit_totals <- tidyr::pivot_longer(metannotate_data_normalized_list$total_normalized_hits, -Dataset,
                                    names_to = "HMM.Family", values_to = "percent_abundance")
  hit_totals$HMM.Family <- factor(hit_totals$HMM.Family, levels = unique(hit_totals$HMM.Family), ordered = TRUE)
  
  # Detect the taxonomy that the data has been collapsed to
  plotting_taxon_colname <- TAXONOMY_NAMING$metannotate_colnames[
    TAXONOMY_NAMING$metannotate_colnames %in% colnames(metannotate_data)] %>%
    tail(n = 1)
  plotting_taxon <- TAXONOMY_NAMING$taxonomy[match(plotting_taxon_colname,
                                                   TAXONOMY_NAMING$metannotate_colnames)]
  flog.debug(paste0("Plotting input dataframe has been collapsed to the '", plotting_taxon, "' level."))
  
  # Subset to the desired top_x cutoff
  # TODO - longer-term move subsetting out of this script for clarity
  metannotate_data <- subset_normalized_metannotate_data(metannotate_data, top_x, percent_mode = percent_mode)

  # Determine normalizing HMM for labelling on the plot
  if (normalizing_HMM == "auto") {
    normalizing_HMM  <- dplyr::group_by(hit_totals, HMM.Family) %>%
      dplyr::summarise(mean_abund = mean(percent_abundance)) %>%
      dplyr::filter(mean_abund == 100) %>%
      dplyr::pull(HMM.Family)

    if (length(normalizing_HMM) != 1) {
      stop("Auto-detection of the normalizing_HMM failed; you'll have to specify it manually.")
    }
  }

  # Remove normalizing_HMM from the hit totals to avoid plotting an extraneous bar at 100%
  hit_totals <- dplyr::filter(hit_totals, HMM.Family != normalizing_HMM)

  # Optionally remove the normalizing HMM from the final plot altogether
  if (plot_normalizing_HMM == FALSE) {
    metannotate_data <- dplyr::filter(metannotate_data, HMM.Family != normalizing_HMM)

    if (nrow(metannotate_data) == 0) {
      stop(paste0("Looks like you have no data left after removing the normalizing_HMM; ",
                  "did you only include one HMM in your data? ",
                  "Probably best that you leave 'plot_normalizing_HMM = TRUE'."))
    }

  } else if (plot_normalizing_HMM != TRUE) {
    stop(paste0("'plot_normalizing_HMM' must be either TRUE or FALSE; you specified '",
                      plot_normalizing_HMM, "'."))
  }

  # Make or read in a plotting colour table; or generate auto-colours
  plotting_colour_data <- process_plotting_colours(metannotate_data, colouring_template_filename)

  # Make the plotting column into an ordered factor based on the plotting_colours order
  metannotate_data[,plotting_taxon_colname] <- factor(dplyr::pull(metannotate_data, plotting_taxon_colname),
                                                      levels = unique(dplyr::pull(plotting_colour_data,
                                                                                  plotting_taxon_colname)),
                                                      ordered = TRUE)

  if (dump_raw_data == TRUE) {
    # TODO - this is bad design; function should not output something totally different given a flag...
    flog.info("Dumping raw data in lieu of a ggplot")
    metannotate_plot <- metannotate_data

  } else {
    flog.info("Creating the ggplot")
    metannotate_plot <- metannotate_ggplot(metannotate_data_normalized = metannotate_data,
                                           hit_totals = hit_totals,
                                           plotting_colour_data = plotting_colour_data,
                                           plotting_taxon = plotting_taxon,
                                           normalizing_HMM = normalizing_HMM,
                                           ...)
  }

  return(metannotate_plot)
}

#' Explore MetAnnotate data
#'
#' @description high-level exploration function for examining MetAnnotate data
#' @param metannotate_data_mapped The tibble output by map_naming_information()
#' @param evalue E-value cutoff for HMM hits
#' @param taxon Character vector (length 1) giving the taxon name to collapse to
#' Can be: domain, phylum, class, order, family, genus, species (case insensitive)
#' @param normalizing_HMM Name of the normalizing HMM (e.g., 'rpoB')]; specify 'auto' to attempt auto-detection
#' @param top_x Numeric vector (length 1) giving the subsetting amount you desire.
#' If top_x >=1, the script will return the "top_x most abundant taxa" for each Dataset/HMM.Family
#' If top_x <1, the script will return "all taxa of (top_x * 100%) abundance or greater for each Dataset/HMM.Family - but see below.
#' @param percent_mode If top_x <1, there are two different methods for keeping the most abundant organisms:
#' - "within_sample" -- the normalized % abundance relative to rpoB is used
#' - "within_HMM" -- the percent abundance of that taxon within the specific HMM gene hits is used.
#' You won't notice much of a different between these modes unless one of your HMMs has very few hits and you want to
#' show some of the taxa that were hit. This would be a good time to use 'within_HMM'.
#' @param colouring_template_filename Filename of the colouring template you want to load
#' If the file does not exist, then this function will write a template to that file
#' If 'NA' is entered, then the function will auto-generate colours and continue on
#' #' @param ... Other fine-tuned plotting options controlled by metannotate_plotter()
#' @return A ggplot of MetAnnotate data
#' @export
explore_metannotate_data <- function(metannotate_data_mapped, evalue = 1e-10, taxon = "Family",
                                     normalizing_HMM = "rpoB", top_x = 0.02, percent_mode = "within_sample",
                                     colouring_template_filename = NA, ...) {
  
  # Filter by e-value cutoff and report stats to user
  flog.info(paste0("Filtering by e-value cutoff of ", evalue))
  metannotate_data_filtered <- filter_by_evalue(metannotate_data_mapped, evalue = evalue)
  metannotate_data <- metannotate_data_filtered$metannotate_data
  flog.info("Percent change from e-value filtration:")
  print(metannotate_data_filtered$read_counts$percent_change)
  # TODO - optionally output the info to the user
  
  # Collapse the table to the desired taxonomic rank
  flog.info(paste0("Collapsing table to taxonomic rank '", taxon, "'"))
  metannotate_data_collapsed <- collapse_metannotate_table_by_taxon(metannotate_data, taxon = taxon)
  
  # Normalize the data by HMM length
  flog.info("Normalizing data")
  metannotate_data_normalized_list <- normalize_collapsed_metannotate_data(metannotate_data_collapsed, 
                                                                           normalizing_HMM = normalizing_HMM)
  flog.info("Total normalized % abundance of analyzed genes compared to the marker gene:")
  print(metannotate_data_normalized_list$total_normalized_hits)
  
  # Make plots
  flog.info("Plotting data")
  metannotate_plot <- metannotate_plotter(metannotate_data_normalized_list = metannotate_data_normalized_list,
                                          colouring_template_filename = colouring_template_filename,
                                          top_x = top_x,
                                          percent_mode = percent_mode,
                                          normalizing_HMM = normalizing_HMM,
                                          ...)

  return(metannotate_plot)
}
