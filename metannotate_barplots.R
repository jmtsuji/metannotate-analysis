# metannotate_barplots.R
# MetAnnotate barplot generator
# Copyright Jackson M. Tsuji, 2018 (Neufeld Lab)

##### Load libraries ######
library(argparser)
library(parallel)
library(futile.logger)
library(roxygen2)
library(tools)
library(glue)
library(tibble)
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(scales)
#######

#' Loads MetAnnotate data table
#' 
#' @param metannotate_table_filename filepath to the metannotate table (TSV format; can be compressed)
#' @return tibble of metannotate data
#' @export
read_metannotate_tibble <- function(metannotate_table_filename) {
  
  # Load the data with care
  metannotate_data <- tibble::as.tibble(read.table(metannotate_table_filename, sep = "\t", header = TRUE, 
                                  comment.char = "", quote = "", stringsAsFactors = FALSE))
  
  return(metannotate_data)
}

#' Creates templates for user-specified HMM.Family and Dataset naming data
#' 
#' @param metannotate_data Tibble of metannotate data
#' @param write_tables Logical; write the tables to files 'hmm_info_template.tsv' and 'dataset_info_template.tsv'?
#' @return List of two tibbles: hmm_info and dataset_info; user can fill these out
#' @export
create_setup_templates <- function(metannotate_data, write_tables = FALSE) {
  # Make template for HMM naming
  # Columns will be: c("HMM.Family", "raw_name", "HMM_length", "notes")
  hmm_info <- tibble::tibble(raw_name = unique(metannotate_data$HMM.Family))
  hmm_info <- tibble::add_column(hmm_info, HMM.Family = "", HMM_length = "", notes = "")
  
  if (write_tables == TRUE) {
    hmm_info_filename <- "hmm_info_template.tsv"
    write.table(hmm_info, file = hmm_info_filename, sep = "\t", row.names = F, col.names = T)
  }
    
  # Make template for sample naming
  dataset_info <- tibble::tibble(raw_name = unique(metannotate_data$Dataset))
  dataset_info <- tibble::add_column(dataset_info, Dataset = "")
  
  if (write_tables == TRUE) {
    dataset_info_filename <- "dataset_info_template.tsv"
    write.table(dataset_info, file = dataset_info_filename, sep = "\t", row.names = F, col.names = T)
  }
  
  output_list <- list(hmm_info, dataset_info)
  names(output_list) <- c("hmm_info", "dataset_info")
  return(output_list)
}

#' Applies user-provided HMM.Family/Dataset naming/length information onto the metannotate data
#' 
#' @description A few interesting quirks about this function:
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
map_naming_information <- function(metannotate_data, hmm_naming_info_filename, dataset_naming_info_filename) {
  
  # Read the data
  hmm_info <- tibble::as.tibble(read.table(hmm_naming_info_filename, sep = "\t", header = TRUE, 
                                           stringsAsFactors = FALSE, comment.char = "", quote = ""))
  dataset_info <- tibble::as.tibble(read.table(dataset_naming_info_filename, sep = "\t", header = TRUE, 
                                               stringsAsFactors = FALSE, comment.char = "", quote = ""))
  
  # Remove elements of the metannotate tablenot in the naming info tables
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
#' @export
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
#' @param format Character vector (length 1) of either 'wide' or 'long' (reshape2 terminology) 
#' for the style of output table
#' @param collapsed Logical (length 1) - has the table already been collapsed to a taxonomy rank?
#' If TRUE, then the function will sum the "hits" column instead of a simple count.
#' @return Tibble of summarized hit counts, wide format
#' @export
summarize_total_reads_all_genes <- function(metannotate_data, format = "wide", collapsed = FALSE) {
  
  # Generate summary for each gene
  gene_summaries <- lapply(unique(metannotate_data$HMM.Family), function(x) {
    summarize_total_reads(metannotate_data, gene = x, collapsed = collapsed)
  })
  names(gene_summaries) <- unique(metannotate_data$HMM.Family)
  
  # Combine the data
  gene_summaries_df <- dplyr::bind_rows(lapply(names(gene_summaries), function(x) {
    gene_summary <- gene_summaries[[x]]
    gene_summary <- tibble::add_column(gene_summary, gene = x)
    return(gene_summary)
  }))
  gene_summaries_df$gene <- factor(gene_summaries_df$gene, levels = levels(metannotate_data$HMM.Family),
                                   ordered = TRUE)
  
  # Make wide format if desired
  if (tolower(format) == "wide") {
    
    flog.debug("Returning wide-format data frame")
    gene_summaries_df <- reshape2::dcast(gene_summaries_df, Dataset ~ gene, value.var = "hits")
    return(gene_summaries_df)
    
  } else if (tolower(format) == "long") {
    
    flog.debug("Returning long-format data frame")
    return(gene_summaries_df)
    
  }
  
  
}

#' Filters the hits to the desired e-value cutoff
#' 
#' @param metannotate_data Tibble of metannotate data
#' @param evalue Numeric vector of the evalue cutoff
#' @return List of two:
#' - Tibble of metannotate data, e-value filtered
#' - List of three read count tables. First and second are from before and after e-value filtration.
#'   See summarize_total_reads_all_genes(). The third is the % change between the above two tables.
#' @export
filter_by_evalue <- function(metannotate_data, evalue = 1e-40) {
  # Check initial stats for the HMMs
  read_counts <- list()
  read_counts$original_data <- summarize_total_reads_all_genes(metannotate_data)
  
  metannotate_data <- dplyr::filter(metannotate_data, HMM.E.val <= evalue)
  
  # Then see how it looks after filtering
  read_counts$evalue_filtered_data <- summarize_total_reads_all_genes(metannotate_data)
  
  # Calculate % change (pseudocount added into the denominator to prevent Inf/NaN values)
  read_counts$percent_change <- round((read_counts$evalue_filtered_data[,c(-1)] - read_counts$original_data[,c(-1)]) /
    (read_counts$original_data[,c(-1)] + 1e-10) * 100, digits = 1)
  read_counts$percent_change <- tibble::add_column(read_counts$percent_change, 
                                                   Dataset = read_counts$original_data$Dataset,
                                                   .before = 1)
  
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
  
  # HARD-CODED table relating canonical taxonomy values to the colnames of the metannotate table
  taxonomy_naming <- tibble::tibble(taxonomy = c("domain", "phylum", "class", "order", "family", "genus", "species"),
                                    metannotate_colnames = c("Closest.Homolog.Superkingdom", "Closest.Homolog.Phylum",
                                                             "Closest.Homolog.Class", "Closest.Homolog.Order", "Closest.Homolog.Family", 
                                                             "Closest.Homolog.Genus", "Closest.Homolog.Species"))
  
  # Convert user input into a metannotate-friendly name. Also get number code.
  taxon_metannotate <- taxonomy_naming$metannotate_colnames[taxonomy_naming$taxonomy %in% tolower(taxon)]
  taxon_number <- match(x = tolower(taxon), table = taxonomy_naming$taxonomy)
  
  # Determine the taxonomy columsn to keep while summarizing (should be everything above the chosen rank)
  # We know that the metannotate taxonomy columns run in reverse order from columns 8-14 (see above)
  taxon_colnum_metannotate <- (7 - taxon_number) + 8
  taxon_colnum_series <- rev(taxon_colnum_metannotate:14)
  
  # Summarize the table, keeping all taxonomic ranks above the desired cutoff
  metannotate_summ <- dplyr::group_by_at(metannotate_data, c(1,2,15,taxon_colnum_series))
  metannotate_summ <- dplyr::summarise(metannotate_summ, hits = n())
  
  return(metannotate_summ)
  
}

#' Double-normlizes metannotate data by HMM length and a single-copy gene marker
#' 
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
  total_hits_join <- dplyr::select(total_hits, Dataset, normalizing_HMM)
  total_hits_join <- dplyr::rename(total_hits_join, normalizing_HMM_total_hits = normalizing_HMM)
  metannotate_data_collapsed <- dplyr::left_join(metannotate_data_collapsed, total_hits_join, by = "Dataset")
  total_hits_join <- NULL
  
  # Normalize metannotate table
  metannotate_data_collapsed$hits <- metannotate_data_collapsed$hits / 
    metannotate_data_collapsed$normalizing_HMM_total_hits * 100
  
  ### Clean up the output
  metannotate_data_collapsed <- dplyr::ungroup(metannotate_data_collapsed)
  metannotate_data_collapsed <- dplyr::select(metannotate_data_collapsed, -HMM_length, -normalizing_HMM_total_hits)
  metannotate_data_collapsed <- dplyr::rename(metannotate_data_collapsed, percent_abundance = hits)
  
  ### Collect overall summary stats for the normalization
  total_hits_normalized <- dplyr::select(total_hits, -Dataset) / dplyr::pull(total_hits, normalizing_HMM) * 100
  total_hits_normalized <- tibble::add_column(total_hits_normalized, Dataset = total_hits$Dataset, .before = 1)
  
  output_list <- list(metannotate_data_collapsed, total_hits_normalized)
  names(output_list) <- c("metannotate_data_normalized", "total_normalized_hits")
  return(output_list)
  
}

#' Subsets normalized metannotate data to some desired threshold of most abundance organisms
#' 
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
    flog.info(glue::glue("Subsetting to the top ", top_x, " most abundant taxa per sample"))
    
    metannotate_data_summ <- dplyr::group_by(metannotate_data_normalized, Dataset, HMM.Family)
    metannotate_data_summ <- dplyr::top_n(metannotate_data_summ, n = top_x, wt = percent_abundance)
    
  } else if (top_x < 1 && top_x > 0) {
    if (percent_mode == "within_sample") {
      flog.info(glue::glue("Subsetting all taxa of at least ", top_x * 100, "% normalized relative abundance"))
      
      metannotate_data_summ <- dplyr::filter(metannotate_data_normalized, percent_abundance >= (top_x * 100))
      
    } else if (percent_mode == "within_HMM") {
      flog.info(glue::glue("Subsetting all taxa of at least ", top_x * 100, "% relative abundance per HMM"))
      
      # Re-calculate the total hits for each HMM (scripts have to be coerced a bit to work here)
      metannotate_data_normalized <- dplyr::rename(metannotate_data_normalized, hits = percent_abundance)
      total_hits <- summarize_total_reads_all_genes(metannotate_data_normalized, format = "long", collapsed = TRUE)
      metannotate_data_normalized <- dplyr::rename(metannotate_data_normalized, percent_abundance = hits)
      total_hits <- dplyr::rename(total_hits, HMM.Family = gene, total_percent_abundance_within_HMM = hits)
      
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
  
  flog.debug(glue::glue("Generating ", number_of_entries, " colours"))
  
  # Choose the best colour scale based on the number of entries to be plotted
  if ( number_of_entries == 1 ) {
    colour_palette <- "#000000"
  } else if ( number_of_entries == 2 ) {
    colour_palette <- RColorBrewer::brewer.pal(n = 3, name = "Dark2")[c(1:2)]
  } else if ( number_of_entries <= 8 ) {
    colour_palette <- RColorBrewer::brewer.pal(n = number_of_entries, name = "Dark2")
  } else if ( number_of_entries <= 12 ) {
    colour_palette <- RColorBrewer::brewer.pal(n = number_of_entries, name = "Set3")
  } else if ( number_of_entries > 12 ) {
    colour_palette <- scales::hue_pal(h = c(20,290))(number_of_entries)
  } else {
    stop("Something is wrong with the number_of_entries ('", number_of_entries, "'). Is it non-numeric? Exiting...")
  }
  
  return(colour_palette)
  
}

#' Generates colour for a metannotate barplot
#' 
#' @param metannotate_data_normalized Tibble of normalized metannotate data - see normalize_collapsed_metannotate_data()
#' @return a tibble of unique taxa with HTML colour codes in the 'colour' column
#' @export
generate_plotting_colours <- function(metannotate_data_normalized) {
  
  metannotate_data_normalized <- dplyr::ungroup(metannotate_data_normalized)
  
  plotting_colours <- unique(dplyr::select(metannotate_data_normalized, -percent_abundance, -Dataset, -HMM.Family))
  plotting_colours <- dplyr::arrange_all(plotting_colours) # Sort by domain, phylum, and so on (in that order)
  
  unique_taxa <- dplyr::pull(plotting_colours, ncol(plotting_colours))
  # TODO confirm somehow that the lowest taxon rank is actually the last column as this assumes.
  
  flog.info(glue::glue("Generating automatic colour scheme for ", length(unique_taxa), " unique taxa"))
  
  plotting_colours$colour <- choose_discrete_colour_scale(length(unique_taxa))
  
  return(plotting_colours)
}

#' Generates colour for a metannotate barplot or loads user-defined colours
#' 
#' @description A wrapper for generate_plotting_colours() to handle bigger-picture decision making
#' @param metannotate_data_normalized Tibble of normalized metannotate data - see normalize_collapsed_metannotate_data()
#' @param colouring_template_filename Filename of the colouring template you want to load
#' If the file does not exist, then this function will write a template to that file
#' If 'NA' is entered, then the function will auto-generate colours and continue on
#' @return a tibble of unique taxa with HTML colour codes in the 'colour' column
#' @export
process_plotting_colours <- function(metannotate_data_normalized, colouring_template_filename = NA) {
  
  if (is.na(colouring_template_filename) == TRUE) {
    
    plotting_colours <- generate_plotting_colours(metannotate_data_normalized)
    
  } else if (file.exists(colouring_template_filename) == FALSE) {
    
    plotting_colours <- generate_plotting_colours(metannotate_data_normalized)
    
    flog.info(glue::glue("Saving plot colour template to '", colouring_template_filename, "'"))
    write.table(plotting_colours, colouring_template_filename, sep = "\t", row.names = FALSE,
                col.names = TRUE)
    
  } else if (file.exists(colouring_template_filename) == TRUE) {
    flog.info(glue::glue("Loading plot colour template from '", colouring_template_filename, "'"))
    
    plotting_colours <- tibble::as.tibble(read.table(colouring_template_filename, sep = "\t", header = TRUE, 
                                                     comment.char = "", stringsAsFactors = FALSE))
    
  }
  
  return(plotting_colours)
  
}


metannotate_ggplot <- function(metannotate_data_normalized, type = "bar", plotting_colours, 
                               collapse_taxonomy_metannotate, hit_totals, normalizing_HMM) {
  
  # Make the base plot
  metannotate_plot <- ggplot(metannotate_data_normalized) +
    facet_grid(HMM.Family ~ ., scales = "free") +
    theme_bw() +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10, colour = "black"),
          strip.text = element_text(size = 11, face = "italic"), strip.background = element_rect(fill = "#e6e6e6"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.ticks = element_line(size = 0.5, colour = "black"), axis.line = element_line(size = 0.5),
          legend.text = element_text(size = 10, face = "italic"), legend.title = element_text(size = 8),
          legend.key = element_rect(colour = "transparent"), legend.key.size = unit(6, "mm"),
          legend.spacing = unit(1, "mm"), legend.box.just = "left") +
    xlab("Sample")
  
  if (type == "bar") {
    flog.debug("Generating barplot")
    
    metannotate_plot <- metannotate_plot +
      geom_bar(data = hit_totals, aes(x = Dataset, weight = percent_abundance), fill = "#808080") +
      geom_bar(aes_string(x = "Dataset", weight = "percent_abundance", 
                                                   fill = collapse_taxonomy_metannotate)) +
      scale_fill_manual(values = plotting_colours$colour) +
      theme(panel.grid = element_blank(), panel.border = element_rect(colour = "transparent"), 
            panel.spacing.y = unit(3, "mm")) +
      guides(fill = guide_legend(ncol = 1, title = element_blank())) +
      ylab(paste("Gene hits relative to ", normalizing_HMM, " (%; normalized)", sep = ""))
    
  } else if (type == "bubble") {
    flog.debug("Generating bubble plot")
    
    metannotate_plot <- metannotate_plot +
      geom_point(aes_string(x = "Dataset", y = collapse_taxonomy_metannotate,
                            size = "percent_abundance", fill = "Closest.Homolog.Phylum"), shape = 21,
                 alpha = 0.8) +
      scale_size_continuous(range = c(1,10)) +
      scale_fill_manual(values = 
                          choose_discrete_colour_scale(length(unique(
                            dplyr::pull(metannotate_data_normalized, Closest.Homolog.Phylum))))) +
      theme(axis.text.y = element_text(size = 5, face = "italic")) +
      guides(size = guide_legend(title = paste("Gene hits relative to \n", normalizing_HMM, 
                                               " (%; normalized)", sep = ""), 
                                 override.aes = list(fill = "#4d4d4d")),
             fill = guide_legend(title = "Phylum")) +
      ylab(paste(strsplit(collapse_taxonomy_metannotate, split = ".", fixed = TRUE)[[1]][3],
                 " of closest homologue", sep = ""))
    
  }
  
  return(metannotate_plot)
  
}

metannotate_plotter <- function(metannotate_data_normalized_list, type = "bar", top_x = NA, highlight_taxon = NA,
                                colouring_template_filename = NA, percent_mode = "within_sample") {
  # # Example column names of the plotting table, if collapsed to family
  # [1] "Dataset"                      "HMM.Family"                   "Closest.Homolog.Superkingdom"
  # [4] "Closest.Homolog.Phylum"       "Closest.Homolog.Class"        "Closest.Homolog.Order"       
  # [7] "Closest.Homolog.Family"       "percent_abundance" 
  
  # HARD-CODED table relating canonical taxonomy values to the colnames of the metannotate table
  taxonomy_naming <- tibble::tibble(taxonomy = c("domain", "phylum", "class", "order", "family", "genus", "species"),
                                    metannotate_colnames = c("Closest.Homolog.Superkingdom", "Closest.Homolog.Phylum",
                                                             "Closest.Homolog.Class", "Closest.Homolog.Order", "Closest.Homolog.Family", 
                                                             "Closest.Homolog.Genus", "Closest.Homolog.Species"))
  
  # Extract list components
  metannotate_data <- metannotate_data_normalized_list$metannotate_data_normalized
  hit_totals_wide <- metannotate_data_normalized_list$total_normalized_hits
  hit_totals <- reshape2::melt(metannotate_data_normalized_list$total_normalized_hits, id.vars = "Dataset",
                               variable.name = "HMM.Family", value.name = "percent_abundance")
  hit_totals$HMM.Family <- factor(hit_totals$HMM.Family, levels = unique(hit_totals$HMM.Family), ordered = TRUE)
  
  # Detect the taxonomy that the data has been collapsed to
  collapse_taxonomy_metannotate <- tail(taxonomy_naming$metannotate_colnames[
    taxonomy_naming$metannotate_colnames %in% colnames(metannotate_data)], n = 1)
  collapse_taxonomy <- taxonomy_naming$taxonomy[
    match(x = collapse_taxonomy_metannotate, table = taxonomy_naming$metannotate_colnames)]
  flog.debug(glue::glue("Plotting input dataframe has been collapsed to the '", collapse_taxonomy, "' level."))
  
  # Subset to the desired top_x cutoff
  metannotate_data <- subset_normalized_metannotate_data(metannotate_data, top_x, percent_mode = percent_mode)
  
  # Optionally subset to a specific taxon
  ########
  
  # Make or read in a plotting colour table; or generate auto-colours
  plotting_colours <- process_plotting_colours(metannotate_data, colouring_template_filename)
  
  # Make the plotting column into an ordered factor based on the plotting_colours order
  metannotate_data[,collapse_taxonomy_metannotate] <- factor(dplyr::pull(metannotate_data, collapse_taxonomy_metannotate),
                                                             levels = unique(dplyr::pull(plotting_colours, collapse_taxonomy_metannotate)),
                                                             ordered = TRUE)
  
  # Determine normalizing HMM for labelling on the plot
  normalizing_HMM <- colnames(hit_totals_wide)[
    match(100, unlist(lapply(2:ncol(hit_totals_wide), function(x) { mean(hit_totals_wide[,x]) }))) + 1 ]
  hit_totals <- dplyr::filter(hit_totals, HMM.Family != normalizing_HMM)
  
  flog.info("Creating the ggplot")
  metannotate_plot <- metannotate_ggplot(metannotate_data, type = type, plotting_colours, 
                                         collapse_taxonomy_metannotate, hit_totals, normalizing_HMM)
  
  return(metannotate_plot)
  
}


explore_metannotate_data <- function(metannotate_data_renamed, evalue = 1e-40, taxon = "Family",
                                     normalizing_HMM = "rpoB", plot_type = "bar", top_x = 0.02,
                                     colouring_template_filename = NA, percent_mode = "within_sample") {
  
  # Filter by e-value cutoff and report stats to user
  flog.info(glue::glue("Filtering by e-value cutoff of ", evalue))
  metannotate_data_filtered <- filter_by_evalue(metannotate_data_renamed, evalue = evalue)
  metannotate_data <- metannotate_data_filtered$metannotate_data
  flog.info(glue::glue("Percent change from e-value filtration:"))
  print(metannotate_data_filtered$read_counts$percent_change)
  # TODO - output the info to the user as a file
  
  # Collapse the table to the desired taxonomic rank
  flog.info(glue::glue("Collapsing table to taxonomic rank '", taxon, "'"))
  metannotate_data_collapsed <- collapse_metannotate_table_by_taxon(metannotate_data, taxon = taxon)
  
  # Normalize the data by HMM length
  flog.info("Normalizing data")
  metannotate_data_normalized_list <- normalize_collapsed_metannotate_data(metannotate_data_collapsed, 
                                                                           normalizing_HMM = normalizing_HMM)
  flog.info("Total normalized % abundance of analyzed genes compared to the marker gene:")
  print(metannotate_data_normalized_list$total_normalized_hits)
  
  # Make plots
  flog.info("Plotting data")
  metannotate_plot <- metannotate_plotter(metannotate_data_normalized_list, type = plot_type, top_x = top_x,
                                          colouring_template_filename = colouring_template_filename)
  
  return(metannotate_plot)
  
}

main <- function(params) {
  
  # Read data
  flog.info("Loading metannotate table (can take time)...")
  metannotate_data <- read_metannotate_tibble(params$metannotate_table_filename)
  
  # Print raw HMM and sample names in template if desired
  if (params$setup == TRUE) {
    setup_templates <- create_setup_templates(metannotate_data, write_tables = TRUE)
    flog.info(glue::glue("Wrote setup templates to 'hmm_info_template.tsv' and 'dataset_info_template.tsv. ",
                         "Please fill these out and then use them as script inputs. Exiting..."))
    quit(save = "no", status = 0)
  }
  
  # Apply user-provided setup templates if desired (required because of HMM length info)
  # TODO - make it possible to somehow run this in 'auto' mode to avoid the iterative usage for quick use.
  if (is.na(params$hmm_naming_info) == FALSE && is.na(params$dataset_naming_info) == FALSE) {
    flog.info("Loading user-provided HMM and dataset naming information")
    metannotate_data <- map_naming_information(metannotate_data, params$hmm_naming_info, params$dataset_naming_info)
  } else {
    flog.fatal(glue::glue("Must provde hmm_naming_info and dataset_naming_info. Exiting..."))
    quit(save = "no", status = 1)
  }
  
  metannotate_plot <- explore_metannotate_data(metannotate_data, evalue = params$evalue, taxon = params$taxon,
                         normalizing_HMM = params$normalizing_HMM, plot_type = params$plot_type, 
                         top_x = params$top_x, colouring_template_filename = params$colouring_template_filename,
                         percent_mode = params$percent_mode)
  print(metannotate_plot)
  ggsave(file = params$output_filename, width = params$plot_dimensions_mm[1], 
         height = params$plot_dimensions_mm[2], units = "mm")
  
}

if (interactive() == FALSE) {
  # TODO - Use a parser. Examples are left here for now.
  
  params <- list()
  setwd("/home/jmtsuji/Downloads/met_raw_reads_data_ELA111314/")
  params$metannotate_table_filename <- "all_annotations_O24C2t28808773.tsv.gz"
  params$setup <- FALSE
  params$hmm_naming_info <- "hmm_info_filled_vs2.tsv"
  params$dataset_naming_info <- "dataset_info_filled_vs1.tsv"
  params$evalue <- 1e-10
  params$taxon <- "Genus"
  params$normalizing_HMM <- "rpoB"
  params$top_x <- 0.01
  params$colouring_template_filename <- NA
  params$plot_type <- "bubble"
  params$percent_mode <- "within_HMM"
  params$output_filename <- "test_vs2.pdf"
  params$plot_dimensions_mm <- c(200,200)
  
  main(params)
  
} 

