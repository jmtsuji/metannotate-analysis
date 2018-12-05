# MetAnnotate barplot generator
# By Jackson M. Tsuji (Neufeld Lab PhD student)
# Under the BSD 3 License: https://opensource.org/licenses/BSD-3-Clause
# Copyright (c) 2017, Jackson M. Tsuji

##### Load libraries ######
library(plyr)
library(dplyr)
library(ggplot2)
library(glue)

#######################################
##### User variables ##################

### Global variables
setwd("/home/jmtsuji/Research_General/PhD/04e_Biogeography/05_JGI_metagenome_analysis/02b_171116_prelim_2/04_results_phase3")
write_tables <- TRUE # Print off summary tables?
printPDF <- TRUE # Print PDF plots to the folder?
PDF_dimension_scaling <- c(1,1.2) # If you want to adjust the relative c(width, height) of the plots compared to the defaults
                                # Default is c(1,1)
output_name_general <- "02a_Fe_S/180110_vs1" # general name of output files to append to

### Inputs
input_filename <- "00b_script_input/all_annotations_phase1and2and3_COMB.tsv"

script_setup <- FALSE # Prints off raw HMM and sample names in template for setting up sample data files, then exits early.
                      # MUST run this the first time you use this script on a given dataset

### Required supplemental data (user needs to modify the template printed in script_setup)
dataset_info_filename <- "00b_script_input/dataset_info_template_FILLED.tsv"  # Includes sample raw names and corrected names for plotting
                                                            # ***NOTE: order of datasets in this file dictates their order 
                                                            # in the final plot

hmm_info_filename <- "00b_script_input/hmm_info_template_FILLED.tsv" # Includes HMM raw names, corrected names for plotting, and HMM lengths

### Basic plot settings
HMMs_to_plot <- c("cyc2-PV1-GSB", "dsrA")
datasets_to_plot <- "all" # or give the names of the datasets in the order you want them plotted
normalizing_HMM <- "rpoB"
tax_rank_to_plot <- "Family"
top_number_to_plot <- 0.005  # If < 1, then plot all taxa with at least this relative abundance in the community.
                            # If > 1 (e.g., 10), then plot the top ___ (e.g., 10) taxa for each gene
percent_sort_method <- "by_dataset" # If top_number_to_plot is < 1, you need to provide guidance for which type of percentage-based sorting to use. Options are either:
                                            # If "by_dataset", gives the taxa above x% relative to the normalizing_HMM for each dataset.
                                            # If "by_HMM", gives the taxa above x% abundance relative to the total hits for each specific HMM.
                                            # If you're sorting by top x taxa (i.e., top_number_to_plot > 1), then it doesn't matter what this variable is set to.


### Advanced feature: making custom plots with user-provided taxa order and colours
# ** NOTE: the basic plot settings (above) MUST match those from when the template file was created!
print_custom_plot <- FALSE          # MUST run this script first to generate custom plot template file, then load in the template below.
custom_plot_template_filename <- "171011_barplot_05_custom_plot_template_Family_0.01_FILLED.tsv" 


### Advanced feature: make a subset plot of specific taxa (names within the tax_rank_to_plot rank)
make_taxa_subset <- TRUE
taxa_names_to_subset <- c("Chlorobiaceae")

#######################################
#######################################

#############################################
### Part 1: read/clean the data #############
#############################################

### Read data
print("Loading data file...")
hits_all <- read.table(input_filename, sep = "\t", header = TRUE, comment.char = "", quote = "", stringsAsFactors = FALSE)

### Print raw HMM and sample names in template if desired, to help build sample info file:
if (script_setup == TRUE) {
  
  # Make template for HMM naming
  hmm_info_cols <- c("HMM.Family", "raw_name", "HMM_length", "notes")
  hmm_info_raw <- lapply(1:length(hmm_info_cols), function(num) { character(length = length(unique(hits_all$HMM.Family))) } )
  names(hmm_info_raw) <- hmm_info_cols
  hmm_info_template <- dplyr::bind_rows(hmm_info_raw)
  hmm_info_template$raw_name <- unique(hits_all$HMM.Family)
  hmm_info_template_filename <- "hmm_info_template.tsv"
  write.table(hmm_info_template, file = hmm_info_template_filename, sep = "\t", row.names = F, col.names = T)
  
  # Make template for sample naming
  dataset_info_cols <- c("Dataset", "raw_name")
  dataset_info_raw <- lapply(1:length(dataset_info_cols), function(num) { character(length = length(unique(hits_all$Dataset))) } )
  names(dataset_info_raw) <- dataset_info_cols
  dataset_info_template <- dplyr::bind_rows(dataset_info_raw)
  dataset_info_template$raw_name <- unique(hits_all$Dataset)
  dataset_info_template_filename <- "dataset_info_template.tsv"
  write.table(dataset_info_template, file = dataset_info_template_filename, sep = "\t", row.names = F, col.names = T)
  
  stop(paste("Wrote raw names of HMMs and Datasets to files ", hmm_info_template_filename, " and ", dataset_info_template_filename, ".\nFill out these templates and then use for the supplemental info file required to run the script further.\nExiting... set script_setup = FALSE to run rest of script when ready.", sep = ""), call. = FALSE)
}

### Fix HMM names
hmm_info <- read.table(hmm_info_filename, sep = "\t", header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
hits_all$HMM.Family <- mapvalues(hits_all$HMM.Family, hmm_info$raw_name, as.character(hmm_info$HMM.Family))

### Fix sample names
dataset_info <- read.table(dataset_info_filename, sep = "\t", header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
hits_all$Dataset <- mapvalues(hits_all$Dataset, dataset_info$raw_name, as.character(dataset_info$Dataset))

### Map HMM lengths and gene category onto table
hmm_info_merge <- hmm_info[,c("HMM.Family", "HMM_length")]
hits_all <- dplyr::left_join(hits_all, hmm_info_merge, by = "HMM.Family")
hmm_info_merge <- NULL

## Subset HMMs of interest and the HMM to normalize to, to save memory for the rest of the script
hits_all <- dplyr::filter(hits_all, HMM.Family %in% c(HMMs_to_plot, normalizing_HMM))


#############################################
### Part 2: normalize data ##################
#############################################
print("Normalizing data...")

### Create normalized column by HMM length (between-HMM normalization)
hits_all$normalized_by_hmm <- (1 / hits_all$HMM_length)
# Can now just sum this column and get the total number of a given group
# This ONLY works for within sample comparisons. For between sample comparison, 
# have to also normalize by rpoB (single-copy taxonomic marker) hits and show as a 
# relative abundance. The current normalization acts kind of like a non-rarefied OTU table.


### Also normalize by rpoB hits for relative abundance comparisons (between samples)
# Sum rpoB (or other normalizing HMM) hit counts pre-normalized by HMM_length for each sample and map to 
# data frame, as a divider tool
normalizing_HMM_subset <- dplyr::filter(hits_all, HMM.Family == normalizing_HMM)
normalizing_HMM_subset_grouped <- dplyr::group_by(normalizing_HMM_subset, Dataset)
normalizing_HMM_subset_summarized <- dplyr::summarise(normalizing_HMM_subset_grouped, 
                                                      normalized_rpob_hits = sum(normalized_by_hmm))
normalizing_HMM_subset <- NULL
normalizing_HMM_subset_grouped <- NULL

# Map to main data frame
hits_all <- dplyr::left_join(hits_all, normalizing_HMM_subset_summarized, by = "Dataset")

# Get normalized hit counts for each sample to rpoB and HMM_length
hits_all$normalized_count_to_rpoB <- hits_all$normalized_by_hmm / hits_all$normalized_rpob_hits

normalizing_HMM_hits_summary <- dplyr::summarise(dplyr::group_by(hits_all, Dataset, HMM.Family), normalized_count_to_rpoB = sum(normalized_count_to_rpoB))  

# Check out the relative importance of all genes, if interested...
if (write_tables == TRUE) {
  normalizing_HMM_hits_summary_filename <- paste(output_name_general, "_01_total_normalized_hits", ".tsv", sep = "")
  write.table(normalizing_HMM_hits_summary, file = normalizing_HMM_hits_summary_filename, sep = "\t", col.names = T, row.names = F)
}

#############################################
### Part 3: Functions for step 4 ###
#############################################

# Functions
make_settings_log <- function() {
  # Relies on global variables set to create a log for the run
  # using cat instead of print to avoid index numbers: see https://stackoverflow.com/a/19943171 (accessed 171117)
  
  cat("Logfile for metannnotate_barplots.R\n")
  cat(paste("Script run on ", date(), " in directory ", getwd(), "\n", sep = ""))
  cat("\n")
  cat(paste("Input file: ", input_filename, "\n", sep = ""))
  cat(paste("Dataset template file: ", dataset_info_filename, "\n", sep = ""))
  cat(paste("HMM template file: ", hmm_info_filename, "\n", sep = ""))
  cat(paste("General output name: ", output_name_general, "\n", sep = ""))
  cat("\n")
  cat(paste("Requested HMMs in plot: ", glue::collapse(HMMs_to_plot, sep = ", "), "\n", sep = ""))
  cat(paste("Normalizing counts to HMM: ", normalizing_HMM, "\n", sep = ""))
  cat(paste("Summarizing at rank: ", tax_rank_to_plot, "\n", sep = ""))
  if (top_number_to_plot >= 1) {
    cat(paste("Abundance filter: plotting top ", top_number_to_plot, " taxa", "\n", sep = ""))
  } else if (top_number_to_plot < 1) {
    if (percent_sort_method == "by_dataset") {
      cat(paste("Abundance filter: plotting taxa over ", top_number_to_plot, "% abundance relative to ", normalizing_HMM, "\n", sep = ""))
    } else if (percent_sort_method == "by_HMM") {
      cat(paste("Abundance filter: plotting taxa over ", top_number_to_plot, "% abundance within each gene", "\n", sep = ""))
    }
  }
  cat("\n")
  if (write_tables == TRUE) {
    cat("write_tables = TRUE; will print output tables\n")
  } else if (write_tables == FALSE) {
    cat("write_tables = FALSE; will not print output tables\n")
  }
  if (printPDF == TRUE) {
    cat(paste("printPDF = TRUE; will print output PDFs with width/height scaling of ", PDF_dimension_scaling[1], "/", PDF_dimension_scaling[2], "\n", sep = ""))
  } else if (printPDF == FALSE) {
    cat("printPDF = FALSE; will not print output PDFs\n")
  }
  cat("\n")
  
  if (print_custom_plot == TRUE) {
    cat(paste("print_custom_plot = TRUE; will create custom plot from provided template ", custom_plot_template_filename, "\n", sep = ""))
  } else if (print_custom_plot == FALSE) {
    cat("print_custom_plot = FALSE; will not create custom plot\n")
  }
  cat("\n")
  
  if (make_taxa_subset == TRUE) {
    cat(paste("make_taxa_subset = TRUE; will create custom plot subsetting these (rank ", tax_rank_to_plot, ") taxa: ", glue::collapse(taxa_names_to_subset, sep = ", "), "\n", sep = ""))
  } else if (make_taxa_subset == FALSE) {
    cat("make_taxa_subset = FALSE; will not create subset plot")
  }
  cat("\n")

}


collapse_table <- function(table, tax_rank) {
  # Test vars
  # table <- hits_all
  # tax_rank <- "Family"
  
  # Find taxonomic rank in table
  tax_rank_name <- paste("Closest.Homolog.", tax_rank, sep = "")
  if (is.na(match(tax_rank_name, colnames(table)))) {
    stop("ERROR: Something is not right with your tax_rank_to_plot. Should be a standard taxonomic rank in the MetAnnotate output table. Exiting...")
  } else {
    tax_rank_col_num <- match(tax_rank_name, colnames(table))
  }
  
  # Collpase the table by the given taxonomic rank
  grouped <- dplyr::group_by(table, Dataset, HMM.Family, table[,tax_rank_col_num])
  collapsed <- dplyr::summarise(grouped, normalized_count_to_rpoB = sum(normalized_count_to_rpoB))
  colnames(collapsed)[3] <- colnames(table)[tax_rank_col_num]
  
  # Add in higher taxonomic rankings data
  ranks <- hits_all[!duplicated(hits_all[,tax_rank_col_num]),tax_rank_col_num:14]
  
  # From http://stackoverflow.com/a/9945116, accessed 170312
  collapsed <- as.data.frame(dplyr::left_join(collapsed, ranks, by = colnames(table)[tax_rank_col_num]))
  
  # Sort by rpoB count for ease of reading
  collapsed <- dplyr::arrange(collapsed, Dataset, HMM.Family, desc(normalized_count_to_rpoB))
  
  return(collapsed)
} 


subset_collapsed_table <- function(collapsed_table, top_num) {
  # test vars
  # collapsed_table <- collapse_table(hits_all, "Family")
  # top_num <- 10
  
  # Determine sort method based on top_num provided
  if (top_num >= 1) {
    sort_method <- "top_num"
  } else if (top_num < 1 && percent_sort_method == "by_dataset") {
    sort_method <- "top_perc_by_dataset"
  } else if (top_num < 1 && percent_sort_method == "by_HMM") {
    sort_method <- "top_perc_by_HMM"
  } else {
    stop("ERROR: top_number_to_plot or percent_sort_method was somehow out of bounds. top_number_to_plot must be greater than 0. If <1, percent_sort_method must be 'by_dataset' or 'by_HMM'.")
  }
  
  if (sort_method == "top_num") {
    # Pick out the most abundant taxa
    collapsed_table_grouped <- dplyr::group_by(collapsed_table, HMM.Family, Dataset)
    collapsed_table_top_values <- dplyr::top_n(collapsed_table_grouped, top_num, wt = normalized_count_to_rpoB)
  } else if (sort_method == "top_perc_by_dataset") {
    # Subset all taxa above the desired percent abundance threshold (relative to rpoB)
    collapsed_table_top_values <- dplyr::filter(collapsed_table, normalized_count_to_rpoB >= top_num)
  } else if (sort_method == "top_perc_by_HMM") {
    # First, get the total sum of HMM hits for each metagenome dataset, per HMM
    collapsed_table_grouped <- dplyr::group_by(collapsed_table, HMM.Family, Dataset)
    collapsed_table_summed_counts <- dplyr::summarise(collapsed_table_grouped, normalized_count_total = sum(normalized_count_to_rpoB))
    
    # Join this to the main data frame
    collapsed_table_joined <- dplyr::left_join(collapsed_table, collapsed_table_summed_counts, by = c("Dataset", "HMM.Family"))
    
    # Use this to determine the relative abundance of each hit within the sample
    collapsed_table_joined$per_sample_rel_abund <- collapsed_table_joined$normalized_count_to_rpoB / collapsed_table_joined$normalized_count_total
    
    # Sort by within-sample relative abundance
    collapsed_table_top_values <- dplyr::filter(collapsed_table_joined, per_sample_rel_abund >= top_num)
    
    # Delete temporary columns and temporary table
    collapsed_table_top_values[,c("normalized_count_total", "per_sample_rel_abund")] <- NULL
    collapsed_table_joined <- NULL
  }
  
  return(collapsed_table_top_values)
}


make_automated_barplot <- function(subsetted_table) {
  # test vars
  # subsetted_table <- subset_collapsed_table(collapse_table(hits_all, "Family"), 0.01)
  
  # Assign HMM levels
  subsetted_table$HMM.Family <- factor(subsetted_table$HMM.Family, levels = HMMs_to_plot, ordered = TRUE)
  
  # Assign Dataset levels based on dataset info filename
  subsetted_table$Dataset <- factor(subsetted_table$Dataset, levels = dataset_info$Dataset, ordered = TRUE)
  
  # Get gene totals and filter (and assign factors) similarly
  normalizing_HMM_hits_summary <- dplyr::summarise(dplyr::group_by(hits_all, Dataset, HMM.Family), normalized_count_to_rpoB = sum(normalized_count_to_rpoB))  
   
  # Similarly assign HMM and Dataset levels
  normalizing_HMM_hits_summary$Dataset <- factor(normalizing_HMM_hits_summary$Dataset, levels = dataset_info$Dataset, ordered = TRUE)
  normalizing_HMM_hits_summary$HMM.Family <- factor(normalizing_HMM_hits_summary$HMM.Family, levels = HMMs_to_plot, ordered = TRUE)
  
  # Remove normalizing HMM from the gene totals, if still present; not helpful to have a grey box at 100% for that gene...
  normalizing_HMM_hits_summary <- dplyr::filter(normalizing_HMM_hits_summary, HMM.Family != normalizing_HMM)
  
  # Convert proportions to percents for plot
  normalizing_HMM_hits_summary$normalized_count_to_rpoB <- normalizing_HMM_hits_summary$normalized_count_to_rpoB * 100
  subsetted_table$normalized_count_to_rpoB <- subsetted_table$normalized_count_to_rpoB * 100
  
  # Build plot
  bar_plot_overlay <- ggplot() +
    geom_bar(data = normalizing_HMM_hits_summary, aes(x = Dataset, weight = normalized_count_to_rpoB), fill = "#808080") +
    geom_bar(data = subsetted_table, aes_string(x = "Dataset", weight = colnames(subsetted_table)[4], 
                                            fill = colnames(subsetted_table)[3]), position = "stack") +
    facet_grid(HMM.Family ~ ., scales = "free") +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.title = element_text(size = 12), 
          strip.text = element_text(size = 11, face = "italic"), strip.background = element_rect(fill = "#e6e6e6"),
          panel.border = element_rect(colour = "transparent"), panel.spacing.y = unit(3, "mm"),
          axis.text = element_text(size = 10, colour = "black"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.ticks = element_line(size = 0.5, colour = "black"), axis.line = element_line(size = 0.5),
          legend.text = element_text(size = 10, face = "italic"), legend.title = element_blank(),
          legend.key = element_rect(colour = "transparent"), legend.key.size = unit(6, "mm"),
          legend.spacing = unit(1, "mm"), legend.box.just = "left", plot.margin = unit(c(2,2,2,2), "mm")) +
    guides(fill = guide_legend(ncol = 1)) +
    xlab("Sample") +
    ylab("Gene hits relative to rpoB (%; normalized)")
  
  # aes_strings idea from http://stackoverflow.com/q/15458526 (accessed 170311)
  
  # print(bar_plot_overlay)
  
  return(bar_plot_overlay)
}


make_barplot_template <- function(subsetted_table) {
  # test vars
  # subsetted_table <- subset_collapsed_table(collapse_table(hits_all, "Family"), 0.01)
  
  # Add on extra columns for desired info 
  plotting_template_table_TEMP <- data.frame("HTML_colour_code" = character(length = nrow(subsetted_table)), 
                                             "Notes" = character(length = nrow(subsetted_table)), stringsAsFactors = FALSE)
  plotting_template_table <- dplyr::bind_cols(subsetted_table, plotting_template_table_TEMP)
  
  # Remove extraneous info and reduce to unique elements in main plotting column
  cols_to_remove <- match(c("Dataset", "HMM.Family", "normalized_count_to_rpoB"), colnames(plotting_template_table))
  plotting_template_table <- plotting_template_table[,-cols_to_remove]
  plotting_template_table <- plotting_template_table[!duplicated(plotting_template_table[,1]),]
  
  return(plotting_template_table)
}


subset_barplot <- function(plotting_table, normalized_hits_table, hmms_to_plot, taxa_to_plot) {
  # Function to highlight specific community members of interest ONLY
  
  # # test vars for developing function
  # plotting_table <- hits_collapsed
  # normalized_hits_table <- normalizing_HMM_hits_summary
  # hmms_to_plot <- HMMs_to_plot
  # taxa_to_plot <- taxa_subset
  
  # Plotting table must be COLLAPSED but can be either before OR after subset step
  # *** REQUIRES global variable 'normalizing_HMM'
  
  plotting_taxon_rank <- colnames(plotting_table)[3]
  
  # Filter input tables
  plotting_table_filtered <- dplyr::filter(plotting_table, HMM.Family %in% hmms_to_plot , get(plotting_taxon_rank) %in% taxa_to_plot) # See https://stackoverflow.com/a/34220245, accessed 171117
  # Also remove normalizing HMM if present
  normalized_hits_filtered <- dplyr::filter(normalized_hits_table, HMM.Family %in% hmms_to_plot & HMM.Family != normalizing_HMM)
  colnames(normalized_hits_filtered)[3] <- "normalized_TOTAL_count_to_rpoB"
  
  # Switch from proportion to percent
  normalized_hits_filtered$normalized_TOTAL_count_to_rpoB <- normalized_hits_filtered$normalized_TOTAL_count_to_rpoB * 100
  plotting_table_filtered$normalized_count_to_rpoB <- plotting_table_filtered$normalized_count_to_rpoB * 100
  
  # Assign HMM and Dataset levels
  normalized_hits_filtered$Dataset <- factor(normalized_hits_filtered$Dataset, levels = dataset_info$Dataset, ordered = TRUE)
  normalized_hits_filtered$HMM.Family <- factor(normalized_hits_filtered$HMM.Family, levels = HMMs_to_plot, ordered = TRUE)
  plotting_table_filtered$Dataset <- factor(plotting_table_filtered$Dataset, levels = dataset_info$Dataset, ordered = TRUE)
  plotting_table_filtered$HMM.Family <- factor(plotting_table_filtered$HMM.Family, levels = HMMs_to_plot, ordered = TRUE)
  
  # plotting_table_filtered <- dplyr::left_join(plotting_table_filtered, normalized_hits_filtered, by = c("Dataset", "HMM.Family"))
  # plotting_table_filtered$rel_abund <- plotting_table_filtered$normalized_count_to_rpoB / plotting_table_filtered$normalized_TOTAL_count_to_rpoB
  
  # # Decide on colours
  # plotting_colours <- c(hue_pal()(length(taxa)), "#808080")
  
  # Make plot
  bar_plot_subset <- ggplot() +
    geom_bar(data = normalized_hits_filtered, aes(x = Dataset, weight = normalized_TOTAL_count_to_rpoB), fill = "#808080") +
    geom_bar(data = plotting_table_filtered, aes_string(x = "Dataset", weight = colnames(plotting_table_filtered)[4], 
                                                        fill = colnames(plotting_table_filtered)[3]), position = "stack") +
    facet_grid(HMM.Family ~ ., scales = "fixed") +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.title = element_text(size = 12), 
          strip.text = element_text(size = 11, face = "italic"), strip.background = element_rect(fill = "#e6e6e6"),
          panel.border = element_rect(colour = "transparent"), panel.spacing.y = unit(3, "mm"),
          axis.text = element_text(size = 10, colour = "black"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.ticks = element_line(size = 0.5, colour = "black"), axis.line = element_line(size = 0.5),
          legend.text = element_text(size = 10, face = "italic"), legend.title = element_blank(),
          legend.key = element_rect(colour = "transparent"), legend.key.size = unit(6, "mm"),
          legend.spacing = unit(1, "mm"), legend.box.just = "left", plot.margin = unit(c(2,2,2,2), "mm")) +
    guides(fill = guide_legend(ncol = 1)) +
    xlab("Sample") +
    ylab("Gene hits relative to rpoB (%; normalized)")
  
  return(bar_plot_subset)
}


#############################################
### Part 4: collapse table to taxonomic rank, subset, and make auto plots ###
#############################################
print("Collapsing to desired rank...")

# Subset to HMMs of interest
hits_all <- dplyr::filter(hits_all, HMM.Family %in% HMMs_to_plot)

# Run functions
# 1. Collapse the table to the given taxonomic rank. Write table if desired.
hits_collapsed <- collapse_table(hits_all, tax_rank_to_plot)
if (write_tables == TRUE) {
  hits_collapsed_filename <- paste(output_name_general, "_02_collapsed_to_", tax_rank_to_plot, ".tsv", sep = "")
  write.table(hits_collapsed, file = hits_collapsed_filename, sep = "\t", col.names = T, row.names = F)
}

# 2. Subset the top x entries for each HMM/dataset. Write table if desired.
hits_collapsed_subset <- subset_collapsed_table(hits_collapsed, top_number_to_plot)

if (write_tables == TRUE) {
  hits_collapsed_subset_filename <- paste(output_name_general, "_03_plotting_table_", tax_rank_to_plot, "_", top_number_to_plot, ".tsv", sep = "")
  write.table(hits_collapsed_subset, file = hits_collapsed_subset_filename, sep = "\t", col.names = T, row.names = F)
}

# 3. Generate an automated barplot of the subsetted table. Print if desired
print("Making plots...")
auto_barplot <- make_automated_barplot(hits_collapsed_subset)

print(auto_barplot)

if (printPDF == TRUE) {
  auto_barplot_filename <- paste(output_name_general, "_04_auto_barplot_", tax_rank_to_plot, "_", top_number_to_plot, ".pdf", sep = "")
  
  # Scale plot dimensions relative to input data size
  # Width: consider the number of datasets (x axis) and the longest legend entry
  auto_barplot_width <- (length(unique(hits_collapsed_subset$Dataset)) * 13 + max(nchar(unique(hits_collapsed_subset[,3]))) * 2) * PDF_dimension_scaling[1]
  # Height: consider the number of HMMs (y panels) and the longest dataset name (x axis label)
  auto_barplot_height <- (length(unique(hits_collapsed_subset$HMM.Family)) * 40 + max(nchar(unique(hits_collapsed_subset$Dataset))) * 2) * PDF_dimension_scaling[2]
  
  ggsave(auto_barplot_filename, width = auto_barplot_width, height = auto_barplot_height, units = c("mm"))
}

# 4. Print custom plotting template
barplot_template <- make_barplot_template(hits_collapsed_subset)
plotting_template_table_filename <- paste(output_name_general, "_05_custom_plot_template_", tax_rank_to_plot, "_", top_number_to_plot, ".tsv", sep = "")
write.table(barplot_template, file = plotting_template_table_filename, sep = "\t", row.names = F, col.names = T)

# 5. Optionally make subset barplot if desired
if (make_taxa_subset == TRUE) {
  subset_plot <- subset_barplot(hits_collapsed, normalizing_HMM_hits_summary, HMMs_to_plot, taxa_names_to_subset)
  
  if (printPDF == TRUE) {
    sub_barplot_filename <- paste(output_name_general, "_06_subset_barplot_", tax_rank_to_plot, "_", top_number_to_plot, ".pdf", sep = "")
    
    # Scale plot dimensions relative to input data size
    # Width: consider the number of datasets (x axis) and the longest legend entry
    sub_barplot_width <- (length(unique(hits_collapsed_subset$Dataset)) * 13 + max(nchar(taxa_names_to_subset)) * 2) * PDF_dimension_scaling[1]
    # Height: consider the number of HMMs (y panels) and the longest dataset name (x axis label)
    sub_barplot_height <- (length(HMMs_to_plot) * 40 + max(nchar(unique(hits_collapsed_subset$Dataset))) * 2) * PDF_dimension_scaling[2]
    
    ggsave(sub_barplot_filename, width = sub_barplot_width, height = sub_barplot_height, units = c("mm"))
  }
  
}


#############################################
### Part 4: Make custom features plots (advanced) ###
#############################################

### Function to make plot with custom features
make_custom_plot <- function(subsetted_table, plot_template) {
  # test vars
  # subsetted_table <- subset_collapsed_table(collapse_table(hits_all, "Family"), 0.01)
  # plot_template <- custom_plot_template
  
  # Apply taxon order to plotting table
  subsetted_table[,3] <- factor(subsetted_table[,3], levels = plot_template[,1], ordered = TRUE)
  
  # Assign HMM levels
  subsetted_table$HMM.Family <- factor(subsetted_table$HMM.Family, levels = HMMs_to_plot, ordered = TRUE)
  
  # Assign Dataset levels based on dataset info filename
  subsetted_table$Dataset <- factor(subsetted_table$Dataset, levels = dataset_info$Dataset, ordered = TRUE)
  
  # Get gene totals and filter (and assign factors) similarly
  normalizing_HMM_hits_summary <- dplyr::summarise(dplyr::group_by(hits_all, Dataset, HMM.Family), normalized_count_to_rpoB = sum(normalized_count_to_rpoB))  
  
  # Similarly assign HMM and Dataset levels
  normalizing_HMM_hits_summary$Dataset <- factor(normalizing_HMM_hits_summary$Dataset, levels = dataset_info$Dataset, ordered = TRUE)
  normalizing_HMM_hits_summary$HMM.Family <- factor(normalizing_HMM_hits_summary$HMM.Family, levels = HMMs_to_plot, ordered = TRUE)
  
  # Remove normalizing HMM from the gene totals, if still present; not helpful to have a grey box at 100% for that gene...
  normalizing_HMM_hits_summary <- dplyr::filter(normalizing_HMM_hits_summary, HMM.Family != normalizing_HMM)
  
  # Convert proportions to percents for plot
  normalizing_HMM_hits_summary$normalized_count_to_rpoB <- normalizing_HMM_hits_summary$normalized_count_to_rpoB * 100
  subsetted_table$normalized_count_to_rpoB <- subsetted_table$normalized_count_to_rpoB * 100
  
  # Make new plot
  custom_barplot_overlay <- ggplot() +
    geom_bar(data = normalizing_HMM_hits_summary, aes(x = Dataset, weight = normalized_count_to_rpoB), fill = "#808080") +
    geom_bar(data = subsetted_table, aes_string(x = "Dataset", weight = colnames(subsetted_table)[4], 
                                                fill = colnames(subsetted_table)[3]), position = "stack") +
    facet_grid(HMM.Family ~ ., scales = "free") +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.title = element_text(size = 12), 
          strip.text = element_text(size = 11, face = "italic"), strip.background = element_rect(fill = "#e6e6e6"),
          panel.border = element_rect(colour = "transparent"), panel.spacing.y = unit(3, "mm"),
          axis.text = element_text(size = 10, colour = "black"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.ticks = element_line(size = 0.5, colour = "black"), axis.line = element_line(size = 0.5),
          legend.text = element_text(size = 10, face = "italic"), legend.title = element_blank(),
          legend.key = element_rect(colour = "transparent"), legend.key.size = unit(6, "mm"),
          legend.spacing = unit(1, "mm"), legend.box.just = "left", plot.margin = unit(c(2,2,2,2), "mm")) +
    guides(fill = guide_legend(ncol = 1)) +
    scale_fill_manual(values = plot_template$HTML_colour_code) +
    xlab("Sample") +
    ylab("Gene hits relative to rpoB (%; normalized)")
  
  return(custom_barplot_overlay)
}

if (print_custom_plot == TRUE) {
  # Read in the filled in template
  custom_plot_template <- read.table(custom_plot_template_filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
  
  # Check if the internally generated table appears to match the provided template file (user might have changed the settings after generating the template, in which case the mapping will have problems)
  barplot_template_test <- make_barplot_template(hits_collapsed_subset)
  if (length(setdiff(barplot_template_test[,1], custom_plot_template[,1])) != 0) {
    stop("ERROR: based on first column, plotting template does not appear to match the current basic plot settings. Please check the settings and try again. Exiting...")
  }
  
  # Check the style of HTML colour code
  for (i in 1:nrow(custom_plot_template)) {
    check <- custom_plot_template$HTML_colour_code[i]
    if (substr(check, start = 1, stop = 1) != "#") {
      stop("Something appears to be wrong with your HTML colour codes provided with the plotting template. Should all start with hash symbols ('#'), e.g., #FFFFFF. Exiting...")
    } else if (nchar(check) != 7) {
      stop("Something appears to be wrong with your HTML colour codes provided with the plotting template. Should all start with hash symbols ('#') and then have six letters/numbers, e.g., #FFFFFF. Exiting...")
    }
  }
  
  # Run function
  custom_barplot <- make_custom_plot(hits_collapsed_subset, custom_plot_template)
  
  print(custom_barplot)
  if (printPDF == TRUE) {
    custom_barplot_filename <- paste(output_name_general, "_07_custom_barplot_", tax_rank_to_plot, "_", top_number_to_plot, ".pdf", sep = "")
    
        # Scale plot dimensions relative to input data size
    # Width: consider the number of datasets (x axis) and the longest legend entry
    custom_barplot_width <- (length(unique(hits_collapsed_subset$Dataset)) * 13 + max(nchar(unique(hits_collapsed_subset[,3]))) * 2) * PDF_dimension_scaling[1]
    # Height: consider the number of HMMs (y panels) and the longest dataset name (x axis label)
    custom_barplot_height <- (length(unique(hits_collapsed_subset$HMM.Family)) * 40 + max(nchar(unique(hits_collapsed_subset$Dataset))) * 2)  * PDF_dimension_scaling[2]
    
    ggsave(custom_barplot_filename, width = custom_barplot_width, height = custom_barplot_height, units = c("mm"))
  }
}

### Lastly, if run was successful, make a logfile of settings used:
custom_barplot_filename <- paste(output_name_general, "_00_LOG_", tax_rank_to_plot, "_", top_number_to_plot, ".log", sep = "")
sink(custom_barplot_filename, type = c("output"))
make_settings_log()
sink()

print("metannotate_barplots.R: finished running successfully.")
