# MetAnnotate barplot generator
# By Jackson M. Tsuji (Neufeld Lab PhD student)
# Under the BSD 3 License: https://opensource.org/licenses/BSD-3-Clause
# Copyright (c) 2017, Jackson M. Tsuji

#######################################
##### User variables ##################

setwd("/Users/")
input_filename <- "frag_0_all_annotations_bLx8GY742162434.tsv"

script_setup <- FALSE # Prints off raw HMM and sample names in template for setting up sample data files, then exits early.
                      # MUST run this the first time you use this script on a given dataset

### Required supplemental data (user needs to modify the template printed in script_setup)
dataset_info_filename <- "dataset_info_template_FILLED.tsv"  # Includes sample raw names and corrected names for plotting
                                                            # ***NOTE: order of datasets in this file dictates their order 
                                                            # in the final plot

hmm_info_filename <- "hmm_info_template_FILLED.tsv" # Includes HMM raw names, corrected names for plotting, and HMM lengths

### Basic plot settings
HMMs_to_plot <- c("pmoA", "dsrA", "cyc2-PV1")
normalizing_HMM <- "rpoB"
tax_rank_to_plot <- "Family"
top_number_to_plot <- 0.01  # If < 1, then plot all taxa with at least this relative abundance in the community.
                            # If > 1 (e.g., 10), then plot the top ___ (e.g., 10) taxa for each gene

### Options for output:
write_tables <- TRUE # Print off summary tables?
printPDF <- TRUE
output_name_general <- "171011_barplot" # general name of output files to append to


### Advanced features for making custom plots

# ** NOTE: the basic plot settings (above) MUST match those from when the template file was created!
print_custom_plot_template <- FALSE # Prints a template for the user to generate a plot with their own colours, then EXITS
print_custom_plot <- TRUE           # ONLY works if you have already run print_custom_plot_template and filled in the table

# Required for custom_plot_template_filename
custom_plot_template_filename <- "171011_barplot_05_custom_plot_template_Family_0.01_FILLED.tsv" 


#######################################
#######################################

##### Load libraries ######
library(plyr)
library(dplyr)
library(ggplot2)


#############################################
### Part 1: read/clean the data #############
#############################################

### Read data
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
### Part 3: collapse table to taxonomic rank, subset, and make auto plots ###
#############################################

# Functions
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
  
  return(collapsed)
} 


subset_collapsed_table <- function(collapsed_table, top_num) {
  # test vars
  # collapsed_table <- collapse_table(hits_all, "Family")
  # top_num <- 10
  
  # Determine sort method based on top_num provided
  if (top_num < 1) {
    sort_method <- "top_perc"
  } else if (top_num >= 1) {
    sort_method <- "top_num"
  } else {
    stop("ERROR: top_number_to_plot was somehow out of bounds. Must be greater than 0.")
  }
  
  # Split into individual list elements by sample/hmm to be able to sort effectively
  hits_by_sample <- split(collapsed_table, collapsed_table$Dataset)
  hits_by_hmm <- lapply(unique(collapsed_table$Dataset), function(dataset_name) { split(hits_by_sample[[dataset_name]], hits_by_sample[[dataset_name]]$HMM.Family) } )
  names(hits_by_hmm) <- unique(collapsed_table$Dataset)
  
  # Combine two functions together (nested) to sort within the nested loop
  # 1: main function; sorts for each individual data frame within the list
  sort_top_num1 <- function(hmm_name, dataset_table) {
    # test vars
    # dataset_table <- hits_by_hmm[[1]]
    # hmm_name <- "dsrA"
    
    hmm_table <- dataset_table[[hmm_name]]
    
    # Sort from most to least abundant
    hmm_table <- hmm_table[order(hmm_table$normalized_count_to_rpoB, decreasing = T),]
    
    if (sort_method == "top_num") {
      if (nrow(hmm_table) < top_num) {
        # If there are fewer entries than desired top #, just return the table as is
        hmm_table <- hmm_table
      } else {
        # Take top portion, because table is already sorted
        hmm_table <- hmm_table[1:top_num,]
      }
      
    } else if (sort_method == "top_perc") {
      
      hmm_table <- dplyr::filter(hmm_table, normalized_count_to_rpoB >= top_num)
      
    } else if (sort_method == "top_perc_alt") {
      # Make additional normalized count within sample. Can then sort based on rel. abund within sample...
      temp_sum <- sum(hmm_table$normalized_count_to_rpoB)
      hmm_table$normalized_count_within_sample_TEMP <- hmm_table$normalized_count_to_rpoB / temp_sum
      
      hmm_table <- dplyr::filter(hmm_table, normalized_count_within_sample_TEMP >= top_num)
      
      # Remove temp column
      hmm_table$normalized_count_within_sample_TEMP <- NULL
      
    } else {
      stop("Invalid sort method. Should be 'top_num' or 'top_perc'. Exiting...")
    }
    
    return(hmm_table)
  }
  
  # 2: helper function; guides function 1 to each list
  sort_top_num2 <- function(dataset_name) {
    dataset_table <- hits_by_hmm[[dataset_name]]
    dataset_table_sort <- lapply(names(dataset_table), sort_top_num1, dataset_table = dataset_table)
    # For working with a multi-parameter function in lapply, see https://stackoverflow.com/a/6827519 (accessed 170606)
    return(dataset_table_sort)
  }
  # Run helper function 2 using lapply (and function 2 calls function 1...)
  sorted <- lapply(names(hits_by_hmm), sort_top_num2)
  names(sorted) <- names(hits_by_hmm)
  
  # Put the separate data frames back together
  merge_df <- function(dataset_name) {
    merged_dataset <- dplyr::bind_rows(sorted[[dataset_name]])
    return(merged_dataset)
  }
  sorted_final <- dplyr::bind_rows(lapply(names(sorted), merge_df))
  
  return(sorted_final)
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
  normalizing_HMM_hits_summary_subset <- dplyr::filter(normalizing_HMM_hits_summary, HMM.Family %in% HMMs_to_plot)
  normalizing_HMM_hits_summary_subset$HMM.Family <- factor(normalizing_HMM_hits_summary_subset$HMM.Family, levels = HMMs_to_plot, ordered = TRUE)
  
  # Remove normalizing HMM from the gene totals, if still present; not helpful to have a grey box at 100% for that gene...
  normalizing_HMM_hits_summary_subset <- dplyr::filter(normalizing_HMM_hits_summary, HMM.Family != normalizing_HMM)
  
  # Convert proportions to percents for plot
  normalizing_HMM_hits_summary_subset$normalized_count_to_rpoB <- normalizing_HMM_hits_summary_subset$normalized_count_to_rpoB * 100
  subsetted_table$normalized_count_to_rpoB <- subsetted_table$normalized_count_to_rpoB * 100
  
  # Build plot
  bar_plot_overlay <- ggplot() +
    geom_bar(data = normalizing_HMM_hits_summary_subset, aes(x = Dataset, weight = normalized_count_to_rpoB), fill = "#808080") +
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
auto_barplot <- make_automated_barplot(hits_collapsed_subset)

print(auto_barplot)

if (printPDF == TRUE) {
  auto_barplot_filename <- paste(output_name_general, "_04_auto_barplot_", tax_rank_to_plot, "_", top_number_to_plot, ".pdf", sep = "")
  
  # Scale plot dimensions relative to input data size
  # Width: consider the number of datasets (x axis) and the longest legend entry
  auto_barplot_width <- length(unique(hits_collapsed_subset$Dataset)) * 30 + max(nchar(unique(hits_collapsed_subset[,3]))) * 3
  # Height: consider the number of HMMs (y panels) and the longest dataset name (x axis label)
  auto_barplot_height <- length(unique(hits_collapsed_subset$HMM.Family)) * 40 + max(nchar(unique(hits_collapsed_subset$Dataset))) * 2
  
  ggsave(auto_barplot_filename, width = auto_barplot_width, height = auto_barplot_height, units = c("mm"))
}

# 4. Print custom plotting template if desired, then exit
if (print_custom_plot_template == TRUE) {
  barplot_template <- make_barplot_template(hits_collapsed_subset)
  
  # Write to file
  plotting_template_table_filename <- paste(output_name_general, "_05_custom_plot_template_", tax_rank_to_plot, "_", top_number_to_plot, ".tsv", sep = "")
  write.table(barplot_template, file = plotting_template_table_filename, sep = "\t", row.names = F, col.names = T)
  
  stop(paste("Wrote plotting template table to ", plotting_template_table_filename, ".\nPlease fill out desired plotting colours AND put taxa (rows) in the order in which you want them to appear on the plot.\nThen, use this as input for the custom_plot_template variable and set print_custom_plot to TRUE.\nExiting...", sep = ""))
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
  normalizing_HMM_hits_summary_subset <- dplyr::filter(normalizing_HMM_hits_summary, HMM.Family %in% HMMs_to_plot)
  normalizing_HMM_hits_summary_subset$HMM.Family <- factor(normalizing_HMM_hits_summary_subset$HMM.Family, levels = HMMs_to_plot, ordered = TRUE)
  
  # Remove normalizing HMM from the gene totals, if still present; not helpful to have a grey box at 100% for that gene...
  normalizing_HMM_hits_summary_subset <- dplyr::filter(normalizing_HMM_hits_summary, HMM.Family != normalizing_HMM)
  
  # Convert proportions to percents for plot
  normalizing_HMM_hits_summary_subset$normalized_count_to_rpoB <- normalizing_HMM_hits_summary_subset$normalized_count_to_rpoB * 100
  subsetted_table$normalized_count_to_rpoB <- subsetted_table$normalized_count_to_rpoB * 100
  
  # Make new plot
  custom_barplot_overlay <- ggplot() +
    geom_bar(data = normalizing_HMM_hits_summary_subset, aes(x = Dataset, weight = normalized_count_to_rpoB), fill = "#808080") +
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
    custom_barplot_filename <- paste(output_name_general, "_06_custom_barplot_", tax_rank_to_plot, "_", top_number_to_plot, ".pdf", sep = "")
    
    # Scale plot dimensions relative to input data size
    # Width: consider the number of datasets (x axis) and the longest legend entry
    custom_barplot_width <- length(unique(hits_collapsed_subset$Dataset)) * 30 + max(nchar(unique(hits_collapsed_subset[,3]))) * 3
    # Height: consider the number of HMMs (y panels) and the longest dataset name (x axis label)
    custom_barplot_height <- length(unique(hits_collapsed_subset$HMM.Family)) * 40 + max(nchar(unique(hits_collapsed_subset$Dataset))) * 2
    
    ggsave(custom_barplot_filename, width = custom_barplot_width, height = custom_barplot_height, units = c("mm"))
  }
}

