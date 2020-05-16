# colours.R
# colour manipulation for metannotate data
# Copyright Jackson M. Tsuji, 2020 (Neufeld Lab)

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
#' @param metannotate_data_normalized Tibble of normalized metannotate data - see
#' \code{\link{normalize_collapsed_metannotate_data}}
#' @return a tibble of unique taxa with HTML colour codes in the 'colour' column
#' @keywords internal
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
#' @description A wrapper for \code{\link{generate_plotting_colours}} to handle bigger-picture decision making
#' @param metannotate_data_normalized Tibble of normalized metannotate data - see
#' \code{\link{normalize_collapsed_metannotate_data}}
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
