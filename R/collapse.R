# collapse.R
# Collapse metannotate data
# Copyright Jackson M. Tsuji, 2020 (Neufeld Lab)

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
