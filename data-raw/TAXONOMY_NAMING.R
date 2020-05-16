# Table relating canonical taxonomy values to the colnames of the metannotate table
TAXONOMY_NAMING <- tibble::tibble(taxonomy = c("domain", "phylum", "class", "order", "family", "genus", "species"),
                                  metannotate_colnames = c("Closest.Homolog.Superkingdom", "Closest.Homolog.Phylum",
                                                           "Closest.Homolog.Class", "Closest.Homolog.Order",
                                                           "Closest.Homolog.Family", "Closest.Homolog.Genus",
                                                           "Closest.Homolog.Species"))

usethis::use_data(TAXONOMY_NAMING)
