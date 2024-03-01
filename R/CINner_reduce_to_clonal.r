CINner_reduce_to_clonal <- function(simulation) {
    #-----------------------------------Input current mutation & CN data
    segmented_copynumber <- simulation$sample$segmented_copynumber
    mutation_bulk <- simulation$sample_phylogeny$package_cell_phylogeny$mutation_bulk
    #------------------Remove all regions with subclonal CN from CN data
    vec_delete <- which(segmented_copynumber$Clonal_frequency < 1)
    if (length(vec_delete) > 0) segmented_copynumber <- segmented_copynumber[-vec_delete, ]
    #---------------------------------------------Store the updated data
    simulation$sample$segmented_copynumber <- segmented_copynumber
    return(simulation)
}
