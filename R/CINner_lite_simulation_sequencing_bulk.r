simulation_sequencing_bulk <- function(simulation = list(),
                                       n_sample = 0,
                                       ploidy = 1,
                                       purity = 1,
                                       bulk_coverage_model = "",
                                       bulk_coverage_variables = c(),
                                       bulk_min_alt_readcounts = 0) {
    bulk_coverage_variables$n_sample <- n_sample
    #----------------------------------Find the counts for all mutations
    mutation_count <- unlist(simulation$sample_mutational_table_truth_node_mutation_counts)
    n_mutation_true <- sum(mutation_count)
    #--------------------------------Find the true VAF for all mutations
    mutation_VAF <- lengths(simulation$sample_mutational_table_truth_node_tips) / n_sample
    mutation_VAF <- mutation_VAF[mutation_VAF != 0]
    mutation_true_VAF <- rep(purity * mutation_VAF / ploidy, mutation_count)
    mutation_node_markers <- rep(unlist(simulation$sample_mutational_table_truth_node_markers), mutation_count)
    #-------------------Simulate the total read counts for all mutations
    mutation_readcount_tot <- rand_coverage(
        bulk_coverage_model = bulk_coverage_model,
        bulk_coverage_variables = bulk_coverage_variables,
        n_mutation_true = n_mutation_true
    )
    #---------------Simulate the alternate read counts for all mutations
    mutation_readcount_alt <- rbinom(
        n = length(mutation_readcount_tot),
        size = mutation_readcount_tot,
        prob = mutation_true_VAF
    )
    #---------------Remove mutations not satisfying readcount conditions
    vec_delete <- which(mutation_readcount_alt < bulk_min_alt_readcounts)
    if (length(vec_delete) > 0) {
        mutation_readcount_tot <- mutation_readcount_tot[-vec_delete]
        mutation_readcount_alt <- mutation_readcount_alt[-vec_delete]
        mutation_node_markers <- mutation_node_markers[-vec_delete]
    }
    #-------------------Find the reference read counts for all mutations
    mutation_readcount_ref <- mutation_readcount_tot - mutation_readcount_alt
    #---------------------------Create readcount table for the mutations
    sample_mutational_table_bulk <- data.frame(mutation_readcount_ref, mutation_readcount_alt, mutation_node_markers)
    colnames(sample_mutational_table_bulk) <- c("Ref_count", "Alt_count", "Marker")
    #---------------------------------------------Prepare output package
    simulation$sample_mutational_table_bulk <- sample_mutational_table_bulk
    return(simulation)
}
rand_coverage <- function(bulk_coverage_model,
                          bulk_coverage_variables,
                          n_mutation_true) {
    n_end <- bulk_coverage_variables$n_sample
    if (bulk_coverage_model == "constant") {
        n_coverage <- bulk_coverage_variables$mean_coverage
        output <- rep(n_coverage, length = n_mutation_true)
    } else if (bulk_coverage_model == "uniform") {
        r_min <- bulk_coverage_variables$min_coverage
        r_max <- bulk_coverage_variables$max_coverage
        output <- sample(
            size = n_mutation_true,
            x = r_min:r_max,
            replace = TRUE
        )
    } else if (bulk_coverage_model == "binomial") {
        mean_coverage <- bulk_coverage_variables$mean_coverage
        output <- rbinom(
            n = n_mutation_true,
            size = n_end,
            prob = mean_coverage / n_end
        )
    } else if (bulk_coverage_model == "custom") {
        coverage <- bulk_coverage_variables$coverage
        output <- sample(
            size = n_mutation_true,
            x = coverage$Read_count,
            prob = coverage$Frequency / sum(coverage$Frequency),
            replace = TRUE
        )
    }
    return(output)
}
