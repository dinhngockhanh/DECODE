obsSFS_generate <- function(subclonal_cell_count = c(),
                            subclonal_mutation_count = c(),
                            subclonal_parent = c(),
                            tail_mutation_count = 0,
                            tail_power = NULL,
                            purity = 1,
                            coverage_distribution = data.frame(
                                Coverage = 1:100,
                                Probability = dbinom(1:100, size = 100, prob = 0.3)
                            ),
                            min_variant_read = 3,
                            output_dir = NULL,
                            plot = TRUE,
                            run_name = "SFS_simulation") {
    nTrials <- 100000
    #---Check the validity of input parameters
    if (length(subclonal_cell_count) != length(subclonal_mutation_count) || length(subclonal_cell_count) != length(subclonal_parent)) stop("Subclonal vectors must have the same length.")
    invalid_parents <- which(subclonal_parent >= seq_along(subclonal_parent))
    if (length(invalid_parents) > 0) {
        stop(paste0(
            "Invalid parent references: subclonal_parent[", paste(invalid_parents, collapse = ", "),
            "] cannot be >= their position in the vector."
        ))
    }
    clones <- length(subclonal_cell_count)
    total_cancer_cell_count <- sum(subclonal_cell_count)
    total_sample_cell_count <- round(total_cancer_cell_count / purity)
    all_parameters <- data.frame(matrix(ncol = 0, nrow = 1))
    #---Simulate SFS tail
    if ((tail_mutation_count > 0) && (!is.null(tail_power))) {
        #   Compute exact SFS pdf
        exact_SFS <- data.frame(
            VAF = (1:total_cancer_cell_count) / (2 * total_sample_cell_count),
            Probability = (1:total_cancer_cell_count)^(-tail_power)
        )
        exact_SFS$Probability <- exact_SFS$Probability / sum(exact_SFS$Probability)
        #   Simulate observed SFS until enough tail mutations
        Alt_counts <- c()
        Ref_counts <- c()
        Nmut_exact <- 0
        while (length(Alt_counts) < tail_mutation_count) {
            exact_VAFs <- sample(
                exact_SFS$VAF,
                size = nTrials,
                replace = TRUE,
                prob = exact_SFS$Probability
            )
            Tot_counts_next <- sample(
                coverage_distribution$Coverage,
                size = nTrials,
                replace = TRUE,
                prob = coverage_distribution$Probability
            )
            Alt_counts_next <- rbinom(nTrials, size = Tot_counts_next, prob = exact_VAFs)
            valid_indices <- which((Alt_counts_next >= min_variant_read) & (Alt_counts_next < Tot_counts_next))
            if ((length(Alt_counts) + length(valid_indices)) < tail_mutation_count) {
                Nmut_exact <- Nmut_exact + nTrials
            } else {
                valid_indices <- valid_indices[1:(tail_mutation_count - length(Alt_counts))]
                Nmut_exact <- Nmut_exact + valid_indices[length(valid_indices)]
            }
            Alt_counts <- c(Alt_counts, Alt_counts_next[valid_indices])
            Ref_counts <- c(Ref_counts, Tot_counts_next[valid_indices] - Alt_counts_next[valid_indices])
        }
        mutation_table <- data.frame(
            Mutation_ID = "Tail",
            Alt_count = Alt_counts,
            Ref_count = Ref_counts
        )
        all_parameters[["Tail_power"]] <- tail_power
        all_parameters[["Tail_Nmut_exact"]] <- Nmut_exact
        all_parameters[["Tail_Nmut_observed"]] <- tail_mutation_count
    } else {
        mutation_table <- data.frame(
            Mutation_ID = character(0),
            Alt_count = integer(0),
            Ref_count = integer(0)
        )
    }
    #---Tabulate cluster information based on clonal structure
    if (clones >= 1) {
        for (cluster in clones:1) {
            parent <- subclonal_parent[cluster]
            if (parent <= 0) next
            subclonal_cell_count[parent] <- subclonal_cell_count[parent] + subclonal_cell_count[cluster]
        }
    }
    #---Simulate SFS clusters
    for (cluster in seq_len(clones)) {
        cluster_frequency <- subclonal_cell_count[cluster] / (2 * total_sample_cell_count)
        cluster_mutation_count <- subclonal_mutation_count[cluster]
        #   Simulate observed SFS until enough cluster mutations
        Alt_counts <- c()
        Ref_counts <- c()
        Nmut_exact <- 0
        while (length(Alt_counts) < cluster_mutation_count) {
            Tot_counts_next <- sample(
                coverage_distribution$Coverage,
                size = nTrials,
                replace = TRUE,
                prob = coverage_distribution$Probability
            )
            Alt_counts_next <- rbinom(nTrials, size = Tot_counts_next, prob = cluster_frequency)
            valid_indices <- which((Alt_counts_next >= min_variant_read) & (Alt_counts_next < Tot_counts_next))
            if ((length(Alt_counts) + length(valid_indices)) < cluster_mutation_count) {
                Nmut_exact <- Nmut_exact + nTrials
            } else {
                valid_indices <- valid_indices[1:(cluster_mutation_count - length(Alt_counts))]
                Nmut_exact <- Nmut_exact + valid_indices[length(valid_indices)]
            }
            Alt_counts <- c(Alt_counts, Alt_counts_next[valid_indices])
            Ref_counts <- c(Ref_counts, Tot_counts_next[valid_indices] - Alt_counts_next[valid_indices])
        }
        mutation_table <- rbind(
            mutation_table,
            data.frame(
                Mutation_ID = paste0("Cluster_", cluster),
                Alt_count = Alt_counts,
                Ref_count = Ref_counts
            )
        )
        all_parameters[[paste0("Cluster_", cluster, "_freq")]] <- cluster_frequency
        all_parameters[[paste0("Cluster_", cluster, "_Nmut_exact")]] <- Nmut_exact
        all_parameters[[paste0("Cluster_", cluster, "_Nmut_observed")]] <- cluster_mutation_count
    }
    #---Output simulated results
    if (!is.null(output_dir)) {
        write.csv(mutation_table, file = file.path(output_dir, paste0(run_name, "_mutation_table.csv")), row.names = FALSE)
        write.csv(all_parameters, file = file.path(output_dir, paste0(run_name, "_parameter_summary.csv")), row.names = FALSE)
        if (plot) {
            library(ggplot2)
            observed_SFS_plot <- ggplot(mutation_table, aes(x = Alt_count / (Alt_count + Ref_count))) +
                geom_histogram(binwidth = 0.01) +
                labs(title = "Observed SFS", x = "Variant Allele Frequency", y = "Number of Mutations") +
                theme_minimal()
            ggsave(filename = file.path(output_dir, paste0(run_name, "_Observed_SFS_plot.jpg")), plot = observed_SFS_plot)
        }
    }
    return(list(
        mutation_table = mutation_table,
        parameters = all_parameters
    ))
}
