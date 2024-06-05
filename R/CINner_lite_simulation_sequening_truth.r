simulation_sequening_truth <- function(simulation = list(),
                                       choice_theta = "",
                                       vec_theta_parameters = c(0.1, 0.1),
                                       truncal_mutations = 0,
                                       n_sample = 0,
                                       n_selective_clones = 1) {
    phylogeny_origin <- simulation$phylogeny_origin
    phylogeny_genotype <- simulation$phylogeny_genotype
    phylogeny_birthtime <- simulation$phylogeny_birthtime
    phylogeny_deathtime <- simulation$phylogeny_deathtime
    hclust_nodes <- simulation$hclust_nodes
    #--------------------------Count the number of progeny for each node
    #-------------------------------------------in sample phylogeny tree
    phylogeny_progeny_count <- rep(0, (2 * n_sample - 1))
    phylogeny_progeny_count[n_sample:(2 * n_sample - 1)] <- 1
    phylogeny_progeny_list <- vector("list", length = (2 * n_sample - 1))
    for (node_leaf in n_sample:(2 * n_sample - 1)) {
        phylogeny_progeny_list[[node_leaf]] <- node_leaf - n_sample + 1
    }
    for (node_child in (2 * n_sample - 1):1) {
        node_mother <- phylogeny_origin[node_child]
        if (node_mother <= 0) {
            next
        }
        phylogeny_progeny_count[node_mother] <-
            phylogeny_progeny_count[node_mother] +
            phylogeny_progeny_count[node_child]
        phylogeny_progeny_list[[node_mother]] <- c(
            phylogeny_progeny_list[[node_mother]],
            phylogeny_progeny_list[[node_child]]
        )
    }
    #--------------------------------Find background nodes of each clone
    #   Find leaves of each clone
    clonal_leaves <- vector("list", length = (n_selective_clones + 1))
    for (clone in 0:n_selective_clones) {
        clonal_leaves[[clone + 1]] <- which(
            phylogeny_genotype[n_sample:(2 * n_sample - 1)] == clone
        )
    }
    #   Find background nodes of each clone
    clonal_background <- vector("list", length = (n_selective_clones + 1))
    for (node in 1:(2 * n_sample - 1)) {
        for (clone in n_selective_clones:0) {
            if (all(is.element(clonal_leaves[[clone + 1]], phylogeny_progeny_list[[node]]))) {
                clonal_background[[clone + 1]] <- c(clonal_background[[clone + 1]], node)
            }
        }
    }
    #--------------------------Find the clonal identity of sampled cells
    sample_genotype <- phylogeny_genotype[n_sample:(2 * n_sample - 1)]
    #-------------Create the true cell-mutation matrix for sampled cells
    #   Find the lifetime of every node in the phylogeny
    phylogeny_lifetime <- phylogeny_deathtime - phylogeny_birthtime
    phylogeny_lifetime[which(phylogeny_lifetime < 0)] <- 0
    #   Find the theta for every node in the phylogeny
    if (choice_theta == "constant") {
        phylogeny_theta <- rep(0, length = (2 * n_sample - 1))
        for (clone in 0:n_selective_clones) {
            vec_loc <- which(phylogeny_genotype == clone)
            phylogeny_theta[vec_loc] <- vec_theta_parameters[clone + 1]
        }
    } else {
        if (choice_theta == "log-normal") {
            phylogeny_theta <- rep(0, length = (2 * n_sample - 1))
            for (clone in 0:n_selective_clones) {
                mu <- vec_theta_parameters[1, clone + 1]
                sd <- vec_theta_parameters[2, clone + 1]
                vec_loc <- which(phylogeny_genotype == clone)
                phylogeny_theta[vec_loc] <- rlnorm(
                    n = length(vec_loc),
                    meanlog = mu, sdlog = sd
                )
            }
        } else {
            if (choice_theta == "gamma") {
                phylogeny_theta <- rep(0, length = (2 * n_sample - 1))
                for (clone in 0:n_selective_clones) {
                    shape <- vec_theta_parameters[1, clone + 1]
                    scale <- vec_theta_parameters[2, clone + 1]
                    vec_loc <- which(phylogeny_genotype == clone)
                    phylogeny_theta[vec_loc] <- rgamma(
                        n = length(vec_loc),
                        shape = shape, scale = scale
                    )
                }
            }
        }
    }
    #   Simulate the mutation count for every node in the phylogeny
    phylogeny_mutation_count <- rpois(
        n = (2 * n_sample - 1),
        lambda = (phylogeny_theta * phylogeny_lifetime)
    )
    #   Add truncal mutations
    phylogeny_mutation_count[1] <- phylogeny_mutation_count[1] + truncal_mutations


    n_mutations <- sum(phylogeny_mutation_count)
    #-----------------Compute the cell-mutation matrix for sampled cells
    sample_mutational_table_truth_node_tips <- vector("list", length = (2 * n_sample - 1))
    sample_mutational_table_truth_node_mutation_counts <- vector("list", length = (2 * n_sample - 1))
    sample_mutational_table_truth_node_markers <- vector("list", length = (2 * n_sample - 1))
    sample_mutational_table_truth_node_hclust_indices <- vector("list", length = (2 * n_sample - 1))
    # sample_mutational_table_list <- vector("list", length = (2 * n_sample - 1))
    vec_colnames <- paste("Cell_", 1:n_sample, sep = "")
    N_rows <- 0
    for (node in 1:(2 * n_sample - 1)) {
        #   Find list of progeny of this node in the sample
        node_progeny_list <- phylogeny_progeny_list[[node]]
        #   Find mutation count of this node
        node_mutation_count <- phylogeny_mutation_count[node]
        if (node_mutation_count == 0) next
        #   Find node index in hclust style
        node_hclust_index <- hclust_nodes[node]
        #   Find node marker
        node_marker <- ""
        if (node %in% Reduce(intersect, clonal_background)) {
            node_marker <- "Truncal"
        } else {
            for (clone in 0:n_selective_clones) {
                if (is.element(node, clonal_background[[clone + 1]])) {
                    if (node_marker == "") {
                        node_marker <- paste0("Background_", clone)
                    } else {
                        node_marker <- paste0(node_marker, "&", clone)
                    }
                }
            }
        }
        if (node_marker == "") {
            node_genotype <- phylogeny_genotype[node]
            node_marker <- paste("Foreground_", node_genotype, sep = "")
        }
        ################################################################
        ################################################################
        ################################################################
        #   Store descendant tips
        sample_mutational_table_truth_node_tips[[node]] <- node_progeny_list
        #   Store mutational count
        sample_mutational_table_truth_node_mutation_counts[[node]] <- node_mutation_count
        #   Store mutational node marker (clonal background/foreground)
        sample_mutational_table_truth_node_markers[[node]] <- node_marker
        #   Store mutational node hclust index
        sample_mutational_table_truth_node_hclust_indices[[node]] <- node_hclust_index
        ################################################################
        ################################################################
        ################################################################
        # #   Compose pure mutational table
        # node_mutational_table <- matrix(0, nrow = node_mutation_count, ncol = n_sample)
        # node_mutational_table[, node_progeny_list] <- 1
        # sample_mutational_table_node <- as.data.frame(node_mutational_table)
        # colnames(sample_mutational_table_node) <- vec_colnames
        # #   Add column for node hclust index
        # sample_mutational_table_node <- cbind(hclust_index = node_hclust_index, sample_mutational_table_node)
        # #   Add column for node marker
        # sample_mutational_table_node <- cbind(Marker = node_marker, sample_mutational_table_node)
        # #   Add column for node row index
        # N_rows_start <- N_rows + 1
        # N_rows_end <- N_rows + node_mutation_count
        # sample_mutational_table_node <- cbind(Row = paste("Mutation_", N_rows_start:N_rows_end, sep = ""), sample_mutational_table_node)
        # N_rows <- N_rows_end
        # #   Store the mutational table for this node
        # sample_mutational_table_list[[node]] <- sample_mutational_table_node
    }
    # sample_mutational_table <- as.data.frame(do.call(rbind, list_node_mutational_table))
    # colnames(sample_mutational_table) <- vec_colnames
    # sample_mutational_table <- cbind(hclust_index = unlist(sample_mutational_table_truth_node_hclust_indices), sample_mutational_table)
    # sample_mutational_table <- cbind(Marker = unlist(sample_mutational_table_truth_node_markers), sample_mutational_table)
    # sample_mutational_table <- cbind(Count = unlist(sample_mutational_table_truth_node_mutation_counts), sample_mutational_table)
    #---------------------------------------------Prepare output package
    # simulation$sample_mutational_table_truth <- sample_mutational_table
    simulation$sample_mutational_table_truth_node_tips <- sample_mutational_table_truth_node_tips
    simulation$sample_mutational_table_truth_node_mutation_counts <- sample_mutational_table_truth_node_mutation_counts
    simulation$sample_mutational_table_truth_node_markers <- sample_mutational_table_truth_node_markers
    simulation$sample_mutational_table_truth_node_hclust_indices <- sample_mutational_table_truth_node_hclust_indices
    simulation$sample_genotype <- sample_genotype
    return(simulation)
}
