simulation_sample_phylogeny <- function(simulation = list(),
                                        t_end_time = 0,
                                        vec_s_mut_time = 0,
                                        vec_s_mut_hierarchy = 0,
                                        n_selective_clones = 1,
                                        n_sample = 0) {
    t_current <- t_end_time
    record_t_previous <- simulation$record_t_previous
    record_t_now <- simulation$record_t_now
    record_vec_populations <- simulation$record_vec_populations
    record_vec_count_division <- simulation$record_vec_count_division
    record_vec_count_death <- simulation$record_vec_count_death
    #---------------------Initialize the phylogeny record - hclust style
    #   Initialize information to build phylogeny in hclust style
    hclust_row <- 0
    hclust_nodes <- rep(0, length = (2 * n_sample - 1))
    hclust_nodes[n_sample:(2 * n_sample - 1)] <- -1:-n_sample
    #   Initialize phylogeny in hclust style
    hclust_merge <- matrix(0, nrow = (n_sample - 1), ncol = 2)
    hclust_height <- rep(0, length = (n_sample - 1))
    #------------------------Initialize the phylogeny record - our style
    phylogeny_origin <- rep(0, length = (2 * n_sample - 1))
    phylogeny_genotype <- rep(0, length = (2 * n_sample - 1))
    phylogeny_birthtime <- rep(0, length = (2 * n_sample - 1))
    phylogeny_deathtime <- rep(0, length = (2 * n_sample - 1))
    #---------------------------Find a random sample of final population
    #   Initialize the current list of nodes in the sample phylogeny
    node_list_current <- n_sample:(2 * n_sample - 1)
    #   Initialize the current list of genoytpes of the nodes
    final_clonal_ID <- 0:n_selective_clones
    final_clonal_population <-
        record_vec_populations[nrow(record_vec_populations), ]
    final_population <- c()
    for (i in 1:length(final_clonal_ID)) {
        final_population <-
            c(
                final_population,
                rep(final_clonal_ID[i], length = final_clonal_population[i])
            )
    }
    node_genotype_current <- sort(sample(x = final_population, size = n_sample, replace = FALSE))
    #   Initialize data for leaves of sample phylogeny
    phylogeny_genotype[node_list_current] <- node_genotype_current
    phylogeny_deathtime[node_list_current] <- t_current
    #------------------------------------Build the sample phylogeny tree
    for (i in seq(length(record_t_now), 1, -1)) {
        #   Get time points
        t_previous <- record_t_previous[i]
        t_now <- record_t_now[i]
        #   Get current total clonal population (after divisions)
        total_clonal_ID <- final_clonal_ID
        total_clonal_population <- record_vec_populations[i, ]
        #   Get current sample clonal population (after divisions)
        sample_clonal_population <- rep(0, length = (n_selective_clones + 1))
        for (clone in 0:n_selective_clones) {
            sample_clonal_population[clone + 1] <- length(which(node_genotype_current == clone))
        }
        #   Get list of eligible nodes of each genotype
        sample_eligible_nodes <- vector("list",
            length = (n_selective_clones + 1)
        )
        for (clone in 0:n_selective_clones) {
            sample_eligible_nodes[[clone + 1]] <- node_list_current[which(node_genotype_current == clone)]
        }
        #---------------------------------------Perform normal divisions
        #----------(mother cell and daughter cells belong to same clone)
        #   Get list of normal divisions
        matrix_division <- record_vec_count_division[i, ]
        #   For each type of divisions...
        for (event_type in 1:(n_selective_clones + 1)) {
            #   Get number of divisions
            no_divisions <- matrix_division[event_type]
            if (no_divisions == 0) {
                next
            }
            #   Get genotype of mother and daughters
            genotype <- event_type - 1
            #   If daughter genotypes are not in current nodes, move on
            if (sample_clonal_population[genotype + 1] <= 0) {
                next
            }
            #   Get list of nodes in total population to merge
            list_nodes_available <- sample_eligible_nodes[[genotype + 1]]
            list_nodes_available <- c(
                list_nodes_available,
                rep(0, length = (total_clonal_population[genotype + 1] -
                    sample_clonal_population[genotype + 1]))
            )
            #   Decide which nodes in total population to merge
            list_nodes_two_daughters <- sample(
                x = list_nodes_available,
                size = 2 * no_divisions,
                replace = FALSE
            )
            list_nodes_daughter_1 <- list_nodes_two_daughters[1:no_divisions]
            list_nodes_daughter_2 <- list_nodes_two_daughters[(no_divisions + 1):
            (2 * no_divisions)]
            #--------------------------------------Perform node mergings
            list_events <- which((list_nodes_daughter_1 > 0) & (list_nodes_daughter_2 > 0))
            if (length(list_events) == 0) {
                next
            }
            #   Get list of nodes to merge
            vec_node_1 <- list_nodes_daughter_1[list_events]
            vec_node_2 <- list_nodes_daughter_2[list_events]
            #   Find list of internal nodes to merge into
            vec_node_mother <-
                (min(node_list_current) - length(list_events)):
                (min(node_list_current) - 1)
            #   Update phylogeny in hclust style
            for (event in 1:length(list_events)) {
                node_1 <- vec_node_1[event]
                node_2 <- vec_node_2[event]
                node_mother <- vec_node_mother[event]
                hclust_row <- hclust_row + 1
                hclust_nodes[node_mother] <- hclust_row
                hclust_merge[hclust_row, ] <- c(
                    hclust_nodes[node_1],
                    hclust_nodes[node_2]
                )
                hclust_height[hclust_row] <- t_end_time - t_previous
            }
            #   Update phylogeny in our style
            phylogeny_origin[c(vec_node_1, vec_node_2)] <- vec_node_mother
            phylogeny_genotype[vec_node_mother] <- genotype
            phylogeny_birthtime[c(vec_node_1, vec_node_2)] <- t_previous
            phylogeny_deathtime[vec_node_mother] <- t_previous
            vec_pos_delete <- which(is.element(
                node_list_current,
                c(vec_node_1, vec_node_2)
            ))
            node_list_current <- node_list_current[-vec_pos_delete]
            node_list_current <- c(vec_node_mother, node_list_current)
            node_genotype_current <- node_genotype_current[-vec_pos_delete]
            node_genotype_current <- c(
                rep(genotype, length = length(list_events)),
                node_genotype_current
            )
        }
        #---------------------------------------Perform selective sweeps
        #--------------------------(mother cell gains selective mutation
        #-----------------------------------to become one daughter cell)
        i_next_selective_mut <- which(t_previous <= vec_s_mut_time)
        if (length(i_next_selective_mut) == 0) {
            next
        }
        i_next_selective_mut <- i_next_selective_mut[1]
        t_next_selective_mut <- vec_s_mut_time[i_next_selective_mut]
        if (t_now > t_next_selective_mut) {
            #   Find which mother clone gives rise to which daughter clone
            genotype_mother <- vec_s_mut_hierarchy[i_next_selective_mut]
            # genotype_mother <- i_next_selective_mut - 1
            genotype_daughter <- i_next_selective_mut
            #   Find which node belongs to the daughter clone
            node_selection <- node_list_current[
                which(node_genotype_current == genotype_daughter)
            ]
            #   Update the node
            node_genotype_current[which(node_list_current ==
                node_selection)] <- genotype_mother
        }
    }
    #----------------------------------------Get MRCA age for each clone
    #   Find all subclonal identities for each clone
    list_subclones <- vector("list", length = (n_selective_clones + 1))
    for (clone in n_selective_clones:0) {
        clone_index <- clone + 1
        list_subclones[[clone_index]] <- sort(c(clone, list_subclones[[clone_index]]))
        if (clone > 0) {
            clone_mother <- vec_s_mut_hierarchy[clone]
            clone_mother_index <- clone_mother + 1
            list_subclones[[clone_mother_index]] <- sort(c(list_subclones[[clone_mother_index]], list_subclones[[clone_index]]))
        }
    }
    #   Find ancestors for every leaf
    list_leaves <- phylogeny_genotype[(length(phylogeny_genotype) - n_sample + 1):length(phylogeny_genotype)]
    list_ancestors <- vector("list", length = length(list_leaves))
    for (leaf in (length(phylogeny_genotype) - n_sample + 1):length(phylogeny_genotype)) {
        vec_ancestors <- c()
        node <- leaf
        while (node > 0) {
            vec_ancestors <- c(node, vec_ancestors)
            node <- phylogeny_origin[node]
        }
        list_ancestors[[leaf - (length(phylogeny_genotype) - n_sample + 1) + 1]] <- vec_ancestors
    }
    #   Find MRCA age for every clone
    MRCA_ages <- rep(0, n_selective_clones + 1)
    for (clone in 0:n_selective_clones) {
        subclones <- list_subclones[[clone + 1]]
        ancestors <- NULL
        for (leaf in 1:length(list_leaves)) {
            if (list_leaves[leaf] %in% subclones) {
                if (is.null(ancestors)) {
                    ancestors <- list_ancestors[[leaf]]
                } else {
                    ancestors <- intersect(ancestors, list_ancestors[[leaf]])
                }
            }
        }
        MRCA_ages[clone + 1] <- max(phylogeny_deathtime[ancestors])
    }
    #--------------------------------------------------Reorder the nodes
    list_roots <- which(phylogeny_origin == 0)
    #---Find an order on all nodes of the phylogeny in our style
    #   Find number of progeny of each node
    progeny_count <- rep(0, length = (2 * n_sample - 1))
    progeny_count[n_sample:(2 * n_sample - 1)] <- 1
    for (node in ((2 * n_sample - 1)):1) {
        mother_node <- phylogeny_origin[node]
        if (mother_node > 0) {
            progeny_count[mother_node] <- progeny_count[mother_node] + progeny_count[node]
        }
    }
    #   Reorder the sample phylogeny tree based on progeny counts
    phylogeny_order <- rep(0, length = (2 * n_sample - 1))
    phylogeny_order[list_roots] <- 1
    for (node in 0:(2 * n_sample - 1)) {
        vec_daughter_nodes <- which(phylogeny_origin == node)
        if (length(vec_daughter_nodes) == 0) {
            next
        }
        vec_progeny_counts <- progeny_count[vec_daughter_nodes]
        tmp <- sort(vec_progeny_counts, index.return = TRUE)
        vec_progeny_counts <- tmp$x
        vec_order <- tmp$ix
        vec_daughter_nodes <- vec_daughter_nodes[vec_order]
        for (i in 1:length(vec_daughter_nodes)) {
            daughter_node <- vec_daughter_nodes[i]
            if (i > 1) {
                progeny_count_extra <- sum(vec_progeny_counts[1:i - 1])
            } else {
                progeny_count_extra <- 0
            }
            if (node == 0) {
                phylogeny_order[daughter_node] <- phylogeny_order[daughter_node] + progeny_count_extra
            } else {
                phylogeny_order[daughter_node] <- phylogeny_order[node] + progeny_count_extra
            }
        }
    }
    #---Extract the order for phylogeny in hclust style
    hclust_order_inverse <- phylogeny_order[n_sample:(2 * n_sample - 1)]
    hclust_order <- rep(0, n_sample)
    for (i_cell in 1:n_sample) {
        loc <- hclust_order_inverse[i_cell]
        hclust_order[loc] <- i_cell
    }
    #--------------------------------------------Create clustering table
    sample_cell_ID <- paste("Cell_", 1:n_sample, sep = "")
    sample_clone_ID <- paste("Clone_", phylogeny_genotype[n_sample:(2 * n_sample - 1)], sep = "")
    cell_phylogeny_clustering <- data.frame(sample_cell_ID, sample_clone_ID)
    names(cell_phylogeny_clustering) <- c("cell_id", "clone_id")
    #----------------------------Create phylogeny object in hclust style
    cell_phylogeny_hclust <- list()
    cell_phylogeny_hclust$merge <- hclust_merge
    cell_phylogeny_hclust$height <- hclust_height
    cell_phylogeny_hclust$order <- hclust_order
    cell_phylogeny_hclust$labels <- sample_cell_ID
    class(cell_phylogeny_hclust) <- "hclust"
    #---------------------------------------------Prepare output package
    simulation$cell_phylogeny_hclust <- cell_phylogeny_hclust
    simulation$cell_phylogeny_clustering <- cell_phylogeny_clustering
    simulation$phylogeny_origin <- phylogeny_origin
    simulation$phylogeny_genotype <- phylogeny_genotype
    simulation$phylogeny_birthtime <- phylogeny_birthtime
    simulation$phylogeny_deathtime <- phylogeny_deathtime
    simulation$hclust_nodes <- hclust_nodes
    simulation$MRCA_ages <- MRCA_ages
    return(simulation)
}
