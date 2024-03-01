CINner_add_passenger_mutations <- function(simulation, passenger_mutation_rate) {
    library(data.table)
    #----------------------------------------Input the general variables
    size_CN_block_DNA <- as.numeric(simulation$parameters$general_variables$Value[which(simulation$parameters$general_variables$Variable == "size_CN_block_DNA")])
    #-----------------------------------------Input the clonal evolution
    genotype_list_ploidy_chrom <- simulation$clonal_evolution$genotype_list_ploidy_chrom
    genotype_list_ploidy_block <- simulation$clonal_evolution$genotype_list_ploidy_block
    evolution_genotype_changes <- simulation$clonal_evolution$evolution_genotype_changes
    #---------------------------------------------------Input the sample
    sample_genotype_unique <- simulation$sample$sample_genotype_unique
    all_sample_genotype <- simulation$sample$all_sample_genotype
    N_sample <- length(simulation$sample$sample_cell_ID)
    #-----------------------------------------Input the sample phylogeny
    phylogeny_origin <- simulation$sample_phylogeny$package_cell_phylogeny$phylogeny_origin
    phylogeny_elapsed_gens <- simulation$sample_phylogeny$package_cell_phylogeny$phylogeny_elapsed_gens
    phylogeny_elapsed_genotypes <- simulation$sample_phylogeny$package_cell_phylogeny$phylogeny_elapsed_genotypes
    phylogeny_genotype <- simulation$sample_phylogeny$package_cell_phylogeny$phylogeny_genotype
    phylogeny_birthtime <- simulation$sample_phylogeny$package_cell_phylogeny$phylogeny_birthtime
    phylogeny_deathtime <- simulation$sample_phylogeny$package_cell_phylogeny$phylogeny_deathtime
    phylogeny_order <- simulation$sample_phylogeny$package_cell_phylogeny$phylogeny_order
    #----------------Get information about genotypes in sample phylogeny
    #   Find all genotypes in the sample phylogeny
    all_genotypes_ID <- unique(unlist(phylogeny_elapsed_genotypes))
    #   Find total length of each chromosome strand in each genotype
    all_genotypes_strand_length <- vector("list", length = length(all_genotypes_ID))
    for (clone in 1:length(all_genotypes_ID)) {
        clone_ID <- all_genotypes_ID[clone]
        ploidy_chrom <- genotype_list_ploidy_chrom[[clone_ID]]
        ploidy_block <- genotype_list_ploidy_block[[clone_ID]]
        chroms <- c()
        strands <- c()
        lengths <- c()
        for (chrom in 1:length(ploidy_chrom)) {
            if (ploidy_chrom[chrom] == 0) next
            for (strand in 1:ploidy_chrom[chrom]) {
                chroms <- c(chroms, chrom)
                strands <- c(strands, strand)
                lengths <- c(lengths, sum(ploidy_block[[chrom]][[strand]]))
            }
        }
        all_genotypes_strand_length[[clone]] <- rbind(chroms, strands, lengths)
    }
    #   Find Poisson probability of passenger mutations in each genotype
    #   = (mutation rate / nucleotide / division) X (genome length in blocks) X (block size in nucleotides)
    all_genotypes_prob_new_passengers <- rep(0, length(all_genotypes_ID))
    for (clone in 1:length(all_genotypes_ID)) {
        all_genotypes_prob_new_passengers[clone] <- passenger_mutation_rate * sum(all_genotypes_strand_length[[clone]][3, ]) * size_CN_block_DNA
    }
    #-----------------------------------Simulate the passenger mutations
    phylogeny_passenger_mutations <- vector("list", length = length(phylogeny_genotype))
    passenger_ID <- 0
    pb <- txtProgressBar(
        min = 1, max = length(phylogeny_origin),
        style = 3, width = 50, char = "="
    )
    for (branch in 1:length(phylogeny_origin)) {
        setTxtProgressBar(pb, branch)
        # ==Get mother genotype with existing mutations
        mother_branch <- phylogeny_origin[branch]
        if (mother_branch <= 0) {
            mother_passenger_mutations <- c()
            mother_genotype <- 0
        } else {
            mother_passenger_mutations <- phylogeny_passenger_mutations[[mother_branch]]
            mother_genotype <- phylogeny_genotype[mother_branch]
        }
        # ==Get list of elapsed genotypes in original phylogeny
        elapsed_genotypes <- phylogeny_elapsed_genotypes[[branch]]
        # ==Update passenger mutations after each elapsed cell generation
        passenger_mutations <- mother_passenger_mutations
        for (elapsed_gen in 1:length(elapsed_genotypes)) {
            #---Find genotype of the previous cell generation
            if (elapsed_gen == 1) {
                previous_genotype <- mother_genotype
            } else {
                previous_genotype <- elapsed_genotypes[elapsed_gen - 1]
            }
            #---Modify existing mutations if CNAs have occurred
            if (elapsed_genotypes[elapsed_gen] != previous_genotype) {
                ploidy_chrom <- genotype_list_ploidy_chrom[[elapsed_genotypes[elapsed_gen]]]
                ploidy_block <- genotype_list_ploidy_block[[elapsed_genotypes[elapsed_gen]]]
                events <- evolution_genotype_changes[[elapsed_genotypes[elapsed_gen]]]
                if (!is.null(events)) {
                    for (j in 1:length(events)) {
                        event <- events[[j]]
                        if (event[1] == "new-driver") {
                            #   If a driver mutation occurs, no passenger mutations are affected...
                        } else if (event[1] == "whole-genome-duplication") {
                            #   If a whole-genome duplication occurs,
                            #   all passenger mutations are duplicated...
                            if (nrow(passenger_mutations) > 0) {
                                new_passenger_mutations <- passenger_mutations
                                if (!is.matrix(new_passenger_mutations)) new_passenger_mutations <- matrix(new_passenger_mutations, nrow = 1)
                                for (mutation in 1:nrow(new_passenger_mutations)) {
                                    chrom <- new_passenger_mutations[mutation, 2]
                                    chrom_ploidy <- ploidy_chrom[chrom]
                                    new_passenger_mutations[mutation, 3] <- new_passenger_mutations[mutation, 3] + chrom_ploidy / 2
                                }
                                passenger_mutations <- rbind(passenger_mutations, new_passenger_mutations)
                            }
                        } else if (event[1] == "missegregation") {
                            #   If a missegregation occurs,
                            #   passenger mutations on the affected chromosome
                            #   are doubled (if gain) or deleted (if lost)...
                            if (nrow(passenger_mutations) > 0) {
                                affected_chrom <- event[2]
                                affected_strand <- event[3]
                                affected_status <- event[4]
                                pos_affected_mutations <- intersect(which((passenger_mutations[, 2] == affected_chrom)), which((passenger_mutations[, 3] == affected_strand)))
                                if (length(pos_affected_mutations) > 0) {
                                    if (affected_status == 1) {
                                        new_passenger_mutations <- passenger_mutations[pos_affected_mutations, ]
                                        if (!is.matrix(new_passenger_mutations)) new_passenger_mutations <- matrix(new_passenger_mutations, nrow = 1)
                                        new_passenger_mutations[, 3] <- ploidy_chrom[affected_chrom]
                                        passenger_mutations <- rbind(passenger_mutations, new_passenger_mutations)
                                    } else if (affected_status == -1) {
                                        passenger_mutations <- passenger_mutations[-pos_affected_mutations, ]
                                    }
                                }
                            }
                        } else if (event[1] == "chromosome-arm-missegregation") {
                            #   If an arm-missegregation occurs,
                            #   passenger mutations on the affected chromosome arm
                            #   are doubled (if gain) or deleted (if lost)...
                            if (nrow(passenger_mutations) > 0) {
                                affected_chrom <- event[2]
                                affected_chrom_bin_count <- simulation$parameters$cn_info$Bin_count[which(simulation$parameters$cn_info$Chromosome == affected_chrom)]
                                affected_chrom_centromere_location <- simulation$parameters$cn_info$Centromere_location[which(simulation$parameters$cn_info$Chromosome == affected_chrom)]
                                affected_strand <- event[3]
                                affected_arm <- event[4]
                                if (affected_arm == 1) {
                                    affected_block_start <- 1
                                    affected_block_end <- affected_chrom_centromere_location
                                } else if (affected_arm == 2) {
                                    affected_block_start <- affected_chrom_centromere_location + 1
                                    affected_block_end <- affected_chrom_bin_count
                                }
                                affected_status <- event[5]
                                pos_affected_mutations <- intersect(intersect(which(passenger_mutations[, 2] == affected_chrom), which(passenger_mutations[, 3] == affected_strand)), intersect(which(passenger_mutations[, 4] >= affected_block_start), which(passenger_mutations[, 4] <= affected_block_end)))
                                if (length(pos_affected_mutations) > 0) {
                                    if (affected_status == 1) {
                                        new_passenger_mutations <- passenger_mutations[pos_affected_mutations, ]
                                        if (!is.matrix(new_passenger_mutations)) new_passenger_mutations <- matrix(new_passenger_mutations, nrow = 1)
                                        new_passenger_mutations[, 3] <- ploidy_chrom[affected_chrom]
                                        passenger_mutations <- rbind(passenger_mutations, new_passenger_mutations)
                                    } else if (affected_status == -1) {
                                        passenger_mutations <- passenger_mutations[-pos_affected_mutations, ]
                                    }
                                }
                            }
                        } else if (event[1] == "focal-amplification") {
                            #   If a focal amplification occurs,
                            #   passenger mutations in the affected region
                            #   are duplicated...
                            if (nrow(passenger_mutations) > 0) {
                                affected_chrom <- event[2]
                                affected_strand <- event[3]
                                affected_block_start <- event[4]
                                affected_block_end <- event[5]
                                pos_affected_mutations <- intersect(intersect(which(passenger_mutations[, 2] == affected_chrom), which(passenger_mutations[, 3] == affected_strand)), intersect(which(passenger_mutations[, 4] >= affected_block_start), which(passenger_mutations[, 4] <= affected_block_end)))
                                if (length(pos_affected_mutations) > 0) {
                                    new_passenger_mutations <- passenger_mutations[pos_affected_mutations, ]
                                    if (!is.matrix(new_passenger_mutations)) new_passenger_mutations <- matrix(new_passenger_mutations, nrow = 1)
                                    for (mutation in 1:nrow(new_passenger_mutations)) {
                                        block <- new_passenger_mutations[mutation, 4]
                                        mutation_ploidy_block <- ploidy_block[affected_chrom][[affected_strand]][block] / 2
                                        new_passenger_mutations[mutation, 5] <- new_passenger_mutations[mutation, 5] + mutation_ploidy_block
                                    }
                                    passenger_mutations <- rbind(passenger_mutations, new_passenger_mutations)
                                }
                            }
                        } else if (event[1] == "focal-deletion") {
                            #   If a focal deletion occurs,
                            #   passenger mutations in the affected region
                            #   are deleted...
                            if (nrow(passenger_mutations) > 0) {
                                affected_chrom <- event[2]
                                affected_strand <- event[3]
                                affected_block_start <- event[4]
                                affected_block_end <- event[5]
                                pos_affected_mutations <- intersect(intersect(which(passenger_mutations[, 2] == affected_chrom), which(passenger_mutations[, 3] == affected_strand)), intersect(which(passenger_mutations[, 4] >= affected_block_start), which(passenger_mutations[, 4] <= affected_block_end)))
                                if (length(pos_affected_mutations) > 0) {
                                    passenger_mutations <- passenger_mutations[-pos_affected_mutations, ]
                                }
                            }
                        } else if (event[1] == "cnloh-interstitial") {
                            #   ...
                            #   ...
                            #   ...
                        } else if (event[1] == "cnloh-terminal") {
                            #   ...
                            #   ...
                            #   ...
                        }
                    }
                }
            }
            if (!is.null(passenger_mutations) & !is.matrix(passenger_mutations)) passenger_mutations <- matrix(passenger_mutations, nrow = 1)
            #---Add new passenger mutations
            genotype <- elapsed_genotypes[elapsed_gen]
            prob_new_passengers <- all_genotypes_prob_new_passengers[which(all_genotypes_ID == genotype)]
            count_new_passengers <- rpois(1, prob_new_passengers)
            for (new_passenger in 1:count_new_passengers) {
                passenger_ID <- passenger_ID + 1
                #   Choose a random chromosome & strand
                lengths <- all_genotypes_strand_length[[which(all_genotypes_ID == genotype)]][3, ]
                tmp <- sample(1:length(lengths), size = 1, prob = lengths / sum(lengths))
                chrom <- all_genotypes_strand_length[[which(all_genotypes_ID == genotype)]][1, tmp]
                strand <- all_genotypes_strand_length[[which(all_genotypes_ID == genotype)]][2, tmp]
                #   Choose a random position in the chromosome strand
                block_sizes <- genotype_list_ploidy_block[[genotype]][[chrom]][[strand]]
                block <- sample(1:length(block_sizes), size = 1, prob = block_sizes / sum(block_sizes))
                unit <- sample(1:block_sizes[block], size = 1)
                location <- size_CN_block_DNA * (block - 1) + sample(1:size_CN_block_DNA, size = 1)
                #   Add the new passenger mutation
                passenger_mutations <- rbind(passenger_mutations, c(passenger_ID, chrom, strand, block, unit, location))
            }
        }
        phylogeny_passenger_mutations[[branch]] <- passenger_mutations
    }
    cat("\n")
    #--------------------------Find sample passenger mutation statistics
    #   Find total copy number of each clone in the sample
    sample_genotype_unique_copy_number <- vector("list", length = length(sample_genotype_unique))
    for (clone in 1:length(sample_genotype_unique)) {
        clone_ID <- sample_genotype_unique[clone]
        ploidy_chrom <- genotype_list_ploidy_chrom[[clone_ID]]
        ploidy_block <- genotype_list_ploidy_block[[clone_ID]]
        copy_number <- list()
        for (chrom in 1:length(ploidy_chrom)) {
            copy_number[[chrom]] <- rep(0, simulation$parameters$cn_info$Bin_count[which(simulation$parameters$cn_info$Chromosome == chrom)])
            if (ploidy_chrom[chrom] == 0) next
            for (strand in 1:ploidy_chrom[chrom]) {
                copy_number[[chrom]] <- copy_number[[chrom]] + ploidy_block[[chrom]][[strand]]
            }
        }
        sample_genotype_unique_copy_number[[clone]] <- copy_number
    }
    #   Find total copy number of all cells in the sample
    total_copy_number <- list()
    for (clone in 1:length(sample_genotype_unique)) {
        clone_cell_count <- length(which(all_sample_genotype == sample_genotype_unique[clone]))
        for (chrom in 1:length(sample_genotype_unique_copy_number[[1]])) {
            if (clone == 1) {
                total_copy_number[[chrom]] <- clone_cell_count * sample_genotype_unique_copy_number[[clone]][[chrom]]
            } else {
                total_copy_number[[chrom]] <- total_copy_number[[chrom]] + clone_cell_count * sample_genotype_unique_copy_number[[clone]][[chrom]]
            }
        }
    }
    #   Find sample passenger mutation statistics
    sample_passenger_mutations <- phylogeny_passenger_mutations[N_sample:(2 * N_sample - 1)]
    sample_passenger_mutations <- rbindlist(lapply(sample_passenger_mutations, as.data.table))
    colnames(sample_passenger_mutations) <- c("ID", "Chromosome", "Strand", "Block", "Unit", "Position")
    mutation_bulk <- unique(sample_passenger_mutations[, c("ID", "Chromosome", "Position", "Block")])

    mutation_bulk <- mutation_bulk %>% arrange(Chromosome, Position)

    id_counts <- table(sample_passenger_mutations$ID)
    mutation_bulk$Variant_Allele_count <- id_counts[match(mutation_bulk$ID, names(id_counts))]
    mutation_bulk$Total_Allele_count <- mapply(
        function(chrom, block) total_copy_number[[chrom]][block],
        mutation_bulk$Chromosome, mutation_bulk$Block
    )
    mutation_bulk$VAF <- mutation_bulk$Variant_Allele_count / mutation_bulk$Total_Allele_count
    #--------------------------------------Store the passenger mutations
    simulation$sample_phylogeny$package_cell_phylogeny$phylogeny_passenger_mutations <- phylogeny_passenger_mutations
    simulation$sample_phylogeny$package_cell_phylogeny$mutation_bulk <- mutation_bulk
    return(simulation)
}
