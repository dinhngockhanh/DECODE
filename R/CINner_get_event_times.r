CINner_get_event_times <- function(simulation) {
    #----------------------------------------Input the general variables
    size_CN_block_DNA <- as.numeric(simulation$parameters$general_variables$Value[which(simulation$parameters$general_variables$Variable == "size_CN_block_DNA")])
    cn_info <- simulation$parameters$cn_info
    #-----------------------------------------Input the clonal evolution
    T_final <- simulation$clonal_evolution$T_current
    evolution_traj_time <- simulation$clonal_evolution$evolution_traj_time
    evolution_traj_clonal_ID <- simulation$clonal_evolution$evolution_traj_clonal_ID
    evolution_origin <- simulation$clonal_evolution$evolution_origin
    evolution_genotype_changes <- simulation$clonal_evolution$evolution_genotype_changes
    genotype_list_ploidy_chrom <- simulation$clonal_evolution$genotype_list_ploidy_chrom
    genotype_list_ploidy_block <- simulation$clonal_evolution$genotype_list_ploidy_block
    #-----------------------------------------Input the sample phylogeny
    phylogeny_origin <- simulation$sample_phylogeny$package_cell_phylogeny$phylogeny_origin
    phylogeny_elapsed_genotypes <- simulation$sample_phylogeny$package_cell_phylogeny$phylogeny_elapsed_genotypes
    phylogeny_genotype <- simulation$sample_phylogeny$package_cell_phylogeny$phylogeny_genotype
    #---------------------------------------------------Input the sample
    N_sample <- length(simulation$sample$sample_cell_ID)
    #--------------Get list of genotypes in history of each sampled cell
    phylogeny_tip_counts <- rep(0, length(phylogeny_origin))
    for (tip in 1:N_sample) {
        node <- tip + N_sample - 1
        while (node > 0) {
            phylogeny_tip_counts[node] <- phylogeny_tip_counts[node] + 1
            node <- phylogeny_origin[node]
        }
    }
    #-----------------------------Function to find birth time of a clone
    find_clonal_birthtime <- function(genotype) {
        birthtime <- NA
        i <- 0
        while (is.na(birthtime)) {
            i <- i + 1
            if (genotype %in% evolution_traj_clonal_ID[[i]]) birthtime <- evolution_traj_time[i]
        }
        return(birthtime)
    }
    #-----Function to find mean CN of a given region in a given genotype
    find_mean_CN <- function(ploidy_chrom, ploidy_block, spec_chrom, spec_bin_start, spec_bin_end) {
        spec_chrom <- as.numeric(spec_chrom)
        vec_CN <- rep(0, spec_bin_end - spec_bin_start + 1)
        ploidy <- ploidy_chrom[spec_chrom]
        if (ploidy > 0) {
            for (strand in 1:ploidy) {
                vec_CN <- vec_CN + ploidy_block[[spec_chrom]][[strand]][spec_bin_start:spec_bin_end]
            }
        }
        mean_CN <- mean(vec_CN)
        return(mean_CN)
    }
    #----------------Find the time of each CNA event in sample phylogeny
    all_event_df <- c()
    for (branch in 1:length(phylogeny_origin)) {
        # ==Find how many sampled cells descend from this branch
        branch_tip_count <- phylogeny_tip_counts[branch]
        Clonal_fraction <- branch_tip_count / N_sample
        if (Clonal_fraction == 1) {
            Clonality <- "clonal"
        } else {
            Clonality <- "subclonal"
        }
        # ==Get mother genotype with existing mutations
        mother_branch <- phylogeny_origin[branch]
        if (mother_branch <= 0) {
            mother_genotype <- 0
        } else {
            mother_genotype <- phylogeny_genotype[mother_branch]
        }
        # ==Get list of elapsed genotypes in original phylogeny
        elapsed_genotypes <- unique(phylogeny_elapsed_genotypes[[branch]])
        # ==Get information about new clonal CNA events
        for (elapsed_gen in 1:length(elapsed_genotypes)) {
            current_genotype <- elapsed_genotypes[elapsed_gen]
            #---Find next genotype
            if (elapsed_gen == 1) {
                previous_genotype <- mother_genotype
            } else {
                previous_genotype <- evolution_origin[current_genotype]
                # elapsed_genotypes[elapsed_gen - 1]
            }
            #---Record information about CNA events causing next genotype
            if (current_genotype != previous_genotype) {
                #   Find the birth time of the current genotype
                current_genotype_birthtime <- find_clonal_birthtime(current_genotype)
                Time <- current_genotype_birthtime / T_final
                #   Record information about new CNAs in the current genotype
                events <- evolution_genotype_changes[[current_genotype]]
                if (!is.null(events)) {
                    for (j in 1:length(events)) {
                        event <- events[[j]]
                        Event <- event[1]
                        if (Event == "new-driver") {
                            Chromosome <- NA
                            Type <- NA
                            Start <- NA
                            End <- NA
                            CN_parent <- NA
                            CN_child <- NA
                        } else if (Event == "whole-genome-duplication") {
                            Chromosome <- NA
                            Type <- "+"
                            Start <- NA
                            End <- NA
                            CN_parent <- NA
                            CN_child <- NA
                        } else if (Event == "missegregation") {
                            Chromosome <- event[2]
                            Type <- ifelse(event[4] == 1, "+", "-")
                            Start <- 1
                            End <- cn_info$Bin_count[which(cn_info$Chromosome == event[2])]
                            Event <- paste0(Event, ":", Type, Chromosome)
                            CN_parent <- find_mean_CN(genotype_list_ploidy_chrom[[previous_genotype]], genotype_list_ploidy_block[[previous_genotype]], Chromosome, Start, End)
                            CN_child <- find_mean_CN(genotype_list_ploidy_chrom[[current_genotype]], genotype_list_ploidy_block[[current_genotype]], Chromosome, Start, End)
                            Start <- size_CN_block_DNA * (Start - 1) + 1
                            End <- size_CN_block_DNA * End
                        } else if (Event == "chromosome-arm-missegregation") {
                            Chromosome <- event[2]
                            Type <- ifelse(event[5] == 1, "+", "-")
                            if (event[4] == 1) {
                                Start <- 1
                                End <- cn_info$Centromere_location[which(cn_info$Chromosome == event[2])]
                                Event <- paste0(Event, ":", Type, Chromosome, "p")
                            } else if (event[4] == 2) {
                                Start <- cn_info$Centromere_location[which(cn_info$Chromosome == event[2])] + 1
                                End <- cn_info$Bin_count[which(cn_info$Chromosome == event[2])]
                                Event <- paste0(Event, ":", Type, Chromosome, "q")
                            }
                            CN_parent <- find_mean_CN(genotype_list_ploidy_chrom[[previous_genotype]], genotype_list_ploidy_block[[previous_genotype]], Chromosome, Start, End)
                            CN_child <- find_mean_CN(genotype_list_ploidy_chrom[[current_genotype]], genotype_list_ploidy_block[[current_genotype]], Chromosome, Start, End)
                            Start <- size_CN_block_DNA * (Start - 1) + 1
                            End <- size_CN_block_DNA * End
                        } else if (event[1] == "focal-amplification") {
                            Chromosome <- event[2]
                            Type <- "+"
                            Start <- event[4]
                            End <- event[5]
                            CN_parent <- find_mean_CN(genotype_list_ploidy_chrom[[previous_genotype]], genotype_list_ploidy_block[[previous_genotype]], Chromosome, Start, End)
                            CN_child <- find_mean_CN(genotype_list_ploidy_chrom[[current_genotype]], genotype_list_ploidy_block[[current_genotype]], Chromosome, Start, End)
                            Start <- size_CN_block_DNA * (Start - 1) + 1
                            End <- size_CN_block_DNA * End
                            Event <- paste0(Event, ":", Type, Chromosome, "(", Start, "-", End, ")")
                        } else if (event[1] == "focal-deletion") {
                            Chromosome <- event[2]
                            Type <- "-"
                            Start <- event[4]
                            End <- event[5]
                            CN_parent <- find_mean_CN(genotype_list_ploidy_chrom[[previous_genotype]], genotype_list_ploidy_block[[previous_genotype]], Chromosome, Start, End)
                            CN_child <- find_mean_CN(genotype_list_ploidy_chrom[[current_genotype]], genotype_list_ploidy_block[[current_genotype]], Chromosome, Start, End)
                            Start <- size_CN_block_DNA * (Start - 1) + 1
                            End <- size_CN_block_DNA * End
                            Event <- paste0(Event, ":", Type, Chromosome, "(", Start, "-", End, ")")
                        } else if (event[1] == "cnloh-interstitial") {
                            #   ...
                            #   ...
                            #   ...
                        } else if (event[1] == "cnloh-terminal") {
                            #   ...
                            #   ...
                            #   ...
                        }
                        event_df <- data.frame(
                            Time = Time,
                            Clonality = Clonality,
                            Clonal_fraction = Clonal_fraction,
                            Event = Event,
                            Chromosome = Chromosome,
                            Start = Start,
                            End = End,
                            Type = Type,
                            CN_parent = CN_parent,
                            CN_child = CN_child
                        )
                        if (is.null(all_event_df)) {
                            all_event_df <- event_df
                        } else {
                            all_event_df <- rbind(all_event_df, event_df)
                        }
                    }
                }
            }
        }
    }
    #-----------------------------------Store the ground-truth CNA times
    simulation$sample$event_details <- all_event_df
    return(simulation)
}
