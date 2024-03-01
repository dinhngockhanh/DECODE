CINner_postprocessing <- function(model_name,
                                  folder_workplace,
                                  n_simulations,
                                  passenger_mutation_rate,
                                  bulk_coverage_model,
                                  bulk_coverage_variables,
                                  bulk_min_alt_readcounts,
                                  compute_parallel = FALSE,
                                  n_cores = NULL) {
    if (compute_parallel == FALSE) {
        library(dplyr)
        library(GenomicRanges)
        #-------------------------Process simulations in sequential mode
        for (i_simulation in 1:n_simulations) {
            CINner_postprocessing_one_simulation(
                i_simulation = i_simulation,
                model_name = model_name,
                folder_workplace = folder_workplace,
                n_simulations = n_simulations,
                passenger_mutation_rate = passenger_mutation_rate,
                bulk_coverage_model = bulk_coverage_model,
                bulk_coverage_variables = bulk_coverage_variables,
                bulk_min_alt_readcounts = bulk_min_alt_readcounts,
                compute_parallel = compute_parallel
            )
        }
    } else {
        library(parallel)
        library(pbapply)
        #---------------------------Process simulations in parallel mode
        #   Start parallel cluster
        if (is.null(n_cores)) {
            numCores <- detectCores()
        } else {
            numCores <- n_cores
        }
        cl <- makePSOCKcluster(numCores - 1)
        #   Prepare input parameters
        CINner_postprocessing_one_simulation <<- CINner_postprocessing_one_simulation
        CINner_get_event_times <<- CINner_get_event_times
        CINner_add_passenger_mutations <<- CINner_add_passenger_mutations
        CINner_get_cn_segments <<- CINner_get_cn_segments
        CINner_reduce_to_clonal <<- CINner_reduce_to_clonal
        rand_coverage <<- rand_coverage
        folder_workplace <<- folder_workplace
        model_name <<- model_name
        passenger_mutation_rate <<- passenger_mutation_rate
        bulk_coverage_model <<- bulk_coverage_model
        bulk_coverage_variables <<- bulk_coverage_variables
        bulk_min_alt_readcounts <<- bulk_min_alt_readcounts
        clusterExport(cl, varlist = c(
            "CINner_postprocessing_one_simulation",
            "CINner_get_event_times",
            "CINner_add_passenger_mutations",
            "CINner_get_cn_segments",
            "CINner_reduce_to_clonal",
            "rand_coverage",
            "folder_workplace",
            "model_name",
            "passenger_mutation_rate",
            "bulk_coverage_model",
            "bulk_coverage_variables",
            "bulk_min_alt_readcounts"
        ))
        clusterEvalQ(cl, library(dplyr))
        clusterEvalQ(cl, library(GenomicRanges))
        #   Process simulations in parallel
        pblapply(cl = cl, X = 1:n_simulations, FUN = function(i_simulation) {
            CINner_postprocessing_one_simulation(
                i_simulation = i_simulation,
                model_name = model_name,
                folder_workplace = folder_workplace,
                n_simulations = n_simulations,
                passenger_mutation_rate = passenger_mutation_rate,
                bulk_coverage_model = bulk_coverage_model,
                bulk_coverage_variables = bulk_coverage_variables,
                bulk_min_alt_readcounts = bulk_min_alt_readcounts,
                compute_parallel = compute_parallel
            )
            return("")
        })
    }
}

CINner_postprocessing_one_simulation <- function(i_simulation,
                                                 model_name,
                                                 folder_workplace,
                                                 n_simulations,
                                                 passenger_mutation_rate,
                                                 bulk_coverage_model,
                                                 bulk_coverage_variables,
                                                 bulk_min_alt_readcounts,
                                                 compute_parallel) {
    cat("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
    cat(paste0("MODIFYING SIMULATION-", i_simulation, "...\n"))
    #   Load the simulation output
    load(paste0(folder_workplace, "/", model_name, "_simulation_", i_simulation, ".rda"))
    #   Get ground-truth CNA timing
    cat("\nGet ground-truth CNA timing...\n")
    simulation <- CINner_get_event_times(simulation)
    #   Add passenger mutations to the simulation
    cat("\nSimulate passenger mutations...\n")
    simulation <- CINner_add_passenger_mutations(simulation, passenger_mutation_rate)
    #   Get clonal & subclonal copy number segments
    cat("\nGet clonal & subclonal segmented CN profiles...\n")
    simulation <- CINner_get_cn_segments(simulation)
    #   Reduce mutation & copy number data to only clonal region
    cat("\nReduce passenger mutations & CN profiles to only clonal regions...\n")
    simulation <- CINner_reduce_to_clonal(simulation)
    # #   Save the simulation output
    # cat("Save simulation package...\n")
    # save(simulation, file = paste0(folder_workplace, "/", model_name, "_simulation_", i_simulation, ".rda"))
    #   Save the ground-truth CNA timing as a CSV file
    event_details <- simulation$sample$event_details
    event_filename <- paste0(folder_workplace, "/", model_name, "_simulation_", i_simulation, "_true_timing.csv")
    write.csv(event_details, file = event_filename, row.names = FALSE)
    #   Save the copy number segments as a CSV file
    cat("Save copy number segmentation...\n")
    segmented_copynumber <- simulation$sample$segmented_copynumber
    cn_filename <- paste0(folder_workplace, "/", model_name, "_simulation_", i_simulation, "_bulk_copynumber_data.csv")
    write.csv(segmented_copynumber, file = cn_filename, row.names = FALSE)
    #   Save the copy number segments as a Granges file
    cn_chr <- segmented_copynumber$chr
    cn_start <- as.integer(segmented_copynumber$start)
    cn_end <- as.integer(segmented_copynumber$end)
    cn_major_cn <- as.integer(segmented_copynumber$Maj)
    cn_minor_cn <- as.integer(segmented_copynumber$Min)
    cn_clonal_frequency <- segmented_copynumber$Clonal_frequency
    cn_granges <- GRanges(
        seqnames = Rle(cn_chr),
        ranges = IRanges(start = cn_start, end = cn_end),
        strand = Rle(strand(rep("*", length(cn_chr)))),
        major_cn = cn_major_cn,
        minor_cn = cn_minor_cn,
        clonal_frequency = cn_clonal_frequency
    )
    cn_filename <- paste0(folder_workplace, "/", model_name, "_simulation_", i_simulation, "_cn.RData")
    save(cn_granges, file = cn_filename)
    #   Simulate mutations from bulk sequencing data
    cat("Simulate mutations from bulk sequencing...\n")
    mutation_bulk <- simulation$sample_phylogeny$package_cell_phylogeny$mutation_bulk
    mutation_bulk$readcount_tot <- rand_coverage(
        bulk_coverage_model = bulk_coverage_model,
        bulk_coverage_variables = bulk_coverage_variables,
        n_mutation_true = nrow(mutation_bulk)
    )
    mutation_bulk$readcount_alt <- rbinom(
        n = nrow(mutation_bulk),
        size = mutation_bulk$readcount_tot,
        prob = mutation_bulk$VAF
    )
    vec_delete <- which(mutation_bulk$readcount_alt <= bulk_min_alt_readcounts)
    if (length(vec_delete) > 0) mutation_bulk <- mutation_bulk[-vec_delete, ]
    mutation_bulk$readcount_ref <- mutation_bulk$readcount_tot - mutation_bulk$readcount_alt
    #   Save the mutations from bulk sequencing data as a CSV file
    cat("Save mutations from bulk sequencing...\n")
    mut_filename <- paste0(folder_workplace, "/", model_name, "_simulation_", i_simulation, "_bulk_mutational_data.csv")
    write.csv(mutation_bulk, file = mut_filename, row.names = FALSE)
    #   Save the mutations from bulk sequencing data as a VCF file
    mut_chr <- mutation_bulk$Chromosome
    mut_pos <- mutation_bulk$Position
    mut_alt_count <- mutation_bulk$readcount_alt
    mut_ref_count <- mutation_bulk$readcount_ref
    vcf <- c(
        "##fileformat=VCFv4.2",
        "##INFO=<ID=t_alt_count,Number=1,Type=Integer,Description=\"Tumor alternate count\">",
        "##INFO=<ID=t_ref_count,Number=1,Type=Integer,Description=\"Tumor reference count\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    )
    vcf <- c(vcf, paste0(mut_chr, "\t", mut_pos, "\t.\tN\tA\tNA\tNA\tt_alt_count=", mut_alt_count, ";t_ref_count=", mut_ref_count))
    mut_filename <- paste0(folder_workplace, "/", model_name, "_simulation_", i_simulation, "_mut.vcf")
    writeLines(vcf, mut_filename)
}
