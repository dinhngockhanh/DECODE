simulator_full <- function(file_prefix = "",
                           n_simulations = 0,
                           n_selective_clones = 1,
                           t_end_time = 0,
                           t_tau_step = 1,
                           vec_s_mut_time = 0,
                           vec_s_mut_hierarchy = 0,
                           vec_propensity = c(0, 0),
                           vec_propensity_sum = c(0, 0),
                           range_clonal_perc = c(40, 60),
                           mindiff_clonal_perc = 0,
                           range_population = c(800, 1200),
                           ploidy = 1,
                           purity = 1,
                           truncal_mutations = 0,
                           choice_theta = "constant",
                           vec_theta_parameters = c(0.6, 0.6),
                           n_sample = 0,
                           save_rda = save_rda,
                           save_true_mutation_table = save_true_mutation_table,
                           output_bulk = FALSE,
                           output_sc = FALSE,
                           compute_parallel = TRUE,
                           n_cores = NULL,
                           sc_rate_false_positive = 0,
                           sc_rate_false_negative = 0,
                           sc_rate_unknown = 0,
                           bulk_coverage_model = "binomial",
                           bulk_coverage_variables = c(0, 100),
                           bulk_min_alt_readcounts = 0,
                           subfolder = "",
                           R_libPaths = NULL) {
    if (compute_parallel == FALSE) {
        #-------------------------Compute simulations in sequential mode
        pb <- txtProgressBar(
            min = 0,
            max = n_simulations,
            style = 3,
            width = 50,
            char = "="
        )
        for (n_simulation in 1:n_simulations) {
            setTxtProgressBar(pb, n_simulation)
            tmp <- simulator_one_simulation(
                n_simulation,
                file_prefix,
                n_simulations,
                n_selective_clones,
                t_end_time,
                t_tau_step,
                vec_s_mut_time,
                vec_s_mut_hierarchy,
                vec_propensity,
                vec_propensity_sum,
                range_clonal_perc,
                mindiff_clonal_perc,
                range_population,
                ploidy,
                purity,
                truncal_mutations,
                choice_theta,
                vec_theta_parameters,
                n_sample,
                save_rda,
                save_true_mutation_table,
                output_bulk,
                output_sc,
                sc_rate_false_positive,
                sc_rate_false_negative,
                sc_rate_unknown,
                bulk_coverage_model,
                bulk_coverage_variables,
                bulk_min_alt_readcounts,
                subfolder
            )
            tmp <- ""
            return(tmp)
        }
    } else {
        #---------------------------Compute simulations in parallel mode
        #   Start parallel cluster
        if (is.null(n_cores)) {
            numCores <- detectCores()
        } else {
            numCores <- n_cores
        }
        cl <- makePSOCKcluster(numCores - 1)
        if (is.null(R_libPaths) == FALSE) {
            R_libPaths <<- R_libPaths
            clusterExport(cl, varlist = c("R_libPaths"))
            clusterEvalQ(cl = cl, .libPaths(R_libPaths))
        }
        #   Prepare input parameters
        simulator_one_simulation <<- simulator_one_simulation
        simulation_clonal_evolution <<- simulation_clonal_evolution
        simulation_sample_phylogeny <<- simulation_sample_phylogeny
        simulation_sequening_truth <<- simulation_sequening_truth
        simulation_sequencing_sc <<- simulation_sequencing_sc
        simulation_sequencing_bulk <<- simulation_sequencing_bulk
        rand_coverage <<- rand_coverage
        file_prefix <<- file_prefix
        n_simulations <<- n_simulations
        n_selective_clones <<- n_selective_clones
        t_end_time <<- t_end_time
        t_tau_step <<- t_tau_step
        vec_s_mut_time <<- vec_s_mut_time
        vec_s_mut_hierarchy <<- vec_s_mut_hierarchy
        vec_propensity <<- vec_propensity
        vec_propensity_sum <<- vec_propensity_sum
        range_clonal_perc <<- range_clonal_perc
        mindiff_clonal_perc <<- mindiff_clonal_perc
        range_population <<- range_population
        ploidy <<- ploidy
        purity <<- purity
        truncal_mutations <<- truncal_mutations
        choice_theta <<- choice_theta
        vec_theta_parameters <<- vec_theta_parameters
        n_sample <<- n_sample
        output_bulk <<- output_bulk
        output_sc <<- output_sc
        sc_rate_false_positive <<- sc_rate_false_positive
        sc_rate_false_negative <<- sc_rate_false_negative
        sc_rate_unknown <<- sc_rate_unknown
        bulk_coverage_model <<- bulk_coverage_model
        bulk_coverage_variables <<- bulk_coverage_variables
        bulk_min_alt_readcounts <<- bulk_min_alt_readcounts
        subfolder <<- subfolder
        clusterExport(cl, varlist = c(
            "simulator_one_simulation",
            "simulation_clonal_evolution",
            "simulation_sample_phylogeny",
            "simulation_sequening_truth",
            "simulation_sequencing_sc",
            "simulation_sequencing_bulk",
            "rand_coverage",
            "n_simulations",
            "file_prefix",
            "n_simulations",
            "n_selective_clones",
            "t_end_time",
            "t_tau_step",
            "vec_s_mut_time",
            "vec_s_mut_hierarchy",
            "vec_propensity",
            "vec_propensity_sum",
            "range_clonal_perc",
            "mindiff_clonal_perc",
            "range_population",
            "ploidy",
            "purity",
            "truncal_mutations",
            "choice_theta",
            "vec_theta_parameters",
            "n_sample",
            "output_bulk",
            "output_sc",
            "sc_rate_false_positive",
            "sc_rate_false_negative",
            "sc_rate_unknown",
            "bulk_coverage_model",
            "bulk_coverage_variables",
            "bulk_min_alt_readcounts",
            "subfolder"
        ))
        clusterEvalQ(cl, library(data.table))
        #   Compute simulations in parallel
        pblapply(cl = cl, X = 1:n_simulations, FUN = function(n_simulation) {
            simulator_one_simulation(
                n_simulation,
                file_prefix,
                n_simulations,
                n_selective_clones,
                t_end_time,
                t_tau_step,
                vec_s_mut_time,
                vec_s_mut_hierarchy,
                vec_propensity,
                vec_propensity_sum,
                range_clonal_perc,
                mindiff_clonal_perc,
                range_population,
                ploidy,
                purity,
                truncal_mutations,
                choice_theta,
                vec_theta_parameters,
                n_sample,
                save_rda,
                save_true_mutation_table,
                output_bulk,
                output_sc,
                sc_rate_false_positive,
                sc_rate_false_negative,
                sc_rate_unknown,
                bulk_coverage_model,
                bulk_coverage_variables,
                bulk_min_alt_readcounts,
                subfolder
            )
            return("")
        })
    }
}

simulator_one_simulation <- function(n_simulation,
                                     file_prefix,
                                     n_simulations,
                                     n_selective_clones,
                                     t_end_time,
                                     t_tau_step,
                                     vec_s_mut_time,
                                     vec_s_mut_hierarchy,
                                     vec_propensity,
                                     vec_propensity_sum,
                                     range_clonal_perc,
                                     mindiff_clonal_perc,
                                     range_population,
                                     ploidy,
                                     purity,
                                     truncal_mutations,
                                     choice_theta,
                                     vec_theta_parameters,
                                     n_sample,
                                     save_rda,
                                     save_true_mutation_table,
                                     output_bulk,
                                     output_sc,
                                     sc_rate_false_positive,
                                     sc_rate_false_negative,
                                     sc_rate_unknown,
                                     bulk_coverage_model,
                                     bulk_coverage_variables,
                                     bulk_min_alt_readcounts,
                                     subfolder) {
    #--------------------------------------Simulate the clonal evolution
    conditioning <- 0
    print("")
    print("Simulating clonal evolution")
    while (min(conditioning) == 0) {
        output <- simulation_clonal_evolution(
            t_end_time = t_end_time,
            n_selective_clones = n_selective_clones,
            vec_s_mut_time = vec_s_mut_time,
            vec_s_mut_hierarchy = vec_s_mut_hierarchy,
            vec_propensity = vec_propensity,
            vec_propensity_sum = vec_propensity_sum,
            range_clonal_perc = range_clonal_perc,
            mindiff_clonal_perc = mindiff_clonal_perc,
            range_population = range_population,
            t_tau_step = t_tau_step
        )
        conditioning <- output$conditioning
    }
    simulation <- output$simulation
    #------------Simulate the sample phylogeny from the clonal evolution
    print("Simulating sample phylogeny")
    start_time <- Sys.time()
    simulation <- simulation_sample_phylogeny(
        simulation = simulation,
        t_end_time = t_end_time,
        vec_s_mut_time = vec_s_mut_time,
        vec_s_mut_hierarchy = vec_s_mut_hierarchy,
        n_selective_clones = n_selective_clones,
        n_sample = n_sample
    )
    end_time <- Sys.time()
    print(end_time - start_time)
    start_time <- Sys.time()
    #   Output the cell phylogeny in hclust style
    filename <- paste(subfolder, file_prefix, n_simulation, "_phylogeny_hclust.rda", sep = "")
    cell_phylogeny_hclust <- simulation$cell_phylogeny_hclust
    save(cell_phylogeny_hclust, file = filename)
    filename <- paste(subfolder, file_prefix, n_simulation, "_phylogeny_clustering.rda", sep = "")
    cell_phylogeny_clustering <- simulation$cell_phylogeny_clustering
    save(cell_phylogeny_clustering, file = filename)
    end_time <- Sys.time()
    print(end_time - start_time)
    #---------------------------------Simulate the true mutational table
    print("Simulating true mutational table")
    start_time <- Sys.time()
    simulation <- simulation_sequening_truth(
        simulation = simulation,
        choice_theta = choice_theta,
        truncal_mutations = truncal_mutations,
        vec_theta_parameters = vec_theta_parameters,
        n_sample = n_sample,
        n_selective_clones = n_selective_clones
    )
    end_time <- Sys.time()
    print(end_time - start_time)
    #   Get ingredients to build true mutational table
    if (save_true_mutation_table == TRUE) {
        start_time <- Sys.time()
        sample_mutational_table_truth <-
            simulation$sample_mutational_table_truth
        sample_genotype <- simulation$sample_genotype
        #   Prepare first row in data frame for cell clones
        sample_mutational_table_cell_clone <- data.frame(Row = "Cell_Clone", Marker = NaN, hclust_index = NaN)
        for (node in 1:n_sample) {
            sample_mutational_table_cell_clone[paste("Cell_", node, sep = "")] <- sample_genotype[node]
        }
        #   Attach two data frames together
        sample_mutational_table_true <- rbind(sample_mutational_table_cell_clone, sample_mutational_table_truth)
        #   Save true mutational table to file
        filename <- paste(subfolder, file_prefix,
            n_simulation, "_mutational_data_TRUTH.csv",
            sep = ""
        )
        write.csv(sample_mutational_table_true, file = filename)
        end_time <- Sys.time()
        print(end_time - start_time)
    }
    #---------------------------Simulate the single-cell mutational data
    if (output_sc == TRUE) {
        print("Simulating single-cell mutational data")
        simulation <- simulation_sequencing_sc(
            simulation = simulation,
            sc_rate_false_positive = sc_rate_false_positive,
            sc_rate_false_negative = sc_rate_false_negative,
            sc_rate_unknown = sc_rate_unknown
        )
        #   Get ingredients to build single-cell mutational table
        sample_mutational_table_sc <- simulation$sample_mutational_table_sc
        sample_genotype <- simulation$sample_genotype
        #   Prepare first row in data frame for cell clones
        sample_mutational_table_cell_clone <- data.frame(Row = "Cell_Clone", Marker = NaN, hclust_index = NaN)
        for (node in 1:n_sample) {
            sample_mutational_table_cell_clone[paste("Cell_", node, sep = "")] <- sample_genotype[node]
        }
        #   Attach two data frames together
        sample_mutational_table_sc <- rbind(sample_mutational_table_cell_clone, sample_mutational_table_sc)
        #   Save single-cell mutational table to file
        filename <- paste(
            subfolder, file_prefix, n_simulation,
            "_mutational_data_SINGLE_CELL_FP=", sc_rate_false_positive,
            "_FN=", sc_rate_false_negative,
            "_NA=", sc_rate_unknown, ".csv",
            sep = ""
        )
        write.csv(sample_mutational_table_sc, file = filename)
    }
    #----------------------------------Simulate the bulk sequencing data
    if (output_bulk == TRUE) {
        print("Simulating bulk sequencing data")
        start_time <- Sys.time()
        simulation <- simulation_sequencing_bulk(
            simulation = simulation,
            n_sample = n_sample,
            ploidy = ploidy,
            purity = purity,
            bulk_coverage_model = bulk_coverage_model,
            bulk_coverage_variables = bulk_coverage_variables,
            bulk_min_alt_readcounts = bulk_min_alt_readcounts
        )
        end_time <- Sys.time()
        print(end_time - start_time)
        start_time <- Sys.time()
        #   Get bulk mutational table
        sample_mutational_table_bulk <-
            simulation$sample_mutational_table_bulk
        #   Save bulk mutational table to file
        filename <- paste(subfolder, file_prefix,
            n_simulation, "_mutational_data_BULK.csv",
            sep = ""
        )
        write.csv(sample_mutational_table_bulk,
            file = filename,
            row.names = FALSE
        )
        end_time <- Sys.time()
        print(end_time - start_time)
    }
    if (save_rda == TRUE) {
        filename <- paste0(subfolder, file_prefix, "simulation_", n_simulation, ".rda")
        save(simulation, file = filename)
    }
    #   Output the simulation variables
    Count_in_sample <- rep(0, n_selective_clones + 1)
    for (clone in 0:n_selective_clones) Count_in_sample[clone + 1] <- sum(simulation$sample_genotype == clone)
    simulation_variables <- data.frame(
        Clone_ID = paste0("Clone_", 0:n_selective_clones),
        MRCA_ages = simulation$MRCA_ages,
        Count_in_population = simulation$record_vec_populations[nrow(simulation$record_vec_populations), ],
        Count_in_sample = Count_in_sample
    )
    filename <- paste(subfolder, file_prefix, n_simulation, "_simulation_variables.csv", sep = "")
    write.csv(simulation_variables, file = filename)
    return(simulation)
}
