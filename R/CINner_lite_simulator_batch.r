simulator_batch <- function(n_simulations = 0,
                            t_end_time = 0,
                            cell_lifespan = cell_lifespan,
                            t_tau_step = NA,
                            n_selective_clones = 1,
                            vec_time_points_s_mut = c(0),
                            vec_hierarchy_s_mut = c(0),
                            expected_end_population = 10^3,
                            vec_expected_percent_select = c(0.5, 0.5),
                            n_sample = 10^2,
                            range_population = c(0, Inf),
                            range_clonal_perc = c(0, Inf),
                            mindiff_clonal_perc = 0,
                            ploidy = 1,
                            purity = 1,
                            truncal_mutations = 0,
                            choice_theta = "constant",
                            vec_theta_parameters = c(0.6, 0.6),
                            vec_theta_mean = c(0.6, 0.6),
                            save_rda = FALSE,
                            save_true_mutation_table = TRUE,
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
                            file_prefix = NA,
                            R_libPaths = NULL) {
    if (is.na(t_tau_step)) t_tau_step <- cell_lifespan / 2
    #------------------------------Make table of time point combinations
    if (length(vec_time_points_s_mut) != length(vec_hierarchy_s_mut)) {
        print("ERROR: lengths of vec_time_points_s_mut and vec_hierarchy_s_mut are not correct")
    }
    if (any(vec_hierarchy_s_mut >= (1:length(vec_hierarchy_s_mut)))) {
        print("ERROR: vec_hierarchy_s_mut is impossible")
    }
    table_timepoint_index <- t(combn(
        1:length(vec_time_points_s_mut),
        n_selective_clones
    ))
    n_batches <- nrow(table_timepoint_index)
    #-----------------Make batches of simulations for each parameter set
    #   Initialize table of parameters for each batch of simulations
    table_parameters <- data.frame(
        Batch_ID = rep(NA, n_batches)
    )
    for (clone in 1:(n_selective_clones + 1)) {
        table_parameters[[paste("Clone ", clone - 1, " birth time (days)", sep = "")]] <- rep(NA, n_batches)
        table_parameters[[paste("Clone ", clone - 1, " mean lifespan (days)", sep = "")]] <- rep(NA, n_batches)
        table_parameters[[paste("Clone ", clone - 1, " growth rate (1/day)", sep = "")]] <- rep(NA, n_batches)
        table_parameters[[paste("Clone ", clone - 1, " division probability", sep = "")]] <- rep(NA, n_batches)
        table_parameters[[paste("Clone ", clone - 1, " doubling time (weeks)", sep = "")]] <- rep(NA, n_batches)
    }
    #   Make one batch of simulations for each time point combination
    for (i_batch in 1:n_batches) {
        vec_timepoint_index <- table_timepoint_index[i_batch, ]
        #   Set up parameters for the simulations
        vec_s_mut_time <- vec_time_points_s_mut[vec_timepoint_index]
        vec_s_mut_hierarchy <- vec_hierarchy_s_mut
        vec_s_mut_time_tmp <- c(0, vec_s_mut_time)

        vec_cell_lifespan <- rep(cell_lifespan, length = (n_selective_clones + 1))
        vec_division_probability <- rep(0, length = (n_selective_clones + 1))

        for (clone in 1:(n_selective_clones + 1)) {
            codename <- "ClonalTimes="
            if (length(vec_s_mut_time) > 0) codename <- paste(codename, vec_s_mut_time[1], sep = "")
            if (n_selective_clones > 1) {
                for (i in 2:length(vec_s_mut_time)) {
                    codename <- paste(codename, ",", vec_s_mut_time[i], sep = "")
                }
            }
            codename <- paste(codename, "_ClonalHierarchy=", vec_s_mut_hierarchy[1], sep = "")
            if (n_selective_clones > 1) {
                for (i in 2:length(vec_s_mut_hierarchy)) {
                    codename <- paste(codename, ",", vec_s_mut_hierarchy[i], sep = "")
                }
            }
            table_parameters[["Batch_ID"]][i_batch] <- codename
            table_parameters[[paste("Clone ", clone - 1, " birth time (days)", sep = "")]][i_batch] <- vec_s_mut_time_tmp[clone]
            table_parameters[[paste("Clone ", clone - 1, " mean lifespan (days)", sep = "")]][i_batch] <- cell_lifespan
            growth_rate <- log(vec_expected_percent_select[clone] * expected_end_population) / (t_end_time - vec_s_mut_time_tmp[clone])
            table_parameters[[paste("Clone ", clone - 1, " growth rate (1/day)", sep = "")]][i_batch] <- growth_rate
            division_probability <- (growth_rate * cell_lifespan + 1) / 2
            vec_division_probability[clone] <- division_probability
            if (division_probability > 1) stop("Error: division_probability cannot be greater than 1")
            table_parameters[[paste("Clone ", clone - 1, " division probability", sep = "")]][i_batch] <- division_probability
            doubling_time_in_weeks <- log(2) / growth_rate / 7
            table_parameters[[paste("Clone ", clone - 1, " doubling time (weeks)", sep = "")]][i_batch] <- doubling_time_in_weeks
        }
        #   Make the batch of simulations for the corresponding parameter set
        if (is.na(file_prefix)) {
            file_prefix <- paste(codename, "_simulated_SFS_", sep = "")
        }
        cat(paste("\nMaking simulations for batch ", i_batch, "/", n_batches, "...\n", sep = ""))
        simulator_full(
            file_prefix = file_prefix,
            n_simulations = n_simulations,
            n_selective_clones = n_selective_clones,
            t_end_time = t_end_time,
            vec_s_mut_time = vec_s_mut_time,
            vec_s_mut_hierarchy = vec_s_mut_hierarchy,
            vec_cell_lifespan = vec_cell_lifespan,
            vec_division_probability = vec_division_probability,
            range_clonal_perc = range_clonal_perc,
            mindiff_clonal_perc = mindiff_clonal_perc,
            range_population = range_population,
            ploidy = ploidy,
            purity = purity,
            truncal_mutations = truncal_mutations,
            choice_theta = choice_theta,
            vec_theta_parameters = vec_theta_parameters,
            t_tau_step = t_tau_step,
            n_sample = n_sample,
            save_rda = save_rda,
            save_true_mutation_table = save_true_mutation_table,
            output_bulk = output_bulk,
            output_sc = output_sc,
            compute_parallel = compute_parallel,
            n_cores = n_cores,
            sc_rate_false_positive = sc_rate_false_positive,
            sc_rate_false_negative = sc_rate_false_negative,
            sc_rate_unknown = sc_rate_unknown,
            bulk_coverage_model = bulk_coverage_model,
            bulk_coverage_variables = bulk_coverage_variables,
            bulk_min_alt_readcounts = bulk_min_alt_readcounts,
            subfolder = subfolder,
            R_libPaths = R_libPaths
        )
    }
    write.csv(table_parameters, paste(subfolder, "table_parameters_",
        n_selective_clones, ".csv",
        sep = ""
    ), row.names = FALSE)
    return(table_parameters)
}

simulator_batch_old <- function(n_simulations = 0,
                                t_end_time = 0,
                                t_tau_step = 1,
                                n_selective_clones = 1,
                                vec_time_points_s_mut = c(0),
                                vec_hierarchy_s_mut = c(0),
                                expected_end_population = 10^3,
                                vec_expected_percent_select = c(0.5, 0.5),
                                death_rate = 0.01,
                                n_sample = 10^2,
                                range_population = c(0, Inf),
                                range_clonal_perc = c(0, Inf),
                                mindiff_clonal_perc = 0,
                                ploidy = 1,
                                purity = 1,
                                truncal_mutations = 0,
                                choice_theta = "constant",
                                vec_theta_parameters = c(0.6, 0.6),
                                vec_theta_mean = c(0.6, 0.6),
                                save_rda = FALSE,
                                save_true_mutation_table = TRUE,
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
                                file_prefix = NA,
                                R_libPaths = NULL) {
    #------------------------------Make table of time point combinations
    if (length(vec_time_points_s_mut) != length(vec_hierarchy_s_mut)) {
        print("ERROR: lengths of vec_time_points_s_mut and vec_hierarchy_s_mut are not correct")
    }
    if (any(vec_hierarchy_s_mut >= (1:length(vec_hierarchy_s_mut)))) {
        print("ERROR: vec_hierarchy_s_mut is impossible")
    }
    table_timepoint_index <- t(combn(
        1:length(vec_time_points_s_mut),
        n_selective_clones
    ))
    n_batches <- nrow(table_timepoint_index)
    #-----------------Make batches of simulations for each parameter set
    #   Initialize table of parameters for each batch of simulations
    table_parameters_codename <- rep("", length = n_batches)

    table_parameters_birthrate <- matrix(0,
        ncol = (n_selective_clones + 1), nrow = n_batches
    )
    table_parameters_deathrate <- matrix(0,
        ncol = (n_selective_clones + 1), nrow = n_batches
    )
    table_parameters_growthrate <- matrix(0,
        ncol = (n_selective_clones + 1), nrow = n_batches
    )
    table_parameters_doublingtime <- matrix(0,
        ncol = (n_selective_clones + 1), nrow = n_batches
    )
    table_parameters_mutationcount <- matrix(0,
        ncol = (n_selective_clones + 1), nrow = n_batches
    )
    table_parameters_birthtime <- matrix(0,
        ncol = n_selective_clones, nrow = n_batches
    )
    #   Make one batch of simulations for each time point combination
    for (i_batch in 1:n_batches) {
        vec_timepoint_index <- table_timepoint_index[i_batch, ]
        #   Set up parameters for the simulations
        vec_s_mut_time <- vec_time_points_s_mut[vec_timepoint_index]
        # vec_s_mut_time <- c(0, vec_time_points_s_mut[vec_timepoint_index])
        vec_s_mut_hierarchy <- vec_hierarchy_s_mut
        vec_d_rate <- rep(death_rate, length = (n_selective_clones + 1))
        vec_b_rate <- rep(0, length = (n_selective_clones + 1))
        vec_r_rate <- rep(0, length = (n_selective_clones + 1))
        vec_propensity <- matrix(0, nrow = 2, ncol = (n_selective_clones + 1))
        vec_propensity_sum <- rep(0, length = (n_selective_clones + 1))
        vec_s_mut_time_tmp <- c(0, vec_s_mut_time)
        for (i in 1:(n_selective_clones + 1)) {
            vec_r_rate[i] <- log(vec_expected_percent_select[i] *
                expected_end_population) / (t_end_time - vec_s_mut_time_tmp[i])
            vec_b_rate[i] <- vec_d_rate[i] + vec_r_rate[i]
            vec_propensity[1, i] <- vec_b_rate[i]
            vec_propensity[2, i] <- vec_d_rate[i]
            vec_propensity_sum[i] <- vec_b_rate[i] + vec_d_rate[i]
        }
        #   Update table of parameters
        codename <- "ClonalTimes="
        if (length(vec_s_mut_time) > 0) codename <- paste(codename, vec_s_mut_time[1], sep = "")
        if (n_selective_clones > 1) {
            for (i in 2:length(vec_s_mut_time)) {
                codename <- paste(codename, ",", vec_s_mut_time[i], sep = "")
            }
        }
        codename <- paste(codename, "_ClonalHierarchy=", vec_s_mut_hierarchy[1], sep = "")
        if (n_selective_clones > 1) {
            for (i in 2:length(vec_s_mut_hierarchy)) {
                codename <- paste(codename, ",", vec_s_mut_hierarchy[i], sep = "")
            }
        }

        table_parameters_codename[i_batch] <- codename

        table_parameters_birthrate[i_batch, ] <- vec_b_rate
        table_parameters_deathrate[i_batch, ] <- vec_d_rate
        table_parameters_growthrate[i_batch, ] <- vec_r_rate
        table_parameters_doublingtime[i_batch, ] <- log(2) / vec_r_rate
        table_parameters_mutationcount[i_batch, ] <- vec_theta_mean / vec_b_rate

        table_parameters_birthtime[i_batch, ] <-
            vec_time_points_s_mut[vec_timepoint_index]
        #   Make the batch of simulations for the corresponding parameter set
        if (is.na(file_prefix)) {
            file_prefix <- paste(codename, "_simulated_SFS_", sep = "")
        }
        cat(paste("\nMaking simulations for batch ", i_batch, "/", n_batches, "...\n", sep = ""))
        simulator_full_old(
            file_prefix = file_prefix,
            n_simulations = n_simulations,
            n_selective_clones = n_selective_clones,
            t_end_time = t_end_time,
            vec_s_mut_time = vec_s_mut_time,
            vec_s_mut_hierarchy = vec_s_mut_hierarchy,
            vec_propensity = vec_propensity,
            vec_propensity_sum = vec_propensity_sum,
            range_clonal_perc = range_clonal_perc,
            mindiff_clonal_perc = mindiff_clonal_perc,
            range_population = range_population,
            ploidy = ploidy,
            purity = purity,
            truncal_mutations = truncal_mutations,
            choice_theta = choice_theta,
            vec_theta_parameters = vec_theta_parameters,
            t_tau_step = t_tau_step,
            n_sample = n_sample,
            save_rda = save_rda,
            save_true_mutation_table = save_true_mutation_table,
            output_bulk = output_bulk,
            output_sc = output_sc,
            compute_parallel = compute_parallel,
            n_cores = n_cores,
            sc_rate_false_positive = sc_rate_false_positive,
            sc_rate_false_negative = sc_rate_false_negative,
            sc_rate_unknown = sc_rate_unknown,
            bulk_coverage_model = bulk_coverage_model,
            bulk_coverage_variables = bulk_coverage_variables,
            bulk_min_alt_readcounts = bulk_min_alt_readcounts,
            subfolder = subfolder,
            R_libPaths = R_libPaths
        )
    }
    #-----------------------------------------Output table of parameters
    table_parameters <- data.frame(Code_Name = table_parameters_codename)
    for (clone in 1:(n_selective_clones + 1)) {
        table_parameters[paste("Clone_", clone - 1, "_Birth_rate", sep = "")] <-
            table_parameters_birthrate[, clone]
        table_parameters[paste("Clone_", clone - 1, "_Death_rate", sep = "")] <-
            table_parameters_deathrate[, clone]
        table_parameters[paste("Clone_", clone - 1, "_Growth_rate", sep = "")] <-
            table_parameters_growthrate[, clone]
        table_parameters[paste("Clone_", clone - 1, "_Doubling_time_weeks", sep = "")] <-
            table_parameters_doublingtime[, clone]
        table_parameters[paste("Clone_", clone - 1, "_Average_mutation_after_division", sep = "")] <-
            table_parameters_mutationcount[, clone]
        if (clone > 1) {
            table_parameters[paste("Clone_", clone - 1, "_Birth_time", sep = "")] <-
                table_parameters_birthtime[, (clone - 1)]
        }
    }
    write.csv(table_parameters, paste(subfolder, "table_parameters_",
        n_selective_clones, ".csv",
        sep = ""
    ), row.names = FALSE)
}
