DECODE_ABC <- function(sample_id = "",
                       mutation_table,
                       criterion = "BIC",
                       criterion_ratio = 0.999, # <<<<<<<<<<<<<<<<<<<<<<
                       neutral_power_min = 0.5,
                       neutral_power_max = 5,
                       cluster_frequency_min = 0.01,
                       cluster_frequency_max = 1,
                       max_total_read = NULL,
                       sample_size = 1000,
                       matrix_binomial_sample_size = 1000, # <<<<<<<<<<<
                       matrix_binomial_ploidy = 1, # <<<<<<<<<<<<<<<<<<<
                       sfs_bincount = 100, # <<<<<<<<<<<<<<<<<<<<<<<<<<<
                       inference_retained_freq = 75,
                       validation_mutation_count = 5000,
                       validation_N_trials = 100, # <<<<<<<<<<<<<<<<<<<<
                       coverage_distribution = "sample-specific",
                       coverage_variables = NULL,
                       N_trials = 10000,
                       N_trials_kept = 100,
                       neutral_tail = NA,
                       min_N_humps = 1,
                       max_N_humps = Inf,
                       pi_cutoff = 0.02,
                       zero_cutoff = 1e-50,
                       compute_parallel = TRUE,
                       n_cores = NULL) {
    suppressPackageStartupMessages(library(crayon))
    cat(paste0("\n\n\n", bold(red("PERFORMING DECODE FOR SAMPLE ")), bold(yellow(sample_id)), bold(red("...")), "\n"))
    mutation_table$Tot_count <- mutation_table$Ref_count + mutation_table$Alt_count
    mutation_table$VAF <- mutation_table$Alt_count / mutation_table$Tot_count
    if (is.null(max_total_read)) max_total_read <- max(mutation_table$Tot_count)
    #---Choose mutation thresholds, get resulting SFS from data
    SFS_data_frequencies <- seq(1, sfs_bincount) / sfs_bincount
    threshold_results <- choose_mutation_thresholds_ABC(
        mutation_table = mutation_table,
        max_total_read = max_total_read,
        SFS_data_frequencies = SFS_data_frequencies,
        inference_retained_freq = inference_retained_freq,
        validation_mutation_count = validation_mutation_count,
        validation_N_trials = validation_N_trials
    )
    #---Prepare the SFS convolution matrix
    SFS_convolution_inference_A <- build_convolution_matrix_ABC(
        sfs_bincount = sfs_bincount,
        mode = "inference A",
        sample_size = matrix_binomial_sample_size,
        min_variant_read = threshold_results$min_variant_read_inference_A,
        min_total_read = threshold_results$min_total_read_inference_A,
        max_total_read = max_total_read,
        matrix_binomial_ploidy = matrix_binomial_ploidy,
        coverage_distribution = coverage_distribution,
        coverage_variables = coverage_variables,
        sample_coverage = threshold_results$sample_coverage_inference_A
    )
    SFS_convolution_inference_B <- build_convolution_matrix_ABC(
        sfs_bincount = sfs_bincount,
        mode = "inference B",
        sample_size = matrix_binomial_sample_size,
        min_variant_read = threshold_results$min_variant_read_inference_B,
        min_total_read = threshold_results$min_total_read_inference_B,
        max_total_read = max_total_read,
        matrix_binomial_ploidy = matrix_binomial_ploidy,
        coverage_distribution = coverage_distribution,
        coverage_variables = coverage_variables,
        sample_coverage = threshold_results$sample_coverage_inference_B
    )
    SFS_convolution_validation <- build_convolution_matrix_ABC(
        sfs_bincount = sfs_bincount,
        mode = "validation",
        sample_size = matrix_binomial_sample_size,
        min_variant_read = threshold_results$min_variant_read_validation,
        min_total_read = threshold_results$min_total_read_validation,
        max_total_read = max_total_read,
        matrix_binomial_ploidy = matrix_binomial_ploidy,
        coverage_distribution = coverage_distribution,
        coverage_variables = coverage_variables,
        sample_coverage = threshold_results$sample_coverage_validation
    )
    #---DECODE
    DECODE_result <- list()
    DECODE_result$sample_id <- sample_id
    DECODE_result$criterion <- criterion
    DECODE_result$criterion_ratio <- criterion_ratio
    DECODE_result$sfs_bincount <- sfs_bincount
    DECODE_result$mutational_table <- mutation_table
    DECODE_result$readcount_distribution <- threshold_results$readcount_distribution
    DECODE_result$min_variant_read_inference_A <- threshold_results$min_variant_read_inference_A
    DECODE_result$min_total_read_inference_A <- threshold_results$min_total_read_inference_A
    DECODE_result$min_variant_read_inference_B <- threshold_results$min_variant_read_inference_B
    DECODE_result$min_total_read_inference_B <- threshold_results$min_total_read_inference_B
    DECODE_result$min_variant_read_validation <- threshold_results$min_variant_read_validation
    DECODE_result$min_total_read_validation <- threshold_results$min_total_read_validation
    DECODE_result$max_total_read <- max_total_read
    DECODE_result$SFS_frequencies <- SFS_data_frequencies
    DECODE_result$SFS_data_inference_A <- threshold_results$SFS_data_inference_A
    DECODE_result$SFS_data_inference_B <- threshold_results$SFS_data_inference_B
    DECODE_result$SFS_data_validation <- threshold_results$SFS_data_validation
    if (is.na(neutral_tail)) {
        result_with_tail <- DECODE_given_tail_status_ABC(
            SFS_data_inference_A = threshold_results$SFS_data_inference_A,
            SFS_data_inference_B = threshold_results$SFS_data_inference_B,
            SFS_data_validation = threshold_results$SFS_data_validation,
            sfs_bincount = sfs_bincount,
            with_tail = TRUE,
            criterion = criterion,
            criterion_ratio = criterion_ratio,
            min_N_humps = min_N_humps,
            max_N_humps = max_N_humps,
            N_trials = N_trials,
            N_trials_kept = N_trials_kept,
            SFS_convolution_inference_A = SFS_convolution_inference_A,
            SFS_convolution_inference_B = SFS_convolution_inference_B,
            SFS_convolution_validation = SFS_convolution_validation,
            neutral_power_min = neutral_power_min,
            neutral_power_max = neutral_power_max,
            cluster_frequency_min = cluster_frequency_min,
            cluster_frequency_max = cluster_frequency_max,
            pi_cutoff = pi_cutoff,
            zero_cutoff = zero_cutoff,
            compute_parallel = compute_parallel,
            n_cores = n_cores
        )
        DECODE_result$fits_with_tail <- result_with_tail
        result_without_tail <- DECODE_given_tail_status_ABC(
            SFS_data_inference_A = threshold_results$SFS_data_inference_A,
            SFS_data_inference_B = threshold_results$SFS_data_inference_B,
            SFS_data_validation = threshold_results$SFS_data_validation,
            sfs_bincount = sfs_bincount,
            with_tail = FALSE,
            criterion = criterion,
            criterion_ratio = criterion_ratio,
            min_N_humps = min_N_humps,
            max_N_humps = max_N_humps,
            N_trials = N_trials,
            N_trials_kept = N_trials_kept,
            SFS_convolution_inference_A = SFS_convolution_inference_A,
            SFS_convolution_inference_B = SFS_convolution_inference_B,
            SFS_convolution_validation = SFS_convolution_validation,
            neutral_power_min = neutral_power_min,
            neutral_power_max = neutral_power_max,
            cluster_frequency_min = cluster_frequency_min,
            cluster_frequency_max = cluster_frequency_max,
            pi_cutoff = pi_cutoff,
            zero_cutoff = zero_cutoff,
            compute_parallel = compute_parallel,
            n_cores = n_cores
        )
        DECODE_result$fits_without_tail <- result_without_tail
        if (result_with_tail$best_fit$selected_criterion_value < result_without_tail$best_fit$selected_criterion_value) {
            N_humps <- min_N_humps
            while (!is.null(DECODE_result$fits_with_tail$all_fits[[paste0(N_humps, "_clusters")]])) {
                if (DECODE_result$fits_with_tail$all_fits[[paste0(N_humps, "_clusters")]]$note == "best given tail status") {
                    DECODE_result$fits_with_tail$all_fits[[paste0(N_humps, "_clusters")]]$note <- "best"
                }
                N_humps <- N_humps + 1
            }
            final_result <- result_with_tail
        } else {
            N_humps <- min_N_humps
            while (!is.null(DECODE_result$fits_without_tail$all_fits[[paste0(N_humps, "_clusters")]])) {
                if (DECODE_result$fits_without_tail$all_fits[[paste0(N_humps, "_clusters")]]$note == "best given tail status") {
                    DECODE_result$fits_without_tail$all_fits[[paste0(N_humps, "_clusters")]]$note <- "best"
                }
                N_humps <- N_humps + 1
            }
            final_result <- result_without_tail
        }
    } else if (neutral_tail == TRUE) {
        final_result <- DECODE_given_tail_status_ABC(
            SFS_data_inference_A = threshold_results$SFS_data_inference_A,
            SFS_data_inference_B = threshold_results$SFS_data_inference_B,
            SFS_data_validation = threshold_results$SFS_data_validation,
            sfs_bincount = sfs_bincount,
            with_tail = TRUE,
            criterion = criterion,
            criterion_ratio = criterion_ratio,
            min_N_humps = min_N_humps,
            max_N_humps = max_N_humps,
            N_trials = N_trials,
            N_trials_kept = N_trials_kept,
            SFS_convolution_inference_A = SFS_convolution_inference_A,
            SFS_convolution_inference_B = SFS_convolution_inference_B,
            SFS_convolution_validation = SFS_convolution_validation,
            neutral_power_min = neutral_power_min,
            neutral_power_max = neutral_power_max,
            cluster_frequency_min = cluster_frequency_min,
            cluster_frequency_max = cluster_frequency_max,
            pi_cutoff = pi_cutoff,
            zero_cutoff = zero_cutoff,
            compute_parallel = compute_parallel,
            n_cores = n_cores
        )
        DECODE_result$fits_with_tail <- final_result
        N_humps <- min_N_humps
        while (!is.null(DECODE_result$fits_with_tail$all_fits[[paste0(N_humps, "_clusters")]])) {
            if (DECODE_result$fits_with_tail$all_fits[[paste0(N_humps, "_clusters")]]$note == "best given tail status") {
                DECODE_result$fits_with_tail$all_fits[[paste0(N_humps, "_clusters")]]$note <- "best"
            }
            N_humps <- N_humps + 1
        }
    } else if (neutral_tail == FALSE) {
        final_result <- DECODE_given_tail_status_ABC(
            SFS_data_inference_A = threshold_results$SFS_data_inference_A,
            SFS_data_inference_B = threshold_results$SFS_data_inference_B,
            SFS_data_validation = threshold_results$SFS_data_validation,
            sfs_bincount = sfs_bincount,
            with_tail = FALSE,
            criterion = criterion,
            criterion_ratio = criterion_ratio,
            min_N_humps = min_N_humps,
            max_N_humps = max_N_humps,
            N_trials = N_trials,
            N_trials_kept = N_trials_kept,
            SFS_convolution_inference_A = SFS_convolution_inference_A,
            SFS_convolution_inference_B = SFS_convolution_inference_B,
            SFS_convolution_validation = SFS_convolution_validation,
            neutral_power_min = neutral_power_min,
            neutral_power_max = neutral_power_max,
            cluster_frequency_min = cluster_frequency_min,
            cluster_frequency_max = cluster_frequency_max,
            pi_cutoff = pi_cutoff,
            zero_cutoff = zero_cutoff,
            compute_parallel = compute_parallel,
            n_cores = n_cores
        )
        DECODE_result$fits_without_tail <- final_result
        N_humps <- min_N_humps
        while (!is.null(DECODE_result$fits_without_tail$all_fits[[paste0(N_humps, "_clusters")]])) {
            if (DECODE_result$fits_without_tail$all_fits[[paste0(N_humps, "_clusters")]]$note == "best given tail status") {
                DECODE_result$fits_without_tail$all_fits[[paste0(N_humps, "_clusters")]]$note <- "best"
            }
            N_humps <- N_humps + 1
        }
    }
    DECODE_result$final_fit <- final_result
    #---Report the best fit
    tail_status_final_result <- final_result$best_fit$tail_status
    parameters_final_result <- final_result$best_fit$parameters
    criterion_value_final_result <- final_result$best_fit$selected_criterion_value
    N_humps_final_result <- parameters_final_result[["Cluster_count"]]
    if (tail_status_final_result) {
        report <- bold(underline(red(paste0("Best fit = neutral tail + ", N_humps_final_result, " clusters:\n"))))
        report <- paste0(report, red("Score            : "), yellow(paste0(criterion, " = ", format(round(criterion_value_final_result, 3), nsmall = 3))), "\n")
        report <- paste0(report, red("Neutral component: "), yellow(paste0("\u03B1   = ", format(round(parameters_final_result[["Tail_power"]], 3), nsmall = 3))), red(", "))
        report <- paste0(report, yellow(paste0("\u03C0 = ", format(round(parameters_final_result[["Tail_proportion_inference_A"]], 3), nsmall = 3), " [A], ", format(round(parameters_final_result[["Tail_proportion_inference_B"]], 3), nsmall = 3), " [B]")), "\n")
        ii <- 0
    } else {
        report <- bold(underline(red(paste0("Best fit = no neutral tail + ", N_humps_final_result, " clusters:\n"))))
        report <- paste0(report, red("Score            : "), yellow(paste0(criterion, " = ", format(round(criterion_value_final_result, 3), nsmall = 3))), "\n")
        ii <- -1
    }
    if (N_humps_final_result > 0) {
        for (i in 1:N_humps_final_result) {
            report <- paste0(report, red(paste0("Cluster ", i, "        : ")), yellow(paste0("VAF = ", format(round(parameters_final_result[[paste0("Cluster_VAF_", i)]], 3), nsmall = 3))), red(", "))
            report <- paste0(report, yellow(paste0("\u03C0 = ", format(round(parameters_final_result[[paste0("Cluster_proportion_inference_A_", i)]], 3), nsmall = 3), " [A], ", format(round(parameters_final_result[[paste0("Cluster_proportion_inference_B_", i)]], 3), nsmall = 3), " [B]")), "\n")
        }
    }
    cat(report)
    # #---Translation to parameters of cancer evolution in the sample
    # tmp <- parameter_conversion_ABC(
    #     result = final_result,
    #     mutation_count_for_fitting = sum(threshold_results$SFS_data_inference_A),
    #     sample_size = sample_size,
    #     matrix_binomial_sample_size = matrix_binomial_sample_size,
    #     matrix_binomial_ploidy = matrix_binomial_ploidy
    # )
    # DECODE_result$final_fit$parameters_df <- tmp$parameters_df
    #---Return the SFS deconvolution results
    return(DECODE_result)
}

DECODE_given_tail_status_ABC <- function(SFS_data_inference_A,
                                         SFS_data_inference_B,
                                         SFS_data_validation,
                                         criterion,
                                         criterion_ratio,
                                         min_N_humps,
                                         max_N_humps,
                                         with_tail,
                                         N_trials,
                                         N_trials_kept,
                                         sfs_bincount,
                                         SFS_convolution_inference_A,
                                         SFS_convolution_inference_B,
                                         SFS_convolution_validation,
                                         neutral_power_min,
                                         neutral_power_max,
                                         cluster_frequency_min,
                                         cluster_frequency_max,
                                         pi_cutoff,
                                         zero_cutoff,
                                         compute_parallel,
                                         n_cores) {
    N_humps <- min_N_humps
    criterion_best_final <- Inf
    all_fits <- list()
    while (TRUE) {
        #---Find best parameter set, given the number of humps
        cat(bold(blue(paste0("Inference for ", N_humps, " clusters ", ifelse(with_tail, "with", "without"), " neutral tail component...\n"))))
        fit_results <- DECODE_given_tail_status_and_Ncluster_ABC(
            SFS_data_inference_A = SFS_data_inference_A,
            SFS_data_inference_B = SFS_data_inference_B,
            SFS_data_validation = SFS_data_validation,
            N_humps = N_humps,
            with_tail = with_tail,
            N_trials = N_trials,
            N_trials_kept = N_trials_kept,
            sfs_bincount = sfs_bincount,
            SFS_convolution_inference_A = SFS_convolution_inference_A,
            SFS_convolution_inference_B = SFS_convolution_inference_B,
            SFS_convolution_validation = SFS_convolution_validation,
            neutral_power_min = neutral_power_min,
            neutral_power_max = neutral_power_max,
            cluster_frequency_min = cluster_frequency_min,
            cluster_frequency_max = cluster_frequency_max,
            zero_cutoff = zero_cutoff,
            compute_parallel = compute_parallel,
            n_cores = n_cores
        )
        fit_results$note <- "none"
        all_fits[[paste0(N_humps, "_clusters")]] <- fit_results
        parameters_best_current <- fit_results$summary$parameters
        component_distributions_inference_A_best_current <- fit_results$summary$component_distributions_inference_A
        component_distributions_inference_B_best_current <- fit_results$summary$component_distributions_inference_B
        component_distributions_validation_best_current <- fit_results$summary$component_distributions_validation
        criterion_all_best_current <- fit_results$summary$criteria
        criterion_best_current <- fit_results$summary$criteria[[criterion]]
        #   Report the best fit for the current hump count
        cluster_pis_inference_A <- Inf
        cluster_pis_inference_B <- Inf
        report <- paste0(blue("Score            : "), cyan(paste0(criterion, " = ", format(round(criterion_best_current, 3), nsmall = 3), " (+/-", format(round(fit_results$summary$criteria[[paste0(criterion, "_sd")]], 3), nsmall = 3), ")")), "\n")
        if (with_tail) {
            report <- paste0(report, blue("Neutral component: "), cyan(paste0("\u03B1   = ", format(round(parameters_best_current[["Tail_power"]], 3), nsmall = 3))), blue(", "))
            report <- paste0(report, cyan(paste0("\u03C0 = ", format(round(parameters_best_current[["Tail_proportion_inference_A"]], 3), nsmall = 3), " [A], ", format(round(parameters_best_current[["Tail_proportion_inference_B"]], 3), nsmall = 3), " [B]")), "\n")
            ii <- 0
        } else {
            ii <- -1
        }
        if (N_humps > 0) {
            for (i in 1:N_humps) {
                report <- paste0(report, blue(paste0("Cluster ", i, "        : ")), cyan(paste0("VAF = ", format(round(parameters_best_current[[paste0("Cluster_VAF_", i)]], 3), nsmall = 3))), blue(", "))
                report <- paste0(report, cyan(paste0("\u03C0 = ", format(round(parameters_best_current[[paste0("Cluster_proportion_inference_A_", i)]], 3), nsmall = 3), " [A], ", format(round(parameters_best_current[[paste0("Cluster_proportion_inference_B_", i)]], 3), nsmall = 3), " [B]")), "\n")
                cluster_pis_inference_A <- c(cluster_pis_inference_A, parameters_best_current[[paste0("Cluster_proportion_inference_A_", i)]])
                cluster_pis_inference_B <- c(cluster_pis_inference_B, parameters_best_current[[paste0("Cluster_proportion_inference_B_", i)]])
            }
        }
        cat(report)
        #   Check if the increased hump count leads to lower criterion score without tiny selective components...
        if ((N_humps == min_N_humps) | ((criterion_best_current < criterion_ratio * criterion_best_final) & (min(pmax(cluster_pis_inference_A, cluster_pis_inference_B)) >= pi_cutoff))) {
            #   ... if yes, then update the best fit and continue with 1 more hump
            N_humps_best_final <- N_humps
            fit_results_best_final <- fit_results
            criterion_best_final <- criterion_best_current
            criterion_all_final <- criterion_all_best_current
            parameters_best_final <- parameters_best_current
            component_distributions_inference_A_best_final <- component_distributions_inference_A_best_current
            component_distributions_inference_B_best_final <- component_distributions_inference_B_best_current
            component_distributions_validation_best_final <- component_distributions_validation_best_current
            N_humps <- N_humps + 1
            #   ... except if exceeding maximum number of clusters
            if (N_humps > max_N_humps) {
                break
            }
        } else {
            #   ... if no, then stop
            break
        }
    }
    #---Report the best fit
    all_fits[[paste0(N_humps_best_final, "_clusters")]]$note <- "best given tail status"
    result <- list()
    result$all_fits <- all_fits
    result$best_fit <- list()
    result$best_fit$parameters <- parameters_best_final
    result$best_fit$component_distributions_inference_A <- component_distributions_inference_A_best_final
    result$best_fit$component_distributions_inference_B <- component_distributions_inference_B_best_final
    result$best_fit$component_distributions_validation <- component_distributions_validation_best_final
    result$best_fit$selected_criterion <- criterion
    result$best_fit$all_criteria <- criterion_all_final
    result$best_fit$selected_criterion_value <- criterion_best_final
    result$best_fit$tail_status <- with_tail
    return(result)
}

DECODE_given_tail_status_and_Ncluster_ABC <- function(SFS_data_inference_A,
                                                      SFS_data_inference_B,
                                                      SFS_data_validation,
                                                      N_humps,
                                                      with_tail,
                                                      N_trials,
                                                      N_trials_kept,
                                                      sfs_bincount,
                                                      SFS_convolution_inference_A,
                                                      SFS_convolution_inference_B,
                                                      SFS_convolution_validation,
                                                      neutral_power_min,
                                                      neutral_power_max,
                                                      cluster_frequency_min,
                                                      cluster_frequency_max,
                                                      zero_cutoff,
                                                      compute_parallel,
                                                      n_cores,
                                                      progress_bar = TRUE) {
    N_end <- SFS_convolution_inference_A$N_end
    SFS_convolution_matrix_inference_A <- SFS_convolution_inference_A$convolution_matrix
    SFS_convolution_matrix_inference_B <- SFS_convolution_inference_B$convolution_matrix
    SFS_convolution_matrix_validation <- SFS_convolution_validation$convolution_matrix
    #---Function to perform one trial to find A & K's
    pi_conversion <- function(from_pis,
                              from_component_distributions,
                              to_component_distributions,
                              with_tail) {
        to_pis <- from_pis
        ii <- ifelse(with_tail, 0, 1)
        for (i in 1:length(to_pis)) {
            to_pis[i] <- to_pis[i] *
                sum(to_component_distributions$SFS_expected[i + ii, ]) /
                sum(from_component_distributions$SFS_expected[i + ii, ])
        }
        to_pis <- to_pis / sum(to_pis)
        return(to_pis)
    }
    func_ABC_trial <- function(with_tail,
                               N_humps,
                               neutral_power_min,
                               neutral_power_max,
                               cluster_frequency_min,
                               cluster_frequency_max,
                               SFS_data_inference_A,
                               SFS_data_inference_B,
                               #    SFS_data_validation,
                               sfs_bincount,
                               N_end,
                               SFS_convolution_matrix_inference_A,
                               SFS_convolution_matrix_inference_B,
                               #    SFS_convolution_matrix_validation,
                               zero_cutoff) {
        #   Sample neutral component power
        neutral_power <- ifelse(with_tail, runif(1, neutral_power_min, neutral_power_max), NA)
        #   Sample cluster frequencies
        cluster_frequencies <- sort(runif(N_humps, cluster_frequency_min, cluster_frequency_max), decreasing = TRUE)
        #   Build the SFS component libraries
        component_distributions_inference_A <- build_SFS_library_ABC(
            neutral_power = neutral_power,
            cluster_frequencies = cluster_frequencies,
            sfs_bincount = sfs_bincount,
            SFS_convolution_matrix = SFS_convolution_matrix_inference_A,
            N_end = N_end
        )
        component_distributions_inference_B <- build_SFS_library_ABC(
            neutral_power = neutral_power,
            cluster_frequencies = cluster_frequencies,
            sfs_bincount = sfs_bincount,
            SFS_convolution_matrix = SFS_convolution_matrix_inference_B,
            N_end = N_end
        )
        #   Sample component proportions for inference A
        library(DirichletReg)
        cluster_proportions_inference_A <- rdirichlet(n = 1, alpha = rep(1, N_humps + with_tail))
        # cluster_proportions_inference_A <- runif(N_humps + with_tail, 0, 1)
        # cluster_proportions_inference_A <- cluster_proportions_inference_A / sum(cluster_proportions_inference_A)
        #   Sample component proportions for inference B
        cluster_proportions_inference_B <- pi_conversion(
            from_pis = cluster_proportions_inference_A,
            from_component_distributions = component_distributions_inference_A,
            to_component_distributions = component_distributions_inference_B,
            with_tail = with_tail
        )
        #   Compute log-likelihoods
        logLikelihood_inference_A <- compute_loglikelihood_ABC(
            A = ifelse(with_tail, cluster_proportions_inference_A[1], NA),
            vec_K = ifelse(with_tail, cluster_proportions_inference_A[-1], cluster_proportions_inference_A),
            component_distributions = component_distributions_inference_A,
            vec_SFS_real = SFS_data_inference_A,
            zero_cutoff = zero_cutoff
        )
        logLikelihood_inference_B <- compute_loglikelihood_ABC(
            A = ifelse(with_tail, cluster_proportions_inference_B[1], NA),
            vec_K = ifelse(with_tail, cluster_proportions_inference_B[-1], cluster_proportions_inference_B),
            component_distributions = component_distributions_inference_B,
            vec_SFS_real = SFS_data_inference_B,
            zero_cutoff = zero_cutoff
        )
        logLikelihood <- logLikelihood_inference_A + logLikelihood_inference_B
        # 	Prepare the results to be returned
        result <- data.frame(
            logLikelihood = logLikelihood,
            logLikelihood_inference_A = logLikelihood_inference_A,
            logLikelihood_inference_B = logLikelihood_inference_B,
            Nmut_inference_A = sum(SFS_data_inference_A),
            Nmut_inference_B = sum(SFS_data_inference_B),
            Tail = with_tail,
            Tail_proportion_inference_A = ifelse(with_tail, cluster_proportions_inference_A[1], NA),
            Tail_Nmut_inference_A = ifelse(with_tail, sum(SFS_data_inference_A) * cluster_proportions_inference_A[1], NA),
            Tail_proportion_inference_B = ifelse(with_tail, cluster_proportions_inference_B[1], NA),
            Tail_Nmut_inference_B = ifelse(with_tail, sum(SFS_data_inference_B) * cluster_proportions_inference_B[1], NA),
            Tail_power = neutral_power,
            Cluster_count = N_humps
        )
        if (N_humps > 0) {
            for (i in 1:N_humps) {
                result[[paste0("Cluster_proportion_inference_A_", i)]] <- ifelse(with_tail, cluster_proportions_inference_A[i + 1], cluster_proportions_inference_A[i])
                result[[paste0("Cluster_Nmut_inference_A_", i)]] <- sum(SFS_data_inference_A) * result[[paste0("Cluster_proportion_inference_A_", i)]]
                result[[paste0("Cluster_proportion_inference_B_", i)]] <- ifelse(with_tail, cluster_proportions_inference_B[i + 1], cluster_proportions_inference_B[i])
                result[[paste0("Cluster_Nmut_inference_B_", i)]] <- sum(SFS_data_inference_B) * result[[paste0("Cluster_proportion_inference_B_", i)]]
                result[[paste0("Cluster_VAF_", i)]] <- cluster_frequencies[i]
            }
        }
        #   Return results from this ABC trial
        output <- list()
        output$result <- result
        output$component_distributions_inference_A <- component_distributions_inference_A
        output$component_distributions_inference_B <- component_distributions_inference_B
        return(output)
    }
    #---Find log-likelihoods for ABC trials
    if (compute_parallel == FALSE) {
        prior_results <- c()
        prior_component_distributions_inference_A <- list()
        prior_component_distributions_inference_B <- list()
        if (progress_bar) {
            pb <- txtProgressBar(
                min = 0,
                max = N_trials,
                style = 3,
                width = 50,
                char = "+"
            )
        }
        for (i in 1:N_trials) {
            if (progress_bar) setTxtProgressBar(pb, i)
            trial_result <- func_ABC_trial(
                with_tail = with_tail,
                N_humps = N_humps,
                neutral_power_min = neutral_power_min,
                neutral_power_max = neutral_power_max,
                cluster_frequency_min = cluster_frequency_min,
                cluster_frequency_max = cluster_frequency_max,
                SFS_data_inference_A = SFS_data_inference_A,
                SFS_data_inference_B = SFS_data_inference_B,
                sfs_bincount = sfs_bincount,
                N_end = N_end,
                SFS_convolution_matrix_inference_A = SFS_convolution_matrix_inference_A,
                SFS_convolution_matrix_inference_B = SFS_convolution_matrix_inference_B,
                zero_cutoff = zero_cutoff
            )
            prior_results <- rbind(prior_results, trial_result$result)
            prior_component_distributions_inference_A[[i]] <- trial_result$component_distributions_inference_A
            prior_component_distributions_inference_B[[i]] <- trial_result$component_distributions_inference_B
        }
        if (progress_bar) cat("\n")
    } else {
        suppressPackageStartupMessages(library(parallel))
        suppressPackageStartupMessages(library(pbapply))
        #   Start parallel cluster
        n_cores <- ifelse(is.null(n_cores), detectCores() - 1, n_cores)
        cl <- makePSOCKcluster(n_cores)
        #   Prepare input parameters
        clusterExport(cl, varlist = c(
            "with_tail", "N_humps", "sfs_bincount", "N_end", "zero_cutoff",
            "SFS_data_inference_A", "SFS_data_inference_B",
            "SFS_convolution_matrix_inference_A", "SFS_convolution_matrix_inference_B",
            "neutral_power_min", "neutral_power_max",
            "cluster_frequency_min", "cluster_frequency_max",
            "build_SFS_library_ABC", "build_SFS_library_Griffiths_Tavare_ABC",
            "compute_loglikelihood_ABC", "compute_SFS_ABC"
        ), envir = environment())
        #   Find log-likelihoods for ABC trials in parallel mode
        if (progress_bar) {
            pboptions(type = "timer", nout = 2)
            output <- pblapply(cl = cl, X = 1:N_trials, FUN = function(i) {
                func_ABC_trial(
                    with_tail = with_tail,
                    N_humps = N_humps,
                    neutral_power_min = neutral_power_min,
                    neutral_power_max = neutral_power_max,
                    cluster_frequency_min = cluster_frequency_min,
                    cluster_frequency_max = cluster_frequency_max,
                    SFS_data_inference_A = SFS_data_inference_A,
                    SFS_data_inference_B = SFS_data_inference_B,
                    sfs_bincount = sfs_bincount,
                    N_end = N_end,
                    SFS_convolution_matrix_inference_A = SFS_convolution_matrix_inference_A,
                    SFS_convolution_matrix_inference_B = SFS_convolution_matrix_inference_B,
                    zero_cutoff = zero_cutoff
                )
            })
        } else {
            output <- parLapply(cl = cl, X = 1:N_trials, fun = function(i) {
                func_ABC_trial(
                    with_tail = with_tail,
                    N_humps = N_humps,
                    neutral_power_min = neutral_power_min,
                    neutral_power_max = neutral_power_max,
                    cluster_frequency_min = cluster_frequency_min,
                    cluster_frequency_max = cluster_frequency_max,
                    SFS_data_inference_A = SFS_data_inference_A,
                    SFS_data_inference_B = SFS_data_inference_B,
                    sfs_bincount = sfs_bincount,
                    N_end = N_end,
                    SFS_convolution_matrix_inference_A = SFS_convolution_matrix_inference_A,
                    SFS_convolution_matrix_inference_B = SFS_convolution_matrix_inference_B,
                    zero_cutoff = zero_cutoff
                )
            })
        }
        stopCluster(cl)
        #   Extract the results
        prior_results <- do.call(rbind, lapply(output, function(x) x$result))
        prior_component_distributions_inference_A <- lapply(output, function(x) x$component_distributions_inference_A)
        prior_component_distributions_inference_B <- lapply(output, function(x) x$component_distributions_inference_B)
    }
    #---Find posterior distributions for DECODE parameters
    posterior_indices <- order(prior_results$logLikelihood, decreasing = TRUE)[1:N_trials_kept]
    posterior_results <- prior_results[posterior_indices, ]
    posterior_component_distributions_inference_A <- prior_component_distributions_inference_A[posterior_indices]
    posterior_component_distributions_inference_B <- prior_component_distributions_inference_B[posterior_indices]
    #---Compute posterior results and SFS component libraries for validation
    posterior_component_distributions_validation <- list()
    posterior_results[["Tail_proportion_validation"]] <- NA
    posterior_results[["Tail_Nmut_validation"]] <- NA
    if (N_humps > 0) {
        for (i in 1:N_humps) {
            posterior_results[[paste0("Cluster_proportion_validation_", i)]] <- NA
            posterior_results[[paste0("Cluster_Nmut_validation_", i)]] <- NA
        }
    }
    for (posterior_trial in 1:N_trials_kept) {
        neutral_power <- posterior_results[["Tail_power"]][posterior_trial]
        if (with_tail) {
            cluster_proportions_inference_A <- posterior_results[["Tail_proportion_inference_A"]][posterior_trial]
        } else {
            cluster_proportions_inference_A <- c()
        }
        cluster_frequencies <- c()
        for (i in 1:N_humps) {
            cluster_proportions_inference_A <- c(cluster_proportions_inference_A, posterior_results[[paste0("Cluster_proportion_inference_A_", i)]][posterior_trial])
            cluster_frequencies <- c(cluster_frequencies, posterior_results[[paste0("Cluster_VAF_", i)]][posterior_trial])
        }
        #   Build the SFS component libraries for validation
        component_distributions_validation <- build_SFS_library_ABC(
            neutral_power = neutral_power,
            cluster_frequencies = cluster_frequencies,
            sfs_bincount = sfs_bincount,
            SFS_convolution_matrix = SFS_convolution_matrix_validation,
            N_end = N_end
        )
        posterior_component_distributions_validation[[posterior_trial]] <- component_distributions_validation
        #   Sample component proportions for validation
        cluster_proportions_validation <- pi_conversion(
            from_pis = cluster_proportions_inference_A,
            from_component_distributions = posterior_component_distributions_inference_A[[posterior_trial]],
            to_component_distributions = component_distributions_validation,
            with_tail = with_tail
        )
        #   Update posterior results
        posterior_results[["Tail_proportion_validation"]][posterior_trial] <- ifelse(with_tail, cluster_proportions_validation[1], NA)
        posterior_results[["Tail_Nmut_validation"]][posterior_trial] <- ifelse(with_tail, sum(SFS_data_validation[1, ]) * cluster_proportions_validation[1], NA)
        if (N_humps > 0) {
            for (i in 1:N_humps) {
                posterior_results[[paste0("Cluster_proportion_validation_", i)]][posterior_trial] <- ifelse(with_tail, cluster_proportions_validation[i + 1], cluster_proportions_validation[i])
                posterior_results[[paste0("Cluster_Nmut_validation_", i)]][posterior_trial] <- sum(SFS_data_validation[1, ]) * posterior_results[[paste0("Cluster_proportion_validation_", i)]][posterior_trial]
            }
        }
    }
    #---Compute stopping criteria from posterior distribution
    num_parameters <- ifelse(with_tail, 2 * N_humps + 1, 2 * N_humps - 1)
    if (neutral_power_max > neutral_power_min) num_parameters <- num_parameters + 1
    criteria_results <- cluster_count_criteria_ABC(
        num_parameters = num_parameters,
        SFS_data_validation = SFS_data_validation,
        parameters = posterior_results,
        component_distributions_validation = posterior_component_distributions_validation,
        with_tail = with_tail,
        zero_cutoff = zero_cutoff
    )
    posterior_results <- criteria_results$parameters
    summary_criteria <- criteria_results$criteria
    #---Get the single best fit
    best_result <- posterior_results[1, ]
    best_component_distributions_inference_A <- posterior_component_distributions_inference_A[[1]]
    best_component_distributions_inference_B <- posterior_component_distributions_inference_B[[1]]
    best_component_distributions_validation <- posterior_component_distributions_validation[[1]]
    #---Find the best fit
    fit_results <- list()
    fit_results$prior <- list()
    fit_results$prior$parameters <- prior_results
    fit_results$posterior <- list()
    fit_results$posterior$parameters <- posterior_results
    # fit_results$posterior$component_distributions_inference_A <- posterior_component_distributions_inference_A
    # fit_results$posterior$component_distributions_inference_B <- posterior_component_distributions_inference_B
    # fit_results$posterior$component_distributions_validation <- posterior_component_distributions_validation
    fit_results$summary <- list()
    fit_results$summary$parameters <- best_result
    fit_results$summary$criteria <- summary_criteria
    fit_results$summary$component_distributions_inference_A <- best_component_distributions_inference_A
    fit_results$summary$component_distributions_inference_B <- best_component_distributions_inference_B
    fit_results$summary$component_distributions_validation <- best_component_distributions_validation
    return(fit_results)
}

choose_mutation_thresholds_ABC <- function(mutation_table,
                                           max_total_read,
                                           SFS_data_frequencies,
                                           inference_retained_freq,
                                           validation_mutation_count,
                                           validation_N_trials) {
    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(data.table))
    cat(bold(blue("Choose mutation thresholds for inference and validation...\n")))
    #---Find joint distribution of variant and total readcounts
    mutation_table_tmp <- mutation_table %>%
        group_by(Alt_count, Tot_count) %>%
        summarise(count = n(), .groups = "drop") %>%
        ungroup()
    vec_min_variant_read <- min(mutation_table_tmp$Alt_count):max(mutation_table_tmp$Alt_count)
    vec_min_total_read <- min(mutation_table_tmp$Tot_count):max_total_read
    readcount_distribution <- expand.grid(min_variant_read = vec_min_variant_read, min_total_read = vec_min_total_read, mutation_count = 0)
    pb <- txtProgressBar(
        min = 1,
        max = nrow(mutation_table_tmp),
        style = 3,
        width = 50,
        char = "+"
    )
    for (row in 1:nrow(mutation_table_tmp)) {
        setTxtProgressBar(pb, row)
        min_variant_read <- mutation_table_tmp$Alt_count[row]
        min_total_read <- mutation_table_tmp$Tot_count[row]
        rows <- which(readcount_distribution$min_variant_read <= min_variant_read & readcount_distribution$min_total_read <= min_total_read)
        readcount_distribution$mutation_count[rows] <- readcount_distribution$mutation_count[rows] + mutation_table_tmp$count[row]
    }
    cat("\n")
    readcount_distribution$freq <- 100 * readcount_distribution$mutation_count / sum(mutation_table_tmp$count)
    #---Find the mutation thresholds for inference and validation
    filtered_df <- readcount_distribution[readcount_distribution$min_total_read == min(readcount_distribution$min_total_read), ]
    min_variant_read_tmp <- min(filtered_df$min_variant_read[filtered_df$freq < 100])
    #   Choose mutation thresholds for inference A set
    min_variant_read_inference_A <- min_variant_read_tmp + 2
    filtered_df <- readcount_distribution[readcount_distribution$min_variant_read == min_variant_read_inference_A & readcount_distribution$freq >= inference_retained_freq, ]
    min_total_read_inference_A <- ifelse(nrow(filtered_df) > 0, max(filtered_df$min_total_read), min(readcount_distribution$min_total_read))
    Nmut_inference_A <- readcount_distribution$mutation_count[readcount_distribution$min_variant_read == min_variant_read_inference_A & readcount_distribution$min_total_read == min_total_read_inference_A]
    freq_inference_A <- readcount_distribution$freq[readcount_distribution$min_variant_read == min_variant_read_inference_A & readcount_distribution$min_total_read == min_total_read_inference_A]
    #   Choose mutation thresholds for inference B set
    min_variant_read_inference_B <- min_variant_read_tmp + 6
    filtered_df <- readcount_distribution[readcount_distribution$min_variant_read == min_variant_read_inference_B & readcount_distribution$freq >= inference_retained_freq, ]
    min_total_read_inference_B <- ifelse(nrow(filtered_df) > 0, max(filtered_df$min_total_read), min(readcount_distribution$min_total_read))
    Nmut_inference_B <- readcount_distribution$mutation_count[readcount_distribution$min_variant_read == min_variant_read_inference_B & readcount_distribution$min_total_read == min_total_read_inference_B]
    freq_inference_B <- readcount_distribution$freq[readcount_distribution$min_variant_read == min_variant_read_inference_B & readcount_distribution$min_total_read == min_total_read_inference_B]
    #   Choose mutation thresholds for validation set
    min_variant_read_validation <- min_variant_read_tmp + 4
    min_total_read_validation <- min(readcount_distribution$min_total_read)
    Nmut_validation <- readcount_distribution$mutation_count[readcount_distribution$min_variant_read == min_variant_read_validation & readcount_distribution$min_total_read == min_total_read_validation]
    freq_validation <- readcount_distribution$freq[readcount_distribution$min_variant_read == min_variant_read_validation & readcount_distribution$min_total_read == min_total_read_validation]
    #   Report the chosen mutation thresholds
    report <- paste0(blue("Complete data    : "), cyan(paste0(min(mutation_table$Alt_count), " \u2264 variant reads, ", min(mutation_table$Tot_count), " \u2264 total reads \u2264 ", max(mutation_table$Tot_count), "; ", nrow(mutation_table), " mutations\n")))
    report <- paste0(report, blue("Inference A      : "), cyan(paste0(min_variant_read_inference_A, " \u2264 variant reads, ", min_total_read_inference_A, " \u2264 total reads \u2264 ", max_total_read, "; ", Nmut_inference_A, " mutations (", format(round(freq_inference_A, 3), nsmall = 3), "%)")), "\n")
    report <- paste0(report, blue("Inference B      : "), cyan(paste0(min_variant_read_inference_B, " \u2264 variant reads, ", min_total_read_inference_B, " \u2264 total reads \u2264 ", max_total_read, "; ", Nmut_inference_B, " mutations (", format(round(freq_inference_B, 3), nsmall = 3), "%)")), "\n")
    report <- paste0(report, blue("Validation       : "), cyan(paste0(min_variant_read_validation, " \u2264 variant reads, ", min_total_read_validation, " \u2264 total reads \u2264 ", max_total_read, "; ", Nmut_validation, " mutations (", format(round(freq_validation, 3), nsmall = 3), "%)")), "\n")
    cat(report)
    #---Prepare the real SFS data
    mutation_table <- as.data.table(mutation_table)
    mutation_table_inference_A <- mutation_table[Alt_count >= min_variant_read_inference_A & Tot_count >= min_total_read_inference_A & Tot_count <= max_total_read]
    mutation_table_inference_B <- mutation_table[Alt_count >= min_variant_read_inference_B & Tot_count >= min_total_read_inference_B & Tot_count <= max_total_read]
    mutation_table_validation <- mutation_table[Alt_count >= min_variant_read_validation & Tot_count >= min_total_read_validation & Tot_count <= max_total_read]
    #---Get coverage distribution
    func_coverage <- function(mutation_table) {
        sample_coverage <- mutation_table[, .(pdf = .N), by = .(total_readcount = Tot_count)]
        sample_coverage[, pdf := pdf / sum(pdf)]
        return(sample_coverage)
    }
    sample_coverage_inference_A <- func_coverage(mutation_table_inference_A)
    sample_coverage_inference_B <- func_coverage(mutation_table_inference_B)
    sample_coverage_validation <- func_coverage(mutation_table_validation)
    #---Prepare the real SFS data for the inference sets
    func_SFS <- function(mutation_table) {
        VAF <- mutation_table$VAF
        lower_bound <- c(0, SFS_data_frequencies[-length(SFS_data_frequencies)])
        upper_bound <- SFS_data_frequencies
        SFS_data <- numeric(length(SFS_data_frequencies))
        for (i in seq_along(SFS_data_frequencies)) SFS_data[i] <- sum(VAF > lower_bound[i] & VAF <= upper_bound[i])
        return(SFS_data)
    }
    SFS_data_inference_A <- func_SFS(mutation_table_inference_A)
    SFS_data_inference_B <- func_SFS(mutation_table_inference_B)
    #---Prepare the real SFS data for the validation set
    SFS_data_validation <- lapply(1:validation_N_trials, function(i) {
        mutation_table_validation_tmp <- mutation_table_validation[sample(nrow(mutation_table_validation), validation_mutation_count, replace = TRUE), ]
        func_SFS(mutation_table_validation_tmp)
    })
    SFS_data_validation <- do.call(rbind, SFS_data_validation)
    #---Report the chosen thresholds
    results <- list()
    results$readcount_distribution <- readcount_distribution
    results$SFS_data_inference_A <- SFS_data_inference_A
    results$SFS_data_inference_B <- SFS_data_inference_B
    results$SFS_data_validation <- SFS_data_validation
    results$sample_coverage_inference_A <- sample_coverage_inference_A
    results$sample_coverage_inference_B <- sample_coverage_inference_B
    results$sample_coverage_validation <- sample_coverage_validation
    results$min_variant_read_inference_A <- min_variant_read_inference_A
    results$min_total_read_inference_A <- min_total_read_inference_A
    results$min_variant_read_inference_B <- min_variant_read_inference_B
    results$min_total_read_inference_B <- min_total_read_inference_B
    results$min_variant_read_validation <- min_variant_read_validation
    results$min_total_read_validation <- min_total_read_validation
    return(results)
}

cluster_count_criteria_ABC <- function(num_parameters,
                                       SFS_data_validation,
                                       parameters,
                                       component_distributions_validation,
                                       with_tail,
                                       zero_cutoff) {
    #---Functions for stopping criteria
    compute_AIC <- function(logLikelihood, num_parameters) {
        AIC <- 2 * num_parameters - 2 * logLikelihood
        return(AIC)
    }
    compute_BIC <- function(logLikelihood, num_parameters, num_samples) {
        BIC <- num_parameters * log(num_samples) - 2 * logLikelihood
        return(BIC)
    }
    # compute_ICL <- function() {
    #     #   Compute the latent variable distributions
    #     latent_variable_distributions <- component_distributions_validation$SFS_expected_normalized
    #     for (row in 1:nrow(latent_variable_distributions)) {
    #         latent_variable_distributions[row, ] <- parameters_validation[2 * row - 1] * latent_variable_distributions[row, ]
    #     }
    #     for (col in 1:ncol(latent_variable_distributions)) {
    #         latent_variable_distributions[, col] <- latent_variable_distributions[, col] / sum(latent_variable_distributions[, col])
    #     }
    #     latent_variable_distributions[which(latent_variable_distributions <= zero_cutoff | is.na(latent_variable_distributions))] <- zero_cutoff
    #     #   Compute the entropy
    #     entropy <- sum(SFS_data_validation * colSums(latent_variable_distributions * log(latent_variable_distributions)))
    #     #   Compute the Bayesian Information Criterion
    #     BIC <- log_L - 0.5 * num_parameters * log(num_samples_validation)
    #     #   Compute the Integrated Completed Log-Likelihood
    #     ICL <- -BIC - entropy
    #     return(ICL)
    # }
    # compute_ICL_MAP <- function() {
    #     #   Compute the latent variable distributions
    #     latent_variable_distributions <- component_distributions_validation$SFS_expected_normalized
    #     for (row in 1:nrow(latent_variable_distributions)) {
    #         latent_variable_distributions[row, ] <- parameters_validation[2 * row - 1] * latent_variable_distributions[row, ]
    #     }
    #     for (col in 1:ncol(latent_variable_distributions)) {
    #         latent_variable_distributions[, col] <- latent_variable_distributions[, col] / sum(latent_variable_distributions[, col])
    #     }
    #     latent_variable_distributions[which(latent_variable_distributions <= zero_cutoff | is.na(latent_variable_distributions))] <- zero_cutoff
    #     #   Compute the MAP allocations for mutations to clusters
    #     indicator_latent_variable_distributions <- matrix(0, nrow = nrow(latent_variable_distributions), ncol = ncol(latent_variable_distributions))
    #     for (col in 1:ncol(latent_variable_distributions)) {
    #         max_p <- which(latent_variable_distributions[, col] == max(latent_variable_distributions[, col]))[1]
    #         indicator_latent_variable_distributions[max_p, col] <- 1
    #     }
    #     #   Compute the entropy
    #     entropy_MAP <- sum(SFS_data_validation * colSums(indicator_latent_variable_distributions * log(latent_variable_distributions)))
    #     #   Compute the Bayesian Information Criterion
    #     BIC <- log_L - 0.5 * num_parameters * log(num_samples_validation)
    #     #   Compute the Integrated Completed Log-Likelihood
    #     ICL_MAP <- -BIC - entropy_MAP
    #     return(ICL_MAP)
    # }
    #---Compute stopping criteria for each parameter & validation set
    all_criteria <- data.frame(matrix(NA, nrow = 0, ncol = 3), stringsAsFactors = FALSE)
    colnames(all_criteria) <- c("parameter_id", "validation_id", "logLikelihood")
    for (i_parameter in 1:nrow(parameters)) {
        A <- parameters["Tail_proportion_validation"][i_parameter, ]
        vec_K <- c()
        for (i in 1:parameters[["Cluster_count"]][i_parameter]) {
            vec_K <- c(vec_K, parameters[[paste0("Cluster_proportion_validation_", i)]][i_parameter])
        }
        for (i_validation in 1:nrow(SFS_data_validation)) {
            logLikelihood <- compute_loglikelihood_ABC(
                A = A,
                vec_K = vec_K,
                component_distributions = component_distributions_validation[[i_parameter]],
                vec_SFS_real = SFS_data_validation[i_validation, ],
                zero_cutoff = zero_cutoff
            )
            all_criteria[nrow(all_criteria) + 1, ] <- c(i_parameter, i_validation, logLikelihood)
        }
    }
    all_criteria$AIC <- compute_AIC(
        logLikelihood = all_criteria$logLikelihood,
        num_parameters = num_parameters
    )
    all_criteria$BIC <- compute_BIC(
        logLikelihood = all_criteria$logLikelihood,
        num_parameters = num_parameters,
        num_samples = sum(SFS_data_validation[1, ])
    )
    #---Compute stopping criteria mean for each parameter set
    for (col_name in setdiff(colnames(all_criteria), c("parameter_id", "validation_id", "logLikelihood"))) {
        for (parameter_id in 1:nrow(parameters)) {
            parameters[[col_name]][parameter_id] <- mean(all_criteria[[col_name]][all_criteria$parameter_id == parameter_id])
            parameters[[paste0(col_name, "_sd")]][parameter_id] <- sd(all_criteria[[col_name]][all_criteria$parameter_id == parameter_id])
        }
    }
    #---Compute stopping criteria mean and standard deviation
    criteria <- data.frame(matrix(ncol = 0, nrow = 1))
    for (col_name in setdiff(colnames(all_criteria), c("parameter_id", "validation_id", "logLikelihood"))) {
        criteria[[col_name]] <- mean(all_criteria[[col_name]])
        criteria[[paste0(col_name, "_sd")]] <- sd(all_criteria[[col_name]])
    }
    #---Return the stopping criteria results
    output <- list()
    output$parameters <- parameters
    output$criteria <- criteria
    return(output)
}

parameter_conversion_ABC <- function(result,
                                     output_parameters_df = TRUE,
                                     mutation_count_for_fitting,
                                     sample_size,
                                     matrix_binomial_sample_size,
                                     matrix_binomial_ploidy) {
    tail_status <- result$best_fit$tail_status
    parameters_inference_A <- result$best_fit$parameters_inference_A
    parameters_inference_B <- result$best_fit$parameters_inference_B
    parameters_validation <- result$best_fit$parameters_validation
    component_distributions_inference_A <- result$best_fit$component_distributions_inference_A
    component_distributions_inference_B <- result$best_fit$component_distributions_inference_B
    component_distributions_validation <- result$best_fit$component_distributions_validation
    func_get_parameters <- function(parameters, tail_status) {
        if (tail_status) {
            vec_A <- parameters[1:2]
            N_humps <- length(parameters) / 2 - 1
            ii <- 0
        } else {
            vec_A <- c(NA, NA)
            N_humps <- length(parameters) / 2
            ii <- -1
        }
        if (N_humps == 0) {
            vec_p <- c()
            vec_K <- c()
        } else {
            vec_p <- parameters[seq(4 + 2 * ii, length(parameters), by = 2)]
            sorted_indices <- order(vec_p, decreasing = TRUE)
            vec_p <- vec_p[sorted_indices]
            vec_K <- parameters[seq(3 + 2 * ii, length(parameters), by = 2)]
            vec_K <- vec_K[sorted_indices]
        }
        output <- list()
        output$vec_A <- vec_A
        output$vec_p <- vec_p
        output$vec_K <- vec_K
        output$N_humps <- N_humps
        output$tail_status <- tail_status
        return(output)
    }
    output <- list()
    output$inference_A <- func_get_parameters(parameters_inference_A, tail_status)
    output$inference_B <- func_get_parameters(parameters_inference_B, tail_status)
    output$validation <- func_get_parameters(parameters_validation, tail_status)
    if (output_parameters_df) {
        parameters_df <- data.frame()
        parameters_df[1, "Mutation_count_for_fitting"] <- mutation_count_for_fitting
        parameters_df[1, "Tail"] <- tail_status
        if (tail_status) {
            parameters_df[1, "Tail_power"] <- output$inference_A$vec_A[2]
            parameters_df[1, "Tail_mutcount_observed"] <-
                output$inference_A$vec_A[1] * mutation_count_for_fitting
            parameters_df[1, "Tail_mutcount_predicted"] <-
                output$inference_A$vec_A[1] * mutation_count_for_fitting *
                    sum(component_distributions_inference_A$SFS_exact[1, ]) /
                    sum(component_distributions_inference_A$SFS_expected[1, ]) *
                    sample_size / matrix_binomial_sample_size
        } else {
            parameters_df[1, "Tail_power"] <- NA
            parameters_df[1, "Tail_mutcount_observed"] <- NA
            parameters_df[1, "Tail_mutcount_predicted"] <- NA
        }
        parameters_df[1, "Cluster_count"] <- output$inference_A$N_humps
        if (output$inference_A$N_humps > 0) {
            for (k in 1:output$inference_A$N_humps) {
                parameters_df[1, paste0("Cluster_frequency_", k)] <- output$inference_A$vec_p[k] / matrix_binomial_ploidy
                parameters_df[1, paste0("Cluster_mutcount_observed_", k)] <-
                    output$inference_A$vec_K[k] * mutation_count_for_fitting
                parameters_df[1, paste0("Cluster_mutcount_predicted_", k)] <-
                    output$inference_A$vec_K[k] * mutation_count_for_fitting /
                        sum(component_distributions_inference_A$SFS_expected[k + 1, ]) /
                        sum(component_distributions_inference_A$SFS_exact[k + 1, ])
            }
        }
        output$parameters_df <- parameters_df
    }
    return(output)
}

compute_loglikelihood_ABC <- function(A, vec_K, component_distributions, vec_SFS_real, zero_cutoff) {
    #----------------Compute the SFS probability distribution from model
    vec_SFS_model <- compute_SFS_ABC(
        A = A,
        vec_K = vec_K,
        component_distributions = component_distributions
    )
    vec_SFS_model[which(vec_SFS_model <= zero_cutoff)] <- zero_cutoff
    vec_SFS_model_normalized <- vec_SFS_model / sum(vec_SFS_model)
    #-----------------------------Compute the log-likelihood(data|model)
    loglikelihood <- sum(log(vec_SFS_model_normalized) * vec_SFS_real)
    return(loglikelihood)
}

compute_SFS_ABC <- function(A, vec_K, component_distributions) {
    # 	Add the neutral component
    vec_SFS_model <- component_distributions$SFS_expected_normalized[1, ]
    # print(component_distributions$SFS_expected_normalized)
    if (!is.na(A)) vec_SFS_model <- A * vec_SFS_model
    # 	Add the binomial humps
    for (i_hump in seq_along(vec_K)) {
        vec_SFS_model <- vec_SFS_model + vec_K[i_hump] * component_distributions$SFS_expected_normalized[i_hump + 1, ]
    }
    # 	Return the full SFS
    return(vec_SFS_model)
}

build_SFS_library_ABC <- function(neutral_power, cluster_frequencies, sfs_bincount, SFS_convolution_matrix, N_end) {
    SFS_exact <- c()
    SFS_expected <- c()
    SFS_expected_normalized <- c()
    #   Build the neutral component
    if (is.na(neutral_power)) {
        vec_SFS_GT <- numeric(N_end)
        vec_SFS_expected <- rep(0, sfs_bincount)
        vec_SFS_expected_normalized <- rep(0, sfs_bincount)
    } else {
        vec_para <- c(1, neutral_power)
        vec_SFS_GT <- build_SFS_library_Griffiths_Tavare_ABC(vec_para = vec_para, N_end = N_end)
        vec_SFS_expected <- rep(0, sfs_bincount)
        for (j in 1:sfs_bincount) {
            vec_SFS_expected[j] <- sum(vec_SFS_GT * SFS_convolution_matrix[, j])
        }
        vec_SFS_expected_normalized <- vec_SFS_expected / sum(vec_SFS_expected)
    }
    SFS_exact <- rbind(SFS_exact, vec_SFS_GT)
    SFS_expected <- rbind(SFS_expected, vec_SFS_expected)
    SFS_expected_normalized <- rbind(SFS_expected_normalized, vec_SFS_expected_normalized)
    #   Build the cluster components
    if (length(cluster_frequencies) > 0) {
        for (i in 1:length(cluster_frequencies)) {
            vec_para <- c(0, 0, 1, cluster_frequencies[i])
            vec_SFS_GT <- build_SFS_library_Griffiths_Tavare_ABC(vec_para = vec_para, N_end = N_end)
            vec_SFS_expected <- rep(0, sfs_bincount)
            for (j in 1:sfs_bincount) {
                vec_SFS_expected[j] <- sum(vec_SFS_GT * SFS_convolution_matrix[, j])
            }
            SFS_exact <- rbind(SFS_exact, vec_SFS_GT)
            SFS_expected <- rbind(SFS_expected, vec_SFS_expected)
            SFS_expected_normalized <- rbind(SFS_expected_normalized, vec_SFS_expected / sum(vec_SFS_expected))
        }
    }
    #   Ensure library is in matrix format, where each row = 1 component
    if (is.vector(SFS_exact)) SFS_exact <- matrix(SFS_exact, nrow = 1)
    if (is.vector(SFS_expected)) SFS_expected <- matrix(SFS_expected, nrow = 1)
    if (is.vector(SFS_expected_normalized)) SFS_expected_normalized <- matrix(SFS_expected_normalized, nrow = 1)
    #   Return SFS component library
    component_distributions <- list()
    component_distributions$SFS_exact <- SFS_exact
    component_distributions$SFS_expected <- SFS_expected
    component_distributions$SFS_expected_normalized <- SFS_expected_normalized
    return(component_distributions)
}

build_SFS_library_Griffiths_Tavare_ABC <- function(vec_para, N_end) {
    #---Get the parameters
    no_hump <- (length(vec_para) - 2) / 2
    para_A <- vec_para[1]
    para_alpha <- vec_para[2]
    if (no_hump > 0) {
        para_K <- numeric(no_hump)
        para_P <- numeric(no_hump)
        for (i in 1:no_hump) {
            para_K[i] <- vec_para[2 * i + 1]
            para_P[i] <- vec_para[2 * i + 2]
        }
    } else {
        para_K <- NULL
        para_P <- NULL
    }
    #---Compute the Griffiths-Tavare SFS
    vec_SFS_GT <- numeric(N_end)
    # for (m in 2:N_end) {
    for (m in 1:N_end) {
        # if (para_alpha > 0) vec_SFS_GT[m] <- para_A * N_end / (m^para_alpha)
        if (para_alpha > 0) vec_SFS_GT[m] <- para_A / (m^para_alpha)
        # vec_SFS_GT[m] <- para_A * N_end / (m * (m - 1))
        if (no_hump > 0) {
            for (i in 1:no_hump) {
                K <- para_K[i]
                P <- para_P[i]
                vec_SFS_GT[m] <- vec_SFS_GT[m] + K * dbinom(m, N_end, P)
            }
        }
    }
    return(vec_SFS_GT)
}

build_convolution_matrix_ABC <- function(sfs_bincount,
                                         mode = NULL,
                                         sample_size,
                                         min_variant_read,
                                         min_total_read,
                                         max_total_read,
                                         matrix_binomial_ploidy,
                                         coverage_distribution,
                                         coverage_variables,
                                         sample_coverage,
                                         compute_parallel = FALSE,
                                         n_cores = NULL) {
    suppressPackageStartupMessages(library(progress))
    report <- "Prepare the SFS convolution matrix"
    if (!is.null(mode)) report <- paste0(report, " for ", mode)
    cat(bold(blue(paste0(report, "...\n"))))
    N_end <- sample_size
    SFS_totalsteps_base <- sfs_bincount
    #---Compute the total readcount PDF
    if (coverage_distribution == "uniform") {
        sample_coverage_distribution <- data.frame(total_readcount = min_total_read:max_total_read)
        sample_coverage_distribution$pdf <- 1 / nrow(sample_coverage_distribution)
    } else if (coverage_distribution == "binomial") {
        sample_coverage_distribution <- data.frame(total_readcount = min_total_read:max_total_read)
        sample_coverage_distribution$pdf <- dbinom(sample_coverage_distribution$total_readcount, size = N_end, prob = coverage_variables$mean / N_end)
    } else if (coverage_distribution == "sample-specific") {
        sample_coverage_distribution <- sample_coverage
    }
    #---Function to compute the convolution vector for a given SFS bin
    func_convolution_vector <- function(N_end,
                                        SFS_totalsteps_base,
                                        min_variant_read,
                                        matrix_binomial_ploidy,
                                        sample_coverage_distribution,
                                        sfs_bincount,
                                        i) {
        r <- sample_coverage_distribution$total_readcount
        r1 <- r * (i - 1) / SFS_totalsteps_base
        r1 <- ifelse(r1 %% 1 == 0, r1 + 1, ceiling(r1))
        r1 <- pmax(r1, min_variant_read, 1)
        r2 <- floor(r * i / SFS_totalsteps_base)
        prob <- (1:N_end) / (N_end * matrix_binomial_ploidy)
        vec_convolution <- rep(0, N_end)
        for (m in 1:N_end) {
            B_values <- ifelse(r1 > r2, 0,
                pbinom(r2, size = r, prob = prob[m]) - pbinom(r1 - 1, size = r, prob = prob[m])
            )
            vec_convolution[m] <- sum(sample_coverage_distribution$pdf * B_values)
        }
        return(vec_convolution)
    }
    #---Build convolution matrix to transform Griffiths-Tavare SFS to expected SFS
    if (compute_parallel == FALSE) {
        mat_convolution <- matrix(0, nrow = N_end, ncol = sfs_bincount)
        pb <- txtProgressBar(
            min = 0,
            max = sfs_bincount,
            style = 3,
            width = 50,
            char = "+"
        )
        for (i in 1:sfs_bincount) {
            setTxtProgressBar(pb, i)
            mat_convolution[, i] <- func_convolution_vector(
                N_end = N_end,
                SFS_totalsteps_base = SFS_totalsteps_base,
                min_variant_read = min_variant_read,
                matrix_binomial_ploidy = matrix_binomial_ploidy,
                sample_coverage_distribution = sample_coverage_distribution,
                sfs_bincount = sfs_bincount,
                i = i
            )
        }
        cat("\n")
    } else {
        suppressPackageStartupMessages(library(parallel))
        suppressPackageStartupMessages(library(pbapply))
        #   Start parallel cluster
        numCores <- ifelse(is.null(n_cores), detectCores(), n_cores)
        cl <- makePSOCKcluster(numCores - 1)
        clusterExport(cl, varlist = c(
            "N_end", "SFS_totalsteps_base", "min_variant_read",
            "matrix_binomial_ploidy", "sample_coverage_distribution", "sfs_bincount", "func_convolution_vector"
        ), envir = environment())
        #   Compute each convolution vector
        output <- pblapply(cl = cl, X = 1:sfs_bincount, FUN = function(i) {
            return(func_convolution_vector(
                N_end = N_end,
                SFS_totalsteps_base = SFS_totalsteps_base,
                min_variant_read = min_variant_read,
                matrix_binomial_ploidy = matrix_binomial_ploidy,
                sample_coverage_distribution = sample_coverage_distribution,
                sfs_bincount = sfs_bincount,
                i = i
            ))
        })
        stopCluster(cl)
        mat_convolution <- do.call(cbind, output)
    }
    #---Create the convolution matrix output
    output <- list()
    output$convolution_matrix <- mat_convolution
    output$N_end <- N_end
    output$SFS_totalsteps_base <- SFS_totalsteps_base
    output$min_total_read <- min_total_read
    output$max_total_read <- max_total_read
    output$coverage_distribution <- coverage_distribution
    output$coverage_variables <- coverage_variables
    return(output)
}
