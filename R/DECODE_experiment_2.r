DECODE_experiment_2 <- function(sample_id = "",
                                mutation_table,
                                criterion = "BIC",
                                criterion_ratio = 0.999, # <<<<<<<<<<<<<<
                                neutral_power_min = 0.5,
                                neutral_power_max = 5,
                                cluster_frequency_min = 0.01,
                                cluster_frequency_max = 1,
                                max_total_read = NULL,
                                sample_size = 1000,
                                sfs_bincount,
                                coverage_distribution = "sample-specific",
                                coverage_variables = NULL,
                                N_trials = 10000,
                                neutral_tail = NA,
                                min_N_humps = 1,
                                max_N_humps = Inf,
                                pi_cutoff = 0.02,
                                zero_cutoff = 1e-50,
                                compute_parallel = TRUE,
                                n_cores = NULL) {
    suppressWarnings(library(crayon))
    cat(paste0("\n\n\n", bold(red("PERFORMING DECODE FOR SAMPLE ")), bold(yellow(sample_id)), bold(red("...")), "\n"))
    mutation_table$Tot_count <- mutation_table$Ref_count + mutation_table$Alt_count
    mutation_table$VAF <- mutation_table$Alt_count / mutation_table$Tot_count
    if (is.null(max_total_read)) max_total_read <- max(mutation_table$Tot_count)
    ####################################################################
    ####################################################################
    ####################################################################
    matrix_binomial_sample_size <- 1000
    matrix_binomial_sfs_bincount <- sfs_bincount
    matrix_binomial_ploidy <- 2
    #---Choose mutation thresholds, get resulting SFS from data
    SFS_data_frequencies <- seq(1, sfs_bincount) / sfs_bincount
    results <- choose_mutation_thresholds_experiment_2(
        mutation_table = mutation_table,
        max_total_read = max_total_read,
        SFS_data_frequencies = SFS_data_frequencies
    )
    readcount_distribution <- results$readcount_distribution
    min_variant_read_inference_A <- results$min_variant_read_inference_A
    min_total_read_inference_A <- results$min_total_read_inference_A
    SFS_data_inference_A <- results$SFS_data_inference_A
    sample_coverage_inference_A <- results$sample_coverage_inference_A
    min_variant_read_inference_B <- results$min_variant_read_inference_B
    min_total_read_inference_B <- results$min_total_read_inference_B
    SFS_data_inference_B <- results$SFS_data_inference_B
    sample_coverage_inference_B <- results$sample_coverage_inference_B
    min_variant_read_validation <- results$min_variant_read_validation
    min_total_read_validation <- results$min_total_read_validation
    SFS_data_validation <- results$SFS_data_validation
    sample_coverage_validation <- results$sample_coverage_validation
    ####################################################################
    ####################################################################
    ####################################################################
    #---Prepare the SFS convolution matrix
    SFS_convolution_inference_A <- build_convolution_matrix_experiment_2(
        sfs_bincount = matrix_binomial_sfs_bincount,
        mode = "inference A",
        sample_size = matrix_binomial_sample_size,
        min_variant_read = min_variant_read_inference_A,
        min_total_read = min_total_read_inference_A,
        max_total_read = max_total_read,
        matrix_binomial_ploidy = matrix_binomial_ploidy,
        coverage_distribution = coverage_distribution,
        coverage_variables = coverage_variables,
        sample_coverage = sample_coverage_inference_A
    )
    SFS_convolution_inference_B <- build_convolution_matrix_experiment_2(
        sfs_bincount = matrix_binomial_sfs_bincount,
        mode = "inference B",
        sample_size = matrix_binomial_sample_size,
        min_variant_read = min_variant_read_inference_B,
        min_total_read = min_total_read_inference_B,
        max_total_read = max_total_read,
        matrix_binomial_ploidy = matrix_binomial_ploidy,
        coverage_distribution = coverage_distribution,
        coverage_variables = coverage_variables,
        sample_coverage = sample_coverage_inference_B
    )
    SFS_convolution_validation <- build_convolution_matrix_experiment_2(
        sfs_bincount = matrix_binomial_sfs_bincount,
        mode = "validation",
        sample_size = matrix_binomial_sample_size,
        min_variant_read = min_variant_read_validation,
        min_total_read = min_total_read_validation,
        max_total_read = max_total_read,
        matrix_binomial_ploidy = matrix_binomial_ploidy,
        coverage_distribution = coverage_distribution,
        coverage_variables = coverage_variables,
        sample_coverage = sample_coverage_validation
    )
    #---DECODE
    DECODE_result <- list()
    DECODE_result$sample_id <- sample_id
    DECODE_result$criterion <- criterion
    DECODE_result$criterion_ratio <- criterion_ratio
    DECODE_result$sfs_bincount <- sfs_bincount
    DECODE_result$mutational_table <- mutation_table
    DECODE_result$readcount_distribution <- readcount_distribution
    DECODE_result$min_variant_read_inference_A <- min_variant_read_inference_A
    DECODE_result$min_total_read_inference_A <- min_total_read_inference_A
    DECODE_result$min_variant_read_inference_B <- min_variant_read_inference_B
    DECODE_result$min_total_read_inference_B <- min_total_read_inference_B
    DECODE_result$min_variant_read_validation <- min_variant_read_validation
    DECODE_result$min_total_read_validation <- min_total_read_validation
    DECODE_result$max_total_read <- max_total_read
    DECODE_result$SFS_frequencies <- SFS_data_frequencies
    DECODE_result$SFS_data_inference_A <- SFS_data_inference_A
    DECODE_result$SFS_data_inference_B <- SFS_data_inference_B
    DECODE_result$SFS_data_validation <- SFS_data_validation
    if (is.na(neutral_tail)) {
        # result_with_tail <- DECODE_given_tail_status(
        #     vec_SFS_real = vec_SFS_real,
        #     sfs_bincount = sfs_bincount,
        #     with_tail = TRUE,
        #     criterion = criterion,
        #     criterion_ratio = criterion_ratio,
        #     min_N_humps = min_N_humps,
        #     max_N_humps = max_N_humps,
        #     N_trials = N_trials,
        #     SFS_convolution = SFS_convolution,
        #     neutral_power_min = neutral_power_min,
        #     neutral_power_max = neutral_power_max,
        #     cluster_frequency_min = cluster_frequency_min,
        #     cluster_frequency_max = cluster_frequency_max,
        #     pi_cutoff = pi_cutoff,
        #     zero_cutoff = zero_cutoff,
        #     compute_parallel = compute_parallel,
        #     n_cores = n_cores
        # )
        # DECODE_result$fits_with_tail <- result_with_tail
        # result_without_tail <- DECODE_given_tail_status(
        #     vec_SFS_real = vec_SFS_real,
        #     sfs_bincount = sfs_bincount,
        #     with_tail = FALSE,
        #     criterion = criterion,
        #     criterion_ratio = criterion_ratio,
        #     min_N_humps = min_N_humps,
        #     max_N_humps = max_N_humps,
        #     N_trials = N_trials,
        #     SFS_convolution = SFS_convolution,
        #     neutral_power_min = neutral_power_min,
        #     neutral_power_max = neutral_power_max,
        #     cluster_frequency_min = cluster_frequency_min,
        #     cluster_frequency_max = cluster_frequency_max,
        #     pi_cutoff = pi_cutoff,
        #     zero_cutoff = zero_cutoff,
        #     compute_parallel = compute_parallel,
        #     n_cores = n_cores
        # )
        # DECODE_result$fits_without_tail <- result_without_tail
        # if (result_with_tail$best_fit$selected_criterion_value < result_without_tail$best_fit$selected_criterion_value) {
        #     final_result <- result_with_tail
        # } else {
        #     final_result <- result_without_tail
        # }
    } else if (neutral_tail == TRUE) {
        final_result <- DECODE_given_tail_status_experiment_2(
            SFS_data_inference_A = SFS_data_inference_A,
            SFS_data_inference_B = SFS_data_inference_B,
            SFS_data_validation = SFS_data_validation,
            sfs_bincount = sfs_bincount,
            with_tail = TRUE,
            criterion = criterion,
            criterion_ratio = criterion_ratio,
            min_N_humps = min_N_humps,
            max_N_humps = max_N_humps,
            N_trials = N_trials,
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
    } else if (neutral_tail == FALSE) {
        # final_result <- DECODE_given_tail_status(
        #     vec_SFS_real = vec_SFS_real,
        #     sfs_bincount = sfs_bincount,
        #     with_tail = FALSE,
        #     criterion = criterion,
        #     criterion_ratio = criterion_ratio,
        #     min_N_humps = min_N_humps,
        #     max_N_humps = max_N_humps,
        #     N_trials = N_trials,
        #     SFS_convolution = SFS_convolution,
        #     neutral_power_min = neutral_power_min,
        #     neutral_power_max = neutral_power_max,
        #     cluster_frequency_min = cluster_frequency_min,
        #     cluster_frequency_max = cluster_frequency_max,
        #     pi_cutoff = pi_cutoff,
        #     zero_cutoff = zero_cutoff,
        #     compute_parallel = compute_parallel,
        #     n_cores = n_cores
        # )
        # DECODE_result$fits_without_tail <- final_result
    }
    DECODE_result$final_fit <- final_result
    #---Report the best fit
    tail_status_final_result <- final_result$best_fit$tail_status
    parameters_inference_A_final_result <- final_result$best_fit$parameters_inference_A
    parameters_inference_B_final_result <- final_result$best_fit$parameters_inference_B
    criterion_final_result <- final_result$best_fit$selected_criterion_value
    if (tail_status_final_result) {
        N_humps_final_result <- length(parameters_inference_A_final_result) / 2 - 1
        report <- bold(underline(red(paste0("Best fit = neutral tail + ", N_humps_final_result, " clusters:\n"))))
        report <- paste0(report, red("Score            : "), yellow(paste0(criterion, " = ", format(round(criterion_final_result, 3), nsmall = 3))), "\n")
        report <- paste0(report, red("Neutral component: "), yellow(paste0("power     = ", format(round(parameters_inference_A_final_result[2], 3), nsmall = 3))), red(", "))
        report <- paste0(report, yellow(paste0("\u03C0 = ", format(round(parameters_inference_A_final_result[1], 3), nsmall = 3), " [A], ", format(round(parameters_inference_B_final_result[1], 3), nsmall = 3), " [B]")), "\n")
        ii <- 0
    } else {
        N_humps_final_result <- length(parameters_inference_A_final_result) / 2
        report <- bold(underline(red(paste0("Best fit = no neutral tail + ", N_humps_final_result, " clusters:\n"))))
        report <- paste0(report, red("Score            : "), yellow(paste0(criterion, " = ", format(round(criterion_final_result, 3), nsmall = 3))), "\n")
        ii <- -1
    }
    if (N_humps_final_result > 0) {
        for (i in 1:N_humps_final_result) {
            report <- paste0(report, red(paste0("Cluster ", i, "        : ")), yellow(paste0("frequency = ", format(round(parameters_inference_A_final_result[2 * (i + ii) + 2], 3), nsmall = 3))), red(", "))
            report <- paste0(report, yellow(paste0("\u03C0 = ", format(round(parameters_inference_A_final_result[2 * (i + ii) + 1], 3), nsmall = 3), " [A], ", format(round(parameters_inference_B_final_result[2 * (i + ii) + 1], 3), nsmall = 3), " [B]")), "\n")
        }
    }
    cat(report)
    #---Translation to parameters of cancer evolution in the sample
    tmp <- parameter_conversion_experiment_2(
        result = final_result,
        mutation_count_for_fitting = sum(SFS_data_inference_A),
        sample_size = sample_size,
        matrix_binomial_sample_size = matrix_binomial_sample_size,
        matrix_binomial_ploidy = matrix_binomial_ploidy
    )
    DECODE_result$final_fit$parameters_df <- tmp$parameters_df
    #---Return the SFS deconvolution results
    return(DECODE_result)
}

DECODE_given_tail_status_experiment_2 <- function(SFS_data_inference_A,
                                                  SFS_data_inference_B,
                                                  SFS_data_validation,
                                                  criterion,
                                                  criterion_ratio,
                                                  min_N_humps,
                                                  max_N_humps,
                                                  with_tail,
                                                  N_trials,
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
        fit_results <- DECODE_given_tail_status_and_Ncluster_experiment_2(
            SFS_data_inference_A = SFS_data_inference_A,
            SFS_data_inference_B = SFS_data_inference_B,
            SFS_data_validation = SFS_data_validation,
            N_humps = N_humps,
            with_tail = with_tail,
            N_trials = N_trials,
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
        all_fits[[paste0(N_humps, "_clusters")]] <- fit_results
        parameters_inference_A_best_current <- fit_results$best$parameters_inference_A
        parameters_inference_B_best_current <- fit_results$best$parameters_inference_B
        parameters_validation_best_current <- fit_results$best$parameters_validation
        component_distributions_inference_A_best_current <- fit_results$best$component_distributions_inference_A
        component_distributions_inference_B_best_current <- fit_results$best$component_distributions_inference_B
        component_distributions_validation_best_current <- fit_results$best$component_distributions_validation
        criterion_all_best_current <- fit_results$best$criteria
        criterion_best_current <- criterion_all_best_current[[criterion]]
        #   Report the best fit for the current hump count
        cluster_pis_inference_A <- Inf
        cluster_pis_inference_B <- Inf
        report <- paste0(blue("Score            : "), cyan(paste0(criterion, " = ", format(round(criterion_best_current, 3), nsmall = 3))), "\n")
        if (with_tail) {
            N_humps <- length(parameters_inference_A_best_current) / 2 - 1
            report <- paste0(report, blue("Neutral component: "), cyan(paste0("power     = ", format(round(parameters_inference_A_best_current[2], 3), nsmall = 3))), blue(", "))
            report <- paste0(report, cyan(paste0("\u03C0 = ", format(round(parameters_inference_A_best_current[1], 3), nsmall = 3), " [A], ", format(round(parameters_inference_B_best_current[1], 3), nsmall = 3), " [B]")), "\n")
            ii <- 0
        } else {
            N_humps <- length(parameters_inference_A_best_current) / 2
            ii <- -1
        }
        if (N_humps > 0) {
            for (i in 1:N_humps) {
                report <- paste0(report, blue(paste0("Cluster ", i, "        : ")), cyan(paste0("frequency = ", format(round(parameters_inference_A_best_current[2 * (i + ii) + 2], 3), nsmall = 3))), blue(", "))
                report <- paste0(report, cyan(paste0("\u03C0 = ", format(round(parameters_inference_A_best_current[2 * (i + ii) + 1], 3), nsmall = 3), " [A], ", format(round(parameters_inference_B_best_current[2 * (i + ii) + 1], 3), nsmall = 3), " [B]")), "\n")
                cluster_pis_inference_A <- c(cluster_pis_inference_A, parameters_inference_A_best_current[2 * (i + ii) + 1])
                cluster_pis_inference_B <- c(cluster_pis_inference_B, parameters_inference_B_best_current[2 * (i + ii) + 1])
            }
        }
        cat(report)
        #   Check if the increased hump count leads to lower criterion score without tiny selective components...
        if ((N_humps == min_N_humps) | ((criterion_best_current < criterion_ratio * criterion_best_final) & (min(pmax(cluster_pis_inference_A, cluster_pis_inference_B)) >= pi_cutoff))) {
            #   ... if yes, then update the best fit and continue with 1 more hump
            fit_results_best_final <- fit_results
            criterion_best_final <- criterion_best_current
            criterion_all_final <- criterion_all_best_current
            parameters_inference_A_best_final <- parameters_inference_A_best_current
            parameters_inference_B_best_final <- parameters_inference_B_best_current
            parameters_validation_best_final <- parameters_validation_best_current
            component_distributions_inference_A_best_final <- component_distributions_inference_A_best_current
            component_distributions_inference_B_best_final <- component_distributions_inference_B_best_current
            component_distributions_validation_best_final <- component_distributions_validation_best_current
            N_humps <- N_humps + 1
            #   ... except if exceeding maximum number of clusters
            if (N_humps > max_N_humps) break
        } else {
            #   ... if no, then stop
            break
        }
    }
    #---Check if the neutral tail component is too tiny
    if (with_tail == TRUE & (max(parameters_inference_A_best_final[1], parameters_inference_B_best_final[1]) < pi_cutoff)) {
        with_tail <- FALSE
        parameters_inference_A_best_final[seq(3, length(parameters_inference_A_best_final), by = 2)] <- parameters_inference_A_best_final[seq(3, length(parameters_inference_A_best_final), by = 2)] / sum(parameters_inference_A_best_final[seq(3, length(parameters_inference_A_best_final), by = 2)])
        parameters_inference_A_best_final <- parameters_inference_A_best_final[-c(1, 2)]
        parameters_inference_B_best_final[seq(3, length(parameters_inference_B_best_final), by = 2)] <- parameters_inference_B_best_final[seq(3, length(parameters_inference_B_best_final), by = 2)] / sum(parameters_inference_B_best_final[seq(3, length(parameters_inference_B_best_final), by = 2)])
        parameters_inference_B_best_final <- parameters_inference_B_best_final[-c(1, 2)]
        parameters_validation_best_final[seq(3, length(parameters_validation_best_final), by = 2)] <- parameters_validation_best_final[seq(3, length(parameters_validation_best_final), by = 2)] / sum(parameters_validation_best_final[seq(3, length(parameters_validation_best_final), by = 2)])
        parameters_validation_best_final <- parameters_validation_best_final[-c(1, 2)]
        component_distributions_inference_A_best_final$SFS_exact[1, ] <- rep(0, length(component_distributions_inference_A_best_final$SFS_exact[1, ]))
        component_distributions_inference_A_best_final$SFS_expected[1, ] <- rep(0, length(component_distributions_inference_A_best_final$SFS_expected[1, ]))
        component_distributions_inference_A_best_final$SFS_expected_normalized[1, ] <- rep(0, length(component_distributions_inference_A_best_final$SFS_expected_normalized[1, ]))
        component_distributions_inference_B_best_final$SFS_exact[1, ] <- rep(0, length(component_distributions_inference_B_best_final$SFS_exact[1, ]))
        component_distributions_inference_B_best_final$SFS_expected[1, ] <- rep(0, length(component_distributions_inference_B_best_final$SFS_expected[1, ]))
        component_distributions_inference_B_best_final$SFS_expected_normalized[1, ] <- rep(0, length(component_distributions_inference_B_best_final$SFS_expected_normalized[1, ]))
        component_distributions_validation_best_final$SFS_exact[1, ] <- rep(0, length(component_distributions_validation_best_final$SFS_exact[1, ]))
        component_distributions_validation_best_final$SFS_expected[1, ] <- rep(0, length(component_distributions_validation_best_final$SFS_expected[1, ]))
        component_distributions_validation_best_final$SFS_expected_normalized[1, ] <- rep(0, length(component_distributions_validation_best_final$SFS_expected_normalized[1, ]))
    }
    #---Report the best fit
    result <- list()
    result$all_fits <- all_fits
    result$best_fit <- list()
    result$best_fit$parameters_inference_A <- parameters_inference_A_best_final
    result$best_fit$parameters_inference_B <- parameters_inference_B_best_final
    result$best_fit$parameters_validation <- parameters_validation_best_final
    result$best_fit$component_distributions_inference_A <- component_distributions_inference_A_best_final
    result$best_fit$component_distributions_inference_B <- component_distributions_inference_B_best_final
    result$best_fit$component_distributions_validation <- component_distributions_validation_best_final
    result$best_fit$selected_criterion <- criterion
    result$best_fit$all_criteria <- criterion_all_final
    result$best_fit$selected_criterion_value <- criterion_best_final
    result$best_fit$tail_status <- with_tail
    return(result)
}

DECODE_given_tail_status_and_Ncluster_experiment_2 <- function(SFS_data_inference_A,
                                                               SFS_data_inference_B,
                                                               SFS_data_validation,
                                                               N_humps,
                                                               with_tail,
                                                               N_trials,
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
    func_one_trial <- function(with_tail,
                               N_humps,
                               neutral_power_min,
                               neutral_power_max,
                               cluster_frequency_min,
                               cluster_frequency_max,
                               SFS_data_inference_A,
                               SFS_data_inference_B,
                               SFS_data_validation,
                               sfs_bincount,
                               N_end,
                               SFS_convolution_matrix_inference_A,
                               SFS_convolution_matrix_inference_B,
                               SFS_convolution_matrix_validation,
                               zero_cutoff) {
        #   Sample neutral component power and cluster frequencies
        neutral_power <- ifelse(with_tail, runif(1, neutral_power_min, neutral_power_max), NA)
        cluster_frequencies <- sort(runif(N_humps, cluster_frequency_min, cluster_frequency_max), decreasing = TRUE)
        #   Build the SFS component library
        component_distributions_inference_A <- build_SFS_library(
            neutral_power = neutral_power,
            cluster_frequencies = cluster_frequencies,
            sfs_bincount = sfs_bincount,
            SFS_convolution_matrix = SFS_convolution_matrix_inference_A,
            N_end = N_end
        )
        component_distributions_inference_B <- build_SFS_library(
            neutral_power = neutral_power,
            cluster_frequencies = cluster_frequencies,
            sfs_bincount = sfs_bincount,
            SFS_convolution_matrix = SFS_convolution_matrix_inference_B,
            N_end = N_end
        )
        component_distributions_validation <- build_SFS_library(
            neutral_power = neutral_power,
            cluster_frequencies = cluster_frequencies,
            sfs_bincount = sfs_bincount,
            SFS_convolution_matrix = SFS_convolution_matrix_validation,
            N_end = N_end
        )
        #   Find component proportions
        results <- DECODE_for_pis_experiment_2(
            SFS_data_inference_A = SFS_data_inference_A,
            SFS_data_inference_B = SFS_data_inference_B,
            SFS_data_validation = SFS_data_validation,
            neutral_power = neutral_power,
            cluster_frequencies = cluster_frequencies,
            component_distributions_inference_A = component_distributions_inference_A,
            component_distributions_inference_B = component_distributions_inference_B,
            component_distributions_validation = component_distributions_validation,
            zero_cutoff = zero_cutoff
        )
        output <- list()
        output$logLikelihood <- results$log_L
        output$parameters_inference_A <- results$parameters_A
        output$parameters_inference_B <- results$parameters_B
        output$parameters_validation <- results$parameters_validation
        output$component_distributions_inference_A <- component_distributions_inference_A
        output$component_distributions_inference_B <- component_distributions_inference_B
        output$component_distributions_validation <- component_distributions_validation
        return(output)
    }
    #---Find best variable parameters (A & K's) for each fixed parameter set from many trials
    if (compute_parallel == FALSE) {
        all_logLikelihood <- c()
        all_parameters_inference_A <- c()
        all_parameters_inference_B <- c()
        all_parameters_validation <- c()
        all_component_distributions_inference_A <- list()
        all_component_distributions_inference_B <- list()
        all_component_distributions_validation <- list()
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
            trial_result <- func_one_trial(
                with_tail = with_tail,
                N_humps = N_humps,
                neutral_power_min = neutral_power_min,
                neutral_power_max = neutral_power_max,
                cluster_frequency_min = cluster_frequency_min,
                cluster_frequency_max = cluster_frequency_max,
                SFS_data_inference_A = SFS_data_inference_A,
                SFS_data_inference_B = SFS_data_inference_B,
                SFS_data_validation = SFS_data_validation,
                sfs_bincount = sfs_bincount,
                N_end = N_end,
                SFS_convolution_matrix_inference_A = SFS_convolution_matrix_inference_A,
                SFS_convolution_matrix_inference_B = SFS_convolution_matrix_inference_B,
                SFS_convolution_matrix_validation = SFS_convolution_matrix_validation,
                zero_cutoff = zero_cutoff
            )
            all_parameters_inference_A <- rbind(all_para, trial_result$parameters_inference_A)
            all_parameters_inference_B <- rbind(all_para, trial_result$parameters_inference_B)
            all_parameters_validation <- rbind(all_para, trial_result$parameters_validation)
            all_component_distributions_inference_A[[i]] <- trial_result$component_distributions_inference_A
            all_component_distributions_inference_B[[i]] <- trial_result$component_distributions_inference_B
            all_component_distributions_validation[[i]] <- trial_result$component_distributions_validation
            all_logLikelihood <- c(all_logLikelihood, trial_result$logLikelihood)
        }
        if (progress_bar) cat("\n")
    } else {
        suppressWarnings(library(parallel))
        suppressWarnings(library(pbapply))
        #   Start parallel cluster
        numCores <- ifelse(is.null(n_cores), detectCores(), n_cores)
        cl <- makePSOCKcluster(numCores - 1)
        #   Prepare input parameters
        clusterExport(cl, varlist = c(
            "with_tail", "N_humps", "sfs_bincount", "N_end", "zero_cutoff",
            "SFS_data_inference_A", "SFS_data_inference_B", "SFS_data_validation",
            "SFS_convolution_matrix_inference_A", "SFS_convolution_matrix_inference_B", "SFS_convolution_matrix_validation",
            "neutral_power_min", "neutral_power_max",
            "cluster_frequency_min", "cluster_frequency_max",
            "build_SFS_library", "build_SFS_library_Griffiths_Tavare",
            "DECODE_for_pis_experiment_2", "compute_loglikelihood", "compute_SFS"
        ), envir = environment())
        #   Find best variable parameters in parallel mode
        if (progress_bar) {
            output <- pblapply(cl = cl, X = 1:N_trials, FUN = function(i) {
                func_one_trial(
                    with_tail = with_tail,
                    N_humps = N_humps,
                    neutral_power_min = neutral_power_min,
                    neutral_power_max = neutral_power_max,
                    cluster_frequency_min = cluster_frequency_min,
                    cluster_frequency_max = cluster_frequency_max,
                    SFS_data_inference_A = SFS_data_inference_A,
                    SFS_data_inference_B = SFS_data_inference_B,
                    SFS_data_validation = SFS_data_validation,
                    sfs_bincount = sfs_bincount,
                    N_end = N_end,
                    SFS_convolution_matrix_inference_A = SFS_convolution_matrix_inference_A,
                    SFS_convolution_matrix_inference_B = SFS_convolution_matrix_inference_B,
                    SFS_convolution_matrix_validation = SFS_convolution_matrix_validation,
                    zero_cutoff = zero_cutoff
                )
            })
        } else {
            output <- parLapply(cl = cl, X = 1:N_trials, fun = function(i) {
                func_one_trial(
                    with_tail = with_tail,
                    N_humps = N_humps,
                    neutral_power_min = neutral_power_min,
                    neutral_power_max = neutral_power_max,
                    cluster_frequency_min = cluster_frequency_min,
                    cluster_frequency_max = cluster_frequency_max,
                    SFS_data_inference_A = SFS_data_inference_A,
                    SFS_data_inference_B = SFS_data_inference_B,
                    SFS_data_validation = SFS_data_validation,
                    sfs_bincount = sfs_bincount,
                    N_end = N_end,
                    SFS_convolution_matrix_inference_A = SFS_convolution_matrix_inference_A,
                    SFS_convolution_matrix_inference_B = SFS_convolution_matrix_inference_B,
                    SFS_convolution_matrix_validation = SFS_convolution_matrix_validation,
                    zero_cutoff = zero_cutoff
                )
            })
        }
        stopCluster(cl)
        #   Extract the results
        all_parameters_inference_A <- do.call(rbind, lapply(output, function(x) x$parameters_inference_A))
        all_parameters_inference_B <- do.call(rbind, lapply(output, function(x) x$parameters_inference_B))
        all_parameters_validation <- do.call(rbind, lapply(output, function(x) x$parameters_validation))
        all_component_distributions_inference_A <- lapply(output, function(x) x$component_distributions_inference_A)
        all_component_distributions_inference_B <- lapply(output, function(x) x$component_distributions_inference_B)
        all_component_distributions_validation <- lapply(output, function(x) x$component_distributions_validation)
        all_logLikelihood <- sapply(output, function(x) x$logLikelihood)
    }
    #---Compute stopping criteria for best fit
    best_index <- which.max(all_logLikelihood)
    num_parameters <- ifelse(with_tail, 2 * N_humps + 1, 2 * N_humps - 1)
    if (neutral_power_max > neutral_power_min) num_parameters <- num_parameters + 1
    criteria <- cluster_count_criteria_experiment_2(
        num_parameters = num_parameters,
        SFS_data_validation = SFS_data_validation,
        parameters_validation = all_parameters_validation[best_index, ],
        component_distributions_validation = all_component_distributions_validation[[best_index]],
        with_tail = with_tail,
        zero_cutoff = zero_cutoff
    )
    #---Find the best fit
    fit_results <- list()
    fit_results$all <- list()
    fit_results$all$parameters_inference_A <- all_parameters_inference_A
    fit_results$all$parameters_inference_B <- all_parameters_inference_B
    fit_results$all$parameters_validation <- all_parameters_validation
    fit_results$all$logLikelihood <- all_logLikelihood
    fit_results$best <- list()
    fit_results$best$parameters_inference_A <- all_parameters_inference_A[best_index, ]
    fit_results$best$parameters_inference_B <- all_parameters_inference_B[best_index, ]
    fit_results$best$parameters_validation <- all_parameters_validation[best_index, ]
    fit_results$best$logLikelihood <- all_logLikelihood[best_index]
    fit_results$best$criteria <- criteria
    fit_results$best$component_distributions_inference_A <- all_component_distributions_inference_A[[best_index]]
    fit_results$best$component_distributions_inference_B <- all_component_distributions_inference_B[[best_index]]
    fit_results$best$component_distributions_validation <- all_component_distributions_validation[[best_index]]
    return(fit_results)
}

DECODE_for_pis_experiment_2 <- function(SFS_data_inference_A,
                                        SFS_data_inference_B,
                                        SFS_data_validation,
                                        neutral_power,
                                        cluster_frequencies,
                                        component_distributions_inference_A,
                                        component_distributions_inference_B,
                                        component_distributions_validation,
                                        zero_cutoff) {
    #   Function for conversion of pi's between different thresholds (L,M)
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
    # 	Function for parameter transformation
    parameter_transform <- function(parameters) {
        vec_pi_A <- exp(parameters)
        vec_pi_A <- vec_pi_A / sum(vec_pi_A)
        vec_pi_B <- pi_conversion(
            from_pis = vec_pi_A,
            from_component_distributions = component_distributions_inference_A,
            to_component_distributions = component_distributions_inference_B,
            with_tail = !is.na(neutral_power)
        )
        vec_pi_validation <- pi_conversion(
            from_pis = vec_pi_A,
            from_component_distributions = component_distributions_inference_A,
            to_component_distributions = component_distributions_validation,
            with_tail = !is.na(neutral_power)
        )
        output <- list()
        output$vec_pi_A <- vec_pi_A
        output$vec_pi_B <- vec_pi_B
        output$vec_pi_validation <- vec_pi_validation
        return(output)
    }
    # 	Function for optimization
    func_fit <- function(parameters) {
        transformed_parameters <- parameter_transform(parameters)
        if (is.na(neutral_power)) {
            A_A <- NA
            A_B <- NA
            vec_K_A <- transformed_parameters$vec_pi_A
            vec_K_B <- transformed_parameters$vec_pi_B
        } else {
            A_A <- transformed_parameters$vec_pi_A[1]
            A_B <- transformed_parameters$vec_pi_B[1]
            vec_K_A <- transformed_parameters$vec_pi_A[-1]
            vec_K_B <- transformed_parameters$vec_pi_B[-1]
        }
        loglikelihood <-
            compute_loglikelihood(
                A = A_A,
                vec_K = vec_K_A,
                component_distributions = component_distributions_inference_A,
                vec_SFS_real = SFS_data_inference_A,
                zero_cutoff = zero_cutoff
            ) +
            compute_loglikelihood(
                A = A_B,
                vec_K = vec_K_B,
                component_distributions = component_distributions_inference_B,
                vec_SFS_real = SFS_data_inference_B,
                zero_cutoff = zero_cutoff
            )
        return(loglikelihood)
    }
    N_humps <- length(cluster_frequencies)
    if (N_humps == 0) {
        parameters <- c(1)
    } else {
        # 	Initial values for parameters
        if (is.na(neutral_power)) {
            parameters_initial <- rep(0, N_humps)
        } else {
            parameters_initial <- rep(0, N_humps + 1)
        }
        # 	Optimization using optim function
        optim_results <- optim(
            par = parameters_initial,
            fn = func_fit,
            method = "Nelder-Mead",
            control = list(fnscale = -1)
        )
        # 	Extract the results
        parameters <- optim_results$par
    }
    transformed_parameters <- parameter_transform(parameters)
    vec_pi_A <- transformed_parameters$vec_pi_A
    vec_pi_B <- transformed_parameters$vec_pi_B
    vec_pi_validation <- transformed_parameters$vec_pi_validation
    log_L <- func_fit(parameters)
    # 	Prepare the parameters to be returned
    if (is.na(neutral_power)) {
        parameters_A <- c()
        parameters_B <- c()
        parameters_validation <- c()
    } else {
        parameters_A <- c(vec_pi_A[1], neutral_power)
        parameters_B <- c(vec_pi_B[1], neutral_power)
        parameters_validation <- c(vec_pi_validation[1], neutral_power)
    }
    if (N_humps > 0) {
        for (i in 1:N_humps) {
            if (is.na(neutral_power)) {
                parameters_A <- c(parameters_A, vec_pi_A[i], cluster_frequencies[i])
                parameters_B <- c(parameters_B, vec_pi_B[i], cluster_frequencies[i])
                parameters_validation <- c(parameters_validation, vec_pi_validation[i], cluster_frequencies[i])
            } else {
                parameters_A <- c(parameters_A, vec_pi_A[i + 1], cluster_frequencies[i])
                parameters_B <- c(parameters_B, vec_pi_B[i + 1], cluster_frequencies[i])
                parameters_validation <- c(parameters_validation, vec_pi_validation[i + 1], cluster_frequencies[i])
            }
        }
    }
    output <- list()
    output$log_L <- log_L
    output$parameters_A <- parameters_A
    output$parameters_B <- parameters_B
    output$parameters_validation <- parameters_validation
    return(output)
}

choose_mutation_thresholds_experiment_2 <- function(mutation_table,
                                                    max_total_read,
                                                    SFS_data_frequencies,
                                                    compute_parallel = TRUE,
                                                    n_cores = NULL) {
    suppressWarnings(library(data.table))



    mutation_table <- as.data.table(mutation_table)



    #---Find joint distribution of variant and total readcounts
    cat(bold(blue("Choose mutation thresholds for inference and validation...\n")))
    report <- paste0(blue("Complete data    : "), cyan(paste0(min(mutation_table$Alt_count), " \u2264 variant reads, ", min(mutation_table$Tot_count), " \u2264 total reads \u2264 ", max(mutation_table$Tot_count), "; ", nrow(mutation_table), " mutations\n")))
    cat(report)
    vec_min_variant_read <- min(mutation_table$Alt_count):max(mutation_table$Alt_count)
    vec_min_total_read <- min(mutation_table$Tot_count):max_total_read
    mutation_table <- mutation_table[Tot_count <= max_total_read]
    func_dist_given_min_variant_read <- function(min_variant_read, mutation_table) {
        mutation_table <- mutation_table[Alt_count >= min_variant_read]
        df <- data.frame()
        for (min_total_read in vec_min_total_read) {
            df <- rbind(
                df,
                data.frame(
                    min_total_read = min_total_read,
                    min_variant_read = min_variant_read,
                    mutation_count = mutation_table[Tot_count >= min_total_read, .N]
                    # mutation_count = mutation_table[Alt_count >= min_variant_read & Tot_count >= min_total_read & Tot_count <= max_total_read, .N]
                    # mutation_count = sum(mutation_table$Alt_count >= min_variant_read & mutation_table$Tot_count >= min_total_read & mutation_table$Tot_count <= max_total_read)
                )
            )
        }
        return(df)
    }
    if (compute_parallel == FALSE) {
        readcount_distribution <- data.frame()
        pb <- txtProgressBar(
            min = 1,
            max = length(vec_min_variant_read),
            style = 3,
            width = 50,
            char = "+"
        )
        for (min_variant_read in vec_min_variant_read) {
            readcount_distribution <- rbind(
                readcount_distribution,
                func_dist_given_min_variant_read(min_variant_read, mutation_table)
            )
            setTxtProgressBar(pb, min_variant_read)
        }
    } else {
        suppressWarnings(library(parallel))
        suppressWarnings(library(pbapply))
        #   Start parallel cluster
        numCores <- ifelse(is.null(n_cores), detectCores(), n_cores)
        cl <- makePSOCKcluster(numCores - 1)
        #   Prepare input parameters
        clusterExport(cl, varlist = c("mutation_table"))
        clusterEvalQ(cl, library(data.table))
        #   Compute each sub-dataframe
        output <- pblapply(cl = cl, X = vec_min_variant_read, FUN = function(min_variant_read) {
            return(func_dist_given_min_variant_read(min_variant_read, mutation_table))
        })
        stopCluster(cl)
        readcount_distribution <- do.call(rbind, output)
    }
    readcount_distribution$freq <- 100 * readcount_distribution$mutation_count / mutation_table[, .N]
    #---Choose mutation thresholds for inference and validation
    target_inference_freq <- 75
    target_validation_Nmut <- 1000


    min_total_read_inference_A <- min(readcount_distribution$min_total_read)
    filtered_df <- readcount_distribution[readcount_distribution$min_total_read == min_total_read_inference_A & readcount_distribution$freq >= target_inference_freq, ]
    min_variant_read_inference_A <- ifelse(nrow(filtered_df) > 0, max(filtered_df$min_variant_read), min(readcount_distribution$min_variant_read))
    Nmut_inference_A <- readcount_distribution$mutation_count[readcount_distribution$min_variant_read == min_variant_read_inference_A & readcount_distribution$min_total_read == min_total_read_inference_A]
    freq_inference_A <- readcount_distribution$freq[readcount_distribution$min_variant_read == min_variant_read_inference_A & readcount_distribution$min_total_read == min_total_read_inference_A]



    min_variant_read_inference_B <- min(readcount_distribution$min_variant_read) + 1
    filtered_df <- readcount_distribution[readcount_distribution$min_variant_read == min_variant_read_inference_B & readcount_distribution$freq >= target_inference_freq, ]
    min_total_read_inference_B <- ifelse(nrow(filtered_df) > 0, max(filtered_df$min_total_read), min(readcount_distribution$min_total_read))
    Nmut_inference_B <- readcount_distribution$mutation_count[readcount_distribution$min_variant_read == min_variant_read_inference_B & readcount_distribution$min_total_read == min_total_read_inference_B]
    freq_inference_B <- readcount_distribution$freq[readcount_distribution$min_variant_read == min_variant_read_inference_B & readcount_distribution$min_total_read == min_total_read_inference_B]



    min_variant_read_validation <- floor(0.5 * (min_variant_read_inference_A + min_variant_read_inference_B))
    min_total_read_validation <- floor(0.5 * (min_total_read_inference_A + min_total_read_inference_B))
    Nmut_validation <- readcount_distribution$mutation_count[readcount_distribution$min_variant_read == min_variant_read_validation & readcount_distribution$min_total_read == min_total_read_validation]
    freq_validation <- readcount_distribution$freq[readcount_distribution$min_variant_read == min_variant_read_validation & readcount_distribution$min_total_read == min_total_read_validation]



    report <- paste0(blue("Inference A      : "), cyan(paste0(min_variant_read_inference_A, " \u2264 variant reads, ", min_total_read_inference_A, " \u2264 total reads \u2264 ", max_total_read, "; ", Nmut_inference_A, " mutations (", format(round(freq_inference_A, 3), nsmall = 3), "%)")), "\n")
    report <- paste0(report, blue("Inference B      : "), cyan(paste0(min_variant_read_inference_B, " \u2264 variant reads, ", min_total_read_inference_B, " \u2264 total reads \u2264 ", max_total_read, "; ", Nmut_inference_B, " mutations (", format(round(freq_inference_B, 3), nsmall = 3), "%)")), "\n")
    report <- paste0(report, blue("Validation       : "), cyan(paste0(min_variant_read_validation, " \u2264 variant reads, ", min_total_read_validation, " \u2264 total reads \u2264 ", max_total_read, "; ", Nmut_validation, " mutations (", format(round(freq_validation, 3), nsmall = 3), "%)")), "\n")
    cat(report)


    #---Prepare the real SFS data
    mutation_table_inference_A <- mutation_table[Alt_count >= min_variant_read_inference_A & Tot_count >= min_total_read_inference_A & Tot_count <= max_total_read]
    mutation_table_inference_B <- mutation_table[Alt_count >= min_variant_read_inference_B & Tot_count >= min_total_read_inference_B & Tot_count <= max_total_read]
    mutation_table_validation <- mutation_table[Alt_count >= min_variant_read_validation & Tot_count >= min_total_read_validation & Tot_count <= max_total_read]



    mutation_table_validation <- mutation_table_validation[sample(.N, target_validation_Nmut, replace = TRUE)]



    func_SFS <- function(mutation_table) {
        SFS_data <- rep(0, length(SFS_data_frequencies))
        for (j in 1:length(SFS_data_frequencies)) {
            SFS_data[j] <- length(which(
                mutation_table$VAF >= ifelse(j == 1, 0, SFS_data_frequencies[j - 1]) &
                    mutation_table$VAF < SFS_data_frequencies[j]
            ))
        }
        return(SFS_data)
    }
    SFS_data_inference_A <- func_SFS(mutation_table_inference_A)
    SFS_data_inference_B <- func_SFS(mutation_table_inference_B)
    SFS_data_validation <- func_SFS(mutation_table_validation)
    #---Get coverage distribution
    func_coverage <- function(mutation_table) {
        sample_coverage <- data.frame(
            total_readcount = sort(unique(mutation_table$Tot_count)),
            pdf = 0
        )
        for (i in 1:nrow(sample_coverage)) {
            sample_coverage$pdf[i] <- sum(mutation_table$Tot_count == sample_coverage$total_readcount[i])
        }
        sample_coverage$pdf <- sample_coverage$pdf / sum(sample_coverage$pdf)
        return(sample_coverage)
    }
    sample_coverage_inference_A <- func_coverage(mutation_table_inference_A)
    sample_coverage_inference_B <- func_coverage(mutation_table_inference_B)
    sample_coverage_validation <- func_coverage(mutation_table_validation)
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

cluster_count_criteria_experiment_2 <- function(num_parameters,
                                                SFS_data_validation,
                                                parameters_validation,
                                                component_distributions_validation,
                                                with_tail,
                                                zero_cutoff) {
    num_samples_validation <- sum(SFS_data_validation)
    #   Compute the log-likelihood
    if (with_tail) {
        A <- parameters_validation[1]
        vec_K <- parameters_validation[seq(3, length(parameters_validation), by = 2)]
    } else {
        A <- NA
        vec_K <- parameters_validation[seq(1, length(parameters_validation), by = 2)]
    }
    log_L <- compute_loglikelihood(
        A = A,
        vec_K = vec_K,
        component_distributions = component_distributions_validation,
        vec_SFS_real = SFS_data_validation,
        zero_cutoff = zero_cutoff
    )
    #   Compute the stopping criteria
    compute_AIC <- function() {
        AIC <- 2 * num_parameters - 2 * log_L
        return(AIC)
    }
    compute_BIC <- function() {
        BIC <- num_parameters * log(num_samples_validation) - 2 * log_L
        return(BIC)
    }
    compute_ICL <- function() {
        #   Compute the latent variable distributions
        latent_variable_distributions <- component_distributions_validation$SFS_expected_normalized
        for (row in 1:nrow(latent_variable_distributions)) {
            latent_variable_distributions[row, ] <- parameters_validation[2 * row - 1] * latent_variable_distributions[row, ]
        }
        for (col in 1:ncol(latent_variable_distributions)) {
            latent_variable_distributions[, col] <- latent_variable_distributions[, col] / sum(latent_variable_distributions[, col])
        }
        latent_variable_distributions[which(latent_variable_distributions <= zero_cutoff | is.na(latent_variable_distributions))] <- zero_cutoff
        #   Compute the entropy
        entropy <- sum(SFS_data_validation * colSums(latent_variable_distributions * log(latent_variable_distributions)))
        #   Compute the Bayesian Information Criterion
        BIC <- log_L - 0.5 * num_parameters * log(num_samples_validation)
        #   Compute the Integrated Completed Log-Likelihood
        ICL <- -BIC - entropy
        return(ICL)
    }
    compute_ICL_MAP <- function() {
        #   Compute the latent variable distributions
        latent_variable_distributions <- component_distributions_validation$SFS_expected_normalized
        for (row in 1:nrow(latent_variable_distributions)) {
            latent_variable_distributions[row, ] <- parameters_validation[2 * row - 1] * latent_variable_distributions[row, ]
        }
        for (col in 1:ncol(latent_variable_distributions)) {
            latent_variable_distributions[, col] <- latent_variable_distributions[, col] / sum(latent_variable_distributions[, col])
        }
        latent_variable_distributions[which(latent_variable_distributions <= zero_cutoff | is.na(latent_variable_distributions))] <- zero_cutoff
        #   Compute the MAP allocations for mutations to clusters
        indicator_latent_variable_distributions <- matrix(0, nrow = nrow(latent_variable_distributions), ncol = ncol(latent_variable_distributions))
        for (col in 1:ncol(latent_variable_distributions)) {
            max_p <- which(latent_variable_distributions[, col] == max(latent_variable_distributions[, col]))[1]
            indicator_latent_variable_distributions[max_p, col] <- 1
        }
        #   Compute the entropy
        entropy_MAP <- sum(SFS_data_validation * colSums(indicator_latent_variable_distributions * log(latent_variable_distributions)))
        #   Compute the Bayesian Information Criterion
        BIC <- log_L - 0.5 * num_parameters * log(num_samples_validation)
        #   Compute the Integrated Completed Log-Likelihood
        ICL_MAP <- -BIC - entropy_MAP
        return(ICL_MAP)
    }
    criteria <- data.frame(
        AIC = compute_AIC(),
        BIC = compute_BIC(),
        ICL = compute_ICL(),
        ICL_MAP = compute_ICL_MAP()
    )
    return(criteria)
}

parameter_conversion_experiment_2 <- function(result,
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

compute_loglikelihood <- function(A, vec_K, component_distributions, vec_SFS_real, zero_cutoff) {
    #----------------Compute the SFS probability distribution from model
    vec_SFS_model <- compute_SFS(
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

compute_SFS <- function(A, vec_K, component_distributions) {
    # 	Add the neutral component
    vec_SFS_model <- component_distributions$SFS_expected_normalized[1, ]
    if (!is.na(A)) vec_SFS_model <- A * vec_SFS_model
    # 	Add the binomial humps
    for (i_hump in seq_along(vec_K)) {
        vec_SFS_model <- vec_SFS_model + vec_K[i_hump] * component_distributions$SFS_expected_normalized[i_hump + 1, ]
    }
    # 	Return the full SFS
    return(vec_SFS_model)
}

build_SFS_library <- function(neutral_power, cluster_frequencies, sfs_bincount, SFS_convolution_matrix, N_end) {
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
        vec_SFS_GT <- build_SFS_library_Griffiths_Tavare(vec_para = vec_para, N_end = N_end)
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
            vec_SFS_GT <- build_SFS_library_Griffiths_Tavare(vec_para = vec_para, N_end = N_end)
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

build_SFS_library_Griffiths_Tavare <- function(vec_para, N_end) {
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

build_convolution_matrix_experiment_2 <- function(sfs_bincount,
                                                  mode = NULL,
                                                  sample_size,
                                                  min_variant_read,
                                                  min_total_read,
                                                  max_total_read,
                                                  matrix_binomial_ploidy,
                                                  coverage_distribution,
                                                  coverage_variables,
                                                  sample_coverage,
                                                  compute_parallel = TRUE,
                                                  n_cores = NULL) {
    suppressWarnings(library(progress))
    report <- "Prepare the SFS convolution matrix"
    if (!is.null(mode)) report <- paste0(report, " for ", mode)
    report <- paste0(report, "...\n")
    cat(bold(blue(report)))
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
        vec_convolution <- rep(0, N_end)
        for (m in 1:N_end) {
            B_values <- sapply(sample_coverage_distribution$total_readcount, function(r) {
                r1 <- r * (i - 1) / SFS_totalsteps_base
                r1 <- ifelse(r1 %% 1 == 0, r1 + 1, ceiling(r1))
                r1 <- max(min_variant_read, r1, 1)
                r2 <- floor(r * i / SFS_totalsteps_base)
                if (r1 > r2) {
                    B <- 0
                } else {
                    B <- pbinom(r2, size = r, prob = m / (N_end * matrix_binomial_ploidy)) -
                        pbinom(r1 - 1, size = r, prob = m / (N_end * matrix_binomial_ploidy))
                }
                return(B)
            })
            vec_convolution[m] <- sum(sample_coverage_distribution$pdf * B_values)
        }
        return(vec_convolution)
    }
    #---Build convolution matrix to transform Griffiths-Tavare SFS to expected SFS
    vec_SFS_freq <- seq(0, 1, length.out = sfs_bincount + 1)
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
        suppressWarnings(library(parallel))
        suppressWarnings(library(pbapply))
        #   Start parallel cluster
        numCores <- ifelse(is.null(n_cores), detectCores(), n_cores)
        cl <- makePSOCKcluster(numCores - 1)
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

DECODE_plot_model_selection_experiment_2 <- function(DECODE_result,
                                                     fit = "best",
                                                     mode = "inference_A",
                                                     data_marker_colors = NULL) {
    suppressWarnings(library(gridExtra))
    # ####################################################################
    # ####################################################################
    # ####################################################################
    # ####################################################################
    # ####################################################################
    # suppressWarnings(library(crayon))
    # cat(paste0("\n\n\n", bold(red(paste0(DECODE_result$sample_id, ":"))), "\n"))
    # for (fit_current in names(DECODE_result$fits_with_tail$all_fits)[2:length(names(DECODE_result$fits_with_tail$all_fits))]) {
    #     fit_previous <- paste0(as.numeric(substr(fit_current, 1, 1)) - 1, substr(fit_current, 2, nchar(fit_current)))
    #     cat(paste0(magenta(paste0(fit_current, "==================================")), "\n"))



    #     # zero_cutoff <- 1e-50
    #     # vec_SFS_real_A <- DECODE_result$SFS_data_inference_A
    #     # vec_SFS_real_B <- DECODE_result$SFS_data_inference_B
    #     # vec_SFS_real_validation <- DECODE_result$SFS_data_validation
    #     # parameters_inference_A_current <- DECODE_result$fits_with_tail$all_fits[[fit_current]]$best$parameters_inference_A
    #     # parameters_inference_B_current <- DECODE_result$fits_with_tail$all_fits[[fit_current]]$best$parameters_inference_B
    #     # parameters_validation_current <- DECODE_result$fits_with_tail$all_fits[[fit_current]]$best$parameters_validation
    #     # component_distributions_inference_A_current <- DECODE_result$fits_with_tail$all_fits[[fit_current]]$best$component_distributions_inference_A
    #     # component_distributions_inference_B_current <- DECODE_result$fits_with_tail$all_fits[[fit_current]]$best$component_distributions_inference_B
    #     # component_distributions_validation_current <- DECODE_result$fits_with_tail$all_fits[[fit_current]]$best$component_distributions_validation
    #     # parameters_inference_A_previous <- DECODE_result$fits_with_tail$all_fits[[fit_previous]]$best$parameters_inference_A
    #     # parameters_inference_B_previous <- DECODE_result$fits_with_tail$all_fits[[fit_previous]]$best$parameters_inference_B
    #     # parameters_validation_previous <- DECODE_result$fits_with_tail$all_fits[[fit_previous]]$best$parameters_validation
    #     # component_distributions_inference_A_previous <- DECODE_result$fits_with_tail$all_fits[[fit_previous]]$best$component_distributions_inference_A
    #     # component_distributions_inference_B_previous <- DECODE_result$fits_with_tail$all_fits[[fit_previous]]$best$component_distributions_inference_B
    #     # component_distributions_validation_previous <- DECODE_result$fits_with_tail$all_fits[[fit_previous]]$best$component_distributions_validation



    #     # logLikelihood_A_previous <- compute_loglikelihood(
    #     #     A = parameters_inference_A_previous[1],
    #     #     vec_K = parameters_inference_A_previous[seq(3, length(parameters_inference_A_previous), by = 2)],
    #     #     component_distributions = component_distributions_inference_A_previous,
    #     #     vec_SFS_real = vec_SFS_real_A,
    #     #     zero_cutoff = zero_cutoff
    #     # )
    #     # logLikelihood_B_previous <- compute_loglikelihood(
    #     #     A = parameters_inference_B_previous[1],
    #     #     vec_K = parameters_inference_B_previous[seq(3, length(parameters_inference_B_previous), by = 2)],
    #     #     component_distributions = component_distributions_inference_B_previous,
    #     #     vec_SFS_real = vec_SFS_real_B,
    #     #     zero_cutoff = zero_cutoff
    #     # )
    #     # logLikelihood_validation_previous <- compute_loglikelihood(
    #     #     A = parameters_validation_previous[1],
    #     #     vec_K = parameters_validation_previous[seq(3, length(parameters_validation_previous), by = 2)],
    #     #     component_distributions = component_distributions_validation_previous,
    #     #     vec_SFS_real = vec_SFS_real_validation,
    #     #     zero_cutoff = zero_cutoff
    #     # )
    #     cat(paste0(blue("Previous BIC:      "), green(DECODE_result$fits_with_tail$all_fits[[fit_previous]]$best$criteria$BIC), "\n"))
    #     # cat(paste0(blue("Previous log-likelihood:   "), green(logLikelihood_validation_previous), "\n"))
    #     # logLikelihood_A_current <- compute_loglikelihood(
    #     #     A = parameters_inference_A_current[1],
    #     #     vec_K = parameters_inference_A_current[seq(3, length(parameters_inference_A_current), by = 2)],
    #     #     component_distributions = component_distributions_inference_A_current,
    #     #     vec_SFS_real = vec_SFS_real_A,
    #     #     zero_cutoff = zero_cutoff
    #     # )
    #     # logLikelihood_B_current <- compute_loglikelihood(
    #     #     A = parameters_inference_B_current[1],
    #     #     vec_K = parameters_inference_B_current[seq(3, length(parameters_inference_B_current), by = 2)],
    #     #     component_distributions = component_distributions_inference_B_current,
    #     #     vec_SFS_real = vec_SFS_real_B,
    #     #     zero_cutoff = zero_cutoff
    #     # )
    #     # logLikelihood_validation_current <- compute_loglikelihood(
    #     #     A = parameters_validation_current[1],
    #     #     vec_K = parameters_validation_current[seq(3, length(parameters_validation_current), by = 2)],
    #     #     component_distributions = component_distributions_validation_current,
    #     #     vec_SFS_real = vec_SFS_real_validation,
    #     #     zero_cutoff = zero_cutoff
    #     # )
    #     # cat(paste0(blue("Current log-likelihood:    "), green(logLikelihood_validation_current), "\n"))
    #     cat(paste0(blue("Current BIC:       "), green(DECODE_result$fits_with_tail$all_fits[[fit_current]]$best$criteria$BIC), "\n"))


    #     # Lambda <- 2 * (logLikelihood_validation_current - logLikelihood_validation_previous)
    #     # pval <- pchisq(Lambda, df = 2, lower.tail = FALSE)
    #     # cat(paste0(blue("Lambda:                    "), green(Lambda), "\n"))
    #     # cat(paste0(blue("\nP-value:                   "), green(pval), "\n\n"))
    #     cat(paste0(blue("BIC ration:        "), green(DECODE_result$fits_with_tail$all_fits[[fit_current]]$best$criteria$BIC / DECODE_result$fits_with_tail$all_fits[[fit_previous]]$best$criteria$BIC), "\n"))
    # }
    # ####################################################################
    # ####################################################################
    # ####################################################################
    # ####################################################################
    # ####################################################################
    #---Get information about the best fit
    if (DECODE_result$final_fit$best_fit$tail_status) {
        text_best_fit <- paste0("T+", length(DECODE_result$final_fit$best_fit$parameters_validation) / 2 - 1)
    } else {
        text_best_fit <- paste0("NT+", length(DECODE_result$final_fit$best_fit$parameters_validation) / 2)
    }
    #---Plot each SFS fit
    func_one_fit <- function(p_right, fit, text) {
        if (length(p_right) == 0) {
            text_ylab_inference_A <- "Inference A"
            text_ylab_inference_B <- "Inference B"
            text_ylab_validation <- "Validation"
        } else {
            text_ylab_inference_A <- NULL
            text_ylab_inference_B <- NULL
            text_ylab_validation <- NULL
            text_xlab_validation <- NULL
        }
        p_right_inference_A <- DECODE_plot_SFS_experiment_2(
            DECODE_result = DECODE_result,
            fit = fit,
            mode = "inference_A",
            text_xlab = NULL,
            text_ylab = text_ylab_inference_A,
            text_legend = NULL,
            notation = FALSE,
            data_marker_colors = data_marker_colors
        )
        p_right_inference_B <- DECODE_plot_SFS_experiment_2(
            DECODE_result = DECODE_result,
            fit = fit,
            mode = "inference_B",
            text_xlab = NULL,
            text_ylab = text_ylab_inference_B,
            text_legend = NULL,
            notation = FALSE,
            data_marker_colors = data_marker_colors
        )
        p_right_validation <- DECODE_plot_SFS_experiment_2(
            DECODE_result = DECODE_result,
            fit = fit,
            mode = "validation",
            text_xlab = text,
            text_ylab = text_ylab_validation,
            text_legend = NULL,
            notation = FALSE,
            data_marker_colors = data_marker_colors
        )
        p_right[[length(p_right) + 1]] <-
            arrangeGrob(
                p_right_inference_A,
                p_right_inference_B,
                p_right_validation,
                ncol = 1
            )
        return(p_right)
    }
    p_right <- list()
    criteria <- data.frame()
    if ("fits_with_tail" %in% names(DECODE_result)) {
        for (fit in names(DECODE_result$fits_with_tail$all_fits)) {
            fit_ID <- paste0("T+", substr(fit, 1, 1))
            criteria_new <- DECODE_result$fits_with_tail$all_fits[[fit]]$best$criteria
            criteria_new$fit <- fit_ID
            criteria <- rbind(criteria, criteria_new)
            if (fit_ID == text_best_fit) fit_ID <- paste0(fit_ID, " [best]")
            p_right <- func_one_fit(
                p_right = p_right,
                fit = paste0("fits_with_tail:", fit),
                text = fit_ID
            )
        }
    }
    if ("fits_without_tail" %in% names(DECODE_result)) {
        for (fit in names(DECODE_result$fits_without_tail$all_fits)) {
            fit_ID <- paste0("WT+", substr(fit, 1, 1))
            criteria_new <- DECODE_result$fits_with_tail$all_fits[[fit]]$best$criteria
            criteria_new$fit <- fit_ID
            criteria <- rbind(criteria, criteria_new)
            if (fit_ID == text_best_fit) fit_ID <- paste0(fit_ID, " [best]")
            p_right <- func_one_fit(
                p_right = p_right,
                fit = paste0("fits_without_tail:", fit),
                text = fit_ID
            )
        }
    }
    p_right <- grid.arrange(grobs = p_right, nrow = 1)
    #---Plot mutation threshold selection & fit selection
    p_left_threshold_selection <- DECODE_plot_readcounts_experiment_2(
        DECODE_result = DECODE_result
    )
    p_left_criteria <- DECODE_plot_criteria_experiment_2(
        criteria = criteria,
        criterion = DECODE_result$criterion,
        criterion_ratio = DECODE_result$criterion_ratio
    )
    p_left <- arrangeGrob(
        p_left_threshold_selection,
        p_left_criteria,
        ncol = 1,
        heights = c(2, 1)
    )
    #---Combine the plots
    p <- grid.arrange(grobs = list(p_left, p_right), ncol = 2, widths = c(1, 2.5))
    return(p)
}

DECODE_plot_criteria_experiment_2 <- function(criteria,
                                              criterion,
                                              criterion_ratio) {
    #---Get ratio of each fit compared to previous fit
    criteria[[paste0(criterion, "_target")]] <- NA
    for (row in 1:nrow(criteria)) {
        fit_current <- criteria[row, "fit"]
        fit_previous <- paste0("T+", as.numeric(sub(".*\\+(\\d+).*", "\\1", fit_current)) - 1)
        if (fit_previous %in% criteria$fit) {
            criteria[row, paste0(criterion, "_target")] <- criteria[criteria$fit == fit_previous, criterion] * criterion_ratio
        }
    }
    #---Plot the criterion selection
    p <- ggplot(criteria) +
        geom_point(aes(x = fit, y = !!sym(criterion), shape = "Actual"), color = "#DF536B", size = 10) +
        geom_point(aes(x = fit, y = !!sym(paste0(criterion, "_target")), shape = "Target"), color = "#2297E6", size = 6, stroke = 2, na.rm = TRUE) +
        labs(x = NULL, y = NULL, shape = NULL) +
        scale_shape_manual(values = c("Actual" = 16, "Target" = 6), labels = c(criterion, paste0(criterion, " acceptance threshold"))) +
        theme(
            text = element_text(size = 30),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5),
            legend.key.width = unit(1.5, "cm")
        )
    return(p)
}

DECODE_plot_SFS_experiment_2 <- function(DECODE_result,
                                         fit = "best",
                                         mode = "inference_A",
                                         text_xlab = "Variant Allele Frequency",
                                         text_ylab = "Mutation count",
                                         text_legend = NULL,
                                         notation = TRUE,
                                         data_marker_colors = NULL) {
    suppressWarnings(library(ggplot2))
    if (is.null(data_marker_colors)) data_marker_colors <- c("Data" = "black")
    vec_freq <- DECODE_result$SFS_frequencies
    SFS_totalsteps <- length(vec_freq)
    mutation_table <- DECODE_result$mutational_table
    max_total_read <- DECODE_result$max_total_read
    if (mode == "inference_A") {
        vec_SFS_real <- DECODE_result$SFS_data_inference_A
        min_variant_read <- DECODE_result$min_variant_read_inference_A
        min_total_read <- DECODE_result$min_total_read_inference_A
        if (fit == "best") {
            parameters <- DECODE_result$final_fit$best_fit$parameters_inference_A
            tail_status <- DECODE_result$final_fit$best_fit$tail_status
            component_distributions <- DECODE_result$final_fit$best_fit$component_distributions_inference_A
        } else {
            detail_tail <- sub(":.*", "", fit)
            detail_Ncluster <- sub(".*:", "", fit)
            parameters <- DECODE_result[[detail_tail]]$all_fits[[detail_Ncluster]]$best$parameters_inference_A
            tail_status <- ifelse(detail_tail == "fits_with_tail", TRUE, FALSE)
            component_distributions <- DECODE_result[[detail_tail]]$all_fits[[detail_Ncluster]]$best$component_distributions_inference_A
        }
    } else if (mode == "inference_B") {
        vec_SFS_real <- DECODE_result$SFS_data_inference_B
        min_variant_read <- DECODE_result$min_variant_read_inference_B
        min_total_read <- DECODE_result$min_total_read_inference_B
        if (fit == "best") {
            parameters <- DECODE_result$final_fit$best_fit$parameters_inference_B
            tail_status <- DECODE_result$final_fit$best_fit$tail_status
            component_distributions <- DECODE_result$final_fit$best_fit$component_distributions_inference_B
        } else {
            detail_tail <- sub(":.*", "", fit)
            detail_Ncluster <- sub(".*:", "", fit)
            parameters <- DECODE_result[[detail_tail]]$all_fits[[detail_Ncluster]]$best$parameters_inference_B
            tail_status <- ifelse(detail_tail == "fits_with_tail", TRUE, FALSE)
            component_distributions <- DECODE_result[[detail_tail]]$all_fits[[detail_Ncluster]]$best$component_distributions_inference_B
        }
    } else if (mode == "validation") {
        vec_SFS_real <- DECODE_result$SFS_data_validation
        min_variant_read <- DECODE_result$min_variant_read_validation
        min_total_read <- DECODE_result$min_total_read_validation
        if (fit == "best") {
            parameters <- DECODE_result$final_fit$best_fit$parameters_validation
            tail_status <- DECODE_result$final_fit$best_fit$tail_status
            component_distributions <- DECODE_result$final_fit$best_fit$component_distributions_validation
        } else {
            detail_tail <- sub(":.*", "", fit)
            detail_Ncluster <- sub(".*:", "", fit)
            parameters <- DECODE_result[[detail_tail]]$all_fits[[detail_Ncluster]]$best$parameters_validation
            tail_status <- ifelse(detail_tail == "fits_with_tail", TRUE, FALSE)
            component_distributions <- DECODE_result[[detail_tail]]$all_fits[[detail_Ncluster]]$best$component_distributions_validation
        }
    }
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
    if ("Marker" %in% colnames(mutation_table)) {
        mutation_markers <- mutation_table$Marker
    } else {
        mutation_markers <- c()
    }
    #---Plot the SFS deconvolution
    color_scheme <- c(
        data_marker_colors,
        "Neutral tail" = "#999999",
        "Cluster 1" = "#D55E00",
        "Cluster 2" = "#0072B2",
        "Cluster 3" = "#009E73",
        "Cluster 4" = "#E69F00",
        "Cluster 5" = "#56B4E9",
        "Cluster 6" = "#CC79A7",
        "Cluster 7" = "#F0E442"
    )
    mutation_count <- sum(vec_SFS_real)
    #   Prepare the data for plotting
    if (is.null(mutation_markers)) {
        df_data <- data.frame(frequency = vec_freq, count = vec_SFS_real, fill = "Data")
    } else {
        unique_markers <- unique(mutation_markers)
        df_data <- data.frame(
            ID = rep(1:length(vec_freq), length(unique_markers)),
            frequency = rep(vec_freq, length(unique_markers)),
            count = rep(0, SFS_totalsteps * length(unique_markers)),
            fill = rep(unique_markers, each = SFS_totalsteps)
        )
        for (j in 1:length(mutation_table$Alt_count)) {
            no_variant <- mutation_table$Alt_count[j]
            no_total <- mutation_table$Ref_count[j] + mutation_table$Alt_count[j]
            if (no_variant >= min_variant_read & no_total >= min_total_read & no_total <= max_total_read) {
                VAF <- no_variant / no_total
                pos <- which(vec_freq >= VAF)[1]
                df_data$count[which(df_data$ID == pos & df_data$fill == mutation_markers[j])] <- df_data$count[which(df_data$ID == pos & df_data$fill == mutation_markers[j])] + 1
            }
        }
        df_data$fill <- gsub("_", " ", df_data$fill)
    }
    #   Prepare the deconvolution inference for plotting
    df_fit <- data.frame()
    df_fit_order <- c()
    if (tail_status) {
        SFS_neutral <- compute_SFS(
            A = vec_A[1],
            vec_K = c(),
            component_distributions = component_distributions
        )
        SFS_neutral <- vec_A[1] * mutation_count * SFS_neutral / sum(SFS_neutral)
        df_fit <- rbind(df_fit, data.frame(frequency = vec_freq, count = SFS_neutral, fill = "Neutral tail"))
        df_fit_order <- c(df_fit_order, "Neutral tail")
    }
    if (N_humps > 0) {
        for (i in 1:N_humps) {
            vec_K_tmp <- rep(0, N_humps)
            vec_K_tmp[i] <- vec_K[i]
            SFS_hump <- compute_SFS(
                A = 0,
                vec_K = vec_K_tmp,
                component_distributions = component_distributions
            )
            SFS_hump <- vec_K[i] * mutation_count * SFS_hump / sum(SFS_hump)
            df_fit <- rbind(df_fit, data.frame(frequency = vec_freq, count = SFS_hump, fill = paste0("Cluster ", i)))
            df_fit_order <- c(df_fit_order, paste0("Cluster ", i))
        }
    }
    df_fit$fill <- factor(df_fit$fill, levels = rev(df_fit_order))
    df_data$fill <- factor(df_data$fill, levels = rev(names(data_marker_colors)))
    p <- ggplot() +
        geom_area(data = df_fit, aes(x = frequency, y = count, fill = fill), position = "stack", alpha = 0.8) +
        geom_bar(data = df_data, aes(x = frequency, y = count, fill = fill), stat = "identity", width = 0.5 / SFS_totalsteps) +
        scale_fill_manual(values = color_scheme, name = "") +
        guides(fill = guide_legend(nrow = 1, keywidth = 2, keyheight = 1)) +
        xlab(text_xlab) +
        ylab(text_ylab) +
        labs(title = text_legend) +
        theme(
            text = element_text(size = 30),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5)
        )
    if (!notation) {
        p <- p +
            theme(
                legend.position = "none",
                axis.ticks.x = element_blank(),
                axis.text.x = element_blank()
            )
    }
    if (is.null(text_ylab)) {
        p <- p +
            theme(
                axis.ticks.y = element_blank(),
                axis.text.y = element_blank()
            )
    }
    return(p)
}

DECODE_plot_readcounts_experiment_2 <- function(DECODE_result,
                                                freq_cutoff = 10) {
    suppressWarnings(library(ggplot2))
    suppressWarnings(library(shadowtext))
    suppressWarnings(library(reshape2))
    mutation_table <- DECODE_result$mutational_table
    vec_min_variant_read <- min(mutation_table$Alt_count):max(mutation_table$Alt_count)
    vec_min_total_read <- min(mutation_table$Tot_count):max(mutation_table$Tot_count)
    readcount_distribution <- DECODE_result$readcount_distribution
    min_variant_read_inference_A <- max(DECODE_result$min_variant_read_inference_A, min(mutation_table$Alt_count))
    min_total_read_inference_A <- max(DECODE_result$min_total_read_inference_A, min(mutation_table$Tot_count))
    min_variant_read_inference_B <- max(DECODE_result$min_variant_read_inference_B, min(mutation_table$Alt_count))
    min_total_read_inference_B <- max(DECODE_result$min_total_read_inference_B, min(mutation_table$Tot_count))
    min_variant_read_validation <- max(DECODE_result$min_variant_read_validation, min(mutation_table$Alt_count))
    min_total_read_validation <- max(DECODE_result$min_total_read_validation, min(mutation_table$Tot_count))
    #---Reduce the distribution to region satisfying the frequency cutoff
    max_total_read <- vec_min_total_read[which(readcount_distribution$freq[which(readcount_distribution$min_variant_read == vec_min_variant_read[1])] < freq_cutoff)[1]]
    max_variant_read <- vec_min_variant_read[which(readcount_distribution$freq[which(readcount_distribution$min_total_read == vec_min_total_read[1])] < freq_cutoff)[1]]
    readcount_distribution <- readcount_distribution[which(readcount_distribution$min_total_read <= max_total_read & readcount_distribution$min_variant_read <= max_variant_read), ]
    #---Plot the readcount distribution
    freq_inference_A <- round(readcount_distribution$freq[readcount_distribution$min_total_read == min_total_read_inference_A & readcount_distribution$min_variant_read == min_variant_read_inference_A], 2)
    freq_inference_B <- round(readcount_distribution$freq[readcount_distribution$min_total_read == min_total_read_inference_B & readcount_distribution$min_variant_read == min_variant_read_inference_B], 2)
    freq_validation <- round(readcount_distribution$freq[readcount_distribution$min_total_read == min_total_read_validation & readcount_distribution$min_variant_read == min_variant_read_validation], 2)
    text_inference_A <- paste0("Inference A (", freq_inference_A, "%)")
    text_inference_B <- paste0("Inference B (", freq_inference_B, "%)")
    text_validation <- paste0("Validation (", freq_validation, "%)")
    color_inference_A <- "#DF536B"
    color_inference_B <- "#DF536B"
    color_validation <- "#2297E6"
    p <- ggplot(readcount_distribution, aes(x = min_total_read, y = min_variant_read, fill = freq)) +
        geom_tile() +
        geom_rect(
            aes(
                xmin = min_total_read_inference_A - 0.5, xmax = min_total_read_inference_A + 0.5,
                ymin = min_variant_read_inference_A - 0.5, ymax = min_variant_read_inference_A + 0.5
            ),
            fill = NA, color = "white", linewidth = 1
        ) +
        geom_shadowtext(
            aes(
                x = min_total_read_inference_A + 2,
                y = min_variant_read_inference_A,
                label = text_inference_A
            ),
            angle = 45,
            hjust = 0,
            vjust = 0,
            size = 6,
            color = color_inference_A,
            bg.color = "white",
            fontface = "bold"
        ) +
        geom_rect(
            aes(
                xmin = min_total_read_inference_B - 0.5, xmax = min_total_read_inference_B + 0.5,
                ymin = min_variant_read_inference_B - 0.5, ymax = min_variant_read_inference_B + 0.5
            ),
            fill = NA, color = "white", linewidth = 1
        ) +
        geom_shadowtext(
            aes(
                x = min_total_read_inference_B + 2,
                y = min_variant_read_inference_B,
                label = text_inference_B
            ),
            angle = 45,
            hjust = 0,
            vjust = 0,
            size = 6,
            color = color_inference_B,
            bg.color = "white",
            fontface = "bold"
        ) +
        geom_rect(
            aes(
                xmin = min_total_read_validation - 0.5, xmax = min_total_read_validation + 0.5,
                ymin = min_variant_read_validation - 0.5, ymax = min_variant_read_validation + 0.5
            ),
            fill = NA, color = "white", linewidth = 1
        ) +
        geom_shadowtext(
            aes(
                x = min_total_read_validation + 2,
                y = min_variant_read_validation,
                label = text_validation
            ),
            angle = 45,
            hjust = 0,
            vjust = 0,
            size = 6,
            color = color_validation,
            bg.color = "white",
            fontface = "bold"
        ) +
        scale_fill_gradientn(
            colors = c("#0072B2", "#56B4E9", "#009E73", "#E69F00", "#D55E00"),
            name = "% mutations retained"
        ) +
        theme_minimal() +
        labs(title = "", x = "Minimum total read count", y = "Minimum variant read count") +
        theme(
            text = element_text(size = 30),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5),
            legend.key.width = unit(1.5, "cm"),
        )
    return(p)
}
