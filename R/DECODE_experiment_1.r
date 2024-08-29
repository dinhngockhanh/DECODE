DECODE_experiment_1 <- function(sample_id = "",
                                mutation_table,
                                criterion = "BIC",
                                criterion_ratio = 1,
                                neutral_power_min = 0.5,
                                neutral_power_max = 5,
                                cluster_frequency_min = 0.01,
                                cluster_frequency_max = 1,
                                max_total_read = NULL,
                                sample_size,
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
    library(crayon)
    cat(paste0("\n\n\n", bold(red("PERFORMING DECODE FOR SAMPLE ")), bold(yellow(sample_id)), bold(red("...")), "\n"))
    mutation_table$Tot_count <- mutation_table$Ref_count + mutation_table$Alt_count
    mutation_table$VAF <- mutation_table$Alt_count / mutation_table$Tot_count
    if (is.null(max_total_read)) max_total_read <- max(mutation_table$Tot_count)
    ####################################################################
    ####################################################################
    ####################################################################
    libPaths_binomial_table <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/DECODE_binomial_matrices"
    matrix_binomial_sample_size <- 1000
    matrix_binomial_sfs_bincount <- sfs_bincount
    matrix_binomial_ploidy <- 2
    thresholds <- choose_mutation_thresholds(mutation_table = mutation_table)
    readcount_distribution <- thresholds$readcount_distribution
    min_variant_read_inference_A <- thresholds$min_variant_read_inference_A
    min_total_read_inference_A <- thresholds$min_total_read_inference_A
    min_variant_read_inference_B <- thresholds$min_variant_read_inference_B
    min_total_read_inference_B <- thresholds$min_total_read_inference_B
    min_variant_read_validation <- thresholds$min_variant_read_validation
    min_total_read_validation <- thresholds$min_total_read_validation
    ####################################################################
    ####################################################################
    ####################################################################
    #---Prepare the real SFS
    vec_freq <- seq(1, sfs_bincount) / sfs_bincount
    SFS_data_inference_A <- rep(0, sfs_bincount)
    for (j in 1:sfs_bincount) {
        SFS_data_inference_A[j] <- length(which(
            mutation_table$Alt_count >= min_variant_read_inference_A &
                mutation_table$Tot_count >= min_total_read_inference_A &
                mutation_table$VAF >= ifelse(j == 1, 0, vec_freq[j - 1]) &
                mutation_table$VAF < vec_freq[j]
        ))
    }
    SFS_data_inference_B <- rep(0, sfs_bincount)
    for (j in 1:sfs_bincount) {
        SFS_data_inference_B[j] <- length(which(
            mutation_table$Alt_count >= min_variant_read_inference_B &
                mutation_table$Tot_count >= min_total_read_inference_B &
                mutation_table$VAF >= ifelse(j == 1, 0, vec_freq[j - 1]) &
                mutation_table$VAF < vec_freq[j]
        ))
    }
    SFS_data_validation <- rep(0, sfs_bincount)
    for (j in 1:sfs_bincount) {
        SFS_data_validation[j] <- length(which(
            mutation_table$Alt_count >= min_variant_read_validation &
                mutation_table$Tot_count >= min_total_read_validation &
                mutation_table$VAF >= ifelse(j == 1, 0, vec_freq[j - 1]) &
                mutation_table$VAF < vec_freq[j]
        ))
    }
    #---Prepare the total readcount distribution
    sample_coverage_inference_A <- prep_distribution_patient(
        mutations_total_read =  mutation_table$Tot_count,
        min_total_read = min_total_read_inference_A,
        max_total_read = max_total_read
    )
    sample_coverage_inference_B <- prep_distribution_patient(
        mutations_total_read =  mutation_table$Tot_count,
        min_total_read = min_total_read_inference_B,
        max_total_read = max_total_read
    )
    sample_coverage_validation <- prep_distribution_patient(
        mutations_total_read =  mutation_table$Tot_count,
        min_total_read = min_total_read_validation,
        max_total_read = max_total_read
    )
    #---Get DECODE binomial table
    binomial_matrix_inference_A <- get_binomial_matrix(
        folder = libPaths_binomial_table,
        matrix_binomial_sample_size = matrix_binomial_sample_size,
        matrix_binomial_sfs_bincount = matrix_binomial_sfs_bincount,
        matrix_binomial_ploidy = matrix_binomial_ploidy,
        min_variant_read = min_variant_read_inference_A,
        min_total_read = min_total_read_inference_A,
        max_total_read = max_total_read,
        compute_parallel = compute_parallel,
        n_cores = n_cores
    )
    binomial_matrix_inference_B <- get_binomial_matrix(
        folder = libPaths_binomial_table,
        matrix_binomial_sample_size = matrix_binomial_sample_size,
        matrix_binomial_sfs_bincount = matrix_binomial_sfs_bincount,
        matrix_binomial_ploidy = matrix_binomial_ploidy,
        min_variant_read = min_variant_read_inference_B,
        min_total_read = min_total_read_inference_B,
        max_total_read = max_total_read,
        compute_parallel = compute_parallel,
        n_cores = n_cores
    )
    binomial_matrix_validation <- get_binomial_matrix(
        folder = libPaths_binomial_table,
        matrix_binomial_sample_size = matrix_binomial_sample_size,
        matrix_binomial_sfs_bincount = matrix_binomial_sfs_bincount,
        matrix_binomial_ploidy = matrix_binomial_ploidy,
        min_variant_read = min_variant_read_validation,
        min_total_read = min_total_read_validation,
        max_total_read = max_total_read,
        compute_parallel = compute_parallel,
        n_cores = n_cores
    )
    #---Prepare the SFS convolution matrix
    cat(bold(blue("Prepare the SFS convolution matrix...\n")))
    SFS_convolution_inference_A <- build_convolution_matrix(
        binomial_matrix = binomial_matrix_inference_A,
        sfs_bincount = sfs_bincount,
        coverage_distribution = coverage_distribution,
        coverage_variables = coverage_variables,
        sample_coverage = sample_coverage_inference_A
    )
    SFS_convolution_inference_B <- build_convolution_matrix(
        binomial_matrix = binomial_matrix_inference_B,
        sfs_bincount = sfs_bincount,
        coverage_distribution = coverage_distribution,
        coverage_variables = coverage_variables,
        sample_coverage = sample_coverage_inference_B
    )
    SFS_convolution_validation <- build_convolution_matrix(
        binomial_matrix = binomial_matrix_validation,
        sfs_bincount = sfs_bincount,
        coverage_distribution = coverage_distribution,
        coverage_variables = coverage_variables,
        sample_coverage = sample_coverage_validation
    )
    #---DECODE
    DECODE_result <- list()
    DECODE_result$sample_id <- sample_id
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
    DECODE_result$SFS_frequencies <- vec_freq
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
        final_result <- DECODE_given_tail_status_experiment_1(
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
        report <- paste0(report, yellow(paste0("pi = ", format(round(parameters_inference_A_final_result[1], 3), nsmall = 3), " [A], ", format(round(parameters_inference_B_final_result[1], 3), nsmall = 3), " [B]")), "\n")
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
            report <- paste0(report, yellow(paste0("pi = ", format(round(parameters_inference_A_final_result[2 * (i + ii) + 1], 3), nsmall = 3), " [A], ", format(round(parameters_inference_B_final_result[2 * (i + ii) + 1], 3), nsmall = 3), " [B]")), "\n")
        }
    }
    cat(report)
    #---Translation to parameters of cancer evolution in the sample
    tmp <- parameter_conversion_experiment_1(
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

DECODE_given_tail_status_experiment_1 <- function(SFS_data_inference_A,
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
        fit_results <- DECODE_given_tail_status_and_Ncluster_experiment_1(
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
            report <- paste0(report, cyan(paste0("pi = ", format(round(parameters_inference_A_best_current[1], 3), nsmall = 3), " [A], ", format(round(parameters_inference_B_best_current[1], 3), nsmall = 3), " [B]")), "\n")
            ii <- 0
        } else {
            N_humps <- length(parameters_inference_A_best_current) / 2
            ii <- -1
        }
        if (N_humps > 0) {
            for (i in 1:N_humps) {
                report <- paste0(report, blue(paste0("Cluster ", i, "        : ")), cyan(paste0("frequency = ", format(round(parameters_inference_A_best_current[2 * (i + ii) + 2], 3), nsmall = 3))), blue(", "))
                report <- paste0(report, cyan(paste0("pi = ", format(round(parameters_inference_A_best_current[2 * (i + ii) + 1], 3), nsmall = 3), " [A], ", format(round(parameters_inference_B_best_current[2 * (i + ii) + 1], 3), nsmall = 3), " [B]")), "\n")
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

DECODE_given_tail_status_and_Ncluster_experiment_1 <- function(SFS_data_inference_A,
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
        results <- DECODE_for_pis_experiment_1(
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
        ################################################################
        ################################################################
        ################################################################
        ################################################################
        ################################################################
        #   Compute stopping criteria
        mutation_count <- sum(SFS_data_validation)
        num_parameters <- ifelse(with_tail, 2 * N_humps + 1, 2 * N_humps - 1)
        if (neutral_power_max > neutral_power_min) num_parameters <- num_parameters + 1
        criteria <- cluster_count_criteria_experiment_1(
            num_parameters = num_parameters,
            num_samples = mutation_count,
            SFS_data = SFS_data_validation,
            parameters = results$parameters_validation,
            with_tail = with_tail,
            component_distributions = component_distributions_validation,
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
        output$criteria <- criteria
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
        all_criteria <- data.frame()
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
            all_criteria <- rbind(all_criteria, trial_result$criteria)
        }
        if (progress_bar) cat("\n")
    } else {
        library(parallel)
        library(pbapply)
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
            "DECODE_for_pis_experiment_1", "compute_loglikelihood", "compute_SFS", "cluster_count_criteria_experiment_1"
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
        all_criteria <- do.call(rbind, lapply(output, function(x) x$criteria))
    }
    #---Find the best fit
    best_index <- which.max(all_logLikelihood)
    fit_results <- list()
    fit_results$all <- list()
    fit_results$all$parameters_inference_A <- all_parameters_inference_A
    fit_results$all$parameters_inference_B <- all_parameters_inference_B
    fit_results$all$parameters_validation <- all_parameters_validation
    fit_results$all$logLikelihood <- all_logLikelihood
    fit_results$all$criteria <- all_criteria
    fit_results$best <- list()
    fit_results$best$parameters_inference_A <- all_parameters_inference_A[best_index, ]
    fit_results$best$parameters_inference_B <- all_parameters_inference_B[best_index, ]
    fit_results$best$parameters_validation <- all_parameters_validation[best_index, ]
    fit_results$best$logLikelihood <- all_logLikelihood[best_index]
    fit_results$best$criteria <- all_criteria[best_index, ]
    fit_results$best$component_distributions_inference_A <- all_component_distributions_inference_A[[best_index]]
    fit_results$best$component_distributions_inference_B <- all_component_distributions_inference_B[[best_index]]
    fit_results$best$component_distributions_validation <- all_component_distributions_validation[[best_index]]
    return(fit_results)
}

DECODE_for_pis_experiment_1 <- function(SFS_data_inference_A,
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

choose_mutation_thresholds <- function(mutation_table,
                                       compute_parallel = TRUE,
                                       n_cores = NULL) {
    #---Find joint distribution of variant and total readcounts
    cat(bold(blue("Choose mutation thresholds for inference and validation...\n")))
    vec_min_variant_read <- min(mutation_table$Alt_count):max(mutation_table$Alt_count)
    vec_min_total_read <- min(mutation_table$Tot_count):max(mutation_table$Tot_count)
    func_dist_given_min_variant_read <- function(min_variant_read, mutation_table) {
        df <- data.frame()
        for (min_total_read in vec_min_total_read) {
            df <- rbind(
                df,
                data.frame(
                    min_total_read = min_total_read,
                    min_variant_read = min_variant_read,
                    freq = 100 * sum(mutation_table$Alt_count >= min_variant_read & mutation_table$Tot_count >= min_total_read) / nrow(mutation_table)
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
        library(parallel)
        library(pbapply)
        #   Start parallel cluster
        numCores <- ifelse(is.null(n_cores), detectCores(), n_cores)
        cl <- makePSOCKcluster(numCores - 1)
        #   Prepare input parameters
        clusterExport(cl, varlist = c("mutation_table"))
        #   Compute each sub-dataframe
        output <- pblapply(cl = cl, X = vec_min_variant_read, FUN = function(min_variant_read) {
            return(func_dist_given_min_variant_read(min_variant_read, mutation_table))
        })
        stopCluster(cl)
        readcount_distribution <- do.call(rbind, output)
    }
    #---Choose mutation thresholds for inference and validation
    min_total_read_inference_A <- min(readcount_distribution$min_total_read)
    filtered_df <- readcount_distribution[
        readcount_distribution$min_total_read == min_total_read_inference_A &
            readcount_distribution$freq >= 50,
    ]
    min_variant_read_inference_A <- ifelse(nrow(filtered_df) > 0, max(filtered_df$min_variant_read), min(readcount_distribution$min_variant_read))



    min_variant_read_inference_B <- 4
    filtered_df <- readcount_distribution[
        readcount_distribution$min_variant_read == min_variant_read_inference_B &
            readcount_distribution$freq >= 50,
    ]
    min_total_read_inference_B <- ifelse(nrow(filtered_df) > 0, max(filtered_df$min_total_read), min(readcount_distribution$min_total_read))



    min_variant_read_validation <- 5
    filtered_df <- readcount_distribution[
        readcount_distribution$min_variant_read == min_variant_read_validation &
            readcount_distribution$freq >= 90,
    ]
    min_total_read_validation <- ifelse(nrow(filtered_df) > 0, max(filtered_df$min_total_read), min(readcount_distribution$min_total_read))



    #---Report the chosen thresholds
    thresholds <- list()
    thresholds$readcount_distribution <- readcount_distribution
    thresholds$min_variant_read_inference_A <- min_variant_read_inference_A
    thresholds$min_total_read_inference_A <- min_total_read_inference_A
    thresholds$min_variant_read_inference_B <- min_variant_read_inference_B
    thresholds$min_total_read_inference_B <- min_total_read_inference_B
    thresholds$min_variant_read_validation <- min_variant_read_validation
    thresholds$min_total_read_validation <- min_total_read_validation
    return(thresholds)
}

cluster_count_criteria_experiment_1 <- function(num_parameters, num_samples, SFS_data, parameters, with_tail, component_distributions, zero_cutoff) {
    #   Compute the log-likelihood
    if (with_tail) {
        A <- parameters[1]
        vec_K <- parameters[seq(3, length(parameters), by = 2)]
    } else {
        A <- NA
        vec_K <- parameters[seq(1, length(parameters), by = 2)]
    }
    log_L <-
        compute_loglikelihood(
            A = A,
            vec_K = vec_K,
            component_distributions = component_distributions,
            vec_SFS_real = SFS_data,
            zero_cutoff = zero_cutoff
        )
    #   Compute the stopping criteria
    compute_AIC <- function() {
        AIC <- 2 * num_parameters - 2 * log_L
        return(AIC)
    }
    compute_BIC <- function() {
        BIC <- log(num_samples) * num_parameters - 2 * log_L
        return(BIC)
    }
    compute_ICL <- function() {
        #   Compute the latent variable distributions
        latent_variable_distributions <- component_distributions$SFS_expected_normalized
        for (row in 1:nrow(latent_variable_distributions)) {
            latent_variable_distributions[row, ] <- parameters[2 * row - 1] * latent_variable_distributions[row, ]
        }
        for (col in 1:ncol(latent_variable_distributions)) {
            latent_variable_distributions[, col] <- latent_variable_distributions[, col] / sum(latent_variable_distributions[, col])
        }
        latent_variable_distributions[which(latent_variable_distributions <= zero_cutoff | is.na(latent_variable_distributions))] <- zero_cutoff
        #   Compute the entropy
        entropy <- sum(SFS_data * colSums(latent_variable_distributions * log(latent_variable_distributions)))
        #   Compute the Bayesian Information Criterion
        BIC <- log_L - 0.5 * num_parameters * log(num_samples)
        #   Compute the Integrated Completed Log-Likelihood
        ICL <- -BIC - entropy
        return(ICL)
    }
    compute_ICL_MAP <- function() {
        #   Compute the latent variable distributions
        latent_variable_distributions <- component_distributions$SFS_expected_normalized
        for (row in 1:nrow(latent_variable_distributions)) {
            latent_variable_distributions[row, ] <- parameters[2 * row - 1] * latent_variable_distributions[row, ]
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
        entropy_MAP <- sum(SFS_data * colSums(indicator_latent_variable_distributions * log(latent_variable_distributions)))
        #   Compute the Bayesian Information Criterion
        BIC <- log_L - 0.5 * num_parameters * log(num_samples)
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

parameter_conversion_experiment_1 <- function(result,
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

DECODE_plot_experiment_1 <- function(DECODE_result,
                                     fit = "best",
                                     mode = "inference_A",
                                     data_marker_colors = NULL) {
    library(gridExtra)
    p_left <- DECODE_plot_readcounts_experiment_1(
        DECODE_result = DECODE_result
    )
    p_right_inference_A <- DECODE_plot_SFS_experiment_1(
        DECODE_result = DECODE_result,
        mode = "inference_A",
        text_legend = "Inference A",
        data_marker_colors = data_marker_colors
    )
    p_right_inference_B <- DECODE_plot_SFS_experiment_1(
        DECODE_result = DECODE_result,
        mode = "inference_B",
        text_legend = "Inference B",
        data_marker_colors = data_marker_colors
    )
    p_right_validation <- DECODE_plot_SFS_experiment_1(
        DECODE_result = DECODE_result,
        mode = "validation",
        text_legend = "Validation",
        data_marker_colors = data_marker_colors
    )
    p <- grid.arrange(
        p_left,
        arrangeGrob(
            p_right_inference_A,
            p_right_inference_B,
            p_right_validation,
            ncol = 1
        ),
        ncol = 2
    )
    return(p)
}

DECODE_plot_SFS_experiment_1 <- function(DECODE_result,
                                         fit = "best",
                                         mode = "inference_A",
                                         text_legend = "",
                                         data_marker_colors = NULL) {
    library(ggplot2)
    if (is.null(data_marker_colors)) data_marker_colors <- c("Data" = "black")


    vec_freq <- DECODE_result$SFS_frequencies
    SFS_totalsteps <- length(vec_freq)
    mutation_table <- DECODE_result$mutational_table
    tmp <- parameter_conversion_experiment_1(
        result = DECODE_result$final_fit,
        output_parameters_df = FALSE
    )
    if (mode == "inference_A") {
        vec_SFS_real <- DECODE_result$SFS_data_inference_A
        min_variant_read <- DECODE_result$min_variant_read_inference_A
        min_total_read <- DECODE_result$min_total_read_inference_A
        if (fit == "best") {
            vec_para_best_final <- DECODE_result$final_fit$best_fit$parameters_inference_A
            tail_status_final <- DECODE_result$final_fit$best_fit$tail_status
            component_distributions_best_final <- DECODE_result$final_fit$best_fit$component_distributions_inference_A
        }
        vec_A <- tmp$inference_A$vec_A
        vec_K <- tmp$inference_A$vec_K
        vec_p <- tmp$inference_A$vec_p
        N_humps <- tmp$inference_A$N_humps
        tail_status <- tmp$inference_A$tail_status
    } else if (mode == "inference_B") {
        vec_SFS_real <- DECODE_result$SFS_data_inference_B
        min_variant_read <- DECODE_result$min_variant_read_inference_B
        min_total_read <- DECODE_result$min_total_read_inference_B
        if (fit == "best") {
            vec_para_best_final <- DECODE_result$final_fit$best_fit$parameters_inference_B
            tail_status_final <- DECODE_result$final_fit$best_fit$tail_status
            component_distributions_best_final <- DECODE_result$final_fit$best_fit$component_distributions_inference_B
        }
        vec_A <- tmp$inference_B$vec_A
        vec_K <- tmp$inference_B$vec_K
        vec_p <- tmp$inference_B$vec_p
        N_humps <- tmp$inference_B$N_humps
        tail_status <- tmp$inference_B$tail_status
    } else if (mode == "validation") {
        vec_SFS_real <- DECODE_result$SFS_data_validation
        min_variant_read <- DECODE_result$min_variant_read_validation
        min_total_read <- DECODE_result$min_total_read_validation
        if (fit == "best") {
            vec_para_best_final <- DECODE_result$final_fit$best_fit$parameters_validation
            tail_status_final <- DECODE_result$final_fit$best_fit$tail_status
            component_distributions_best_final <- DECODE_result$final_fit$best_fit$component_distributions_validation
        }
        vec_A <- tmp$validation$vec_A
        vec_K <- tmp$validation$vec_K
        vec_p <- tmp$validation$vec_p
        N_humps <- tmp$validation$N_humps
        tail_status <- tmp$validation$tail_status
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
            if (no_variant >= min_variant_read && no_total >= min_total_read) {
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
    if (tail_status_final) {
        SFS_neutral <- compute_SFS(
            A = vec_A[1],
            vec_K = c(),
            component_distributions = component_distributions_best_final
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
                component_distributions = component_distributions_best_final
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
        scale_fill_manual(values = color_scheme, name = text_legend) +
        guides(fill = guide_legend(nrow = 1, keywidth = 2, keyheight = 1)) +
        xlab("Variant Allele Frequency") +
        ylab("Mutation count") +
        theme(
            text = element_text(size = 30),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5)
        )
    return(p)
}

DECODE_plot_readcounts_experiment_1 <- function(DECODE_result,
                                                freq_cutoff = 10) {
    library(ggplot2)
    library(shadowtext)
    library(reshape2)
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
            fill = NA, color = "white", size = 1
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
            size = 8,
            color = color_inference_A,
            bg.color = "white",
            fontface = "bold"
        ) +
        geom_rect(
            aes(
                xmin = min_total_read_inference_B - 0.5, xmax = min_total_read_inference_B + 0.5,
                ymin = min_variant_read_inference_B - 0.5, ymax = min_variant_read_inference_B + 0.5
            ),
            fill = NA, color = "white", size = 1
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
            size = 8,
            color = color_inference_B,
            bg.color = "white",
            fontface = "bold"
        ) +
        geom_rect(
            aes(
                xmin = min_total_read_validation - 0.5, xmax = min_total_read_validation + 0.5,
                ymin = min_variant_read_validation - 0.5, ymax = min_variant_read_validation + 0.5
            ),
            fill = NA, color = "white", size = 1
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
            size = 8,
            color = color_validation,
            bg.color = "white",
            fontface = "bold"
        ) +
        scale_fill_gradientn(
            colors = c("#0072B2", "#56B4E9", "#009E73", "#E69F00", "#D55E00"),
            name = "% mutations retained"
        ) +
        theme_minimal() +
        labs(title = "", x = "Minimum total readcount", y = "Minimum alternative readcount") +
        theme(
            text = element_text(size = 30),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5),
            legend.key.width = unit(3.5, "cm"),
        )
    return(p)
}
