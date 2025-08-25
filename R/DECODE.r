DECODE <- function(sample_id = "",
                   mutation_table,
                   criterion = "BIC",
                   criterion_penalty_scale = 1.0, # <<<<<<<<<<<<<<<<<<<<
                   criterion_pvalue_threshold = 0.1, # <<<<<<<<<<<<<<<<<
                   neutral_power_min = 0.5,
                   neutral_power_max = 5,
                   cluster_frequency_min = 0.01,
                   cluster_frequency_max = 1,
                   cluster_frequency_mindiff = 0.01,
                   max_total_read = NULL,
                   allele_count = 1000,
                   matrix_binomial_allele_count = 1000, # <<<<<<<<<<<<<<
                   matrix_binomial_ploidy = 1, # <<<<<<<<<<<<<<<<<<<<<<<
                   sfs_bincount = 100, # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                   n_SMCRF_particles = rep(1000, 2), # <<<<<<<<<<<<<<<<
                   min_variant_read_inference_A = NULL,
                   min_variant_read_inference_B = NULL,
                   min_variant_read_validation = NULL,
                   inference_retained_freq = 75,
                   coverage_distribution = "sample-specific",
                   coverage_variables = NULL,
                   neutral_tail = NA,
                   min_N_humps = 1,
                   max_N_humps = Inf,
                   make_readcount_distribution = TRUE,
                   compute_parallel = TRUE,
                   n_cores = NULL) {
    suppressPackageStartupMessages(library(crayon))
    cat(paste0("\n\n\n", bold(red("PERFORMING DECODE FOR SAMPLE ")), bold(yellow(sample_id)), bold(red("...")), "\n"))
    mutation_table$Tot_count <- mutation_table$Ref_count + mutation_table$Alt_count
    mutation_table$VAF <- mutation_table$Alt_count / mutation_table$Tot_count
    if (is.null(max_total_read)) max_total_read <- max(mutation_table$Tot_count)
    #---Choose mutation thresholds, get resulting SFS from data
    SFS_data_frequencies <- seq(1, sfs_bincount) / sfs_bincount
    threshold_results <- choose_mutation_thresholds(
        mutation_table = mutation_table,
        make_readcount_distribution = make_readcount_distribution,
        max_total_read = max_total_read,
        min_variant_read_inference_A = min_variant_read_inference_A,
        min_variant_read_inference_B = min_variant_read_inference_B,
        min_variant_read_validation = min_variant_read_validation,
        SFS_data_frequencies = SFS_data_frequencies,
        inference_retained_freq = inference_retained_freq
    )
    #---Prepare the SFS convolution matrix
    SFS_convolution_inference_A <- build_convolution_matrix(
        sfs_bincount = sfs_bincount,
        mode = "inference A",
        allele_count = matrix_binomial_allele_count,
        min_variant_read = threshold_results$min_variant_read_inference_A,
        min_total_read = threshold_results$min_total_read_inference_A,
        max_total_read = max_total_read,
        matrix_binomial_ploidy = matrix_binomial_ploidy,
        coverage_distribution = coverage_distribution,
        coverage_variables = coverage_variables,
        sample_coverage = threshold_results$sample_coverage_inference_A
    )
    SFS_convolution_inference_B <- build_convolution_matrix(
        sfs_bincount = sfs_bincount,
        mode = "inference B",
        allele_count = matrix_binomial_allele_count,
        min_variant_read = threshold_results$min_variant_read_inference_B,
        min_total_read = threshold_results$min_total_read_inference_B,
        max_total_read = max_total_read,
        matrix_binomial_ploidy = matrix_binomial_ploidy,
        coverage_distribution = coverage_distribution,
        coverage_variables = coverage_variables,
        sample_coverage = threshold_results$sample_coverage_inference_B
    )
    SFS_convolution_validation <- build_convolution_matrix(
        sfs_bincount = sfs_bincount,
        mode = "validation",
        allele_count = matrix_binomial_allele_count,
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
    DECODE_result$criterion_pvalue_threshold <- criterion_pvalue_threshold
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
    # DECODE_result$SFS_convolution_inference_A <- SFS_convolution_inference_A
    # DECODE_result$SFS_convolution_inference_B <- SFS_convolution_inference_B
    # DECODE_result$SFS_convolution_validation <- SFS_convolution_validation
    if (is.na(neutral_tail)) {
        result_with_tail <- DECODE_given_tail_status(
            SFS_data_inference_A = threshold_results$SFS_data_inference_A,
            SFS_data_inference_B = threshold_results$SFS_data_inference_B,
            SFS_data_validation = threshold_results$SFS_data_validation,
            criterion = criterion,
            criterion_penalty_scale = criterion_penalty_scale,
            criterion_pvalue_threshold = criterion_pvalue_threshold,
            min_N_humps = min_N_humps,
            max_N_humps = max_N_humps,
            with_tail = TRUE,
            n_SMCRF_particles = n_SMCRF_particles,
            allele_count = allele_count,
            sfs_bincount = sfs_bincount,
            SFS_convolution_inference_A = SFS_convolution_inference_A,
            SFS_convolution_inference_B = SFS_convolution_inference_B,
            SFS_convolution_validation = SFS_convolution_validation,
            cluster_frequency_mindiff = cluster_frequency_mindiff,
            neutral_power_min = neutral_power_min,
            neutral_power_max = neutral_power_max,
            cluster_frequency_min = cluster_frequency_min,
            cluster_frequency_max = cluster_frequency_max,
            compute_parallel = compute_parallel,
            n_cores = n_cores
        )
        DECODE_result$fits_with_tail <- result_with_tail
        result_without_tail <- DECODE_given_tail_status(
            SFS_data_inference_A = threshold_results$SFS_data_inference_A,
            SFS_data_inference_B = threshold_results$SFS_data_inference_B,
            SFS_data_validation = threshold_results$SFS_data_validation,
            criterion = criterion,
            criterion_penalty_scale = criterion_penalty_scale,
            criterion_pvalue_threshold = criterion_pvalue_threshold,
            min_N_humps = min_N_humps,
            max_N_humps = max_N_humps,
            with_tail = FALSE,
            n_SMCRF_particles = n_SMCRF_particles,
            allele_count = allele_count,
            sfs_bincount = sfs_bincount,
            SFS_convolution_inference_A = SFS_convolution_inference_A,
            SFS_convolution_inference_B = SFS_convolution_inference_B,
            SFS_convolution_validation = SFS_convolution_validation,
            cluster_frequency_mindiff = cluster_frequency_mindiff,
            neutral_power_min = neutral_power_min,
            neutral_power_max = neutral_power_max,
            cluster_frequency_min = cluster_frequency_min,
            cluster_frequency_max = cluster_frequency_max,
            compute_parallel = compute_parallel,
            n_cores = n_cores
        )
        DECODE_result$fits_without_tail <- result_without_tail

        criterion_values_with_tail <- result_with_tail$best_fit$criterion_values
        criterion_values_without_tail <- result_without_tail$best_fit$criterion_values

        tt <- wilcox.test(criterion_values_with_tail, criterion_values_without_tail, alternative = "greater")
        cat("p-value =", signif(tt$p.value, 3), "\n")
        if (tt$p.value > criterion_pvalue_threshold) {
            final_result <- criterion_values_with_tail
        } else {
            final_result <- criterion_values_without_tail
        }



        # if (result_with_tail$best_fit$selected_criterion_value < result_without_tail$best_fit$selected_criterion_value) {
        #     N_humps <- min_N_humps
        #     while (!is.null(DECODE_result$fits_with_tail$all_fits[[paste0(N_humps, "_clusters")]])) {
        #         if (DECODE_result$fits_with_tail$all_fits[[paste0(N_humps, "_clusters")]]$note == "best given tail status") {
        #             DECODE_result$fits_with_tail$all_fits[[paste0(N_humps, "_clusters")]]$note <- "best"
        #         }
        #         N_humps <- N_humps + 1
        #     }
        #     final_result <- result_with_tail
        # } else {
        #     N_humps <- min_N_humps
        #     while (!is.null(DECODE_result$fits_without_tail$all_fits[[paste0(N_humps, "_clusters")]])) {
        #         if (DECODE_result$fits_without_tail$all_fits[[paste0(N_humps, "_clusters")]]$note == "best given tail status") {
        #             DECODE_result$fits_without_tail$all_fits[[paste0(N_humps, "_clusters")]]$note <- "best"
        #         }
        #         N_humps <- N_humps + 1
        #     }
        #     final_result <- result_without_tail
        # }
    } else if (neutral_tail == TRUE) {
        final_result <- DECODE_given_tail_status(
            SFS_data_inference_A = threshold_results$SFS_data_inference_A,
            SFS_data_inference_B = threshold_results$SFS_data_inference_B,
            SFS_data_validation = threshold_results$SFS_data_validation,
            criterion = criterion,
            criterion_penalty_scale = criterion_penalty_scale,
            criterion_pvalue_threshold = criterion_pvalue_threshold,
            min_N_humps = min_N_humps,
            max_N_humps = max_N_humps,
            with_tail = TRUE,
            n_SMCRF_particles = n_SMCRF_particles,
            allele_count = allele_count,
            sfs_bincount = sfs_bincount,
            SFS_convolution_inference_A = SFS_convolution_inference_A,
            SFS_convolution_inference_B = SFS_convolution_inference_B,
            SFS_convolution_validation = SFS_convolution_validation,
            cluster_frequency_mindiff = cluster_frequency_mindiff,
            neutral_power_min = neutral_power_min,
            neutral_power_max = neutral_power_max,
            cluster_frequency_min = cluster_frequency_min,
            cluster_frequency_max = cluster_frequency_max,
            compute_parallel = compute_parallel,
            n_cores = n_cores
        )
        DECODE_result$fits_with_tail <- final_result
        N_humps <- min_N_humps
        # while (!is.null(DECODE_result$fits_with_tail$all_fits[[paste0(N_humps, "_clusters")]])) {
        #     if (DECODE_result$fits_with_tail$all_fits[[paste0(N_humps, "_clusters")]]$note == "best given tail status") {
        #         DECODE_result$fits_with_tail$all_fits[[paste0(N_humps, "_clusters")]]$note <- "best"
        #     }
        #     N_humps <- N_humps + 1
        # }
    } else if (neutral_tail == FALSE) {
        final_result <- DECODE_given_tail_status(
            SFS_data_inference_A = threshold_results$SFS_data_inference_A,
            SFS_data_inference_B = threshold_results$SFS_data_inference_B,
            SFS_data_validation = threshold_results$SFS_data_validation,
            criterion = criterion,
            criterion_penalty_scale = criterion_penalty_scale,
            criterion_pvalue_threshold = criterion_pvalue_threshold,
            min_N_humps = min_N_humps,
            max_N_humps = max_N_humps,
            with_tail = FALSE,
            n_SMCRF_particles = n_SMCRF_particles,
            allele_count = allele_count,
            sfs_bincount = sfs_bincount,
            SFS_convolution_inference_A = SFS_convolution_inference_A,
            SFS_convolution_inference_B = SFS_convolution_inference_B,
            SFS_convolution_validation = SFS_convolution_validation,
            cluster_frequency_mindiff = cluster_frequency_mindiff,
            neutral_power_min = neutral_power_min,
            neutral_power_max = neutral_power_max,
            cluster_frequency_min = cluster_frequency_min,
            cluster_frequency_max = cluster_frequency_max,
            compute_parallel = compute_parallel,
            n_cores = n_cores
        )
        DECODE_result$fits_without_tail <- final_result
        # N_humps <- min_N_humps
        # while (!is.null(DECODE_result$fits_without_tail$all_fits[[paste0(N_humps, "_clusters")]])) {
        #     if (DECODE_result$fits_without_tail$all_fits[[paste0(N_humps, "_clusters")]]$note == "best given tail status") {
        #         DECODE_result$fits_without_tail$all_fits[[paste0(N_humps, "_clusters")]]$note <- "best"
        #     }
        #     N_humps <- N_humps + 1
        # }
    }
    DECODE_result$final_fit <- final_result
    # #---Report the best fit
    # tail_status_final_result <- final_result$best_fit$tail_status
    # parameters_inference_A_final_result <- final_result$best_fit$parameters_inference_A
    # parameters_inference_B_final_result <- final_result$best_fit$parameters_inference_B
    # criterion_final_result <- final_result$best_fit$selected_criterion_value
    # if (tail_status_final_result) {
    #     N_humps_final_result <- length(parameters_inference_A_final_result) / 2 - 1
    #     report <- bold(underline(red(paste0("Best fit = neutral tail + ", N_humps_final_result, " clusters:\n"))))
    #     report <- paste0(report, red("Score            : "), yellow(paste0(criterion, " = ", format(round(criterion_final_result, 3), nsmall = 3))), "\n")
    #     report <- paste0(report, red("Neutral component: "), yellow(paste0("power     = ", format(round(parameters_inference_A_final_result[2], 3), nsmall = 3))), red(", "))
    #     report <- paste0(report, yellow(paste0("\u03C0 = ", format(round(parameters_inference_A_final_result[1], 3), nsmall = 3), " [A], ", format(round(parameters_inference_B_final_result[1], 3), nsmall = 3), " [B]")), "\n")
    #     ii <- 0
    # } else {
    #     N_humps_final_result <- length(parameters_inference_A_final_result) / 2
    #     report <- bold(underline(red(paste0("Best fit = no neutral tail + ", N_humps_final_result, " clusters:\n"))))
    #     report <- paste0(report, red("Score            : "), yellow(paste0(criterion, " = ", format(round(criterion_final_result, 3), nsmall = 3))), "\n")
    #     ii <- -1
    # }
    # if (N_humps_final_result > 0) {
    #     for (i in 1:N_humps_final_result) {
    #         report <- paste0(report, red(paste0("Cluster ", i, "        : ")), yellow(paste0("frequency = ", format(round(parameters_inference_A_final_result[2 * (i + ii) + 2], 3), nsmall = 3))), red(", "))
    #         report <- paste0(report, yellow(paste0("\u03C0 = ", format(round(parameters_inference_A_final_result[2 * (i + ii) + 1], 3), nsmall = 3), " [A], ", format(round(parameters_inference_B_final_result[2 * (i + ii) + 1], 3), nsmall = 3), " [B]")), "\n")
    #     }
    # }
    # cat(report)
    # #---Translation to parameters of cancer evolution in the sample
    # tmp <- parameter_conversion(
    #     result = final_result,
    #     mutation_count_for_fitting = sum(threshold_results$SFS_data_inference_A),
    #     allele_count = allele_count,
    #     matrix_binomial_allele_count = matrix_binomial_allele_count,
    #     matrix_binomial_ploidy = matrix_binomial_ploidy
    # )
    # DECODE_result$final_fit$parameters_df <- tmp$parameters_df
    #---Return the SFS deconvolution results
    return(DECODE_result)
}

DECODE_given_tail_status <- function(SFS_data_inference_A,
                                     SFS_data_inference_B,
                                     SFS_data_validation,
                                     criterion,
                                     criterion_penalty_scale,
                                     criterion_pvalue_threshold,
                                     min_N_humps,
                                     max_N_humps,
                                     with_tail,
                                     n_SMCRF_particles,
                                     allele_count,
                                     sfs_bincount,
                                     SFS_convolution_inference_A,
                                     SFS_convolution_inference_B,
                                     SFS_convolution_validation,
                                     cluster_frequency_mindiff,
                                     neutral_power_min,
                                     neutral_power_max,
                                     cluster_frequency_min,
                                     cluster_frequency_max,
                                     compute_parallel,
                                     n_cores) {
    all_fits <- list()
    criterion_values_list <- list()
    best_criterion_values <- NULL
    best_N_humps <- min_N_humps
    for (N_humps in min_N_humps:max_N_humps) {
        #---Find best parameter set, given the number of humps
        cat(bold(blue(paste0("Inference for ", N_humps, " clusters ", ifelse(with_tail, "with", "without"), " neutral tail component...\n"))))
        fit_results <- DECODE_given_tail_status_and_Ncluster(
            SFS_data_inference_A = SFS_data_inference_A,
            SFS_data_inference_B = SFS_data_inference_B,
            SFS_data_validation = SFS_data_validation,
            criterion = criterion,
            criterion_penalty_scale = criterion_penalty_scale,
            N_humps = N_humps,
            with_tail = with_tail,
            n_SMCRF_particles = n_SMCRF_particles,
            allele_count = allele_count,
            sfs_bincount = sfs_bincount,
            SFS_convolution_inference_A = SFS_convolution_inference_A,
            SFS_convolution_inference_B = SFS_convolution_inference_B,
            SFS_convolution_validation = SFS_convolution_validation,
            cluster_frequency_mindiff = cluster_frequency_mindiff,
            neutral_power_min = neutral_power_min,
            neutral_power_max = neutral_power_max,
            cluster_frequency_min = cluster_frequency_min,
            cluster_frequency_max = cluster_frequency_max,
            compute_parallel = compute_parallel,
            n_cores = n_cores
        )
        ################################################################
        ################################################################
        ################################################################
        ################################################################
        ################################################################
        ################################################################
        ################################################################
        all_fits[[paste0(N_humps, "_clusters")]] <- fit_results

        current_criterion_values <- fit_results$criterion_values$criterion_value
        criterion_values_list[[paste0(N_humps, "_clusters")]] <- current_criterion_values

        if (!is.null(best_criterion_values)) {
            tt <- wilcox.test(best_criterion_values, current_criterion_values, alternative = "greater")
            cat("p-value =", signif(tt$p.value, 3), "\n")
            if (tt$p.value > criterion_pvalue_threshold) {
                cat(sprintf(
                    "→ No significant improvement for %d clusters; best = %d\n",
                    N_humps, N_humps - 1
                ))
                best_N_humps <- N_humps - 1
                break
            }
        }
        best_criterion_values <- current_criterion_values



        # fit_results$note <- "none"
        # all_fits[[paste0(N_humps, "_clusters")]] <- fit_results
        # parameters_inference_A_best_current <- fit_results$best$parameters_inference_A
        # parameters_inference_B_best_current <- fit_results$best$parameters_inference_B
        # parameters_validation_best_current <- fit_results$best$parameters_validation
        # component_distributions_inference_A_best_current <- fit_results$best$component_distributions_inference_A
        # component_distributions_inference_B_best_current <- fit_results$best$component_distributions_inference_B
        # component_distributions_validation_best_current <- fit_results$best$component_distributions_validation
        # criterion_all_best_current <- fit_results$best$criteria
        # criterion_best_current <- criterion_all_best_current[[criterion]]
        # criteria_validation_index_best_current <- fit_results$best$criteria_validation_index
        # #   Report the best fit for the current hump count
        # cluster_pis_inference_A <- Inf
        # cluster_pis_inference_B <- Inf
        # report <- paste0(blue("Score            : "), cyan(paste0(criterion, " = ", format(round(criterion_best_current, 3), nsmall = 3))), "\n")
        # if (with_tail) {
        #     N_humps <- length(parameters_inference_A_best_current) / 2 - 1
        #     report <- paste0(report, blue("Neutral component: "), cyan(paste0("power     = ", format(round(parameters_inference_A_best_current[2], 3), nsmall = 3))), blue(", "))
        #     report <- paste0(report, cyan(paste0("\u03C0 = ", format(round(parameters_inference_A_best_current[1], 3), nsmall = 3), " [A], ", format(round(parameters_inference_B_best_current[1], 3), nsmall = 3), " [B]")), "\n")
        #     ii <- 0
        # } else {
        #     N_humps <- length(parameters_inference_A_best_current) / 2
        #     ii <- -1
        # }
        # if (N_humps > 0) {
        #     for (i in 1:N_humps) {
        #         report <- paste0(report, blue(paste0("Cluster ", i, "        : ")), cyan(paste0("frequency = ", format(round(parameters_inference_A_best_current[2 * (i + ii) + 2], 3), nsmall = 3))), blue(", "))
        #         report <- paste0(report, cyan(paste0("\u03C0 = ", format(round(parameters_inference_A_best_current[2 * (i + ii) + 1], 3), nsmall = 3), " [A], ", format(round(parameters_inference_B_best_current[2 * (i + ii) + 1], 3), nsmall = 3), " [B]")), "\n")
        #         cluster_pis_inference_A <- c(cluster_pis_inference_A, parameters_inference_A_best_current[2 * (i + ii) + 1])
        #         cluster_pis_inference_B <- c(cluster_pis_inference_B, parameters_inference_B_best_current[2 * (i + ii) + 1])
        #     }
        # }
        # cat(report)
        # #   Check if the increased hump count leads to lower criterion score without tiny selective components...
        # if ((N_humps == min_N_humps) | ((criterion_best_current < criterion_pvalue_threshold * criterion_best_final) & (min(pmax(cluster_pis_inference_A, cluster_pis_inference_B)) >= pi_cutoff))) {
        #     #   ... if yes, then update the best fit and continue with 1 more hump
        #     N_humps_best_final <- N_humps
        #     fit_results_best_final <- fit_results
        #     criterion_best_final <- criterion_best_current
        #     criteria_validation_index_best_final <- criteria_validation_index_best_current
        #     criterion_all_final <- criterion_all_best_current
        #     parameters_inference_A_best_final <- parameters_inference_A_best_current
        #     parameters_inference_B_best_final <- parameters_inference_B_best_current
        #     parameters_validation_best_final <- parameters_validation_best_current
        #     component_distributions_inference_A_best_final <- component_distributions_inference_A_best_current
        #     component_distributions_inference_B_best_final <- component_distributions_inference_B_best_current
        #     component_distributions_validation_best_final <- component_distributions_validation_best_current
        #     N_humps <- N_humps + 1
        #     #   ... except if exceeding maximum number of clusters
        #     if (N_humps > max_N_humps) {
        #         break
        #     }
        # } else {
        #     #   ... if no, then stop
        #     break
        # }
    }
    # #---Check if the neutral tail component is too tiny
    # if (with_tail == TRUE & (max(parameters_inference_A_best_final[1], parameters_inference_B_best_final[1]) < pi_cutoff)) {
    #     with_tail <- FALSE
    #     parameters_inference_A_best_final[seq(3, length(parameters_inference_A_best_final), by = 2)] <- parameters_inference_A_best_final[seq(3, length(parameters_inference_A_best_final), by = 2)] / sum(parameters_inference_A_best_final[seq(3, length(parameters_inference_A_best_final), by = 2)])
    #     parameters_inference_A_best_final <- parameters_inference_A_best_final[-c(1, 2)]
    #     parameters_inference_B_best_final[seq(3, length(parameters_inference_B_best_final), by = 2)] <- parameters_inference_B_best_final[seq(3, length(parameters_inference_B_best_final), by = 2)] / sum(parameters_inference_B_best_final[seq(3, length(parameters_inference_B_best_final), by = 2)])
    #     parameters_inference_B_best_final <- parameters_inference_B_best_final[-c(1, 2)]
    #     parameters_validation_best_final[seq(3, length(parameters_validation_best_final), by = 2)] <- parameters_validation_best_final[seq(3, length(parameters_validation_best_final), by = 2)] / sum(parameters_validation_best_final[seq(3, length(parameters_validation_best_final), by = 2)])
    #     parameters_validation_best_final <- parameters_validation_best_final[-c(1, 2)]
    #     component_distributions_inference_A_best_final$SFS_exact[1, ] <- rep(0, length(component_distributions_inference_A_best_final$SFS_exact[1, ]))
    #     component_distributions_inference_A_best_final$SFS_expected[1, ] <- rep(0, length(component_distributions_inference_A_best_final$SFS_expected[1, ]))
    #     component_distributions_inference_A_best_final$SFS_expected_normalized[1, ] <- rep(0, length(component_distributions_inference_A_best_final$SFS_expected_normalized[1, ]))
    #     component_distributions_inference_B_best_final$SFS_exact[1, ] <- rep(0, length(component_distributions_inference_B_best_final$SFS_exact[1, ]))
    #     component_distributions_inference_B_best_final$SFS_expected[1, ] <- rep(0, length(component_distributions_inference_B_best_final$SFS_expected[1, ]))
    #     component_distributions_inference_B_best_final$SFS_expected_normalized[1, ] <- rep(0, length(component_distributions_inference_B_best_final$SFS_expected_normalized[1, ]))
    #     component_distributions_validation_best_final$SFS_exact[1, ] <- rep(0, length(component_distributions_validation_best_final$SFS_exact[1, ]))
    #     component_distributions_validation_best_final$SFS_expected[1, ] <- rep(0, length(component_distributions_validation_best_final$SFS_expected[1, ]))
    #     component_distributions_validation_best_final$SFS_expected_normalized[1, ] <- rep(0, length(component_distributions_validation_best_final$SFS_expected_normalized[1, ]))
    # }
    #---Report the best fit
    # all_fits[[paste0(N_humps_best_final, "_clusters")]]$note <- "best given tail status"
    result <- list()
    result$all_fits <- all_fits
    result$best_fit <- list()
    result$best_fit$selected_criterion <- all_fits[[paste0(best_N_humps, "_clusters")]]$criterion
    result$best_fit$all_criteria <- criterion_values_list
    result$best_fit$selected_criterion_value <- all_fits[[paste0(best_N_humps, "_clusters")]]$criterion_values
    result$best_fit$tail_status <- with_tail
    result$best_fit$criterion_values <- best_criterion_values
    return(result)
}

DECODE_given_tail_status_and_Ncluster <- function(SFS_data_inference_A,
                                                  SFS_data_inference_B,
                                                  SFS_data_validation,
                                                  criterion,
                                                  criterion_penalty_scale,
                                                  N_humps,
                                                  with_tail,
                                                  n_SMCRF_particles,
                                                  allele_count,
                                                  sfs_bincount,
                                                  SFS_convolution_inference_A,
                                                  SFS_convolution_inference_B,
                                                  SFS_convolution_validation,
                                                  cluster_frequency_mindiff,
                                                  neutral_power_min,
                                                  neutral_power_max,
                                                  cluster_frequency_min,
                                                  cluster_frequency_max,
                                                  compute_parallel,
                                                  n_cores,
                                                  progress_bar = TRUE) {
    library(abcsmcrf)
    #---Ingredients for ABC-SMC-DRF
    rprior <- function(Nparameters) {
        parameters <- data.frame(matrix(nrow = Nparameters, ncol = 0))
        if (with_tail) {
            parameters$alpha <- runif(Nparameters, min = neutral_power_min, max = neutral_power_max)
            parameters$omega_inference_A_0 <- runif(Nparameters, min = 0, max = 2 * sum(SFS_data_inference_A))
        }
        for (i in seq_len(N_humps)) {
            parameters[[paste0("p_", i)]] <- runif(Nparameters, min = cluster_frequency_min, max = cluster_frequency_max)
            parameters[[paste0("omega_inference_A_", i)]] <- runif(Nparameters, min = 0, max = 2 * sum(SFS_data_inference_A))
        }
        return(parameters)
    }
    dprior <- function(parameters, parameter_id = "all") {
        probs <- rep(1, nrow(parameters))
        if (with_tail) {
            if (parameter_id %in% c("all", "alpha") && !is.null(parameters$alpha)) {
                probs <- probs * dunif(parameters$alpha, min = neutral_power_min, max = neutral_power_max)
            }
            if (parameter_id %in% c("all", "omega_inference_A_0") && !is.null(parameters$omega_inference_A_0)) {
                probs <- probs * dunif(parameters$omega_inference_A_0, min = 0, max = 2 * sum(SFS_data_inference_A))
            }
        }
        for (i in seq_len(N_humps)) {
            p_col <- paste0("p_", i)
            omega_col <- paste0("omega_inference_A_", i)
            if (parameter_id %in% c("all", p_col) && !is.null(parameters[[p_col]])) {
                probs <- probs * dunif(parameters[[p_col]], min = cluster_frequency_min, max = cluster_frequency_max)
            }
            if (parameter_id %in% c("all", omega_col) && !is.null(parameters[[omega_col]])) {
                probs <- probs * dunif(parameters[[omega_col]], min = 0, max = 2 * sum(SFS_data_inference_A))
            }
        }
        return(probs)
    }
    model <- function(parameters) {
        if (compute_parallel) {
            library(parallel)
            library(pbapply)
            cl <- makePSOCKcluster(ifelse(is.null(n_cores), detectCores() - 1, n_cores))
            clusterExport(cl, varlist = c(
                "one_SFS", "allele_count", "N_humps", "cluster_frequency_mindiff",
                "SFS_convolution_inference_A", "SFS_convolution_inference_B",
                "SFS_data_inference_A", "SFS_data_inference_B",
                "with_tail"
            ), envir = environment())
            stats_list <- parLapply(cl, seq_len(nrow(parameters)), function(i) {
                data.frame(distance = one_SFS(parameters[i, ])$distance)
            })
            stopCluster(cl)
            stats <- do.call(rbind, stats_list)
        } else {
            stats <- do.call(rbind, lapply(seq_len(nrow(parameters)), function(i) {
                data.frame(distance = one_SFS(parameters[i, ])$distance)
            }))
        }
        return(cbind(parameters, stats))
    }
    one_SFS <- function(one_parameter,
                        output_distance = TRUE,
                        output_validation = FALSE,
                        output_SFS = FALSE,
                        output_SFS_components = FALSE,
                        output_all_parameters = FALSE) {
        #   Check that clusters are adequately spaced
        if (N_humps > 1) {
            if (min(diff(sort(as.numeric(one_parameter[paste0("p_", 1:N_humps)])))) < cluster_frequency_mindiff) {
                return(data.frame(distance = NA))
            }
        }
        #   Compute component distributions in the exact SFS
        if (with_tail) {
            dist <- (1:allele_count)^(-as.numeric(one_parameter[["alpha"]]))
            dist <- dist / sum(dist)
            SFS_exact_components <- matrix(dist, nrow = 1)
        } else {
            SFS_exact_components <- matrix(NA, nrow = 0, ncol = allele_count)
        }
        for (i in seq_len(N_humps)) {
            dist <- dbinom(1:allele_count, size = allele_count, prob = as.numeric(one_parameter[paste0("p_", i)])) / (1 - dbinom(0, size = allele_count, prob = as.numeric(one_parameter[paste0("p_", i)])))
            SFS_exact_components <- rbind(SFS_exact_components, dist)
        }
        #   Compute expected SFS for inference A
        SFS_A_components <- SFS_exact_components %*% SFS_convolution_inference_A$convolution_matrix
        if (with_tail) {
            omegas_A <- as.numeric(one_parameter[paste0("omega_inference_A_", 0:N_humps)])
        } else {
            omegas_A <- as.numeric(one_parameter[paste0("omega_inference_A_", 1:N_humps)])
        }
        omegas_exact <- omegas_A / rowSums(SFS_A_components) # omegas_exact <- omegas_A / pmax(rowSums(SFS_A_components), .Machine$double.eps)
        SFS_A_components <- sweep(SFS_A_components, 1, omegas_exact, `*`)
        SFS_A <- colSums(SFS_A_components)
        #   Compute expected SFS for inference B
        SFS_B_components <- SFS_exact_components %*% SFS_convolution_inference_B$convolution_matrix
        SFS_B_components <- sweep(SFS_B_components, 1, omegas_exact, `*`)
        SFS_B <- colSums(SFS_B_components)
        #   Compute expected SFS for validation
        if (output_validation) {
            SFS_validation_components <- SFS_exact_components %*% SFS_convolution_validation$convolution_matrix
            SFS_validation_components <- sweep(SFS_validation_components, 1, omegas_exact, `*`)
            SFS_validation <- colSums(SFS_validation_components)
        }
        #   Statistic = sum of relative L1 errors for inference A & B
        dA <- sum(abs(SFS_A - SFS_data_inference_A)) / sum(abs(SFS_data_inference_A))
        dB <- sum(abs(SFS_B - SFS_data_inference_B)) / sum(abs(SFS_data_inference_B))
        #   Prepare output data
        output <- list()
        if (output_distance) output$distance <- dA + dB
        if (output_SFS) {
            output$SFS_A <- SFS_A
            output$SFS_B <- SFS_B
            if (output_validation) output$SFS_validation <- SFS_validation
        }
        if (output_SFS_components) {
            rownames(SFS_A_components) <- c(if (with_tail) "Tail" else NULL, paste0("Cluster_", 1:N_humps))
            rownames(SFS_B_components) <- c(if (with_tail) "Tail" else NULL, paste0("Cluster_", 1:N_humps))
            if (output_validation) rownames(SFS_validation_components) <- c(if (with_tail) "Tail" else NULL, paste0("Cluster_", 1:N_humps))
            output$SFS_A_components <- SFS_A_components
            output$SFS_B_components <- SFS_B_components
            if (output_validation) output$SFS_validation_components <- SFS_validation_components
        }
        if (output_all_parameters) {
            all_parameters <- data.frame(matrix(ncol = 0, nrow = 1))
            if (with_tail) {
                all_parameters[["Tail_power"]] <- one_parameter[["alpha"]]
                all_parameters[["Tail_Nmut_exact"]] <- omegas_exact[1]
                all_parameters[["Tail_Nmut_A"]] <- sum(SFS_A_components[1, ])
                all_parameters[["Tail_Nmut_B"]] <- sum(SFS_B_components[1, ])
                if (output_validation) all_parameters[["Tail_Nmut_validation"]] <- sum(SFS_validation_components[1, ])
            }
            for (i in seq_len(N_humps)) {
                all_parameters[[paste0("Cluster_", i, "_freq")]] <- one_parameter[[paste0("p_", i)]]
                all_parameters[[paste0("Cluster_", i, "_Nmut_exact")]] <- omegas_exact[i + with_tail]
                all_parameters[[paste0("Cluster_", i, "_Nmut_A")]] <- sum(SFS_A_components[i + with_tail, ])
                all_parameters[[paste0("Cluster_", i, "_Nmut_B")]] <- sum(SFS_B_components[i + with_tail, ])
                if (output_validation) all_parameters[[paste0("Cluster_", i, "_Nmut_validation")]] <- sum(SFS_validation_components[i + with_tail, ])
            }
            output$all_parameters <- all_parameters
        }
        return(output)
    }
    #---ABC-SMC-DRF
    smcrf_result <- smcrf(
        method = "smcrf-multi-param",
        honesty = FALSE, min.node.size = 1,
        statistics_target = data.frame(distance = 0),
        model = model,
        rprior = rprior,
        dprior = dprior,
        nParticles = n_SMCRF_particles,
        model_redo_if_NA = TRUE,
        verbose = FALSE,
        parallel = compute_parallel
    )
    #---Sample parameter sets from ABC-SMC-DRF posterior distribution
    ####################################################################
    validation_count <- 10 ############################################
    ####################################################################
    final_parameters <- smcrf_result[[paste0("Iteration_", length(n_SMCRF_particles))]]$parameters
    final_weights <- smcrf_result[[paste0("Iteration_", length(n_SMCRF_particles))]]$weights
    posterior_parameters <- final_parameters[
        sample(nrow(final_parameters), size = validation_count, prob = final_weights[, 1], replace = T),
    ]
    #---Prepare results based on posterior parameters
    compute_parallel <- FALSE
    if (compute_parallel) {
        library(parallel)
        library(pbapply)
        cl <- makePSOCKcluster(ifelse(is.null(n_cores), detectCores() - 1, n_cores))
        clusterExport(cl, varlist = c(
            "one_SFS", "allele_count", "N_humps", "cluster_frequency_mindiff",
            "SFS_convolution_inference_A", "SFS_convolution_inference_B",
            "SFS_data_inference_A", "SFS_data_inference_B",
            "with_tail"
        ), envir = environment())
        posterior_results_list <- parLapply(cl, seq_len(nrow(posterior_parameters)), function(i) {
            one_SFS(
                posterior_parameters[i, ],
                output_distance = FALSE,
                output_validation = TRUE,
                output_SFS = TRUE,
                output_SFS_components = TRUE,
                output_all_parameters = TRUE
            )
        })
        stopCluster(cl)
    } else {
        posterior_results_list <- lapply(seq_len(nrow(posterior_parameters)), function(i) {
            one_SFS(
                posterior_parameters[i, ],
                output_distance = FALSE,
                output_validation = TRUE,
                output_SFS = TRUE,
                output_SFS_components = TRUE,
                output_all_parameters = TRUE
            )
        })
    }
    SFS_A <- as.data.frame(do.call(rbind, lapply(posterior_results_list, function(x) x$SFS_A)))
    SFS_B <- as.data.frame(do.call(rbind, lapply(posterior_results_list, function(x) x$SFS_B)))
    SFS_validation <- as.data.frame(do.call(rbind, lapply(posterior_results_list, function(x) x$SFS_validation)))
    all_parameters <- as.data.frame(do.call(rbind, lapply(posterior_results_list, function(x) x$all_parameters)))
    if (with_tail) {
        SFS_A_tail <- as.data.frame(do.call(rbind, lapply(posterior_results_list, function(x) x$SFS_A_components["Tail", ])))
        SFS_B_tail <- as.data.frame(do.call(rbind, lapply(posterior_results_list, function(x) x$SFS_B_components["Tail", ])))
        SFS_validation_tail <- as.data.frame(do.call(rbind, lapply(posterior_results_list, function(x) x$SFS_validation_components["Tail", ])))
    }
    SFS_A_clusters <- list()
    SFS_B_clusters <- list()
    SFS_validation_clusters <- list()
    for (i in seq_len(N_humps)) {
        SFS_A_clusters[[i]] <- as.data.frame(do.call(rbind, lapply(posterior_results_list, function(x) x$SFS_A_components[paste0("Cluster_", i), ])))
        SFS_B_clusters[[i]] <- as.data.frame(do.call(rbind, lapply(posterior_results_list, function(x) x$SFS_B_components[paste0("Cluster_", i), ])))
        SFS_validation_clusters[[i]] <- as.data.frame(do.call(rbind, lapply(posterior_results_list, function(x) x$SFS_validation_components[paste0("Cluster_", i), ])))
    }
    #---Compute criterion values
    nParameters <- if (with_tail) 2 * N_humps + 2 else 2 * N_humps
    nData <- sum(SFS_data_validation)
    if (criterion == "BIC") {
        SFS_validation_normalized <- SFS_validation / rowSums(SFS_validation)
        SFS_validation_normalized[which(SFS_validation_normalized <= .Machine$double.eps)] <- .Machine$double.eps
        loglikelihood <- apply(SFS_validation_normalized, 1, function(SFS_pred_validation) sum(log(SFS_pred_validation) * SFS_data_validation))
        criterion_values <- -2 * loglikelihood + criterion_penalty_scale * nParameters * log(nData)
    } else if (criterion == "GIC_L1") {
        err_L1 <- apply(SFS_validation, 1, function(SFS_pred_validation) sum(abs(SFS_pred_validation - SFS_data_validation)))
        criterion_values <- err_L1 + criterion_penalty_scale * nParameters * log(nData)
    } else if (criterion == "GIC_UB") {
        coords <- cbind(0, (seq_len(sfs_bincount) - 0.5) / sfs_bincount)
        target <- transport::wpp(coords, as.numeric(SFS_data_validation))
        UB_Wasserstein_distance <- apply(SFS_validation, 1, function(SFS_pred_validation) {
            source <- transport::wpp(coords, as.numeric(SFS_pred_validation))
            transport::unbalanced(source, target,
                p = 1,
                C = 1,
                method = "revsimplex",
                output = "dist"
            )
        })
        criterion_values <- as.numeric(UB_Wasserstein_distance) + criterion_penalty_scale * nParameters * log(nData)
    }
    criterion_values <- data.frame(criterion_value = criterion_values)
    #---Output results
    fit_results <- list()
    fit_results[["ABCSMCRF_method"]] <- smcrf_result$method
    fit_results[["ABCSMCRF_nParticles"]] <- smcrf_result$n_SMCRF_particles
    fit_results[["criterion"]] <- criterion
    fit_results[["criterion_penalty_scale"]] <- criterion_penalty_scale
    fit_results[["SFS_A"]] <- SFS_A
    fit_results[["SFS_B"]] <- SFS_B
    fit_results[["SFS_validation"]] <- SFS_validation
    fit_results[["parameters"]] <- all_parameters
    if (with_tail) {
        fit_results[["SFS_A_tail"]] <- SFS_A_tail
        fit_results[["SFS_B_tail"]] <- SFS_B_tail
        fit_results[["SFS_validation_tail"]] <- SFS_validation_tail
    }
    for (i in seq_len(N_humps)) {
        fit_results[[paste0("SFS_A_cluster_", i)]] <- SFS_A_clusters[[i]]
        fit_results[[paste0("SFS_B_cluster_", i)]] <- SFS_B_clusters[[i]]
        fit_results[[paste0("SFS_validation_cluster_", i)]] <- SFS_validation_clusters[[i]]
    }
    fit_results[["criterion_values"]] <- criterion_values
    return(fit_results)
}

choose_mutation_thresholds <- function(mutation_table,
                                       make_readcount_distribution,
                                       max_total_read,
                                       min_variant_read_inference_A,
                                       min_variant_read_inference_B,
                                       min_variant_read_validation,
                                       SFS_data_frequencies,
                                       inference_retained_freq) {
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
    if (make_readcount_distribution) {
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
        readcount_distribution$freq <- 100 * readcount_distribution$mutation_count / sum(mutation_table_tmp$count)
    }
    cat("\n")
    #---Find the mutation thresholds for inference and validation
    if (make_readcount_distribution) {
        filtered_df <- readcount_distribution[readcount_distribution$min_total_read == min(readcount_distribution$min_total_read), ]
        min_variant_read_tmp <- min(filtered_df$min_variant_read[filtered_df$freq < 100])
    } else {
        min_variant_read_tmp <- min(mutation_table$Alt_count) + 1
    }
    #   Choose mutation thresholds for inference A set
    if (is.null(min_variant_read_inference_A)) min_variant_read_inference_A <- min_variant_read_tmp + 2
    if (make_readcount_distribution) {
        filtered_df <- readcount_distribution[readcount_distribution$min_variant_read == min_variant_read_inference_A & readcount_distribution$freq >= inference_retained_freq, ]
        min_total_read_inference_A <- ifelse(nrow(filtered_df) > 0, max(filtered_df$min_total_read), min(readcount_distribution$min_total_read))
        Nmut_inference_A <- readcount_distribution$mutation_count[readcount_distribution$min_variant_read == min_variant_read_inference_A & readcount_distribution$min_total_read == min_total_read_inference_A]
        freq_inference_A <- readcount_distribution$freq[readcount_distribution$min_variant_read == min_variant_read_inference_A & readcount_distribution$min_total_read == min_total_read_inference_A]
    } else {
        min_total_read_inference_A <- min(mutation_table$Tot_count)
        Nmut_inference_A <- sum(mutation_table$Alt_count >= min_variant_read_inference_A & mutation_table$Tot_count >= min_total_read_inference_A)
        freq_inference_A <- 100 * Nmut_inference_A / nrow(mutation_table)
        while (freq_inference_A > inference_retained_freq) {
            min_total_read_inference_A <- min_total_read_inference_A + 1
            Nmut_inference_A <- sum(mutation_table$Alt_count >= min_variant_read_inference_A & mutation_table$Tot_count >= min_total_read_inference_A)
            freq_inference_A <- 100 * Nmut_inference_A / nrow(mutation_table)
        }
        min_total_read_inference_A <- min_total_read_inference_A - 1
        Nmut_inference_A <- sum(mutation_table$Alt_count >= min_variant_read_inference_A & mutation_table$Tot_count >= min_total_read_inference_A)
        freq_inference_A <- 100 * Nmut_inference_A / nrow(mutation_table)
    }
    #   Choose mutation thresholds for inference B set
    if (is.null(min_variant_read_inference_B)) min_variant_read_inference_B <- min_variant_read_tmp + 6
    if (make_readcount_distribution) {
        filtered_df <- readcount_distribution[readcount_distribution$min_variant_read == min_variant_read_inference_B & readcount_distribution$freq >= inference_retained_freq, ]
        min_total_read_inference_B <- ifelse(nrow(filtered_df) > 0, max(filtered_df$min_total_read), min(readcount_distribution$min_total_read))
        Nmut_inference_B <- readcount_distribution$mutation_count[readcount_distribution$min_variant_read == min_variant_read_inference_B & readcount_distribution$min_total_read == min_total_read_inference_B]
        freq_inference_B <- readcount_distribution$freq[readcount_distribution$min_variant_read == min_variant_read_inference_B & readcount_distribution$min_total_read == min_total_read_inference_B]
    } else {
        min_total_read_inference_B <- min(mutation_table$Tot_count)
        Nmut_inference_B <- sum(mutation_table$Alt_count >= min_variant_read_inference_B & mutation_table$Tot_count >= min_total_read_inference_B)
        freq_inference_B <- 100 * Nmut_inference_B / nrow(mutation_table)
        while (freq_inference_B > inference_retained_freq) {
            min_total_read_inference_B <- min_total_read_inference_B + 1
            Nmut_inference_B <- sum(mutation_table$Alt_count >= min_variant_read_inference_B & mutation_table$Tot_count >= min_total_read_inference_B)
            freq_inference_B <- 100 * Nmut_inference_B / nrow(mutation_table)
        }
        min_total_read_inference_B <- min_total_read_inference_B - 1
        Nmut_inference_B <- sum(mutation_table$Alt_count >= min_variant_read_inference_B & mutation_table$Tot_count >= min_total_read_inference_B)
        freq_inference_B <- 100 * Nmut_inference_B / nrow(mutation_table)
    }
    #   Choose mutation thresholds for validation set
    if (is.null(min_variant_read_validation)) min_variant_read_validation <- min_variant_read_tmp + 4
    if (make_readcount_distribution) {
        min_total_read_validation <- min(readcount_distribution$min_total_read)
        Nmut_validation <- readcount_distribution$mutation_count[readcount_distribution$min_variant_read == min_variant_read_validation & readcount_distribution$min_total_read == min_total_read_validation]
        freq_validation <- readcount_distribution$freq[readcount_distribution$min_variant_read == min_variant_read_validation & readcount_distribution$min_total_read == min_total_read_validation]
    } else {
        min_total_read_validation <- min(mutation_table$Tot_count)
        Nmut_validation <- sum(mutation_table$Alt_count >= min_variant_read_validation & mutation_table$Tot_count >= min_total_read_validation)
        freq_validation <- 100 * Nmut_validation / nrow(mutation_table)
    }
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
    SFS_data_validation <- func_SFS(mutation_table_validation)
    #---Report the chosen thresholds
    results <- list()
    if (make_readcount_distribution) results$readcount_distribution <- readcount_distribution
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

build_convolution_matrix <- function(sfs_bincount,
                                     mode = NULL,
                                     allele_count,
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
    N_end <- allele_count
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
