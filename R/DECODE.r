DECODE <- function(sample_id = "",
                   mutation_table,
                   criterion = "GIC",
                   criterion_penalty_scale = 0.0015,
                   wilcoxon_pvalue_threshold = 0.01,
                   criterion_Nsamples = 1000,
                   criterion_tail_weight = 0.3, # <<<<<<<<<<<<<<<<<<<<<<
                   neutral_tail = NA,
                   N_clusters = NULL,
                   min_N_clusters = 1,
                   max_N_clusters = 5,
                   neutral_power_min = 0.5, # <<<<<<<<<<<<<<<<<<<<<<<<<<
                   neutral_power_max = 5, # <<<<<<<<<<<<<<<<<<<<<<<<<<<<
                   cluster_frequency_min = 0.01,
                   cluster_frequency_max = 1,
                   max_total_read = NULL,
                   read_distribution_freq_min = 10,
                   allele_count = 1000,
                   matrix_binomial_allele_count = 1000, # <<<<<<<<<<<<<<
                   ploidy = 2, # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                   sfs_bincount = 100, # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                   n_SMCRF_particles = rep(1000, 10), # <<<<<<<<<<<<<<<<
                   min_variant_read_inference_A = NULL,
                   min_variant_read_inference_B = NULL,
                   min_variant_read_validation = NULL,
                   inference_retained_freq = 95,
                   coverage_distribution = "sample-specific",
                   coverage_variables = NULL,
                   compute_parallel = FALSE,
                   n_cores = NULL) {
    if (!is.null(N_clusters)) {
        min_N_clusters <- N_clusters
        max_N_clusters <- N_clusters
    }
    suppressPackageStartupMessages(library(crayon))
    cat(paste0("\n\n\n", bold(red("PERFORMING DECODE FOR SAMPLE ")), bold(yellow(sample_id)), bold(red("...")), "\n"))
    mutation_table$Tot_count <- mutation_table$Ref_count + mutation_table$Alt_count
    mutation_table$VAF <- mutation_table$Alt_count / mutation_table$Tot_count
    if (is.null(max_total_read)) max_total_read <- max(mutation_table$Tot_count)
    #---Choose mutation thresholds, get resulting SFS from data
    SFS_data_frequencies <- seq(1, sfs_bincount) / sfs_bincount
    threshold_results <- choose_mutation_thresholds(
        mutation_table = mutation_table,
        max_total_read = max_total_read,
        min_variant_read_inference_A = min_variant_read_inference_A,
        min_variant_read_inference_B = min_variant_read_inference_B,
        min_variant_read_validation = min_variant_read_validation,
        read_distribution_freq_min = read_distribution_freq_min,
        SFS_data_frequencies = SFS_data_frequencies,
        inference_retained_freq = inference_retained_freq
    )
    if (sum(threshold_results$SFS_data_inference_A) == 0 | (sum(threshold_results$SFS_data_inference_B) == 0) | (sum(threshold_results$SFS_data_validation) == 0)) {
        stop("No mutations passed the filtering criteria. Please adjust the filtering thresholds.")
    }
    #---Prepare the SFS convolution matrix
    SFS_convolution_inference_A <- build_convolution_matrix(
        sfs_bincount = sfs_bincount,
        mode = "inference A",
        allele_count = matrix_binomial_allele_count,
        min_variant_read = threshold_results$min_variant_read_inference_A,
        min_total_read = threshold_results$min_total_read_inference_A,
        max_total_read = max_total_read,
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
        coverage_distribution = coverage_distribution,
        coverage_variables = coverage_variables,
        sample_coverage = threshold_results$sample_coverage_validation
    )
    #---DECODE
    DECODE_result <- list()
    DECODE_result$sample_id <- sample_id
    DECODE_result$mutational_table <- mutation_table
    DECODE_result$criterion <- criterion
    DECODE_result$criterion_penalty_scale <- criterion_penalty_scale
    DECODE_result$wilcoxon_pvalue_threshold <- wilcoxon_pvalue_threshold
    DECODE_result$criterion_Nsamples <- criterion_Nsamples
    DECODE_result$criterion_tail_weight <- criterion_tail_weight
    DECODE_result$neutral_power_min <- neutral_power_min
    DECODE_result$neutral_power_max <- neutral_power_max
    DECODE_result$cluster_frequency_min <- cluster_frequency_min
    DECODE_result$cluster_frequency_max <- cluster_frequency_max
    DECODE_result$max_total_read <- max_total_read
    DECODE_result$allele_count <- allele_count
    DECODE_result$matrix_binomial_allele_count <- matrix_binomial_allele_count
    DECODE_result$ploidy <- ploidy
    DECODE_result$sfs_bincount <- sfs_bincount
    DECODE_result$SFS_frequencies <- SFS_data_frequencies
    DECODE_result$n_SMCRF_particles <- n_SMCRF_particles
    DECODE_result$min_variant_read_inference_A <- threshold_results$min_variant_read_inference_A
    DECODE_result$min_total_read_inference_A <- threshold_results$min_total_read_inference_A
    DECODE_result$min_variant_read_inference_B <- threshold_results$min_variant_read_inference_B
    DECODE_result$min_total_read_inference_B <- threshold_results$min_total_read_inference_B
    DECODE_result$min_variant_read_validation <- threshold_results$min_variant_read_validation
    DECODE_result$min_total_read_validation <- threshold_results$min_total_read_validation
    DECODE_result$inference_retained_freq <- inference_retained_freq
    DECODE_result$coverage_distribution <- coverage_distribution
    DECODE_result$coverage_variables <- coverage_variables
    DECODE_result$neutral_tail <- neutral_tail
    DECODE_result$min_N_clusters <- min_N_clusters
    DECODE_result$max_N_clusters <- max_N_clusters
    DECODE_result$compute_parallel <- compute_parallel
    DECODE_result$n_cores <- n_cores
    DECODE_result$readcount_distribution <- threshold_results$readcount_distribution
    DECODE_result$SFS_data_inference_A <- threshold_results$SFS_data_inference_A
    DECODE_result$SFS_data_inference_B <- threshold_results$SFS_data_inference_B
    DECODE_result$SFS_data_validation <- threshold_results$SFS_data_validation
    if (is.na(neutral_tail)) {
        result_with_tail <- DECODE_given_tail_status(
            SFS_data_inference_A = threshold_results$SFS_data_inference_A,
            SFS_data_inference_B = threshold_results$SFS_data_inference_B,
            SFS_data_validation = threshold_results$SFS_data_validation,
            ploidy = ploidy,
            criterion = criterion,
            criterion_penalty_scale = criterion_penalty_scale,
            wilcoxon_pvalue_threshold = wilcoxon_pvalue_threshold,
            criterion_Nsamples = criterion_Nsamples,
            criterion_tail_weight = criterion_tail_weight,
            min_N_clusters = min_N_clusters,
            max_N_clusters = max_N_clusters,
            with_tail = TRUE,
            n_SMCRF_particles = n_SMCRF_particles,
            allele_count = allele_count,
            sfs_bincount = sfs_bincount,
            SFS_convolution_inference_A = SFS_convolution_inference_A,
            SFS_convolution_inference_B = SFS_convolution_inference_B,
            SFS_convolution_validation = SFS_convolution_validation,
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
            ploidy = ploidy,
            criterion = criterion,
            criterion_penalty_scale = criterion_penalty_scale,
            wilcoxon_pvalue_threshold = wilcoxon_pvalue_threshold,
            criterion_Nsamples = criterion_Nsamples,
            criterion_tail_weight = criterion_tail_weight,
            min_N_clusters = max(1, min_N_clusters),
            max_N_clusters = max(1, min_N_clusters, max_N_clusters),
            with_tail = FALSE,
            n_SMCRF_particles = n_SMCRF_particles,
            allele_count = allele_count,
            sfs_bincount = sfs_bincount,
            SFS_convolution_inference_A = SFS_convolution_inference_A,
            SFS_convolution_inference_B = SFS_convolution_inference_B,
            SFS_convolution_validation = SFS_convolution_validation,
            neutral_power_min = neutral_power_min,
            neutral_power_max = neutral_power_max,
            cluster_frequency_min = cluster_frequency_min,
            cluster_frequency_max = cluster_frequency_max,
            compute_parallel = compute_parallel,
            n_cores = n_cores
        )
        DECODE_result$fits_without_tail <- result_without_tail
        criterion_values_with_tail <- result_with_tail$all_fits[[paste0(result_with_tail$best_N_clusters, "_clusters")]]$criterion_values$criterion_value
        criterion_values_without_tail <- result_without_tail$all_fits[[paste0(result_without_tail$best_N_clusters, "_clusters")]]$criterion_values$criterion_value
        N_clusters_with_tail <- result_with_tail$best_N_clusters
        N_clusters_without_tail <- result_without_tail$best_N_clusters
        tt <- wilcox.test(criterion_values_with_tail, criterion_values_without_tail, alternative = "greater")
        report <- paste0(bold(blue("With vs without tail   : ")), cyan(paste0("Wilcoxon p-value = ", signif(tt$p.value, 3))))
        if (tt$p.value > wilcoxon_pvalue_threshold) {
            report <- paste0(report, yellow(paste0(" → fit with tail is selected")), "\n")
            cat(report)
            DECODE_result$best_with_tail <- TRUE
            DECODE_result$best_N_clusters <- N_clusters_with_tail
        } else {
            report <- paste0(report, yellow(paste0(" → fit without tail is selected")), "\n")
            cat(report)
            DECODE_result$best_with_tail <- FALSE
            DECODE_result$best_N_clusters <- N_clusters_without_tail
        }
    } else if (neutral_tail == TRUE) {
        result_with_tail <- DECODE_given_tail_status(
            SFS_data_inference_A = threshold_results$SFS_data_inference_A,
            SFS_data_inference_B = threshold_results$SFS_data_inference_B,
            SFS_data_validation = threshold_results$SFS_data_validation,
            ploidy = ploidy,
            criterion = criterion,
            criterion_penalty_scale = criterion_penalty_scale,
            wilcoxon_pvalue_threshold = wilcoxon_pvalue_threshold,
            criterion_Nsamples = criterion_Nsamples,
            criterion_tail_weight = criterion_tail_weight,
            min_N_clusters = min_N_clusters,
            max_N_clusters = max_N_clusters,
            with_tail = TRUE,
            n_SMCRF_particles = n_SMCRF_particles,
            allele_count = allele_count,
            sfs_bincount = sfs_bincount,
            SFS_convolution_inference_A = SFS_convolution_inference_A,
            SFS_convolution_inference_B = SFS_convolution_inference_B,
            SFS_convolution_validation = SFS_convolution_validation,
            neutral_power_min = neutral_power_min,
            neutral_power_max = neutral_power_max,
            cluster_frequency_min = cluster_frequency_min,
            cluster_frequency_max = cluster_frequency_max,
            compute_parallel = compute_parallel,
            n_cores = n_cores
        )
        DECODE_result$fits_with_tail <- result_with_tail
        DECODE_result$best_with_tail <- TRUE
        DECODE_result$best_N_clusters <- result_with_tail$best_N_clusters
    } else if (neutral_tail == FALSE) {
        result_without_tail <- DECODE_given_tail_status(
            SFS_data_inference_A = threshold_results$SFS_data_inference_A,
            SFS_data_inference_B = threshold_results$SFS_data_inference_B,
            SFS_data_validation = threshold_results$SFS_data_validation,
            ploidy = ploidy,
            criterion = criterion,
            criterion_penalty_scale = criterion_penalty_scale,
            wilcoxon_pvalue_threshold = wilcoxon_pvalue_threshold,
            criterion_Nsamples = criterion_Nsamples,
            criterion_tail_weight = criterion_tail_weight,
            min_N_clusters = max(1, min_N_clusters),
            max_N_clusters = max(1, min_N_clusters, max_N_clusters),
            with_tail = FALSE,
            n_SMCRF_particles = n_SMCRF_particles,
            allele_count = allele_count,
            sfs_bincount = sfs_bincount,
            SFS_convolution_inference_A = SFS_convolution_inference_A,
            SFS_convolution_inference_B = SFS_convolution_inference_B,
            SFS_convolution_validation = SFS_convolution_validation,
            neutral_power_min = neutral_power_min,
            neutral_power_max = neutral_power_max,
            cluster_frequency_min = cluster_frequency_min,
            cluster_frequency_max = cluster_frequency_max,
            compute_parallel = compute_parallel,
            n_cores = n_cores
        )
        DECODE_result$fits_without_tail <- result_without_tail
        DECODE_result$best_with_tail <- FALSE
        DECODE_result$best_N_clusters <- result_without_tail$best_N_clusters
    }
    #---Report the best fit
    with_tail <- DECODE_result$best_with_tail
    if (DECODE_result$best_with_tail) {
        fit_results <- DECODE_result$fits_with_tail$all_fits[[paste0(DECODE_result$best_N_clusters, "_clusters")]]
    } else {
        fit_results <- DECODE_result$fits_without_tail$all_fits[[paste0(DECODE_result$best_N_clusters, "_clusters")]]
    }
    report <- bold(underline(red(paste0("Best fit = ", ifelse(with_tail, "", "without "), "neutral tail + ", DECODE_result$best_N_clusters, " clusters:\n"))))
    report <- paste0(report, red("Validation score       : "), yellow(paste0(criterion, " = ", format(round(mean(fit_results$criterion_values$criterion_value), 2), nsmall = 2), " \u00B1 ", format(round(sd(fit_results$criterion_values$criterion_value), 2)))))
    if (criterion == "GIC") {
        report <- paste0(report, yellow(paste0(" (L1 error = ", format(round(mean(fit_results$criterion_values$L1_error), 2), nsmall = 2), " \u00B1 ", format(round(sd(fit_results$criterion_values$L1_error), 2)), ", penalty = ", format(round(mean(fit_results$criterion_values$penalty), 2), nsmall = 2), " \u00B1 ", format(round(sd(fit_results$criterion_values$penalty), 2)), ")")))
    }
    report <- paste0(report, "\n")
    if (with_tail) {
        report <- paste0(report, red("Neutral component      : "), yellow(paste0("\u03B1   = ", format(round(mean(fit_results$parameters[["Tail_power"]]), 2), nsmall = 2), " \u00B1 ", format(round(sd(fit_results$parameters[["Tail_power"]]), 2), nsmall = 2))), red("; "))
        report <- paste0(report, yellow(paste0("mutation count: ", round(mean(fit_results$parameters[["Tail_Nmut_exact"]])), " \u00B1 ", round(sd(fit_results$parameters[["Tail_Nmut_exact"]])), " [True], ", round(mean(fit_results$parameters[["Tail_Nmut_A"]])), " \u00B1 ", round(sd(fit_results$parameters[["Tail_Nmut_A"]])), " [A], ", round(mean(fit_results$parameters[["Tail_Nmut_B"]])), " \u00B1 ", round(sd(fit_results$parameters[["Tail_Nmut_B"]])), " [B]"), "\n"))
    }
    for (i in seq_len(DECODE_result$best_N_clusters)) {
        report <- paste0(report, red(paste0("Cluster ", i, "              : ")), yellow(paste0("f   = ", format(round(mean(fit_results$parameters[[paste0("Cluster_", i, "_freq")]]), 2), nsmall = 2), " \u00B1 ", format(round(sd(fit_results$parameters[[paste0("Cluster_", i, "_freq")]]), 2), nsmall = 2))), red("; "))
        report <- paste0(report, yellow(paste0("mutation count: ", round(mean(fit_results$parameters[[paste0("Cluster_", i, "_Nmut_exact")]])), " \u00B1 ", round(sd(fit_results$parameters[[paste0("Cluster_", i, "_Nmut_exact")]])), " [True], ", round(mean(fit_results$parameters[[paste0("Cluster_", i, "_Nmut_A")]])), " \u00B1 ", round(sd(fit_results$parameters[[paste0("Cluster_", i, "_Nmut_A")]])), " [A], ", round(mean(fit_results$parameters[[paste0("Cluster_", i, "_Nmut_B")]])), " \u00B1 ", round(sd(fit_results$parameters[[paste0("Cluster_", i, "_Nmut_B")]])), " [B]"), "\n"))
    }
    cat(report)
    #---Return the SFS deconvolution results
    return(DECODE_result)
}

DECODE_given_tail_status <- function(SFS_data_inference_A,
                                     SFS_data_inference_B,
                                     SFS_data_validation,
                                     ploidy,
                                     criterion,
                                     criterion_penalty_scale,
                                     wilcoxon_pvalue_threshold,
                                     criterion_Nsamples,
                                     criterion_tail_weight,
                                     min_N_clusters,
                                     max_N_clusters,
                                     with_tail,
                                     n_SMCRF_particles,
                                     allele_count,
                                     sfs_bincount,
                                     SFS_convolution_inference_A,
                                     SFS_convolution_inference_B,
                                     SFS_convolution_validation,
                                     neutral_power_min,
                                     neutral_power_max,
                                     cluster_frequency_min,
                                     cluster_frequency_max,
                                     compute_parallel,
                                     n_cores) {
    all_fits <- list()
    best_criterion_values <- NULL
    for (N_clusters in min_N_clusters:max_N_clusters) {
        #---Find best parameter set, given the number of humps
        cat(bold(blue(paste0("Inference for ", N_clusters, " clusters ", ifelse(with_tail, "with", "without"), " neutral tail component...\n"))))
        fit_results <- DECODE_given_tail_status_and_Ncluster(
            SFS_data_inference_A = SFS_data_inference_A,
            SFS_data_inference_B = SFS_data_inference_B,
            SFS_data_validation = SFS_data_validation,
            ploidy = ploidy,
            criterion = criterion,
            criterion_penalty_scale = criterion_penalty_scale,
            criterion_Nsamples = criterion_Nsamples,
            criterion_tail_weight = criterion_tail_weight,
            N_clusters = N_clusters,
            with_tail = with_tail,
            n_SMCRF_particles = n_SMCRF_particles,
            allele_count = allele_count,
            sfs_bincount = sfs_bincount,
            SFS_convolution_inference_A = SFS_convolution_inference_A,
            SFS_convolution_inference_B = SFS_convolution_inference_B,
            SFS_convolution_validation = SFS_convolution_validation,
            neutral_power_min = neutral_power_min,
            neutral_power_max = neutral_power_max,
            cluster_frequency_min = cluster_frequency_min,
            cluster_frequency_max = cluster_frequency_max,
            compute_parallel = compute_parallel,
            n_cores = n_cores
        )
        if (N_clusters > min_N_clusters) {
            fit_results_assisted <- DECODE_given_tail_status_and_Ncluster(
                SFS_data_inference_A = SFS_data_inference_A,
                SFS_data_inference_B = SFS_data_inference_B,
                SFS_data_validation = SFS_data_validation,
                previous_fit = all_fits[[paste0(N_clusters - 1, "_clusters")]],
                ploidy = ploidy,
                criterion = criterion,
                criterion_penalty_scale = criterion_penalty_scale,
                criterion_Nsamples = criterion_Nsamples,
                criterion_tail_weight = criterion_tail_weight,
                N_clusters = N_clusters,
                with_tail = with_tail,
                n_SMCRF_particles = n_SMCRF_particles,
                allele_count = allele_count,
                sfs_bincount = sfs_bincount,
                SFS_convolution_inference_A = SFS_convolution_inference_A,
                SFS_convolution_inference_B = SFS_convolution_inference_B,
                SFS_convolution_validation = SFS_convolution_validation,
                neutral_power_min = neutral_power_min,
                neutral_power_max = neutral_power_max,
                cluster_frequency_min = cluster_frequency_min,
                cluster_frequency_max = cluster_frequency_max,
                compute_parallel = compute_parallel,
                n_cores = n_cores
            )
            if (!is.null(fit_results_assisted) && mean(fit_results_assisted[["distances"]]) < mean(fit_results[["distances"]])) fit_results <- fit_results_assisted
        }
        all_fits[[paste0(N_clusters, "_clusters")]] <- fit_results
        #   Report the best fit for the current hump count
        report <- paste0(blue("Validation score       : "), cyan(paste0(criterion, " = ", format(round(mean(fit_results$criterion_values$criterion_value), 2), nsmall = 2), " \u00B1 ", format(round(sd(fit_results$criterion_values$criterion_value), 2)))))
        if (criterion == "GIC") {
            report <- paste0(report, cyan(paste0(" (L1 error = ", format(round(mean(fit_results$criterion_values$L1_error), 2), nsmall = 2), " \u00B1 ", format(round(sd(fit_results$criterion_values$L1_error), 2)), ", penalty = ", format(round(mean(fit_results$criterion_values$penalty), 2), nsmall = 2), " \u00B1 ", format(round(sd(fit_results$criterion_values$penalty), 2)), ")")))
        }
        report <- paste0(report, "\n")
        if (with_tail) {
            report <- paste0(report, blue("Neutral component      : "), cyan(paste0("\u03B1   = ", format(round(mean(fit_results$parameters[["Tail_power"]]), 2), nsmall = 2), " \u00B1 ", format(round(sd(fit_results$parameters[["Tail_power"]]), 2), nsmall = 2))), blue("; "))
            report <- paste0(report, cyan(paste0("mutation count: ", round(mean(fit_results$parameters[["Tail_Nmut_exact"]])), " \u00B1 ", round(sd(fit_results$parameters[["Tail_Nmut_exact"]])), " [True], ", round(mean(fit_results$parameters[["Tail_Nmut_A"]])), " \u00B1 ", round(sd(fit_results$parameters[["Tail_Nmut_A"]])), " [A], ", round(mean(fit_results$parameters[["Tail_Nmut_B"]])), " \u00B1 ", round(sd(fit_results$parameters[["Tail_Nmut_B"]])), " [B]"), "\n"))
        }
        for (i in seq_len(N_clusters)) {
            report <- paste0(report, blue(paste0("Cluster ", i, "              : ")), cyan(paste0("f   = ", format(round(mean(fit_results$parameters[[paste0("Cluster_", i, "_freq")]]), 2), nsmall = 2), " \u00B1 ", format(round(sd(fit_results$parameters[[paste0("Cluster_", i, "_freq")]]), 2), nsmall = 2))), blue("; "))
            report <- paste0(report, cyan(paste0("mutation count: ", round(mean(fit_results$parameters[[paste0("Cluster_", i, "_Nmut_exact")]])), " \u00B1 ", round(sd(fit_results$parameters[[paste0("Cluster_", i, "_Nmut_exact")]])), " [True], ", round(mean(fit_results$parameters[[paste0("Cluster_", i, "_Nmut_A")]])), " \u00B1 ", round(sd(fit_results$parameters[[paste0("Cluster_", i, "_Nmut_A")]])), " [A], ", round(mean(fit_results$parameters[[paste0("Cluster_", i, "_Nmut_B")]])), " \u00B1 ", round(sd(fit_results$parameters[[paste0("Cluster_", i, "_Nmut_B")]])), " [B]"), "\n"))
        }
        cat(report)
        #   Statistical test to determine whether to continue with higher cluster counts
        current_criterion_values <- fit_results$criterion_values$criterion_value
        if (!is.null(best_criterion_values)) {
            tt <- wilcox.test(best_criterion_values, current_criterion_values, alternative = "greater")
            report <- paste0(bold(blue(paste0("H = ", N_clusters, " vs H = ", N_clusters - 1, "         : "))), cyan(paste0("Wilcoxon p-value = ", signif(tt$p.value, 3))))
            if (tt$p.value > wilcoxon_pvalue_threshold) {
                report <- paste0(report, yellow(paste0(" → ", N_clusters - 1, " clusters is the best fit ", ifelse(with_tail, "with", "without"), " tail")), "\n")
                cat(report)
                break
            } else {
                if (N_clusters == max_N_clusters) {
                    report <- paste0(report, yellow(" → maximum cluster count has been reached, increased max_N_clusters may improve fit"), "\n")
                } else {
                    report <- paste0(report, blue(" → continue with higher cluster counts"), "\n")
                }
                cat(report)
            }
        }
        best_criterion_values <- current_criterion_values
        best_N_clusters <- N_clusters
    }
    #---Report the best fit
    result <- list()
    result$all_fits <- all_fits
    result$best_N_clusters <- best_N_clusters
    return(result)
}

DECODE_given_tail_status_and_Ncluster <- function(SFS_data_inference_A,
                                                  SFS_data_inference_B,
                                                  SFS_data_validation,
                                                  previous_fit = NULL,
                                                  ploidy,
                                                  criterion,
                                                  criterion_penalty_scale,
                                                  criterion_Nsamples,
                                                  criterion_tail_weight,
                                                  N_clusters,
                                                  with_tail,
                                                  n_SMCRF_particles,
                                                  allele_count,
                                                  sfs_bincount,
                                                  SFS_convolution_inference_A,
                                                  SFS_convolution_inference_B,
                                                  SFS_convolution_validation,
                                                  neutral_power_min,
                                                  neutral_power_max,
                                                  cluster_frequency_min,
                                                  cluster_frequency_max,
                                                  compute_parallel,
                                                  n_cores,
                                                  progress_bar = TRUE) {
    library(truncnorm)
    #---Ingredients for ABC-SMC-DRF
    if (!is.null(previous_fit)) {
        #---Initial clustering composition based on previous fit
        SFS_data_frequencies <- seq(1, sfs_bincount) / sfs_bincount
        SFS_inference_A_residual <- pmax(0, SFS_data_inference_A - colMeans(previous_fit[["SFS_inference_A"]]))
        if (all(SFS_inference_A_residual == 0)) {
            return(NULL)
        }
        SFS_inference_A_residual <- sapply(1:length(SFS_inference_A_residual), function(i) {
            mean(SFS_inference_A_residual[pmax(1, i - 2):min(length(SFS_inference_A_residual), i + 2)])
        })
        initial_guess <- colMeans(previous_fit[["posterior_parameters"]])
        initial_guess[[paste0("p_", N_clusters)]] <- SFS_data_frequencies[which.max(SFS_inference_A_residual)]
        initial_guess[[paste0("omega_inference_A_", N_clusters)]] <- sum(SFS_inference_A_residual)
        rprior <- function(Nparameters) {
            parameters <- data.frame(matrix(nrow = Nparameters, ncol = 0))
            if (with_tail) {
                parameters[["alpha"]] <- rtruncnorm(
                    n = Nparameters,
                    a = neutral_power_min, b = neutral_power_max,
                    mean = initial_guess[["alpha"]],
                    sd = 0.05 # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                )
                parameters[["omega_inference_A_0"]] <- rtruncnorm(
                    n = Nparameters,
                    a = 0, b = 2 * sum(SFS_data_inference_A),
                    mean = initial_guess[["omega_inference_A_0"]],
                    sd = 0.05 * sum(SFS_data_inference_A) # <<<<<<<<<<<<
                )
            }
            for (i in seq_len(N_clusters)) {
                parameters[[paste0("p_", i)]] <- rtruncnorm(
                    n = Nparameters,
                    a = cluster_frequency_min, b = cluster_frequency_max,
                    mean = initial_guess[[paste0("p_", i)]],
                    sd = 0.05 # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                )
                parameters[[paste0("omega_inference_A_", i)]] <- rtruncnorm(
                    n = Nparameters,
                    a = 0, b = 2 * sum(SFS_data_inference_A),
                    mean = initial_guess[[paste0("omega_inference_A_", i)]],
                    sd = 0.05 * sum(SFS_data_inference_A) # <<<<<<<<<<<<
                )
            }
            #---Order parameter sets in decreasing cluster VAFs
            if (N_clusters > 1) {
                p_colnames <- paste0("p_", 1:N_clusters)
                omega_colnames <- paste0("omega_inference_A_", 1:N_clusters)
                reorder_clusters <- function(row) {
                    p_values <- as.numeric(row[p_colnames])
                    omega_values <- as.numeric(row[omega_colnames])
                    order_idx <- order(p_values, decreasing = TRUE)
                    new_row <- row
                    new_row[p_colnames] <- p_values[order_idx]
                    new_row[omega_colnames] <- omega_values[order_idx]
                    return(new_row)
                }
                parameters <- as.data.frame(t(apply(parameters, 1, reorder_clusters)))
                colnames(parameters) <- colnames(parameters)
            }
            return(parameters)
        }
        dprior <- function(parameters, parameter_id = "all") {
            probs <- rep(1, nrow(parameters))
            if (with_tail) {
                if (parameter_id %in% c("all", "alpha") && !is.null(parameters[["alpha"]])) {
                    probs <- probs * dtruncnorm(
                        parameters[["alpha"]],
                        a = neutral_power_min, b = neutral_power_max,
                        mean = initial_guess[["alpha"]],
                        sd = 0.05 # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                    )
                }
                if (parameter_id %in% c("all", "omega_inference_A_0") && !is.null(parameters[["omega_inference_A_0"]])) {
                    probs <- probs * dtruncnorm(
                        parameters[["omega_inference_A_0"]],
                        a = 0, b = 2 * sum(SFS_data_inference_A),
                        mean = initial_guess[["omega_inference_A_0"]],
                        sd = 0.05 * sum(SFS_data_inference_A) # <<<<<<<<
                    )
                }
            }
            for (i in seq_len(N_clusters)) {
                p_col <- paste0("p_", i)
                omega_col <- paste0("omega_inference_A_", i)
                if (parameter_id %in% c("all", p_col) && !is.null(parameters[[p_col]])) {
                    probs <- probs * dtruncnorm(
                        parameters[[p_col]],
                        a = cluster_frequency_min, b = cluster_frequency_max,
                        mean = initial_guess[[p_col]],
                        sd = 0.05 # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                    )
                }
                if (parameter_id %in% c("all", omega_col) && !is.null(parameters[[omega_col]])) {
                    probs <- probs * dtruncnorm(
                        parameters[[omega_col]],
                        a = 0, b = 2 * sum(SFS_data_inference_A),
                        mean = initial_guess[[omega_col]],
                        sd = 0.05 * sum(SFS_data_inference_A) # <<<<<<<<
                    )
                }
            }
            return(probs)
        }
    } else {
        rprior <- function(Nparameters) {
            parameters <- data.frame(matrix(nrow = Nparameters, ncol = 0))
            if (with_tail) {
                parameters[["alpha"]] <- runif(Nparameters, min = neutral_power_min, max = neutral_power_max)
                parameters[["omega_inference_A_0"]] <- runif(Nparameters, min = 0, max = 2 * sum(SFS_data_inference_A))
            }
            for (i in seq_len(N_clusters)) {
                parameters[[paste0("p_", i)]] <- runif(Nparameters, min = cluster_frequency_min, max = cluster_frequency_max)
                parameters[[paste0("omega_inference_A_", i)]] <- runif(Nparameters, min = 0, max = 2 * sum(SFS_data_inference_A))
            }
            #---Order parameter sets in decreasing cluster VAFs
            if (N_clusters > 1) {
                p_colnames <- paste0("p_", 1:N_clusters)
                omega_colnames <- paste0("omega_inference_A_", 1:N_clusters)
                reorder_clusters <- function(row) {
                    p_values <- as.numeric(row[p_colnames])
                    omega_values <- as.numeric(row[omega_colnames])
                    order_idx <- order(p_values, decreasing = TRUE)
                    new_row <- row
                    new_row[p_colnames] <- p_values[order_idx]
                    new_row[omega_colnames] <- omega_values[order_idx]
                    return(new_row)
                }
                parameters <- as.data.frame(t(apply(parameters, 1, reorder_clusters)))
                colnames(parameters) <- colnames(parameters)
            }
            return(parameters)
        }
        dprior <- function(parameters, parameter_id = "all") {
            probs <- rep(1, nrow(parameters))
            if (with_tail) {
                if (parameter_id %in% c("all", "alpha") && !is.null(parameters[["alpha"]])) {
                    probs <- probs * dunif(parameters[["alpha"]], min = neutral_power_min, max = neutral_power_max)
                }
                if (parameter_id %in% c("all", "omega_inference_A_0") && !is.null(parameters[["omega_inference_A_0"]])) {
                    probs <- probs * dunif(parameters[["omega_inference_A_0"]], min = 0, max = 2 * sum(SFS_data_inference_A))
                }
            }
            for (i in seq_len(N_clusters)) {
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
    }
    rperturb <- function(parameters_unperturbed, parameters_previous_sampled, iteration) {
        Beaumont_variances <- 2 * pmax(sapply(parameters_previous_sampled, var), 1e-10)
        parameters_perturbed <- parameters_unperturbed
        if (with_tail) {
            parameters_perturbed[["alpha"]] <- rtruncnorm(
                n = nrow(parameters_perturbed),
                a = neutral_power_min, b = neutral_power_max,
                mean = parameters_perturbed[["alpha"]],
                sd = sqrt(Beaumont_variances[["alpha"]])
            )
            parameters_perturbed[["omega_inference_A_0"]] <- rtruncnorm(
                n = nrow(parameters_perturbed),
                a = 0, b = 2 * sum(SFS_data_inference_A),
                mean = parameters_perturbed[["omega_inference_A_0"]],
                sd = sqrt(Beaumont_variances[["omega_inference_A_0"]])
            )
        }
        for (i in seq_len(N_clusters)) {
            p_col <- paste0("p_", i)
            omega_col <- paste0("omega_inference_A_", i)
            parameters_perturbed[[p_col]] <- rtruncnorm(
                n = nrow(parameters_perturbed),
                a = cluster_frequency_min, b = cluster_frequency_max,
                mean = parameters_perturbed[[p_col]],
                sd = sqrt(Beaumont_variances[[p_col]])
            )
            parameters_perturbed[[omega_col]] <- rtruncnorm(
                n = nrow(parameters_perturbed),
                a = 0, b = 2 * sum(SFS_data_inference_A),
                mean = parameters_perturbed[[omega_col]],
                sd = sqrt(Beaumont_variances[[omega_col]])
            )
        }
        return(parameters_perturbed)
    }
    dperturb <- function(parameters, parameters_previous, parameters_previous_sampled, iteration, parameter_id = "all") {
        Beaumont_variances <- 2 * pmax(sapply(parameters_previous_sampled, var), 1e-10)
        probs <- rep(1, nrow(parameters))
        if (with_tail) {
            if (parameter_id %in% c("all", "alpha") && !is.null(parameters[["alpha"]])) {
                probs <- probs * dtruncnorm(
                    parameters[["alpha"]],
                    a = neutral_power_min, b = neutral_power_max,
                    mean = parameters_previous[["alpha"]],
                    sd = sqrt(Beaumont_variances[["alpha"]])
                )
            }
            if (parameter_id %in% c("all", "omega_inference_A_0") && !is.null(parameters[["omega_inference_A_0"]])) {
                probs <- probs * dtruncnorm(
                    parameters[["omega_inference_A_0"]],
                    a = 0, b = 2 * sum(SFS_data_inference_A),
                    mean = parameters_previous[["omega_inference_A_0"]],
                    sd = sqrt(Beaumont_variances[["omega_inference_A_0"]])
                )
            }
        }
        for (i in seq_len(N_clusters)) {
            p_col <- paste0("p_", i)
            omega_col <- paste0("omega_inference_A_", i)
            if (parameter_id %in% c("all", p_col) && !is.null(parameters[[p_col]])) {
                probs <- probs * dtruncnorm(
                    parameters[[p_col]],
                    a = cluster_frequency_min, b = cluster_frequency_max,
                    mean = parameters_previous[[p_col]],
                    sd = sqrt(Beaumont_variances[[p_col]])
                )
            }
            if (parameter_id %in% c("all", omega_col) && !is.null(parameters[[omega_col]])) {
                probs <- probs * dtruncnorm(
                    parameters[[omega_col]],
                    a = 0, b = 2 * sum(SFS_data_inference_A),
                    mean = parameters_previous[[omega_col]],
                    sd = sqrt(Beaumont_variances[[omega_col]])
                )
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
                "one_SFS", "allele_count", "N_clusters",
                "SFS_convolution_inference_A", "SFS_convolution_inference_B",
                "SFS_data_inference_A", "SFS_data_inference_B",
                "with_tail", "ploidy"
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
        #   Compute component distributions in the exact SFS
        if (with_tail) {
            dist <- c(1:round(allele_count * max(as.numeric(one_parameter[paste0("p_", 1:N_clusters)]))))^(-as.numeric(one_parameter[["alpha"]]))
            dist <- c(dist, rep(0, allele_count - length(dist)))
            # dist <- c((1:round(allele_count / ploidy))^(-as.numeric(one_parameter[["alpha"]])), rep(0, allele_count - round(allele_count / ploidy)))
            dist <- dist / sum(dist)
            SFS_exact_components <- matrix(dist, nrow = 1)
        } else {
            SFS_exact_components <- matrix(NA, nrow = 0, ncol = allele_count)
        }
        for (i in seq_len(N_clusters)) {
            dist <- dbinom(1:allele_count, size = allele_count, prob = as.numeric(one_parameter[paste0("p_", i)])) / (1 - dbinom(0, size = allele_count, prob = as.numeric(one_parameter[paste0("p_", i)])))
            SFS_exact_components <- rbind(SFS_exact_components, dist)
        }
        #   Compute expected SFS for inference A
        SFS_A_components <- SFS_exact_components %*% SFS_convolution_inference_A$convolution_matrix
        if (with_tail) {
            omegas_A <- as.numeric(one_parameter[paste0("omega_inference_A_", 0:N_clusters)])
        } else {
            omegas_A <- as.numeric(one_parameter[paste0("omega_inference_A_", 1:N_clusters)])
        }
        omegas_exact <- omegas_A / rowSums(SFS_A_components)
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
        #   Prepare output data
        output <- list()
        if (output_distance) {
            dA <- sum(abs(SFS_A - SFS_data_inference_A)) / sum(abs(SFS_data_inference_A))
            dB <- sum(abs(SFS_B - SFS_data_inference_B)) / sum(abs(SFS_data_inference_B))
            output$distance <- dA + dB
        }
        if (output_SFS) {
            output$SFS_A <- SFS_A
            output$SFS_B <- SFS_B
            if (output_validation) output$SFS_validation <- SFS_validation
        }
        if (output_SFS_components) {
            rownames(SFS_A_components) <- c(if (with_tail) "Tail" else NULL, if (N_clusters > 0) paste0("Cluster_", 1:N_clusters) else NULL)
            rownames(SFS_B_components) <- c(if (with_tail) "Tail" else NULL, if (N_clusters > 0) paste0("Cluster_", 1:N_clusters) else NULL)
            if (output_validation) rownames(SFS_validation_components) <- c(if (with_tail) "Tail" else NULL, if (N_clusters > 0) paste0("Cluster_", 1:N_clusters) else NULL)
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
            for (i in seq_len(N_clusters)) {
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
    smcrf_result <- smcrf_multi_param(
        statistics_target = data.frame(distance = 0),
        model = model,
        rprior = rprior,
        dprior = dprior,
        rperturb = rperturb,
        dperturb = dperturb,
        nParticles = n_SMCRF_particles,
        parallel = compute_parallel,
        min.node.size = 1
    )
    #---Sample parameter sets from ABC-SMC-DRF posterior distribution
    final_parameters <- smcrf_result[[paste0("Iteration_", smcrf_result[["nIterations"]])]]$parameters
    final_weights <- smcrf_result[[paste0("Iteration_", smcrf_result[["nIterations"]])]]$weights
    posterior_parameters <- final_parameters[
        sample(nrow(final_parameters), size = criterion_Nsamples, prob = final_weights[, 1], replace = T),
    ]
    #---Prepare results based on posterior parameters
    posterior_results_list <- lapply(seq_len(nrow(posterior_parameters)), function(i) {
        one_SFS(
            posterior_parameters[i, ],
            output_distance = TRUE,
            output_validation = TRUE,
            output_SFS = TRUE,
            output_SFS_components = TRUE,
            output_all_parameters = TRUE
        )
    })
    distances <- sapply(posterior_results_list, function(x) x$distance)
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
    for (i in seq_len(N_clusters)) {
        SFS_A_clusters[[i]] <- as.data.frame(do.call(rbind, lapply(posterior_results_list, function(x) x$SFS_A_components[paste0("Cluster_", i), ])))
        SFS_B_clusters[[i]] <- as.data.frame(do.call(rbind, lapply(posterior_results_list, function(x) x$SFS_B_components[paste0("Cluster_", i), ])))
        SFS_validation_clusters[[i]] <- as.data.frame(do.call(rbind, lapply(posterior_results_list, function(x) x$SFS_validation_components[paste0("Cluster_", i), ])))
    }
    #---Compute criterion values
    nParameters <- if (with_tail) 2 * N_clusters + 2 * criterion_tail_weight else 2 * N_clusters
    nData <- sum(SFS_data_validation)
    if (criterion == "BIC") {
        SFS_validation_normalized <- SFS_validation
        SFS_validation_normalized[SFS_validation_normalized <= .Machine$double.eps] <- .Machine$double.eps
        SFS_validation_normalized <- SFS_validation_normalized / rowSums(SFS_validation)
        loglikelihood <- apply(SFS_validation_normalized, 1, function(SFS_pred_validation) sum(log(SFS_pred_validation) * SFS_data_validation))
        criterion_values <- data.frame(
            criterion_value = -2 * loglikelihood + criterion_penalty_scale * nParameters * log(nData)
        )
    } else if (criterion == "GIC") {
        err_L1 <- apply(SFS_validation, 1, function(SFS_pred_validation) sum(abs(SFS_pred_validation - SFS_data_validation)))
        criterion_values <- data.frame(
            L1_error = err_L1 / sum(SFS_data_validation),
            penalty = criterion_penalty_scale * nParameters * log(nData)
        )
        criterion_values$criterion_value <- criterion_values$L1_error + criterion_values$penalty
    }
    #---Output results
    fit_results <- list()
    fit_results[["N_clusters"]] <- N_clusters
    fit_results[["with_tail"]] <- with_tail
    fit_results[["SFS_inference_A"]] <- SFS_A
    fit_results[["SFS_inference_B"]] <- SFS_B
    fit_results[["SFS_validation"]] <- SFS_validation
    fit_results[["parameters"]] <- all_parameters
    fit_results[["posterior_parameters"]] <- posterior_parameters
    fit_results[["distances"]] <- distances
    if (with_tail) {
        fit_results[["SFS_inference_A_tail"]] <- SFS_A_tail
        fit_results[["SFS_inference_B_tail"]] <- SFS_B_tail
        fit_results[["SFS_validation_tail"]] <- SFS_validation_tail
    }
    for (i in seq_len(N_clusters)) {
        fit_results[[paste0("SFS_inference_A_cluster_", i)]] <- SFS_A_clusters[[i]]
        fit_results[[paste0("SFS_inference_B_cluster_", i)]] <- SFS_B_clusters[[i]]
        fit_results[[paste0("SFS_validation_cluster_", i)]] <- SFS_validation_clusters[[i]]
    }
    fit_results[["criterion_values"]] <- criterion_values
    return(fit_results)
}

choose_mutation_thresholds <- function(mutation_table,
                                       max_total_read,
                                       min_variant_read_inference_A,
                                       min_variant_read_inference_B,
                                       min_variant_read_validation,
                                       read_distribution_freq_min,
                                       SFS_data_frequencies,
                                       inference_retained_freq) {
    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(data.table))
    cat(bold(blue("Choose mutation thresholds for inference and validation...\n")))
    #---Find variant readcount thresholds for inference and validation
    variant_read_cdf <- mutation_table %>%
        group_by(Alt_count) %>%
        summarise(count = n(), .groups = "drop") %>%
        ungroup() %>%
        arrange(Alt_count) %>%
        mutate(freq = rev(cumsum(rev(count))) / nrow(mutation_table) * 100)
    variant_read_limit_min <- min(variant_read_cdf$Alt_count[variant_read_cdf$freq < 100])
    variant_read_limit_max <- max(variant_read_cdf$Alt_count[variant_read_cdf$freq > inference_retained_freq])
    variant_read_limit_gap <- max(1, round((variant_read_limit_max - variant_read_limit_min) / 4))
    if (is.null(min_variant_read_inference_A)) min_variant_read_inference_A <- variant_read_limit_min + 1 * variant_read_limit_gap
    if (is.null(min_variant_read_inference_B)) min_variant_read_inference_B <- variant_read_limit_min + 3 * variant_read_limit_gap
    if (is.null(min_variant_read_validation)) min_variant_read_validation <- min_variant_read_inference_A + 1
    #---Find total readcount thresholds for inference and validation
    choose_min_total_read <- function(min_variant_read) {
        total_read_cdf <- mutation_table %>%
            filter(Alt_count >= min_variant_read) %>%
            group_by(Tot_count) %>%
            summarise(count = n(), .groups = "drop") %>%
            ungroup() %>%
            arrange(Tot_count) %>%
            mutate(cumulated_count = rev(cumsum(rev(count)))) %>%
            mutate(freq = cumulated_count / nrow(mutation_table) * 100)
        if (max(total_read_cdf$freq) < inference_retained_freq) {
            min_total_read <- min(total_read_cdf$Tot_count)
        } else {
            min_total_read <- max(total_read_cdf$Tot_count[total_read_cdf$freq >= inference_retained_freq])
        }
        Nmut <- total_read_cdf$cumulated_count[total_read_cdf$Tot_count == min_total_read]
        freq <- total_read_cdf$freq[total_read_cdf$Tot_count == min_total_read]
        return(list(min_total_read = min_total_read, total_read_cdf = total_read_cdf, Nmut = Nmut, freq = freq))
    }
    tmp <- choose_min_total_read(min_variant_read = min_variant_read_inference_A)
    min_total_read_inference_A <- tmp$min_total_read
    Nmut_inference_A <- tmp$Nmut
    freq_inference_A <- tmp$freq
    tmp <- choose_min_total_read(min_variant_read = min_variant_read_inference_B)
    min_total_read_inference_B <- tmp$min_total_read
    Nmut_inference_B <- tmp$Nmut
    freq_inference_B <- tmp$freq
    tmp <- choose_min_total_read(min_variant_read = min_variant_read_validation)
    min_total_read_validation <- min(tmp$total_read_cdf$Tot_count)
    Nmut_validation <- tmp$total_read_cdf$cumulated_count[tmp$total_read_cdf$Tot_count == min_total_read_validation]
    freq_validation <- tmp$total_read_cdf$freq[tmp$total_read_cdf$Tot_count == min_total_read_validation]
    #---Find joint distribution of variant and total readcounts
    vec_variant_read <- min(mutation_table$Alt_count):max(variant_read_cdf$Alt_count[variant_read_cdf$freq > read_distribution_freq_min])
    total_read_cdf <- choose_min_total_read(min_variant_read = 0)$total_read_cdf
    vec_total_read <- min(mutation_table$Tot_count):min(max_total_read, max(total_read_cdf$Tot_count[total_read_cdf$freq > read_distribution_freq_min]))
    mutation_table_new <- mutation_table %>%
        group_by(Alt_count, Tot_count) %>%
        summarise(count = n(), .groups = "drop") %>%
        ungroup()
    readcount_distribution <- expand.grid(min_variant_read = vec_variant_read, min_total_read = vec_total_read, mutation_count = 0)
    pb <- txtProgressBar(
        min = 1,
        max = nrow(mutation_table_new),
        style = 3,
        width = 50,
        char = "+"
    )
    for (row in 1:nrow(mutation_table_new)) {
        setTxtProgressBar(pb, row)
        min_variant_read <- mutation_table_new$Alt_count[row]
        min_total_read <- mutation_table_new$Tot_count[row]
        rows <- which(readcount_distribution$min_variant_read <= min_variant_read & readcount_distribution$min_total_read <= min_total_read)
        readcount_distribution$mutation_count[rows] <- readcount_distribution$mutation_count[rows] + mutation_table_new$count[row]
    }
    readcount_distribution$freq <- 100 * readcount_distribution$mutation_count / sum(mutation_table_new$count)
    cat("\n")
    #   Report the chosen mutation thresholds
    report <- paste0(blue("Complete data          : "), cyan(paste0(min(mutation_table$Alt_count), " \u2264 variant reads, ", min(mutation_table$Tot_count), " \u2264 total reads \u2264 ", max(mutation_table$Tot_count), "; ", nrow(mutation_table), " mutations\n")))
    report <- paste0(report, blue("Inference A            : "), cyan(paste0(min_variant_read_inference_A, " \u2264 variant reads, ", min_total_read_inference_A, " \u2264 total reads \u2264 ", max_total_read, "; ", Nmut_inference_A, " mutations (", format(round(freq_inference_A, 2), nsmall = 2), "%)")), "\n")
    report <- paste0(report, blue("Inference B            : "), cyan(paste0(min_variant_read_inference_B, " \u2264 variant reads, ", min_total_read_inference_B, " \u2264 total reads \u2264 ", max_total_read, "; ", Nmut_inference_B, " mutations (", format(round(freq_inference_B, 2), nsmall = 2), "%)")), "\n")
    report <- paste0(report, blue("Validation             : "), cyan(paste0(min_variant_read_validation, " \u2264 variant reads, ", min_total_read_validation, " \u2264 total reads \u2264 ", max_total_read, "; ", Nmut_validation, " mutations (", format(round(freq_validation, 2), nsmall = 2), "%)")), "\n")
    cat(report)
    #---Get coverage distribution
    mutation_table <- as.data.table(mutation_table)
    sample_coverage <- mutation_table[, .(pdf = .N), by = .(total_readcount = Tot_count)]
    sample_coverage[, pdf := pdf / sum(pdf)]
    sample_coverage_inference_A <- sample_coverage[which(sample_coverage$total_readcount >= min_total_read_inference_A & sample_coverage$total_readcount <= max_total_read), ]
    sample_coverage_inference_B <- sample_coverage[which(sample_coverage$total_readcount >= min_total_read_inference_B & sample_coverage$total_readcount <= max_total_read), ]
    sample_coverage_validation <- sample_coverage[which(sample_coverage$total_readcount >= min_total_read_validation & sample_coverage$total_readcount <= max_total_read), ]
    #---Prepare the real SFS data for the inference sets
    func_SFS <- function(mutation_table) {
        VAF <- mutation_table$VAF
        lower_bound <- c(0, SFS_data_frequencies[-length(SFS_data_frequencies)])
        upper_bound <- SFS_data_frequencies
        SFS_data <- numeric(length(SFS_data_frequencies))
        for (i in seq_along(SFS_data_frequencies)) SFS_data[i] <- sum(VAF > lower_bound[i] & VAF <= upper_bound[i])
        return(SFS_data)
    }
    SFS_data_inference_A <- func_SFS(mutation_table[Alt_count >= min_variant_read_inference_A & Tot_count >= min_total_read_inference_A & Tot_count <= max_total_read])
    SFS_data_inference_B <- func_SFS(mutation_table[Alt_count >= min_variant_read_inference_B & Tot_count >= min_total_read_inference_B & Tot_count <= max_total_read])
    SFS_data_validation <- func_SFS(mutation_table[Alt_count >= min_variant_read_validation & Tot_count >= min_total_read_validation & Tot_count <= max_total_read])
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

build_convolution_matrix <- function(sfs_bincount,
                                     mode = NULL,
                                     allele_count,
                                     min_variant_read,
                                     min_total_read,
                                     max_total_read,
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
                                        sample_coverage_distribution,
                                        sfs_bincount,
                                        i) {
        r <- sample_coverage_distribution$total_readcount
        r1 <- r * (i - 1) / SFS_totalsteps_base
        r1 <- ifelse(r1 %% 1 == 0, r1 + 1, ceiling(r1))
        r1 <- pmax(r1, min_variant_read, 1)
        r2 <- floor(r * i / SFS_totalsteps_base)
        prob <- (1:N_end) / N_end
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
            "sample_coverage_distribution", "sfs_bincount", "func_convolution_vector"
        ), envir = environment())
        #   Compute each convolution vector
        output <- pblapply(cl = cl, X = 1:sfs_bincount, FUN = function(i) {
            return(func_convolution_vector(
                N_end = N_end,
                SFS_totalsteps_base = SFS_totalsteps_base,
                min_variant_read = min_variant_read,
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

smcrf_multi_param <- function(statistics_target,
                              model,
                              rprior,
                              dprior,
                              rperturb,
                              dperturb,
                              nParticles,
                              final_sample = FALSE,
                              model_redo_if_NA = TRUE,
                              verbose = FALSE,
                              honesty = FALSE,
                              parallel,
                              save_model = FALSE,
                              save_rds = FALSE,
                              filename_rds,
                              splitting.rule = "CART",
                              smcrf_multi_param_results = NULL,
                              ...) {
    suppressPackageStartupMessages(library(drf))
    #---Obtain information from previous SMC-DRF
    if (!is.null(smcrf_multi_param_results)) {
        #   ... If continuing from a previous ABC-SMC-DRF chain:
        SMCDRF <- smcrf_multi_param_results
        parameters_ids <- colnames(SMCDRF[["Iteration_1"]]$parameters)
        statistics_ids <- colnames(SMCDRF[["Iteration_1"]]$statistics)
        nIterations <- length(nParticles) + SMCDRF[["nIterations"]]
        iteration_start <- 1 + SMCDRF[["nIterations"]]
        if (is.null(statistics_target)) statistics_target <- SMCDRF[["statistics_target"]]
        nParticles <- c(SMCDRF[["nParticles"]], nParticles)
        SMCDRF[[paste0("Iteration_", iteration_start)]] <- c()
        DRF_weights <- SMCDRF[[paste0("Iteration_", iteration_start - 1)]]$weights
        parameters <- SMCDRF[[paste0("Iteration_", iteration_start - 1)]]$parameters
    } else {
        #   ... If starting a new ABC-SMC-DRF chain:
        SMCDRF <- list()
        parameters_ids <- colnames(rprior(Nparameters = 1))
        statistics_ids <- colnames(statistics_target)
        nIterations <- length(nParticles)
        iteration_start <- 1
    }
    reference_ids <- colnames(model(parameters = rprior(Nparameters = 1)))
    iteration_end <- ifelse(final_sample, nIterations + 1, nIterations)
    SMCDRF[["method"]] <- "smcrf-multi-param"
    SMCDRF[["nIterations"]] <- nIterations
    SMCDRF[["nParticles"]] <- nParticles
    SMCDRF[["statistics_target"]] <- statistics_target
    SMCDRF[["parameters_labels"]] <- data.frame(parameter = parameters_ids)
    SMCDRF[["statistics_labels"]] <- data.frame(ID = statistics_ids)
    if (!verbose) {
        pb <- txtProgressBar(
            min = 0,
            max = iteration_end - iteration_start + 1,
            style = 3,
            width = 50,
            char = "+"
        )
    }
    for (iteration in iteration_start:iteration_end) {
        if (verbose) {
            if (iteration == (nIterations + 1)) {
                cat(bold(red("ABC-SMC-DRF FOR MULTIPLE PARAMETERS:")), paste0(bold(yellow("final posterior distribution", "\n"))))
            } else {
                cat(bold(red("ABC-SMC-DRF FOR MULTIPLE PARAMETERS:")), paste0(bold(yellow("iteration", iteration, "\n"))))
            }
        }
        #---Prepare sampled particles from the previous iteration
        if (iteration > 1) {
            parameters_previous_sampled <- as.data.frame(parameters[sample(nrow(parameters), size = 10000, prob = DRF_weights[, 1], replace = TRUE), , drop = FALSE])
        }
        #---Create training set
        if (verbose) cat(blue("Sampling parameters and computing model simulations...\n"))
        nrow <- ifelse(iteration == (nIterations + 1), nParticles[nIterations], nParticles[iteration])
        invalid_indices <- 1:nrow
        parameters_unperturbed <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids)))
        colnames(parameters_unperturbed) <- parameters_ids
        parameters_next <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids)))
        colnames(parameters_next) <- parameters_ids
        reference_next <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids) + ncol(statistics_target)))
        colnames(reference_next) <- reference_ids
        while (length(invalid_indices) > 0) {
            #   Generate new particles...
            if (iteration == 1) {
                #   ... for iteration 1:
                #   Sample from the prior distribution
                parameters_unperturbed[invalid_indices, ] <- rprior(Nparameters = length(invalid_indices))
                parameters_next[invalid_indices, ] <- parameters_unperturbed[invalid_indices, ]
            } else {
                #   ... for later iterations:
                #   Sample from the previous posterior distribution
                parameters_previous <- parameters
                weights_previous <- DRF_weights
                parameter_replace <- parameters[sample(nrow(parameters), size = length(invalid_indices), prob = DRF_weights[, 1], replace = T), ]
                parameters_unperturbed[invalid_indices, ] <- parameter_replace
                #   Perturb parameters
                if (iteration < (nIterations + 1)) {
                    parameter_replace <- rperturb(
                        parameters_unperturbed = parameter_replace,
                        parameters_previous_sampled = parameters_previous_sampled,
                        iteration = iteration
                    )
                }
                parameters_next[invalid_indices, ] <- parameter_replace
            }
            invalid_indices_next <- which(dprior(parameters_next, parameter_id = "all") <= 0)
            #   Generate statistics for the new particles
            ids <- setdiff(invalid_indices, invalid_indices_next)
            tmp <- parameters_next[ids, ]
            if (!is.data.frame(tmp)) {
                tmp <- as.data.frame(tmp, ncol = length(parameters_ids))
                colnames(tmp) <- parameters_ids
            }
            if (length(ids) > 0) reference_next[ids, ] <- model(parameters = tmp)
            #   Find particles that need to be regenerated
            invalid_indices <- invalid_indices_next
            if (model_redo_if_NA) invalid_indices <- unique(c(invalid_indices, which(apply(reference_next, 1, function(x) any(is.na(x))))))
        }
        reference <- reference_next
        parameters <- reference[, parameters_ids, drop = FALSE]
        statistics <- reference[, statistics_ids, drop = FALSE]
        #   Finish the last iteration
        if (iteration == (nIterations + 1)) {
            SMCDRF_iteration <- list()
            SMCDRF_iteration$reference <- reference
            SMCDRF_iteration$parameters <- parameters
            SMCDRF_iteration$parameters_unperturbed <- parameters_unperturbed
            SMCDRF_iteration$statistics <- statistics
            SMCDRF[[paste0("Iteration_", iteration)]] <- SMCDRF_iteration
            if (save_rds == TRUE) {
                saveRDS(SMCDRF, file = filename_rds)
            }
            if (!verbose) setTxtProgressBar(pb, iteration - iteration_start + 1)
            if (!verbose) cat("\n")
            return(SMCDRF)
        }
        #---Run DRF for all parameters
        if (verbose) cat(blue("Performing Random Forest prediction...\n"))
        drfmodel <- drf(
            X = statistics,
            Y = parameters,
            splitting.rule = splitting.rule,
            honesty = honesty,
            ...
        )
        def_pred <- predict(
            object = drfmodel,
            newdata = statistics_target
        )
        DRF_weights <- as.vector(get_sample_weights(drfmodel, statistics_target))
        #---Modify DRF weights
        if (verbose) cat(blue("Recalibrating Random Forest weights...\n"))
        if (iteration > 1) {
            #   Compute numerators for weight recalibration
            weight_modifiers_numerator <- dprior(parameters, parameter_id = "all")
            #   Compute denominators for weight recalibration
            weight_modifiers_denominator <- rep(0, nrow(parameters))
            for (i in 1:nrow(parameters_previous)) {
                weight_modifiers_denominator_i <- rep(weights_previous[i, 1], nrow(parameters)) *
                    dperturb(
                        parameters = parameters,
                        parameters_previous = parameters_previous[i, , drop = FALSE],
                        parameters_previous_sampled = parameters_previous_sampled,
                        iteration = iteration,
                        parameter_id = "all"
                    )
                weight_modifiers_denominator <- weight_modifiers_denominator + weight_modifiers_denominator_i
            }
            #   Modify weights for new particles
            DRF_weights <- DRF_weights * weight_modifiers_numerator / weight_modifiers_denominator
            DRF_weights <- DRF_weights / sum(DRF_weights)
        }
        DRF_weights <- data.frame(matrix(rep(DRF_weights, length(parameters_ids)), ncol = length(parameters_ids)))
        colnames(DRF_weights) <- parameters_ids
        ################################################################
        ################################################################
        ################################################################
        #---Order parameter sets in decreasing cluster VAFs
        N_clusters <- max(as.numeric(sub("p_", "", grep("^p_", colnames(parameters), value = TRUE))))
        if (N_clusters > 1) {
            p_colnames <- paste0("p_", 1:N_clusters)
            omega_colnames <- paste0("omega_inference_A_", 1:N_clusters)
            reorder_clusters <- function(row) {
                p_values <- as.numeric(row[p_colnames])
                omega_values <- as.numeric(row[omega_colnames])
                order_idx <- order(p_values, decreasing = TRUE)
                new_row <- row
                new_row[p_colnames] <- p_values[order_idx]
                new_row[omega_colnames] <- omega_values[order_idx]
                return(new_row)
            }
            parameters <- as.data.frame(t(apply(parameters, 1, reorder_clusters)))
            colnames(parameters) <- colnames(parameters)
        }
        ################################################################
        ################################################################
        ################################################################
        #---Save SMC-DRF results from this iteration
        SMCDRF_iteration <- list()
        SMCDRF_iteration$reference <- reference
        SMCDRF_iteration$parameters <- parameters
        SMCDRF_iteration$parameters_unperturbed <- parameters_unperturbed
        SMCDRF_iteration$statistics <- statistics
        SMCDRF_iteration$weights <- DRF_weights
        if (save_model == TRUE) {
            SMCDRF_iteration$rf_model <- drfmodel
            SMCDRF_iteration$rf_predict <- def_pred
        }
        SMCDRF[[paste0("Iteration_", iteration)]] <- SMCDRF_iteration
        if (save_rds == TRUE) {
            saveRDS(SMCDRF, file = filename_rds)
        }
        if (!verbose) setTxtProgressBar(pb, iteration - iteration_start + 1)
    }
    if (save_rds == TRUE) {
        saveRDS(SMCDRF, file = filename_rds)
    }
    if (!verbose) cat("\n")
    return(SMCDRF)
}
