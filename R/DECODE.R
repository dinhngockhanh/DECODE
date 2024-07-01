DECODE <- function(mutation_table,
                   criterion = "BIC",
                   criterion_ratio = 1,
                   neutral_power_min = 0.5,
                   neutral_power_max = 5,
                   cluster_frequency_min = 0.01,
                   cluster_frequency_max = 1,
                   matrix_binomial_PDF,
                   matrix_binomial_sample_size,
                   matrix_binomial_sfs_stepcount,
                   matrix_binomial_ploidy,
                   sample_size,
                   SFS_totalsteps,
                   r_min,
                   r_max,
                   coverage_distribution,
                   coverage_variables = NULL,
                   N_trials = 10000,
                   N_trials_Morris_tail_sensitivity = 300,
                   neutral_tail = NA,
                   min_N_humps = 0,
                   max_N_humps = Inf,
                   pi_cutoff = 0.02,
                   zero_cutoff = 1e-50,
                   compute_parallel_library = TRUE,
                   compute_parallel_fit = TRUE,
                   n_cores = NULL,
                   parameter_filename = NULL) {
    mutation_refcounts <- mutation_table$Ref_count
    mutation_altcounts <- mutation_table$Alt_count
    mutation_totcounts <- mutation_refcounts + mutation_altcounts
    #---Prepare the total readcount distribution
    cat("Prepare the total readcount distribution...\n")
    sample_coverage <- prep_distribution_patient(mutation_totcounts)
    #---Prepare the real SFS
    cat("Prepare the real SFS...\n")
    no_mutations_total <- length(mutation_refcounts)
    vec_freq <- seq(1, SFS_totalsteps) / SFS_totalsteps
    vec_SFS_real <- rep(0, SFS_totalsteps)
    for (j in 1:no_mutations_total) {
        no_variant <- mutation_altcounts[j]
        no_total <- mutation_refcounts[j] + mutation_altcounts[j]
        if (no_variant >= min_variant_read && no_total >= min_total_read) {
            VAF <- no_variant / no_total
            pos <- which(vec_freq >= VAF)[1]
            vec_SFS_real[pos] <- vec_SFS_real[pos] + 1
        }
    }
    #---Prepare the SFS convolution matrix
    cat("Prepare the SFS convolution matrix...\n")
    SFS_convolution_matrix <- build_convolution_matrix(
        N_end = matrix_binomial_sample_size,
        SFS_totalsteps = SFS_totalsteps,
        SFS_totalsteps_base = matrix_binomial_sfs_stepcount,
        r_min = r_min,
        r_max = r_max,
        coverage_distribution = coverage_distribution,
        coverage_variables = coverage_variables,
        sample_coverage = sample_coverage,
        compute_parallel = compute_parallel_library,
        n_cores = n_cores
    )
    #---DECODE
    DECODE_result <- list()
    DECODE_result$mutational_table <- mutation_table
    DECODE_result$SFS_frequencies <- vec_freq
    DECODE_result$SFS_for_fitting <- vec_SFS_real
    if (is.na(neutral_tail)) {
        result_with_tail <- DECODE_given_tail_status(
            vec_SFS_real = vec_SFS_real,
            criterion = criterion,
            criterion_ratio = criterion_ratio,
            min_N_humps = min_N_humps,
            max_N_humps = max_N_humps,
            with_tail = TRUE,
            N_trials = N_trials,
            N_trials_Morris_tail_sensitivity = N_trials_Morris_tail_sensitivity,
            SFS_totalsteps = SFS_totalsteps,
            SFS_convolution_matrix = SFS_convolution_matrix,
            neutral_power_min = neutral_power_min,
            neutral_power_max = neutral_power_max,
            cluster_frequency_min = cluster_frequency_min,
            cluster_frequency_max = cluster_frequency_max,
            pi_cutoff = pi_cutoff,
            zero_cutoff = zero_cutoff,
            compute_parallel_fit = compute_parallel_fit,
            n_cores = n_cores
        )
        DECODE_result$with_tail <- result_with_tail
        result_without_tail <- DECODE_given_tail_status(
            vec_SFS_real = vec_SFS_real,
            criterion = criterion,
            criterion_ratio = criterion_ratio,
            min_N_humps = min_N_humps,
            max_N_humps = max_N_humps,
            with_tail = FALSE,
            N_trials = N_trials,
            N_trials_Morris_tail_sensitivity = N_trials_Morris_tail_sensitivity,
            SFS_totalsteps = SFS_totalsteps,
            SFS_convolution_matrix = SFS_convolution_matrix,
            neutral_power_min = neutral_power_min,
            neutral_power_max = neutral_power_max,
            cluster_frequency_min = cluster_frequency_min,
            cluster_frequency_max = cluster_frequency_max,
            pi_cutoff = pi_cutoff,
            zero_cutoff = zero_cutoff,
            compute_parallel_fit = compute_parallel_fit,
            n_cores = n_cores
        )
        DECODE_result$without_tail <- result_without_tail
        if (result_with_tail$best_fit$selected_criterion_value < result_without_tail$best_fit$selected_criterion_value) {
            final_result <- result_with_tail
        } else {
            final_result <- result_without_tail
        }
    } else if (neutral_tail == TRUE) {
        final_result <- DECODE_given_tail_status(
            vec_SFS_real = vec_SFS_real,
            criterion = criterion,
            criterion_ratio = criterion_ratio,
            min_N_humps = min_N_humps,
            max_N_humps = max_N_humps,
            with_tail = TRUE,
            N_trials = N_trials,
            N_trials_Morris_tail_sensitivity = N_trials_Morris_tail_sensitivity,
            SFS_totalsteps = SFS_totalsteps,
            SFS_convolution_matrix = SFS_convolution_matrix,
            neutral_power_min = neutral_power_min,
            neutral_power_max = neutral_power_max,
            cluster_frequency_min = cluster_frequency_min,
            cluster_frequency_max = cluster_frequency_max,
            pi_cutoff = pi_cutoff,
            zero_cutoff = zero_cutoff,
            compute_parallel_fit = compute_parallel_fit,
            n_cores = n_cores
        )
        DECODE_result$with_tail <- final_result
    } else if (neutral_tail == FALSE) {
        final_result <- DECODE_given_tail_status(
            vec_SFS_real = vec_SFS_real,
            criterion = criterion,
            criterion_ratio = criterion_ratio,
            min_N_humps = min_N_humps,
            max_N_humps = max_N_humps,
            with_tail = FALSE,
            N_trials = N_trials,
            N_trials_Morris_tail_sensitivity = N_trials_Morris_tail_sensitivity,
            SFS_totalsteps = SFS_totalsteps,
            SFS_convolution_matrix = SFS_convolution_matrix,
            neutral_power_min = neutral_power_min,
            neutral_power_max = neutral_power_max,
            cluster_frequency_min = cluster_frequency_min,
            cluster_frequency_max = cluster_frequency_max,
            pi_cutoff = pi_cutoff,
            zero_cutoff = zero_cutoff,
            compute_parallel_fit = compute_parallel_fit,
            n_cores = n_cores
        )
        DECODE_result$without_tail <- final_result
    }
    DECODE_result$best_result <- final_result
    #---Report the best fit
    tail_status_final_result <- final_result$best_fit$tail_status
    parameters_final_result <- final_result$best_fit$parameters
    criterion_final_result <- final_result$best_fit$selected_criterion_value
    if (tail_status_final_result) {
        N_humps_final_result <- length(parameters_final_result) / 2 - 1
        report <- paste0("\n\n\nBEST FIT = neutral tail + ", N_humps_final_result, " clusters: ", criterion, " = ", round(criterion_final_result, 3), "; neutral: pi = ", round(parameters_final_result[1], 3), " with power = ", round(parameters_final_result[2], 3))
        ii <- 0
    } else {
        N_humps_final_result <- length(parameters_final_result) / 2
        report <- paste0("\n\n\nBEST FIT = no neutral tail + ", N_humps_final_result, " clusters: ", criterion, " = ", round(criterion_final_result, 3))
        ii <- -1
    }
    if (N_humps_final_result > 0) {
        for (i in 1:N_humps_final_result) {
            report <- paste0(report, "; pi = ", round(parameters_final_result[2 * (i + ii) + 1], 3), " at freq = ", round(parameters_final_result[2 * (i + ii) + 2], 3))
        }
    }
    report <- paste0(report, "\n\n\n\n")
    cat(report)
    #---Translation to parameters of cancer evolution in the sample
    tmp <- parameter_conversion(
        result = final_result,
        mutation_count_for_fitting = sum(vec_SFS_real),
        sample_size = sample_size,
        matrix_binomial_sample_size = matrix_binomial_sample_size,
        matrix_binomial_ploidy = matrix_binomial_ploidy
    )
    best_fit_parameters <- tmp$parameters_df
    if (!is.null(parameter_filename)) {
        write.table(best_fit_parameters, parameter_filename, sep = "\t", quote = FALSE, row.names = FALSE)
    }
    #---Return the SFS deconvolution results
    DECODE_result$best_result$parameters_df <- best_fit_parameters
    return(DECODE_result)
}

DECODE_given_tail_status <- function(vec_SFS_real,
                                     criterion,
                                     criterion_ratio,
                                     min_N_humps,
                                     max_N_humps,
                                     with_tail,
                                     N_trials,
                                     N_trials_Morris_tail_sensitivity,
                                     SFS_totalsteps,
                                     SFS_convolution_matrix,
                                     neutral_power_min,
                                     neutral_power_max,
                                     cluster_frequency_min,
                                     cluster_frequency_max,
                                     pi_cutoff, # <<<<<<<<<<<<<<<<<<<<<<
                                     zero_cutoff,
                                     compute_parallel_fit,
                                     n_cores) {
    N_humps <- min_N_humps
    criterion_best_final <- Inf
    N_fitting_rounds <- 200
    mutation_count <- sum(vec_SFS_real)
    all_fits <- list()
    while (TRUE) {
        #---Find best parameter set, given the number of humps
        if (with_tail) {
            cat(paste0("Inference for ", N_humps, " clusters with neutral tail component...\n"))
        } else {
            cat(paste0("Inference for ", N_humps, " clusters without neutral tail component...\n"))
        }
        fit_results <- DECODE_given_tail_status_and_Ncluster(
            vec_SFS_real = vec_SFS_real,
            N_humps = N_humps,
            with_tail = with_tail,
            N_trials = N_trials,
            SFS_totalsteps = SFS_totalsteps,
            SFS_convolution_matrix = SFS_convolution_matrix,
            neutral_power_min = neutral_power_min,
            neutral_power_max = neutral_power_max,
            cluster_frequency_min = cluster_frequency_min,
            cluster_frequency_max = cluster_frequency_max,
            zero_cutoff = zero_cutoff,
            compute_parallel = compute_parallel_fit,
            n_cores = n_cores
        )
        all_fits[[paste0(N_humps, "_clusters")]] <- fit_results
        vec_para_best_current <- fit_results$best_parameters
        component_distributions_best_current <- fit_results$best_component_distributions
        criterion_all_best_current <- data.frame(
            AIC = fit_results$best_AIC,
            BIC = fit_results$best_BIC,
            ICL = fit_results$best_ICL,
            ICL_MAP = fit_results$best_ICL_MAP
        )
        if (criterion == "AIC") {
            criterion_best_current <- fit_results$best_AIC
        } else if (criterion == "BIC") {
            criterion_best_current <- fit_results$best_BIC
        } else if (criterion == "ICL") {
            criterion_best_current <- fit_results$best_ICL
        } else if (criterion == "ICL_MAP") {
            criterion_best_current <- fit_results$best_ICL_MAP
        }
        #   Report the best fit for the current hump count
        cluster_pis <- Inf
        report <- paste0(criterion, " = ", round(criterion_best_current, 3))
        if (with_tail) {
            report <- paste0(report, "; neutral: pi = ", round(vec_para_best_current[1], 3), " with power = ", round(vec_para_best_current[2], 3))
            ii <- 0
        } else {
            ii <- -1
        }
        if (N_humps > 0) {
            for (i in 1:N_humps) {
                report <- paste0(report, "; pi = ", round(vec_para_best_current[2 * (i + ii) + 1], 3), " at freq = ", round(vec_para_best_current[2 * (i + ii) + 2], 3))
                cluster_pis <- c(cluster_pis, vec_para_best_current[2 * (i + ii) + 1])
            }
        }
        report <- paste0(report, "\n")
        #   Check if the increased hump count leads to lower criterion score without tiny selective components...
        if ((criterion_best_current < criterion_ratio * criterion_best_final) & (min(cluster_pis) >= pi_cutoff)) {
            #   ... if yes, then update the best fit and continue with 1 more hump
            fit_results_best_final <- fit_results
            criterion_best_final <- criterion_best_current
            criterion_all_final <- criterion_all_best_current
            vec_para_best_final <- vec_para_best_current
            component_distributions_best_final <- component_distributions_best_current
            N_humps <- N_humps + 1
            #   ... except if exceeding maximum number of clusters
            if (N_humps > max_N_humps) break
        } else {
            #   ... if no, then stop
            break
        }
    }
    #---Check if the neutral tail component is too tiny
    if (with_tail == TRUE & vec_para_best_final[1] < pi_cutoff) {
        vec_para_best_final[seq(3, length(vec_para_best_final), by = 2)] <- vec_para_best_final[seq(3, length(vec_para_best_final), by = 2)] / sum(vec_para_best_final[seq(3, length(vec_para_best_final), by = 2)])
        vec_para_best_final <- vec_para_best_final[-c(1, 2)]

        component_distributions_best_final$SFS_exact[1, ] <- rep(0, length(component_distributions_best_final$SFS_exact[1, ]))
        component_distributions_best_final$SFS_expected[1, ] <- rep(0, length(component_distributions_best_final$SFS_expected[1, ]))
        component_distributions_best_final$SFS_expected_normalized[1, ] <- rep(0, length(component_distributions_best_final$SFS_expected_normalized[1, ]))

        with_tail <- FALSE
    }
    # ####################################################################
    # ####################################################################
    # ####################################################################
    # #---Sensitivity analysis for the tail parametrization
    tail_sensitivity_test <- FALSE
    # if (with_tail) {
    #     best_final_N_humps <- length(vec_para_best_final) / 2 - 1
    #     if (best_final_N_humps == 0) break
    #     cat(paste0("Sensitibity of tail parametrization with ", best_final_N_humps, " clusters...\n"))
    #     tail_sensitivity_test <- TRUE
    #     best_final_pi_0 <- vec_para_best_final[1]
    #     test_pi_0_min <- max(0, best_final_pi_0 - 0.1) # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    #     test_pi_0_max <- min(1, best_final_pi_0 + 0.1) # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    #     best_final_neutral_power <- vec_para_best_final[2]
    #     test_neutral_power_min <- max(neutral_power_min, best_final_neutral_power - 0.5) # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    #     test_neutral_power_max <- min(neutral_power_max, best_final_neutral_power + 0.5) # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    #     sensitivity_results <- DECODE_tail_parameter_sensitivity(
    #         vec_SFS_real = vec_SFS_real,
    #         N_humps = best_final_N_humps,
    #         N_Morris_trials = N_trials_Morris_tail_sensitivity,
    #         SFS_totalsteps = SFS_totalsteps,
    #         SFS_convolution_matrix = SFS_convolution_matrix,
    #         fit_results = fit_results_best_final,
    #         pi_0 = best_final_pi_0,
    #         pi_0_min = test_pi_0_min,
    #         pi_0_max = test_pi_0_max,
    #         neutral_power = best_final_neutral_power,
    #         neutral_power_min = test_neutral_power_min,
    #         neutral_power_max = test_neutral_power_max,
    #         cluster_frequency_min = cluster_frequency_min,
    #         cluster_frequency_max = cluster_frequency_max,
    #         zero_cutoff = zero_cutoff,
    #         compute_parallel = compute_parallel_fit,
    #         n_cores = n_cores
    #     )
    # }
    # ####################################################################
    # ####################################################################
    # ####################################################################
    #---Report the best fit
    result <- list()
    result$all_fits <- all_fits
    result$best_fit <- list()
    result$best_fit$parameters <- vec_para_best_final
    result$best_fit$component_distributions <- component_distributions_best_final
    result$best_fit$selected_criterion <- criterion
    result$best_fit$all_criteria <- criterion_all_final
    result$best_fit$selected_criterion_value <- criterion_best_final
    result$best_fit$tail_status <- with_tail
    if (tail_sensitivity_test) result$best_fit$tail_sensitivity <- sensitivity_results
    return(result)
}

DECODE_given_tail_status_and_Ncluster <- function(vec_SFS_real,
                                                  N_humps,
                                                  with_tail,
                                                  N_trials,
                                                  SFS_totalsteps,
                                                  SFS_convolution_matrix,
                                                  neutral_power = NA,
                                                  pi_0 = NA,
                                                  neutral_power_min,
                                                  neutral_power_max,
                                                  cluster_frequency_min,
                                                  cluster_frequency_max,
                                                  zero_cutoff,
                                                  compute_criteria = TRUE,
                                                  compute_parallel,
                                                  n_cores,
                                                  progress_bar = TRUE) {
    if (compute_criteria) {
        mutation_count <- sum(vec_SFS_real)
        if (with_tail) {
            num_parameters <- 2 * N_humps + 1
        } else {
            num_parameters <- 2 * N_humps - 1
        }
        if (neutral_power_max > neutral_power_min) num_parameters <- num_parameters + 1
    }
    #---Find best variable parameters (A & K's) for each fixed parameter set
    if (compute_parallel == FALSE) {
        all_logLikelihood <- c()
        all_para <- c()
        all_component_distributions <- list()
        if (compute_criteria) {
            all_AIC <- c()
            all_BIC <- c()
            all_ICL <- c()
            all_ICL_MAP <- c()
        }
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
            #   Sample neutral component power and cluster frequencies
            if (with_tail & is.na(neutral_power)) {
                neutral_power <- runif(1, neutral_power_min, neutral_power_max)
            } else if (!with_tail) {
                neutral_power <- NA
            }
            cluster_frequencies <- sort(runif(N_humps, cluster_frequency_min, cluster_frequency_max), decreasing = TRUE)
            #   Build the SFS component library
            component_distributions <- build_SFS_library(
                neutral_power = neutral_power,
                cluster_frequencies = cluster_frequencies,
                SFS_totalsteps = SFS_totalsteps,
                SFS_convolution_matrix = SFS_convolution_matrix
            )
            #   Find component proportions
            results <- DECODE_for_pis(
                vec_SFS_real = vec_SFS_real,
                neutral_power = neutral_power,
                pi_0 = pi_0,
                cluster_frequencies = cluster_frequencies,
                component_distributions = component_distributions,
                zero_cutoff = zero_cutoff
            )
            logLikelihood <- results$log_L
            vec_para <- results$parameters
            component_distributions <- results$component_distributions
            all_para <- rbind(all_para, vec_para)
            all_component_distributions[[i]] <- component_distributions
            all_logLikelihood <- c(all_logLikelihood, logLikelihood)
            if (compute_criteria) {
                criteria <- cluster_count_criteria(
                    num_parameters = num_parameters,
                    log_L = logLikelihood,
                    num_samples = mutation_count,
                    vec_SFS_real = vec_SFS_real,
                    vec_para = vec_para,
                    component_distributions = component_distributions
                )
                AIC <- criteria$AIC
                BIC <- criteria$BIC
                ICL <- criteria$ICL
                ICL_MAP <- criteria$ICL_MAP
                all_AIC <- c(all_AIC, AIC)
                all_BIC <- c(all_BIC, BIC)
                all_ICL <- c(all_ICL, ICL)
                all_ICL_MAP <- c(all_ICL_MAP, ICL_MAP)
            }
        }
        if (progress_bar) cat("\n")
    } else {
        library(parallel)
        library(pbapply)
        #   Start parallel cluster
        if (is.null(n_cores)) {
            numCores <- detectCores()
        } else {
            numCores <- n_cores
        }
        cl <- makePSOCKcluster(numCores - 1)
        #   Prepare input parameters
        vec_SFS_real <<- vec_SFS_real
        with_tail <<- with_tail
        N_humps <<- N_humps
        SFS_totalsteps <<- SFS_totalsteps
        N_end <<- N_end
        SFS_convolution_matrix <<- SFS_convolution_matrix
        pi_0 <<- pi_0
        neutral_power <<- neutral_power
        neutral_power_min <<- ifelse(missing(neutral_power_min), NA, neutral_power_min)
        neutral_power_max <<- ifelse(missing(neutral_power_max), NA, neutral_power_max)
        cluster_frequency_min <<- cluster_frequency_min
        cluster_frequency_max <<- cluster_frequency_max
        matrix_binomial_PDF <<- matrix_binomial_PDF
        zero_cutoff <<- zero_cutoff
        compute_criteria <<- compute_criteria
        clusterExport(cl, varlist = c(
            "zero_cutoff",
            "compute_criteria",
            "vec_SFS_real",
            "with_tail",
            "N_humps",
            "SFS_totalsteps",
            "N_end",
            "SFS_convolution_matrix",
            "pi_0",
            "neutral_power",
            "neutral_power_min",
            "neutral_power_max",
            "cluster_frequency_min",
            "cluster_frequency_max",
            "DECODE_for_pis",
            "build_SFS_library",
            "build_SFS_library_Griffiths_Tavare",
            "compute_loglikelihood",
            "compute_SFS",
            "cluster_count_criteria"
        ))
        #   Find best variable parameters in parallel mode
        func_parallel <- function(i) {
            #   Sample neutral component power and cluster frequencies
            if (with_tail & is.na(neutral_power)) {
                neutral_power <- runif(1, neutral_power_min, neutral_power_max)
            } else if (!with_tail) {
                neutral_power <- NA
            }
            cluster_frequencies <- sort(runif(N_humps, cluster_frequency_min, cluster_frequency_max), decreasing = TRUE)
            #   Build the SFS component library
            component_distributions <- build_SFS_library(
                neutral_power = neutral_power,
                cluster_frequencies = cluster_frequencies,
                SFS_totalsteps = SFS_totalsteps,
                SFS_convolution_matrix = SFS_convolution_matrix
            )
            #   Find component proportions
            results <- DECODE_for_pis(
                vec_SFS_real = vec_SFS_real,
                neutral_power = neutral_power,
                pi_0 = pi_0,
                cluster_frequencies = cluster_frequencies,
                component_distributions = component_distributions,
                zero_cutoff = zero_cutoff
            )
            logLikelihood <- results$log_L
            vec_para <- results$parameters
            component_distributions <- results$component_distributions
            if (compute_criteria) {
                criteria <- cluster_count_criteria(
                    num_parameters = num_parameters,
                    log_L = logLikelihood,
                    num_samples = mutation_count,
                    vec_SFS_real = vec_SFS_real,
                    vec_para = vec_para,
                    component_distributions = component_distributions
                )
                AIC <- criteria$AIC
                BIC <- criteria$BIC
                ICL <- criteria$ICL
                ICL_MAP <- criteria$ICL_MAP
                return(
                    list(
                        para = vec_para,
                        component_distributions = component_distributions,
                        logLikelihood = logLikelihood,
                        AIC = AIC,
                        BIC = BIC,
                        ICL = ICL,
                        ICL_MAP = ICL_MAP
                    )
                )
            } else {
                return(
                    list(
                        para = vec_para,
                        component_distributions = component_distributions,
                        logLikelihood = logLikelihood
                    )
                )
            }
        }
        if (progress_bar) {
            output <- pblapply(cl = cl, X = 1:N_trials, FUN = function(i) {
                return(func_parallel(i))
            })
        } else {
            output <- parLapply(cl = cl, X = 1:N_trials, fun = function(i) {
                return(func_parallel(i))
            })
        }
        stopCluster(cl)
        #   Extract the results
        all_para <- do.call(rbind, lapply(output, function(x) x$para))
        all_component_distributions <- lapply(output, function(x) x$component_distributions)
        all_logLikelihood <- sapply(output, function(x) x$logLikelihood)
        if (compute_criteria) {
            all_AIC <- sapply(output, function(x) x$AIC)
            all_BIC <- sapply(output, function(x) x$BIC)
            all_ICL <- sapply(output, function(x) x$ICL)
            all_ICL_MAP <- sapply(output, function(x) x$ICL_MAP)
        }
    }
    #---Find the best fit
    best_index <- which.max(all_logLikelihood)
    fit_results <- list()
    fit_results$all_parameters <- all_para
    fit_results$all_logLikelihood <- all_logLikelihood
    fit_results$best_parameters <- all_para[best_index, ]
    fit_results$best_component_distributions <- all_component_distributions[[best_index]]
    fit_results$best_logLikelihood <- all_logLikelihood[best_index]
    if (compute_criteria) {
        fit_results$all_AIC <- all_AIC
        fit_results$all_BIC <- all_BIC
        fit_results$all_ICL <- all_ICL
        fit_results$all_ICL_MAP <- all_ICL_MAP
        fit_results$best_AIC <- all_AIC[best_index]
        fit_results$best_BIC <- all_BIC[best_index]
        fit_results$best_ICL <- all_ICL[best_index]
        fit_results$best_ICL_MAP <- all_ICL_MAP[best_index]
    }
    return(fit_results)
}

DECODE_for_pis <- function(vec_SFS_real,
                           neutral_power,
                           pi_0 = NA,
                           cluster_frequencies,
                           component_distributions,
                           zero_cutoff) {
    # 	Function for parameter transformation
    parameter_transform <- function(parameters) {
        vec_pi <- exp(parameters)
        if (is.na(pi_0)) {
            vec_pi <- vec_pi / sum(vec_pi)
        } else {
            vec_pi <- c(pi_0, (1 - pi_0) * vec_pi / sum(vec_pi))
        }
        return(vec_pi)
    }
    # 	Function for optimization
    func_fit <- function(parameters) {
        vec_pi <- parameter_transform(parameters)
        if (is.na(neutral_power)) {
            loglikelihood <- compute_loglikelihood(
                A = NA,
                vec_K = vec_pi,
                component_distributions = component_distributions,
                vec_SFS_real = vec_SFS_real,
                zero_cutoff = zero_cutoff
            )
        } else {
            loglikelihood <- compute_loglikelihood(
                A = vec_pi[1],
                vec_K = vec_pi[-1],
                component_distributions = component_distributions,
                vec_SFS_real = vec_SFS_real,
                zero_cutoff = zero_cutoff
            )
        }
        return(loglikelihood)
    }
    N_humps <- length(cluster_frequencies)
    if (N_humps == 0) {
        parameters <- c(1)
    } else {
        # 	Initial values for parameters
        if (is.na(neutral_power)) {
            parameters_initial <- rep(0, N_humps)
        } else if (!is.na(pi_0)) {
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
    vec_pi <- parameter_transform(parameters)
    log_L <- func_fit(parameters)
    # 	Prepare the parameters to be returned
    if (is.na(neutral_power)) {
        vec_para <- c()
    } else {
        vec_para <- c(vec_pi[1], neutral_power)
    }
    if (N_humps > 0) {
        for (i in 1:N_humps) {
            if (is.na(neutral_power)) {
                vec_para <- c(vec_para, vec_pi[i], cluster_frequencies[i])
            } else {
                vec_para <- c(vec_para, vec_pi[i + 1], cluster_frequencies[i])
            }
        }
    }
    output <- list()
    output$log_L <- log_L
    output$parameters <- vec_para
    output$component_distributions <- component_distributions
    return(output)
}

DECODE_tail_parameter_sensitivity <- function(vec_SFS_real,
                                              N_humps,
                                              Bayesian_percentile = 0.01,
                                              N_Morris_trials = 300,
                                              SFS_totalsteps,
                                              SFS_convolution_matrix,
                                              fit_results,
                                              pi_0,
                                              pi_0_min,
                                              pi_0_max,
                                              neutral_power,
                                              neutral_power_min,
                                              neutral_power_max,
                                              cluster_frequency_min,
                                              cluster_frequency_max,
                                              zero_cutoff,
                                              compute_parallel,
                                              n_cores,
                                              progress_bar = TRUE) {
    #----------------------------------------Bayesian posterior analysis
    N_keep <- max(10, round(Bayesian_percentile * length(fit_results$all_logLikelihood)))
    indices_keep <- order(fit_results$all_logLikelihood, decreasing = TRUE)[1:N_keep]
    parameters_keep <- fit_results$all_parameters[indices_keep, ]
    pi_0_std <- sd(parameters_keep[, 1])
    neutral_power_std <- sd(parameters_keep[, 2])
    cat("Bayesian posterior analysis for neutral compartment:\n")
    cat(paste0("Neutral power:      std=", neutral_power_std, "\n"))
    cat(paste0("Neutral proportion: std=", pi_0_std, "\n"))
    #----------------------------------------Morris sensitivity analysis
    library(sensitivity)
    #   Find the maximal log-likelihood
    tmp <- DECODE_given_tail_status_and_Ncluster(
        vec_SFS_real = vec_SFS_real,
        N_humps = N_humps,
        with_tail = TRUE,
        N_trials = N_Morris_trials,
        SFS_totalsteps = SFS_totalsteps,
        SFS_convolution_matrix = SFS_convolution_matrix,
        neutral_power = neutral_power,
        pi_0 = pi_0,
        cluster_frequency_min = cluster_frequency_min,
        cluster_frequency_max = cluster_frequency_max,
        zero_cutoff = zero_cutoff,
        compute_criteria = FALSE,
        compute_parallel = compute_parallel,
        n_cores = n_cores,
        progress_bar = FALSE
    )
    max_logLikelihood <- tmp$best_logLikelihood
    #   Define the model to compute Morris sensitivity
    model <- function(parameters) {
        if (compute_parallel == FALSE) {
            ratio_logLikelihood <- rep(0, nrow(parameters))
            if (progress_bar) {
                pb <- txtProgressBar(
                    min = 0,
                    max = nrow(parameters),
                    style = 3,
                    width = 50,
                    char = "+"
                )
            }
            for (i in 1:nrow(parameters)) {
                if (progress_bar) setTxtProgressBar(pb, i)
                #   Find impact of different neutral parameters on likelihood
                tmp <- DECODE_given_tail_status_and_Ncluster(
                    vec_SFS_real = vec_SFS_real,
                    N_humps = N_humps,
                    with_tail = TRUE,
                    N_trials = N_Morris_trials,
                    SFS_totalsteps = SFS_totalsteps,
                    SFS_convolution_matrix = SFS_convolution_matrix,
                    neutral_power = parameters[i, 2],
                    pi_0 = parameters[i, 1],
                    cluster_frequency_min = cluster_frequency_min,
                    cluster_frequency_max = cluster_frequency_max,
                    zero_cutoff = zero_cutoff,
                    compute_criteria = FALSE,
                    compute_parallel = FALSE,
                    progress_bar = FALSE
                )
                ratio_logLikelihood[i] <- tmp$best_logLikelihood / max_logLikelihood
            }
            if (progress_bar) cat("\n")
        } else {
            library(parallel)
            library(pbapply)
            #   Start parallel cluster
            if (is.null(n_cores)) {
                numCores <- detectCores()
            } else {
                numCores <- n_cores
            }
            cl <- makePSOCKcluster(numCores - 1)
            #   Prepare input parameters
            clusterExport(cl, varlist = c(
                "N_end",
                "build_SFS_library",
                "build_SFS_library_Griffiths_Tavare",
                "compute_loglikelihood",
                "compute_SFS",
                "DECODE_for_pis",
                "DECODE_given_tail_status_and_Ncluster"
            ))
            #   Find impact of different neutral parameters on likelihood
            func_parallel <- function(i) {
                tmp <- DECODE_given_tail_status_and_Ncluster(
                    vec_SFS_real = vec_SFS_real,
                    N_humps = N_humps,
                    with_tail = TRUE,
                    N_trials = N_Morris_trials,
                    SFS_totalsteps = SFS_totalsteps,
                    SFS_convolution_matrix = SFS_convolution_matrix,
                    neutral_power = parameters[i, 2],
                    pi_0 = parameters[i, 1],
                    cluster_frequency_min = cluster_frequency_min,
                    cluster_frequency_max = cluster_frequency_max,
                    zero_cutoff = zero_cutoff,
                    compute_criteria = FALSE,
                    compute_parallel = FALSE,
                    progress_bar = FALSE
                )
                ratio_logLikelihood <- tmp$best_logLikelihood / max_logLikelihood
            }
            if (progress_bar) {
                output <- pblapply(cl = cl, X = 1:nrow(parameters), FUN = function(i) {
                    return(func_parallel(i))
                })
            } else {
                output <- parLapply(cl = cl, X = 1:nrow(parameters), fun = function(i) {
                    return(func_parallel(i))
                })
            }
            stopCluster(cl)
            #   Extract the results
            ratio_logLikelihood <- unlist(output)
        }
        return(ratio_logLikelihood)
    }
    #   Compute the Morris sensitivity
    morris_result <- morris(
        model = model,
        factors = 2,
        r = 100,
        design = list(type = "oat", levels = 5),
        binf = c(pi_0_min, neutral_power_min),
        bsup = c(pi_0_max, neutral_power_max),
        scale = TRUE
    )
    #   Extract sensitivities of neutral tail power and mutation count
    mu <- apply(morris_result$ee, 2, mean)
    mu_star <- apply(morris_result$ee, 2, function(x) mean(abs(x)))
    sigma <- apply(morris_result$ee, 2, sd)
    pi_0_mu <- mu[1]
    pi_0_mu_star <- mu_star[1]
    pi_0_sigma <- sigma[1]
    neutral_power_mu <- mu[2]
    neutral_power_mu_star <- mu_star[2]
    neutral_power_sigma <- sigma[2]
    cat("Morris sensitivity analysis for neutral compartment:\n")
    cat(paste0("Neutral power:      mu_*=", neutral_power_mu_star, "; sigma=", neutral_power_sigma, "\n"))
    cat(paste0("Neutral proportion: mu_*=", pi_0_mu_star, "; sigma=", pi_0_sigma, "\n"))
    result <- list()
    result$bayesian_pi_0_std <- pi_0_std
    result$bayesian_neutral_power_std <- neutral_power_std
    result$morris <- morris_result
    result$morris_pi_0_mean <- pi_0_mu
    result$morris_pi_0_mean_abs <- pi_0_mu_star
    result$morris_pi_0_std <- pi_0_sigma
    result$morris_neutral_power_mean <- neutral_power_mu
    result$morris_neutral_power_mean_abs <- neutral_power_mu_star
    result$morris_neutral_power_std <- neutral_power_sigma
    return(result)
}

cluster_count_criteria <- function(num_parameters, log_L, num_samples, vec_SFS_real, vec_para, component_distributions) {
    compute_AIC <- function(log_L) {
        AIC <- 2 * num_parameters - 2 * log_L
        return(AIC)
    }
    compute_BIC <- function(log_L, num_samples) {
        BIC <- num_parameters * log(num_samples) - 2 * log_L
        return(BIC)
    }
    compute_ICL <- function(log_L, num_samples, vec_SFS_real, vec_para, component_distributions) {
        #   Compute the latent variable distributions
        latent_variable_distributions <- component_distributions$SFS_expected_normalized
        for (row in 1:nrow(latent_variable_distributions)) {
            latent_variable_distributions[row, ] <- vec_para[2 * row - 1] * latent_variable_distributions[row, ]
        }
        for (col in 1:ncol(latent_variable_distributions)) {
            latent_variable_distributions[, col] <- latent_variable_distributions[, col] / sum(latent_variable_distributions[, col])
        }
        latent_variable_distributions[which(latent_variable_distributions <= zero_cutoff | is.na(latent_variable_distributions))] <- zero_cutoff
        #   Compute the entropy
        entropy <- sum(vec_SFS_real * colSums(latent_variable_distributions * log(latent_variable_distributions)))
        #   Compute the Bayesian Information Criterion
        BIC <- log_L - 0.5 * num_parameters * log(num_samples)
        #   Compute the Integrated Completed Log-Likelihood
        ICL <- -BIC - entropy
        return(ICL)
    }
    compute_ICL_MAP <- function(log_L, num_samples, vec_SFS_real, vec_para, component_distributions) {
        #   Compute the latent variable distributions
        latent_variable_distributions <- component_distributions$SFS_expected_normalized
        for (row in 1:nrow(latent_variable_distributions)) {
            latent_variable_distributions[row, ] <- vec_para[2 * row - 1] * latent_variable_distributions[row, ]
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
        entropy_MAP <- sum(vec_SFS_real * colSums(indicator_latent_variable_distributions * log(latent_variable_distributions)))
        #   Compute the Bayesian Information Criterion
        BIC <- log_L - 0.5 * num_parameters * log(num_samples)
        #   Compute the Integrated Completed Log-Likelihood
        ICL_MAP <- -BIC - entropy_MAP
        return(ICL_MAP)
    }
    criteria <- list()
    criteria$AIC <- compute_AIC(log_L)
    criteria$BIC <- compute_BIC(log_L, num_samples)
    criteria$ICL <- compute_ICL(log_L, num_samples, vec_SFS_real, vec_para, component_distributions)
    criteria$ICL_MAP <- compute_ICL_MAP(log_L, num_samples, vec_SFS_real, vec_para, component_distributions)
    return(criteria)
}

parameter_conversion <- function(result,
                                 output_parameters_df = TRUE,
                                 mutation_count_for_fitting,
                                 sample_size,
                                 matrix_binomial_sample_size,
                                 matrix_binomial_ploidy) {
    parameters <- result$best_fit$parameters
    component_distributions <- result$best_fit$component_distributions
    tail_status <- result$best_fit$tail_status
    if (!is.null(result$best_fit$tail_sensitivity)) {
        pi_0_std <- result$best_fit$tail_sensitivity$bayesian_pi_0_std
        neutral_power_std <- result$best_fit$tail_sensitivity$bayesian_neutral_power_std
        pi_0_mu <- result$best_fit$tail_sensitivity$morris_pi_0_mean
        pi_0_mu_star <- result$best_fit$tail_sensitivity$morris_pi_0_mean_abs
        pi_0_sigma <- result$best_fit$tail_sensitivity$morris_pi_0_std
        neutral_power_mu <- result$best_fit$tail_sensitivity$morris_neutral_power_mean
        neutral_power_mu_star <- result$best_fit$tail_sensitivity$morris_neutral_power_mean_abs
        neutral_power_sigma <- result$best_fit$tail_sensitivity$morris_neutral_power_std
    } else {
        pi_0_std <- NA
        neutral_power_std <- NA
        pi_0_mu <- NA
        pi_0_mu_star <- NA
        pi_0_sigma <- NA
        neutral_power_mu <- NA
        neutral_power_mu_star <- NA
        neutral_power_sigma <- NA
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
    output <- list()
    output$vec_A <- vec_A
    output$vec_p <- vec_p
    output$vec_K <- vec_K
    output$N_humps <- N_humps
    output$tail_status <- tail_status
    if (output_parameters_df) {
        parameters_df <- data.frame()
        parameters_df[1, "Mutation_count_for_fitting"] <- mutation_count_for_fitting
        parameters_df[1, "Tail"] <- tail_status
        parameters_df[1, "Tail_sensitivity_Bayesian_pi0_std"] <- pi_0_std
        parameters_df[1, "Tail_sensitivity_Bayesian_alpha_std"] <- neutral_power_std
        parameters_df[1, "Tail_sensitivity_Morris_pi0_mean"] <- pi_0_mu
        parameters_df[1, "Tail_sensitivity_Morris_pi0_mean_abs"] <- pi_0_mu_star
        parameters_df[1, "Tail_sensitivity_Morris_pi0_std"] <- pi_0_sigma
        parameters_df[1, "Tail_sensitivity_Morris_alpha_mean"] <- neutral_power_mu
        parameters_df[1, "Tail_sensitivity_Morris_alpha_mean_abs"] <- neutral_power_mu_star
        parameters_df[1, "Tail_sensitivity_Morris_alpha_std"] <- neutral_power_sigma
        if (tail_status) {
            parameters_df[1, "Tail_power"] <- vec_A[2]
            parameters_df[1, "Tail_mutcount_observed"] <-
                vec_A[1] * mutation_count_for_fitting
            parameters_df[1, "Tail_mutcount_predicted"] <-
                vec_A[1] * mutation_count_for_fitting *
                    sum(component_distributions$SFS_exact[1, ]) / sum(component_distributions$SFS_expected[1, ]) *
                    sample_size / matrix_binomial_sample_size
        } else {
            parameters_df[1, "Tail_power"] <- NA
            parameters_df[1, "Tail_mutcount_observed"] <- NA
            parameters_df[1, "Tail_mutcount_predicted"] <- NA
        }
        parameters_df[1, "Cluster_count"] <- N_humps
        if (N_humps > 0) {
            for (k in 1:N_humps) {
                parameters_df[1, paste0("Cluster_frequency_", k)] <- vec_p[k] / matrix_binomial_ploidy
                parameters_df[1, paste0("Cluster_mutcount_observed_", k)] <-
                    vec_K[k] * mutation_count_for_fitting
                parameters_df[1, paste0("Cluster_mutcount_predicted_", k)] <-
                    vec_K[k] * mutation_count_for_fitting *
                        sum(component_distributions$SFS_exact[k + 1, ]) /
                        sum(component_distributions$SFS_expected[k + 1, ])
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

build_SFS_library <- function(neutral_power, cluster_frequencies, SFS_totalsteps, SFS_convolution_matrix) {
    SFS_exact <- c()
    SFS_expected <- c()
    SFS_expected_normalized <- c()
    #   Build the neutral component
    if (is.na(neutral_power)) {
        vec_SFS_GT <- numeric(N_end)
        vec_SFS_expected <- rep(0, SFS_totalsteps)
        vec_SFS_expected_normalized <- rep(0, SFS_totalsteps)
    } else {
        vec_para <- c(1, neutral_power)
        vec_SFS_GT <- build_SFS_library_Griffiths_Tavare(vec_para)
        vec_SFS_expected <- rep(0, SFS_totalsteps)
        for (j in 1:SFS_totalsteps) {
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
            vec_SFS_GT <- build_SFS_library_Griffiths_Tavare(vec_para)
            vec_SFS_expected <- rep(0, SFS_totalsteps)
            for (j in 1:SFS_totalsteps) {
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

build_convolution_matrix <- function(N_end,
                                     SFS_totalsteps,
                                     SFS_totalsteps_base,
                                     r_min,
                                     r_max,
                                     coverage_distribution,
                                     coverage_variables,
                                     sample_coverage,
                                     compute_parallel,
                                     n_cores) {
    #---Build convolution matrix to transform Griffiths-Tavare SFS to expected SFS
    vec_SFS_freq <- seq(0, 1, length.out = SFS_totalsteps + 1)
    if (compute_parallel == FALSE) {
        pb <- txtProgressBar(
            min = 0,
            max = SFS_totalsteps,
            style = 3,
            width = 50,
            char = "+"
        )
        mat_convolution <- matrix(0, nrow = N_end, ncol = SFS_totalsteps)
        for (i in 1:SFS_totalsteps) {
            setTxtProgressBar(pb, i)
            j_lower <- round(SFS_totalsteps_base * vec_SFS_freq[i]) + 1 # x_1*r
            j_upper <- round(SFS_totalsteps_base * vec_SFS_freq[i + 1]) # x_2*r
            for (m in 1:N_end) {
                value <- 0
                for (r in max(r_min, 1):r_max) {
                    value <- value + pdf_coverage(r, sample_coverage) * sum(matrix_binomial_PDF[r, m, j_lower:j_upper])
                }
                mat_convolution[m, i] <- value
            }
        }
    } else {
        library(parallel)
        library(pbapply)
        #   Start parallel cluster
        if (is.null(n_cores)) {
            numCores <- detectCores()
        } else {
            numCores <- n_cores
        }
        cl <- makePSOCKcluster(numCores - 1)
        #   Prepare input parameters
        SFS_totalsteps_base <<- SFS_totalsteps_base
        N_end <<- N_end
        r_min <<- r_min
        r_max <<- r_max
        pdf_coverage <<- pdf_coverage
        coverage_distribution <<- coverage_distribution
        coverage_variables <<- coverage_variables
        sample_coverage <<- sample_coverage
        matrix_binomial_PDF <<- matrix_binomial_PDF
        clusterExport(cl, varlist = c(
            "SFS_totalsteps_base",
            "N_end",
            "r_min",
            "r_max",
            "pdf_coverage",
            "coverage_distribution",
            "coverage_variables",
            "sample_coverage",
            "sample_coverage",
            "matrix_binomial_PDF"
        ))
        #   Build SFS convolution matrix in parallel mode
        output <- pblapply(cl = cl, X = 1:SFS_totalsteps, FUN = function(i) {
            j_lower <- round(SFS_totalsteps_base * vec_SFS_freq[i]) + 1 # x_1*r
            j_upper <- round(SFS_totalsteps_base * vec_SFS_freq[i + 1]) # x_2*r
            mat_convolution_col_i <- rep(0, N_end)
            r_values <- max(r_min, 1):r_max

            mat_convolution_col_i <- sapply(1:N_end, function(m) {
                value <- sum(sapply(r_values, function(r) {
                    pdf_coverage(r, sample_coverage) * sum(matrix_binomial_PDF[r, m, j_lower:j_upper])
                }))
                return(value)
            })
            return(mat_convolution_col_i)
        })
        stopCluster(cl)
        mat_convolution <- do.call(cbind, output)
    }
    return(mat_convolution)
}

build_SFS_library_Griffiths_Tavare <- function(vec_para) {
    #-----------------------------------------------------Get the parameters
    no_hump <- (length(vec_para) - 2) / 2
    para_A <- vec_para[1]
    para_alpha <- vec_para[2]
    para_K <- if (no_hump > 0) numeric(no_hump) else NULL
    para_P <- if (no_hump > 0) numeric(no_hump) else NULL
    for (i in 1:no_hump) {
        para_K[i] <- vec_para[2 * i + 1]
        para_P[i] <- vec_para[2 * i + 2]
    }
    #---------------------------------------Compute the Griffiths-Tavare SFS
    vec_SFS_GT <- numeric(N_end)
    for (m in 2:N_end) {
        if (para_alpha > 0) vec_SFS_GT[m] <- para_A * N_end / (m^para_alpha)
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
