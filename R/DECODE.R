DECODE <- function(mutation_table,
                   criterion = "BIC",
                   criterion_ratio = 1,
                   list_neutral_powers,
                   list_frequencies,
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
                   max_trials = 10000,
                   neutral_tail = NA, # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                   min_N_humps = 1,
                   max_N_humps = Inf,
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
    #---Prepare the SFS library
    cat("Prepare the SFS library...\n")
    library_SFS_component <- build_SFS_library(
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
    if (is.na(neutral_tail)) {
        result_with_tail <- DECODE_given_tail_status(
            vec_SFS_real = vec_SFS_real,
            criterion = criterion,
            criterion_ratio = criterion_ratio,
            min_N_humps = min_N_humps,
            max_N_humps = max_N_humps,
            with_tail = TRUE,
            max_trials = max_trials,
            library_SFS_component = library_SFS_component,
            list_neutral_powers = list_neutral_powers,
            list_frequencies = list_frequencies,
            zero_cutoff = zero_cutoff,
            compute_parallel_fit = compute_parallel_fit,
            n_cores = n_cores
        )
        vec_para_best_final_with_tail <- result_with_tail$best_fit$parameters
        criterion_best_final_with_tail <- result_with_tail$best_fit$selected_criterion_value
        result_without_tail <- DECODE_given_tail_status(
            vec_SFS_real = vec_SFS_real,
            criterion = criterion,
            criterion_ratio = criterion_ratio,
            min_N_humps = min_N_humps,
            max_N_humps = max_N_humps,
            with_tail = FALSE,
            max_trials = max_trials,
            library_SFS_component = library_SFS_component,
            list_neutral_powers = list_neutral_powers,
            list_frequencies = list_frequencies,
            zero_cutoff = zero_cutoff,
            compute_parallel_fit = compute_parallel_fit,
            n_cores = n_cores
        )
        vec_para_best_final_without_tail <- result_without_tail$best_fit$parameters
        criterion_best_final_without_tail <- result_without_tail$best_fit$selected_criterion_value
        if (criterion_best_final_with_tail < criterion_best_final_without_tail) {
            vec_para_best_final <- vec_para_best_final_with_tail
            criterion_best_final <- criterion_best_final_with_tail
            tail_status_final <- TRUE
        } else {
            vec_para_best_final <- vec_para_best_final_without_tail
            criterion_best_final <- criterion_best_final_without_tail
            tail_status_final <- FALSE
        }
    } else if (neutral_tail == TRUE) {
        result <- DECODE_given_tail_status(
            vec_SFS_real = vec_SFS_real,
            criterion = criterion,
            criterion_ratio = criterion_ratio,
            min_N_humps = min_N_humps,
            max_N_humps = max_N_humps,
            with_tail = TRUE,
            max_trials = max_trials,
            library_SFS_component = library_SFS_component,
            list_neutral_powers = list_neutral_powers,
            list_frequencies = list_frequencies,
            zero_cutoff = zero_cutoff,
            compute_parallel_fit = compute_parallel_fit,
            n_cores = n_cores
        )
        vec_para_best_final <- result$best_fit$parameters
        criterion_best_final <- result$best_fit$selected_criterion_value
        tail_status_final <- TRUE
    } else if (neutral_tail == FALSE) {
        result <- DECODE_given_tail_status(
            vec_SFS_real = vec_SFS_real,
            criterion = criterion,
            criterion_ratio = criterion_ratio,
            min_N_humps = min_N_humps,
            max_N_humps = max_N_humps,
            with_tail = FALSE,
            max_trials = max_trials,
            library_SFS_component = library_SFS_component,
            list_neutral_powers = list_neutral_powers,
            list_frequencies = list_frequencies,
            zero_cutoff = zero_cutoff,
            compute_parallel_fit = compute_parallel_fit,
            n_cores = n_cores
        )
        vec_para_best_final <- result$best_fit$parameters
        criterion_best_final <- result$best_fit$selected_criterion_value
        tail_status_final <- FALSE
    }
    #---Report the best fit
    if (tail_status_final) {
        N_humps_best_final <- length(vec_para_best_final) / 2 - 1
        report <- paste0("\n\n\nBEST FIT:\nNeutral tail + ", N_humps_best_final, " humps: ", criterion, " = ", round(criterion_best_final, 3), "; neutral: pi = ", round(vec_para_best_final[1], 3), " with power = ", round(vec_para_best_final[2], 3))
        ii <- 0
    } else {
        N_humps_best_final <- length(vec_para_best_final) / 2
        report <- paste0("\n\n\nBEST FIT:\nNo neutral tail + ", N_humps_best_final, " humps: ", criterion, " = ", round(criterion_best_final, 3))
        ii <- -1
    }
    if (N_humps_best_final > 0) {
        for (i in 1:N_humps_best_final) {
            report <- paste0(report, "; pi = ", round(vec_para_best_final[2 * (i + ii) + 1], 3), " at freq = ", round(vec_para_best_final[2 * (i + ii) + 2], 3))
        }
    }
    report <- paste0(report, "\n\n\n\n")
    cat(report)
    #---Translation to parameters of cancer evolution in the sample
    tmp <- parameter_conversion(
        parameters = vec_para_best_final,
        tail_status = tail_status_final,
        mutation_count_for_fitting = sum(vec_SFS_real),
        library_SFS_component = library_SFS_component,
        list_neutral_powers = list_neutral_powers,
        list_frequencies = list_frequencies,
        sample_size = sample_size,
        matrix_binomial_sample_size = matrix_binomial_sample_size,
        matrix_binomial_ploidy = matrix_binomial_ploidy
    )
    best_fit_parameters <- tmp$parameters_df



    if (!is.null(parameter_filename)) {
        write.table(best_fit_parameters, parameter_filename, sep = "\t", quote = FALSE, row.names = FALSE)
    }
    #---Return the SFS deconvolution results
    DECODE_result <- list()

    DECODE_result$mutational_table <- mutation_table
    DECODE_result$SFS_frequencies <- vec_freq
    DECODE_result$SFS_for_fitting <- vec_SFS_real

    DECODE_result$best_fit <- list()
    DECODE_result$best_fit$parameters <- vec_para_best_final
    DECODE_result$best_fit$selected_criterion <- criterion
    DECODE_result$best_fit$selected_criterion_value <- criterion_best_final
    DECODE_result$best_fit$tail_status <- tail_status_final
    DECODE_result$best_fit$parameters_df <- best_fit_parameters

    DECODE_result$library_SFS_component <- library_SFS_component
    DECODE_result$list_neutral_powers <- list_neutral_powers
    DECODE_result$list_frequencies <- list_frequencies

    # DECODE_result$criterion_best_final <- criterion_best_final
    return(DECODE_result)
}

DECODE_given_tail_status <- function(vec_SFS_real,
                                     criterion,
                                     criterion_ratio,
                                     min_N_humps,
                                     max_N_humps,
                                     with_tail = NA,
                                     max_trials,
                                     library_SFS_component,
                                     list_neutral_powers,
                                     list_frequencies,
                                     zero_cutoff,
                                     compute_parallel_fit,
                                     n_cores) {
    N_humps <- min_N_humps
    criterion_best_final <- Inf
    N_fitting_rounds <- 200
    mutation_count <- sum(vec_SFS_real)
    while (TRUE) {
        #---Find best parameter set, given the number of humps
        if (with_tail) {
            cat(paste0("Inference for ", N_humps, " clusters with neutral tail component...\n"))
        } else {
            cat(paste0("Inference for ", N_humps, " clusters without neutral tail component...\n"))
        }
        fit_results <- DECODE_one_hump(
            vec_SFS_real = vec_SFS_real,
            N_humps = N_humps,
            with_tail = with_tail,
            max_trials = max_trials,
            library_SFS_component = library_SFS_component,
            list_neutral_powers = list_neutral_powers,
            list_frequencies = list_frequencies,
            zero_cutoff = zero_cutoff,
            compute_parallel = compute_parallel_fit,
            n_cores = n_cores
        )
        vec_para_best_current <- fit_results$best_parameters
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
            }
        }
        report <- paste0(report, "\n")
        cat(report)
        # ################################################################################################################################
        # ################################################################################################################################
        # ################################################################################################################################
        # if (with_tail) {
        #     component_distributions <- matrix(unlist(library_SFS_component$neutral$SFS_expected_normalized[[which(list_neutral_powers == vec_para_best_current[2])]]), nrow = 1)
        #     ii <- 0
        # } else {
        #     component_distributions <- c()
        #     ii <- -1
        # }
        # if (N_humps > 0) {
        #     for (i in 1:N_humps) {
        #         component_distributions <- rbind(
        #             component_distributions,
        #             unlist(library_SFS_component$cluster$SFS_expected_normalized[[which(list_frequencies == vec_para_best_current[2 * (i + ii) + 2])]])
        #         )
        #     }
        # }
        # #   Compute the latent variable distributions
        # latent_variable_distributions <- component_distributions
        # for (row in 1:nrow(latent_variable_distributions)) {
        #     latent_variable_distributions[row, ] <- vec_para_best_current[2 * row - 1] * latent_variable_distributions[row, ]
        # }
        # for (col in 1:ncol(latent_variable_distributions)) {
        #     latent_variable_distributions[, col] <- latent_variable_distributions[, col] / sum(latent_variable_distributions[, col])
        # }
        # latent_variable_distributions[which(latent_variable_distributions <= zero_cutoff | is.na(latent_variable_distributions))] <- zero_cutoff
        # #   Compute the MAP allocations for mutations to clusters
        # indicator_latent_variable_distributions <- matrix(0, nrow = nrow(latent_variable_distributions), ncol = ncol(latent_variable_distributions))
        # for (col in 1:ncol(latent_variable_distributions)) {
        #     max_p <- which(latent_variable_distributions[, col] == max(latent_variable_distributions[, col]))[1]
        #     indicator_latent_variable_distributions[max_p, col] <- 1
        # }
        # #   Compute the entropy
        # entropy <- sum(vec_SFS_real * colSums(latent_variable_distributions * log(latent_variable_distributions)))
        # entropy_MAP <- sum(vec_SFS_real * colSums(indicator_latent_variable_distributions * log(latent_variable_distributions)))



        # logL <- fit_results$best_logLikelihood
        # BIC <- logL - 0.5 * length(vec_para_best_current) * log(mutation_count)
        # cat("\n")
        # cat(paste0("Number of parameters: ", length(vec_para_best_current), "\n"))
        # cat(paste0("Sample size:          ", mutation_count, "\n"))
        # cat(paste0("Log-likelihood:       ", logL, "\n"))
        # cat(paste0("\nBIC:                  ", BIC, "\n"))
        # cat(paste0("\nEntropy:              ", entropy, "\n"))
        # # cat(paste0("\nEntropy:              ", entropy_MAP, "\n"))
        # cat(paste0("\nICL:                  ", -BIC - entropy, "\n"))
        # # cat(paste0("\nICL:                  ", -BIC - entropy_MAP, "\n"))



        # # cat("\nLatent distributions:\n")
        # # print(latent_variable_distributions)
        # # cat("Z matrix:\n")
        # # print(latent_variable_distributions * log(latent_variable_distributions))
        # # cat("Observed SFS:\n")
        # # print(vec_SFS_real)
        # # cat("Multiplicative factors:\n")
        # # print(colSums(latent_variable_distributions * log(latent_variable_distributions)))
        # # cat(paste0("Entropy:              ", entropy, "\n"))
        # # cat(paste0("Likelihood:           ", logL, "\n"))
        # # cat(paste0("Completed logL:       ", logL + entropy, "\n"))
        # # cat(paste0("BIC:                  ", BIC, "\n"))
        # # cat(paste0("ICL:                  ", BIC - entropy, "\n"))
        # ################################################################################################################################
        # ################################################################################################################################
        # ################################################################################################################################
        #   Check if the increased hump count leads to lower criterion score...
        # if (criterion_best_current < criterion_best_final) {
        if (criterion_best_current < criterion_ratio * criterion_best_final) {
            #   ... if yes, then update the best fit and continue with 1 more hump
            criterion_best_final <- criterion_best_current
            criterion_all_final <- criterion_all_best_current
            vec_para_best_final <- vec_para_best_current
            N_humps <- N_humps + 1
            if (N_humps > max_N_humps) break
        } else {
            #   ... if no, then stop
            break
        }
    }
    #---Report the best fit
    result <- list()
    result$best_fit <- list()
    result$best_fit$parameters <- vec_para_best_final
    result$best_fit$selected_criterion <- criterion_all_final
    result$best_fit$selected_criterion_value <- criterion_best_final
    result$best_fit$tail_status <- with_tail
    return(result)
}

DECODE_one_hump <- function(vec_SFS_real,
                            N_humps,
                            with_tail,
                            max_trials,
                            library_SFS_component,
                            list_neutral_powers,
                            list_frequencies,
                            zero_cutoff,
                            compute_parallel,
                            n_cores) {
    mutation_count <- sum(vec_SFS_real)
    if (with_tail) {
        num_parameters <- 2 * N_humps + 1
    } else {
        num_parameters <- 2 * N_humps - 1
    }
    if (length(list_neutral_powers) > 1) num_parameters <- num_parameters + 1
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
        latent_variable_distributions <- component_distributions
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
        latent_variable_distributions <- component_distributions
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
    #---Find number of hump frequency combinations there are
    N_trials_true <- choose(length(list_frequencies), N_humps)
    if (with_tail) N_trials_true <- length(list_neutral_powers) * N_trials_true
    N_trials <- min(N_trials_true, max_trials)
    #---Decide list of fixed parameter sets (alpha & p's)
    if (N_trials_true <= N_trials) {
        #   If there are not too many choices: find the best fit among all combinations
        if (with_tail) {
            list_fixed_para <- expand.grid(
                neutral_power = list_neutral_powers,
                cluster_frequencies = combn(list_frequencies, N_humps, simplify = FALSE)
            )
        } else {
            list_fixed_para <- data.frame(neutral_power = rep(NA, N_trials))
            list_fixed_para$cluster_frequencies <- combn(list_frequencies, N_humps, simplify = FALSE)


            # # cluster_frequencies <- combn(list_frequencies, N_humps, simplify = FALSE)
            # # list_fixed_para <- do.call(rbind, cluster_frequencies)

            # list_fixed_para <- expand.grid(
            #     cluster_frequencies = combn(list_frequencies, N_humps, simplify = FALSE)
            # )
            # # print(list_fixed_para)
        }
    } else {
        #   If there are too many choices: try lots of random combinations and find the best fit
        if (with_tail) {
            list_fixed_para <- data.frame(
                neutral_power = sample(list_neutral_powers, N_trials, replace = TRUE)
            )
        } else {
            list_fixed_para <- data.frame(
                neutral_power = rep(NA, N_trials)
            )
        }
        list_fixed_para$cluster_frequencies <- lapply(1:N_trials, function(x) {
            sort(sample(list_frequencies, N_humps, replace = FALSE), decreasing = TRUE)
        })
    }
    #---Find best variable parameters (A & K's) for each fixed parameter set
    if (compute_parallel == FALSE) {
        all_logLikelihood <- c()
        all_para <- c()
        all_AIC <- c()
        all_BIC <- c()
        all_ICL <- c()
        all_ICL_MAP <- c()
        pb <- txtProgressBar(
            min = 0,
            max = nrow(list_fixed_para),
            style = 3,
            width = 50,
            char = "+"
        )
        for (i in 1:nrow(list_fixed_para)) {
            setTxtProgressBar(pb, i)
            neutral_power <- list_fixed_para$neutral_power[i]
            cluster_frequencies <- list_fixed_para$cluster_frequencies[[i]]
            results <- DECODE_hump_component_sizes(
                vec_SFS_real = vec_SFS_real,
                neutral_power = neutral_power,
                vec_p = cluster_frequencies,
                N_humps = N_humps,
                list_neutral_powers = list_neutral_powers,
                list_frequencies = list_frequencies,
                library_SFS_component = library_SFS_component,
                zero_cutoff = zero_cutoff
            )
            logLikelihood <- results$log_L
            vec_para <- results$parameters
            component_distributions <- results$component_distributions
            AIC <- compute_AIC(logLikelihood)
            BIC <- compute_BIC(logLikelihood, mutation_count)
            ICL <- compute_ICL(logLikelihood, mutation_count, vec_SFS_real, vec_para, component_distributions)
            ICL_MAP <- compute_ICL_MAP(logLikelihood, mutation_count, vec_SFS_real, vec_para, component_distributions)
            all_para <- rbind(all_para, vec_para)
            all_logLikelihood <- c(all_logLikelihood, logLikelihood)
            all_AIC <- c(all_AIC, AIC)
            all_BIC <- c(all_BIC, BIC)
            all_ICL <- c(all_ICL, ICL)
            all_ICL_MAP <- c(all_ICL_MAP, ICL_MAP)
        }
        cat("\n")
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
        N_humps <<- N_humps
        list_neutral_powers <<- list_neutral_powers
        list_frequencies <<- list_frequencies
        library_SFS_component <<- library_SFS_component
        matrix_binomial_PDF <<- matrix_binomial_PDF
        zero_cutoff <<- zero_cutoff
        clusterExport(cl, varlist = c(
            "zero_cutoff",
            "vec_SFS_real",
            "N_humps",
            "list_neutral_powers",
            "list_frequencies",
            "library_SFS_component",
            "DECODE_hump_component_sizes",
            "compute_loglikelihood",
            "compute_SFS"
        ))
        #   Find best variable parameters in parallel mode
        output <- pblapply(cl = cl, X = 1:nrow(list_fixed_para), FUN = function(i) {
            neutral_power <- list_fixed_para$neutral_power[i]
            cluster_frequencies <- list_fixed_para$cluster_frequencies[[i]]
            results <- DECODE_hump_component_sizes(
                vec_SFS_real = vec_SFS_real,
                neutral_power = neutral_power,
                vec_p = cluster_frequencies,
                N_humps = N_humps,
                list_neutral_powers = list_neutral_powers,
                list_frequencies = list_frequencies,
                library_SFS_component = library_SFS_component,
                zero_cutoff = zero_cutoff
            )
            logLikelihood <- results$log_L
            vec_para <- results$parameters
            component_distributions <- results$component_distributions
            AIC <- compute_AIC(logLikelihood)
            BIC <- compute_BIC(logLikelihood, mutation_count)
            ICL <- compute_ICL(logLikelihood, mutation_count, vec_SFS_real, vec_para, component_distributions)
            ICL_MAP <- compute_ICL_MAP(logLikelihood, mutation_count, vec_SFS_real, vec_para, component_distributions)
            return(
                list(
                    para = vec_para,
                    logLikelihood = logLikelihood,
                    AIC = AIC,
                    BIC = BIC,
                    ICL = ICL,
                    ICL_MAP = ICL_MAP
                )
            )
        })
        stopCluster(cl)
        #   Extract the results
        all_para <- do.call(rbind, lapply(output, function(x) x$para))
        all_logLikelihood <- sapply(output, function(x) x$logLikelihood)
        all_AIC <- sapply(output, function(x) x$AIC)
        all_BIC <- sapply(output, function(x) x$BIC)
        all_ICL <- sapply(output, function(x) x$ICL)
        all_ICL_MAP <- sapply(output, function(x) x$ICL_MAP)
    }
    #---Find the best fit
    # best_index <- which.max(all_logLikelihood)
    best_index <- which.max(all_logLikelihood)
    fit_results <- list()
    fit_results$best_parameters <- all_para[best_index, ]
    fit_results$best_logLikelihood <- all_logLikelihood[best_index]
    fit_results$best_AIC <- all_AIC[best_index]
    fit_results$best_BIC <- all_BIC[best_index]
    fit_results$best_ICL <- all_ICL[best_index]
    fit_results$best_ICL_MAP <- all_ICL_MAP[best_index]
    return(fit_results)
}

DECODE_hump_component_sizes <- function(vec_SFS_real, neutral_power = NULL, vec_p = NULL, N_humps = NULL, list_neutral_powers, list_frequencies, library_SFS_component, zero_cutoff) {
    if (is.null(neutral_power)) {
        neutral_power <- sample(list_neutral_powers, 1)
    }
    if (is.null(vec_p)) {
        # 	Choose the hump locations at random
        vec_p <- sort(sample(list_frequencies, N_humps, replace = FALSE), decreasing = TRUE)
    } else {
        vec_p <- sort(vec_p, decreasing = TRUE)
    }
    if (is.null(N_humps)) {
        N_humps <- length(vec_p)
    }
    # 	Function for parameter transformation
    parameter_transform <- function(parameters) {
        vec_pi <- exp(parameters)
        vec_pi <- vec_pi / sum(vec_pi)
        return(vec_pi)
    }
    # 	Function for optimization
    func_fit <- function(parameters) {
        vec_pi <- parameter_transform(parameters)
        if (is.na(neutral_power)) {
            loglikelihood <- compute_loglikelihood(
                vec_A = c(NA, NA),
                vec_K = vec_pi,
                vec_p,
                list_neutral_powers = list_neutral_powers,
                list_frequencies = list_frequencies,
                library_SFS_component = library_SFS_component,
                vec_SFS_real = vec_SFS_real,
                zero_cutoff = zero_cutoff
            )
        } else {
            loglikelihood <- compute_loglikelihood(
                vec_A = c(vec_pi[1], neutral_power),
                vec_K = vec_pi[-1],
                vec_p,
                list_neutral_powers = list_neutral_powers,
                list_frequencies = list_frequencies,
                library_SFS_component = library_SFS_component,
                vec_SFS_real = vec_SFS_real,
                zero_cutoff = zero_cutoff
            )
        }
        return(loglikelihood)
    }
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
    vec_pi <- parameter_transform(parameters)
    log_L <- func_fit(parameters)
    # 	Prepare the parameters to be returned
    if (is.na(neutral_power)) {
        vec_para <- c()
        component_distributions <- c()
    } else {
        vec_para <- c(vec_pi[1], neutral_power)
        component_distributions <- matrix(unlist(library_SFS_component$neutral$SFS_expected_normalized[[which(list_neutral_powers == neutral_power)]]), nrow = 1)
    }
    if (N_humps > 0) {
        for (i in 1:N_humps) {
            if (is.na(neutral_power)) {
                vec_para <- c(vec_para, vec_pi[i], vec_p[i])
            } else {
                vec_para <- c(vec_para, vec_pi[i + 1], vec_p[i])
            }
            component_distributions <- rbind(
                component_distributions,
                unlist(library_SFS_component$cluster$SFS_expected_normalized[[which(list_frequencies == vec_p[i])]])
            )
        }
    }
    output <- list()
    output$log_L <- log_L
    output$parameters <- vec_para
    output$component_distributions <- component_distributions
    return(output)
}

parameter_conversion <- function(parameters,
                                 tail_status,
                                 parameters_df = TRUE,
                                 mutation_count_for_fitting,
                                 library_SFS_component,
                                 list_neutral_powers,
                                 list_frequencies,
                                 sample_size,
                                 matrix_binomial_sample_size,
                                 matrix_binomial_ploidy) {
    if (tail_status) {
        vec_A <- parameters[1:2]
        N_humps <- length(parameters) / 2 - 1
        ii <- 0
    } else {
        vec_A <- c(NA, NA)
        N_humps <- length(parameters) / 2 - 1
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
    if (parameters_df) {
        parameters_df <- data.frame()
        parameters_df[1, "Mutation_count_for_fitting"] <- mutation_count_for_fitting
        parameters_df[1, "Tail"] <- tail_status
        if (tail_status) {
            parameters_df[1, "Tail_power"] <- vec_A[2]
            parameters_df[1, "Tail_mutcount_observed"] <-
                vec_A[1] * mutation_count_for_fitting
            parameters_df[1, "Tail_mutcount_predicted"] <-
                vec_A[1] * mutation_count_for_fitting *
                    sum(library_SFS_component$neutral$SFS_exact[[which(list_neutral_powers == vec_A[2])]]) / sum(library_SFS_component$neutral$SFS_expected[[which(list_neutral_powers == vec_A[2])]]) *
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
                        sum(library_SFS_component$cluster$SFS_exact[[which(list_frequencies == vec_p[k])]]) /
                        sum(library_SFS_component$cluster$SFS_expected[[which(list_frequencies == vec_p[k])]])
            }
        }
        output$parameters_df <- parameters_df
    }
    return(output)
}

compute_loglikelihood <- function(vec_A, vec_K, vec_p, list_neutral_powers, list_frequencies, library_SFS_component, vec_SFS_real, zero_cutoff) {
    #----------------Compute the SFS probability distribution from model
    vec_SFS_model <- compute_SFS(
        vec_A = vec_A,
        vec_K = vec_K,
        vec_p = vec_p,
        list_neutral_powers = list_neutral_powers,
        list_frequencies = list_frequencies,
        library_SFS_component = library_SFS_component
    )
    vec_SFS_model[which(vec_SFS_model <= zero_cutoff)] <- zero_cutoff
    vec_SFS_model_normalized <- vec_SFS_model / sum(vec_SFS_model)
    #-----------------------------Compute the log-likelihood(data|model)
    loglikelihood <- sum(log(vec_SFS_model_normalized) * vec_SFS_real)
    return(loglikelihood)
}

compute_SFS <- function(vec_A, vec_K, vec_p, list_neutral_powers, list_frequencies, library_SFS_component) {
    # 	Add the neutral component
    A <- vec_A[1]
    neutral_power <- vec_A[2]
    if (is.na(A) | is.na(neutral_power)) {
        A <- 0
        neutral_power <- list_neutral_powers[1]
    }
    loc <- which(list_neutral_powers == neutral_power)
    vec_SFS_model <- A * unlist(library_SFS_component$neutral$SFS_expected_normalized[[loc]])
    # 	Add the binomial humps
    for (i_hump in seq_along(vec_p)) {
        K <- vec_K[i_hump]
        p <- vec_p[i_hump]
        loc <- which(list_frequencies == p)
        vec_SFS_model <- vec_SFS_model + K * unlist(library_SFS_component$cluster$SFS_expected_normalized[[loc]])
    }
    # 	Return the full SFS
    return(vec_SFS_model)
}

build_SFS_library <- function(N_end,
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
            ############################################################
            ############################################################
            ############################################################
            mat_convolution_col_i <- rep(0, N_end)
            r_values <- max(r_min, 1):r_max

            mat_convolution_col_i <- sapply(1:N_end, function(m) {
                value <- sum(sapply(r_values, function(r) {
                    pdf_coverage(r, sample_coverage) * sum(matrix_binomial_PDF[r, m, j_lower:j_upper])
                }))
                return(value)
            })
            ############################################################
            ############################################################
            ############################################################
            # mat_convolution_col_i <- rep(0, N_end)
            # for (m in 1:N_end) {
            #     value <- 0
            #     for (r in max(r_min, 1):r_max) {
            #         value <- value + pdf_coverage(r,sample_coverage) * sum(matrix_binomial_PDF[r, m, j_lower:j_upper])
            #     }
            #     mat_convolution_col_i[m] <- value
            # }
            ############################################################
            ############################################################
            ############################################################
            return(mat_convolution_col_i)
        })
        stopCluster(cl)
        mat_convolution <- do.call(cbind, output)
    }
    #---Build SFS library for neutral and cluster components
    list_frequencies_tmp <- unique(c(-sort(list_neutral_powers), sort(list_frequencies)))
    library_SFS_component <- list()
    library_SFS_component$neutral$SFS_exact <- list()
    library_SFS_component$neutral$SFS_expected <- list()
    library_SFS_component$neutral$SFS_expected_normalized <- list()
    library_SFS_component$cluster$SFS_exact <- list()
    library_SFS_component$cluster$SFS_expected <- list()
    library_SFS_component$cluster$SFS_expected_normalized <- list()
    for (i in seq_along(list_frequencies_tmp)) {
        p <- list_frequencies_tmp[i]
        if (p <= 0) {
            vec_para <- c(1, -p)
            vec_SFS_GT <- build_SFS_library_Griffiths_Tavare(vec_para)
            vec_SFS_expected <- rep(0, SFS_totalsteps)
            for (j in 1:SFS_totalsteps) {
                vec_SFS_expected[j] <- sum(vec_SFS_GT * mat_convolution[, j])
            }
            library_SFS_component$neutral$SFS_exact[[i]] <- vec_SFS_GT
            library_SFS_component$neutral$SFS_expected[[i]] <- vec_SFS_expected
            library_SFS_component$neutral$SFS_expected_normalized[[i]] <- vec_SFS_expected / sum(vec_SFS_expected)
        } else {
            vec_para <- c(0, 0, 1, p)
            vec_SFS_GT <- build_SFS_library_Griffiths_Tavare(vec_para)
            vec_SFS_expected <- rep(0, SFS_totalsteps)
            for (j in 1:SFS_totalsteps) {
                vec_SFS_expected[j] <- sum(vec_SFS_GT * mat_convolution[, j])
            }
            library_SFS_component$cluster$SFS_exact[[i - length(which(list_frequencies_tmp <= 0))]] <- vec_SFS_GT
            library_SFS_component$cluster$SFS_expected[[i - length(which(list_frequencies_tmp <= 0))]] <- vec_SFS_expected
            library_SFS_component$cluster$SFS_expected_normalized[[i - length(which(list_frequencies_tmp <= 0))]] <- vec_SFS_expected / sum(vec_SFS_expected)
        }
    }
    return(library_SFS_component)
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
                # dbinom change the binopdf in matlab to calculate the binomial prob.
                vec_SFS_GT[m] <- vec_SFS_GT[m] + K * dbinom(m, N_end, P)
            }
        }
    }
    return(vec_SFS_GT)
}
