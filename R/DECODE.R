DECODE <- function(sample_id = "",
                   mutation_table,
                   criterion = "BIC",
                   criterion_ratio = 1,
                   neutral_power_min = 0.5,
                   neutral_power_max = 5,
                   cluster_frequency_min = 0.01,
                   cluster_frequency_max = 1,
                   libPaths_binomial_table,
                   matrix_binomial_sample_size,
                   matrix_binomial_sfs_bincount,
                   matrix_binomial_ploidy,
                   min_variant_read,
                   min_total_read,
                   max_total_read,
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
    mutation_refcounts <- mutation_table$Ref_count
    mutation_altcounts <- mutation_table$Alt_count
    mutation_totcounts <- mutation_refcounts + mutation_altcounts
    mutation_vaf <- mutation_altcounts / mutation_totcounts
    #---Prepare the total readcount distribution
    sample_coverage <- prep_distribution_patient(
        mutations_total_read = mutation_totcounts,
        min_total_read = min_total_read,
        max_total_read = max_total_read
    )
    #---Prepare the real SFS
    vec_freq <- seq(1, sfs_bincount) / sfs_bincount
    vec_SFS_real <- rep(0, sfs_bincount)
    for (j in 1:sfs_bincount) {
        vec_SFS_real[j] <- length(which(
            mutation_altcounts >= min_variant_read &
                mutation_totcounts >= min_total_read &
                mutation_vaf >= ifelse(j == 1, 0, vec_freq[j - 1]) &
                mutation_vaf < vec_freq[j]
        ))
    }
    #---Get DECODE binomial table
    binomial_matrix <- get_binomial_matrix(
        folder = libPaths_binomial_table,
        matrix_binomial_sample_size = matrix_binomial_sample_size,
        matrix_binomial_sfs_bincount = matrix_binomial_sfs_bincount,
        matrix_binomial_ploidy = matrix_binomial_ploidy,
        min_variant_read = min_variant_read,
        min_total_read = min_total_read,
        max_total_read = max_total_read,
        compute_parallel = compute_parallel,
        n_cores = n_cores
    )
    #---Prepare the SFS convolution matrix
    cat(bold(blue("Prepare the SFS convolution matrix...\n")))
    SFS_convolution <- build_convolution_matrix(
        binomial_matrix = binomial_matrix,
        sfs_bincount = sfs_bincount,
        coverage_distribution = coverage_distribution,
        coverage_variables = coverage_variables,
        sample_coverage = sample_coverage
    )
    #---DECODE
    DECODE_result <- list()
    DECODE_result$sample_id <- sample_id
    DECODE_result$sfs_bincount <- sfs_bincount
    DECODE_result$mutational_table <- mutation_table
    DECODE_result$SFS_frequencies <- vec_freq
    DECODE_result$SFS_for_fitting <- vec_SFS_real
    DECODE_result$min_variant_read <- min_variant_read
    DECODE_result$min_total_read <- min_total_read
    if (is.na(neutral_tail)) {
        result_with_tail <- DECODE_given_tail_status(
            vec_SFS_real = vec_SFS_real,
            sfs_bincount = sfs_bincount,
            with_tail = TRUE,
            criterion = criterion,
            criterion_ratio = criterion_ratio,
            min_N_humps = min_N_humps,
            max_N_humps = max_N_humps,
            N_trials = N_trials,
            SFS_convolution = SFS_convolution,
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
        result_without_tail <- DECODE_given_tail_status(
            vec_SFS_real = vec_SFS_real,
            sfs_bincount = sfs_bincount,
            with_tail = FALSE,
            criterion = criterion,
            criterion_ratio = criterion_ratio,
            min_N_humps = min_N_humps,
            max_N_humps = max_N_humps,
            N_trials = N_trials,
            SFS_convolution = SFS_convolution,
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
            final_result <- result_with_tail
        } else {
            final_result <- result_without_tail
        }
    } else if (neutral_tail == TRUE) {
        final_result <- DECODE_given_tail_status(
            vec_SFS_real = vec_SFS_real,
            sfs_bincount = sfs_bincount,
            with_tail = TRUE,
            criterion = criterion,
            criterion_ratio = criterion_ratio,
            min_N_humps = min_N_humps,
            max_N_humps = max_N_humps,
            N_trials = N_trials,
            SFS_convolution = SFS_convolution,
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
        final_result <- DECODE_given_tail_status(
            vec_SFS_real = vec_SFS_real,
            sfs_bincount = sfs_bincount,
            with_tail = FALSE,
            criterion = criterion,
            criterion_ratio = criterion_ratio,
            min_N_humps = min_N_humps,
            max_N_humps = max_N_humps,
            N_trials = N_trials,
            SFS_convolution = SFS_convolution,
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
    }
    DECODE_result$final_fit <- final_result
    #---Report the best fit
    tail_status_final_result <- final_result$best_fit$tail_status
    parameters_final_result <- final_result$best_fit$parameters
    criterion_final_result <- final_result$best_fit$selected_criterion_value
    if (tail_status_final_result) {
        N_humps_final_result <- length(parameters_final_result) / 2 - 1
        report <- bold(underline(red(paste0("Best fit = neutral tail + ", N_humps_final_result, " clusters:\n"))))
        report <- paste0(report, red("Score            : "), yellow(paste0(criterion, " = ", round(criterion_final_result, 3))), "\n")
        report <- paste0(report, red("Neutral component: "), yellow(paste0("pi = ", round(parameters_final_result[1], 3))), red(", "), yellow(paste0("power = ", round(parameters_final_result[2], 3))), "\n")
        ii <- 0
    } else {
        N_humps_final_result <- length(parameters_final_result) / 2
        report <- bold(underline(red(paste0("Best fit = no neutral tail + ", N_humps_final_result, " clusters:\n"))))
        report <- paste0(report, red("Score            : "), yellow(paste0(criterion, " = ", round(criterion_final_result, 3))), "\n")
        ii <- -1
    }
    if (N_humps_final_result > 0) {
        for (i in 1:N_humps_final_result) {
            report <- paste0(report, red(paste0("Cluster ", i, "        : ")), yellow(paste0("pi = ", round(parameters_final_result[2 * (i + ii) + 1], 3))), red(", "), yellow(paste0("f = ", round(parameters_final_result[2 * (i + ii) + 2], 3))), "\n")
        }
    }
    cat(report)
    #---Translation to parameters of cancer evolution in the sample
    tmp <- parameter_conversion(
        result = final_result,
        mutation_count_for_fitting = sum(vec_SFS_real),
        sample_size = sample_size,
        matrix_binomial_sample_size = matrix_binomial_sample_size,
        matrix_binomial_ploidy = matrix_binomial_ploidy
    )
    DECODE_result$final_fit$parameters_df <- tmp$parameters_df
    #---Return the SFS deconvolution results
    return(DECODE_result)
}

DECODE_given_tail_status <- function(vec_SFS_real,
                                     criterion,
                                     criterion_ratio,
                                     min_N_humps,
                                     max_N_humps,
                                     with_tail,
                                     N_trials,
                                     sfs_bincount,
                                     SFS_convolution,
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
        fit_results <- DECODE_given_tail_status_and_Ncluster(
            vec_SFS_real = vec_SFS_real,
            N_humps = N_humps,
            with_tail = with_tail,
            N_trials = N_trials,
            sfs_bincount = sfs_bincount,
            SFS_convolution = SFS_convolution,
            neutral_power_min = neutral_power_min,
            neutral_power_max = neutral_power_max,
            cluster_frequency_min = cluster_frequency_min,
            cluster_frequency_max = cluster_frequency_max,
            zero_cutoff = zero_cutoff,
            compute_parallel = compute_parallel,
            n_cores = n_cores
        )
        all_fits[[paste0(N_humps, "_clusters")]] <- fit_results
        vec_para_best_current <- fit_results$best$parameters
        component_distributions_best_current <- fit_results$best$component_distributions
        criterion_all_best_current <- fit_results$best$criteria
        criterion_best_current <- criterion_all_best_current[[criterion]]
        #   Report the best fit for the current hump count
        cluster_pis <- Inf
        report <- paste0(blue("Score            : "), cyan(paste0(criterion, " = ", round(criterion_best_current, 3))), "\n")
        if (with_tail) {
            N_humps <- length(vec_para_best_current) / 2 - 1
            report <- paste0(report, blue("Neutral component: "), cyan(paste0("pi = ", round(vec_para_best_current[1], 3))), blue(", "), cyan(paste0("power = ", round(vec_para_best_current[2], 3))), "\n")
            ii <- 0
        } else {
            N_humps <- length(vec_para_best_current) / 2
            ii <- -1
        }
        if (N_humps > 0) {
            for (i in 1:N_humps) {
                report <- paste0(report, blue(paste0("Cluster ", i, "        : ")), cyan(paste0("pi = ", round(vec_para_best_current[2 * (i + ii) + 1], 3))), blue(", "), cyan(paste0("f = ", round(vec_para_best_current[2 * (i + ii) + 2], 3))), "\n")
                cluster_pis <- c(cluster_pis, vec_para_best_current[2 * (i + ii) + 1])
            }
        }
        cat(report)
        #   Check if the increased hump count leads to lower criterion score without tiny selective components...
        if ((N_humps == min_N_humps) | ((criterion_best_current < criterion_ratio * criterion_best_final) & (min(cluster_pis) >= pi_cutoff))) {
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
        with_tail <- FALSE
        vec_para_best_final[seq(3, length(vec_para_best_final), by = 2)] <- vec_para_best_final[seq(3, length(vec_para_best_final), by = 2)] / sum(vec_para_best_final[seq(3, length(vec_para_best_final), by = 2)])
        vec_para_best_final <- vec_para_best_final[-c(1, 2)]
        component_distributions_best_final$SFS_exact[1, ] <- rep(0, length(component_distributions_best_final$SFS_exact[1, ]))
        component_distributions_best_final$SFS_expected[1, ] <- rep(0, length(component_distributions_best_final$SFS_expected[1, ]))
        component_distributions_best_final$SFS_expected_normalized[1, ] <- rep(0, length(component_distributions_best_final$SFS_expected_normalized[1, ]))
    }
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
    return(result)
}

DECODE_given_tail_status_and_Ncluster <- function(vec_SFS_real,
                                                  N_humps,
                                                  with_tail,
                                                  N_trials,
                                                  sfs_bincount,
                                                  SFS_convolution,
                                                  neutral_power_min,
                                                  neutral_power_max,
                                                  cluster_frequency_min,
                                                  cluster_frequency_max,
                                                  zero_cutoff,
                                                  compute_criteria = TRUE,
                                                  compute_parallel,
                                                  n_cores,
                                                  progress_bar = TRUE) {
    N_end <- SFS_convolution$N_end
    SFS_convolution_matrix <- SFS_convolution$convolution_matrix
    #---Function to perform one trial to find A & K's
    func_one_trial <- function(with_tail,
                               compute_criteria,
                               N_humps,
                               neutral_power_min,
                               neutral_power_max,
                               cluster_frequency_min,
                               cluster_frequency_max,
                               vec_SFS_real,
                               sfs_bincount,
                               N_end,
                               SFS_convolution_matrix,
                               zero_cutoff) {
        #   Sample neutral component power and cluster frequencies
        neutral_power <- ifelse(with_tail, runif(1, neutral_power_min, neutral_power_max), NA)
        cluster_frequencies <- sort(runif(N_humps, cluster_frequency_min, cluster_frequency_max), decreasing = TRUE)
        #   Build the SFS component library
        component_distributions <- build_SFS_library(
            neutral_power = neutral_power,
            cluster_frequencies = cluster_frequencies,
            sfs_bincount = sfs_bincount,
            SFS_convolution_matrix = SFS_convolution_matrix,
            N_end = N_end
        )
        #   Find component proportions
        results <- DECODE_for_pis(
            vec_SFS_real = vec_SFS_real,
            neutral_power = neutral_power,
            cluster_frequencies = cluster_frequencies,
            component_distributions = component_distributions,
            zero_cutoff = zero_cutoff
        )
        if (compute_criteria) {
            mutation_count <- sum(vec_SFS_real)
            num_parameters <- ifelse(with_tail, 2 * N_humps + 1, 2 * N_humps - 1)
            if (neutral_power_max > neutral_power_min) num_parameters <- num_parameters + 1
            criteria <- cluster_count_criteria(
                num_parameters = num_parameters,
                log_L = results$log_L,
                num_samples = mutation_count,
                vec_SFS_real = vec_SFS_real,
                parameters = results$parameters,
                component_distributions = component_distributions,
                zero_cutoff = zero_cutoff
            )
        }
        output <- list()
        output$parameters <- results$parameters
        output$logLikelihood <- results$log_L
        output$component_distributions <- results$component_distributions
        if (compute_criteria) output$criteria <- criteria
        return(output)
    }
    #---Find best variable parameters (A & K's) for each fixed parameter set from many trials
    if (compute_parallel == FALSE) {
        all_logLikelihood <- c()
        all_para <- c()
        all_component_distributions <- list()
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
                compute_criteria = compute_criteria,
                N_humps = N_humps,
                neutral_power_min = neutral_power_min,
                neutral_power_max = neutral_power_max,
                cluster_frequency_min = cluster_frequency_min,
                cluster_frequency_max = cluster_frequency_max,
                vec_SFS_real = vec_SFS_real,
                sfs_bincount = sfs_bincount,
                N_end = N_end,
                SFS_convolution_matrix = SFS_convolution_matrix,
                zero_cutoff = zero_cutoff
            )
            all_para <- rbind(all_para, trial_result$parameters)
            all_component_distributions[[i]] <- trial_result$component_distributions
            all_logLikelihood <- c(all_logLikelihood, trial_result$logLikelihood)
            if (compute_criteria) {
                all_criteria <- rbind(all_criteria, trial_result$criteria)
            }
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
            "compute_criteria", "vec_SFS_real", "SFS_convolution_matrix",
            "neutral_power_min", "neutral_power_max",
            "cluster_frequency_min", "cluster_frequency_max",
            "build_SFS_library", "build_SFS_library_Griffiths_Tavare",
            "DECODE_for_pis", "compute_loglikelihood", "compute_SFS", "cluster_count_criteria"
        ), envir = environment())
        #   Find best variable parameters in parallel mode
        if (progress_bar) {
            output <- pblapply(cl = cl, X = 1:N_trials, FUN = function(i) {
                func_one_trial(
                    with_tail = with_tail,
                    compute_criteria = compute_criteria,
                    N_humps = N_humps,
                    neutral_power_min = neutral_power_min,
                    neutral_power_max = neutral_power_max,
                    cluster_frequency_min = cluster_frequency_min,
                    cluster_frequency_max = cluster_frequency_max,
                    vec_SFS_real = vec_SFS_real,
                    sfs_bincount = sfs_bincount,
                    N_end = N_end,
                    SFS_convolution_matrix = SFS_convolution_matrix,
                    zero_cutoff = zero_cutoff
                )
            })
        } else {
            output <- parLapply(cl = cl, X = 1:N_trials, fun = function(i) {
                func_one_trial(
                    with_tail = with_tail,
                    compute_criteria = compute_criteria,
                    N_humps = N_humps,
                    neutral_power_min = neutral_power_min,
                    neutral_power_max = neutral_power_max,
                    cluster_frequency_min = cluster_frequency_min,
                    cluster_frequency_max = cluster_frequency_max,
                    vec_SFS_real = vec_SFS_real,
                    sfs_bincount = sfs_bincount,
                    N_end = N_end,
                    SFS_convolution_matrix = SFS_convolution_matrix,
                    zero_cutoff = zero_cutoff
                )
            })
        }
        stopCluster(cl)
        #   Extract the results
        all_para <- do.call(rbind, lapply(output, function(x) x$parameters))
        all_component_distributions <- lapply(output, function(x) x$component_distributions)
        all_logLikelihood <- sapply(output, function(x) x$logLikelihood)
        all_criteria <- do.call(rbind, lapply(output, function(x) x$criteria))
    }
    #---Find the best fit
    best_index <- which.max(all_logLikelihood)
    fit_results <- list()
    fit_results$all <- list()
    fit_results$all$parameters <- all_para
    fit_results$all$logLikelihood <- all_logLikelihood
    fit_results$all$criteria <- all_criteria
    fit_results$best <- list()
    fit_results$best$parameters <- all_para[best_index, ]
    fit_results$best$logLikelihood <- all_logLikelihood[best_index]
    fit_results$best$criteria <- all_criteria[best_index, ]
    fit_results$best$component_distributions <- all_component_distributions[[best_index]]
    return(fit_results)
}

DECODE_for_pis <- function(vec_SFS_real,
                           neutral_power,
                           cluster_frequencies,
                           component_distributions,
                           zero_cutoff) {
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
    return(output)
}

cluster_count_criteria <- function(num_parameters, log_L, num_samples, vec_SFS_real, parameters, component_distributions, zero_cutoff) {
    compute_AIC <- function(log_L) {
        AIC <- 2 * num_parameters - 2 * log_L
        return(AIC)
    }
    compute_BIC <- function(log_L, num_samples) {
        BIC <- num_parameters * log(num_samples) - 2 * log_L
        return(BIC)
    }
    compute_ICL <- function(log_L, num_samples, vec_SFS_real, parameters, component_distributions) {
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
        entropy <- sum(vec_SFS_real * colSums(latent_variable_distributions * log(latent_variable_distributions)))
        #   Compute the Bayesian Information Criterion
        BIC <- log_L - 0.5 * num_parameters * log(num_samples)
        #   Compute the Integrated Completed Log-Likelihood
        ICL <- -BIC - entropy
        return(ICL)
    }
    compute_ICL_MAP <- function(log_L, num_samples, vec_SFS_real, parameters, component_distributions) {
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
        entropy_MAP <- sum(vec_SFS_real * colSums(indicator_latent_variable_distributions * log(latent_variable_distributions)))
        #   Compute the Bayesian Information Criterion
        BIC <- log_L - 0.5 * num_parameters * log(num_samples)
        #   Compute the Integrated Completed Log-Likelihood
        ICL_MAP <- -BIC - entropy_MAP
        return(ICL_MAP)
    }
    criteria <- data.frame(
        AIC = compute_AIC(log_L),
        BIC = compute_BIC(log_L, num_samples),
        ICL = compute_ICL(log_L, num_samples, vec_SFS_real, parameters, component_distributions),
        ICL_MAP = compute_ICL_MAP(log_L, num_samples, vec_SFS_real, parameters, component_distributions)
    )
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
        if (tail_status) {
            parameters_df[1, "Tail_power"] <- vec_A[2]
            parameters_df[1, "Tail_mutcount_observed"] <-
                vec_A[1] * mutation_count_for_fitting
            parameters_df[1, "Tail_mutcount_predicted"] <-
                vec_A[1] * mutation_count_for_fitting *
                    sum(component_distributions$SFS_exact[1, ]) /
                    sum(component_distributions$SFS_expected[1, ]) *
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
                    vec_K[k] * mutation_count_for_fitting /
                        sum(component_distributions$SFS_expected[k + 1, ]) /
                        sum(component_distributions$SFS_exact[k + 1, ])
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

get_binomial_matrix <- function(folder,
                                matrix_binomial_sample_size,
                                matrix_binomial_sfs_bincount,
                                matrix_binomial_ploidy,
                                min_variant_read,
                                min_total_read,
                                max_total_read,
                                compute_parallel = TRUE,
                                n_cores = NULL) {
    #---Get the binomial PDF table if already produced before
    file_list <- list.files(path = folder, pattern = "\\.rds$", full.names = TRUE)
    for (file in file_list) {
        output <- readRDS(file)
        if (identical(class(output), "DECODE_binomial_matrix") &
            output$sample_size == matrix_binomial_sample_size &
            output$sfs_bincount == matrix_binomial_sfs_bincount &
            output$ploidy == matrix_binomial_ploidy &
            output$min_variant_read == min_variant_read &
            output$min_total_read == min_total_read &
            output$max_total_read == max_total_read) {
            cat(bold(blue("Retrieve binomial table...\n")))
            return(output)
        }
    }
    output <- NULL
    #---Function to compute the binomial PDF table for a given r
    func_submatrix_binomial_PDF <- function(r,
                                            matrix_binomial_sample_size,
                                            matrix_binomial_sfs_bincount,
                                            matrix_binomial_ploidy,
                                            min_variant_read,
                                            min_total_read) {
        submatrix_binomial_PDF <- array(0, dim = c(matrix_binomial_sample_size, matrix_binomial_sfs_bincount))
        for (m in 1:matrix_binomial_sample_size) {
            for (i in 1:matrix_binomial_sfs_bincount) {
                if (r < min(min_variant_read, min_total_read)) next
                #   Find boundaries for s = ( (i-1)*r/matrix_binomial_sfs_bincount,i*r/matrix_binomial_sfs_bincount ]
                r1 <- r * (i - 1) / matrix_binomial_sfs_bincount
                r1 <- ifelse(r1 %% 1 == 0, r1 + 1, ceiling(r1))
                r1 <- max(min_variant_read, r1, 1)
                r2 <- floor(r * i / matrix_binomial_sfs_bincount)
                if (r1 > r2) next
                #   Find P{s/r in ((i-1)/matrix_binomial_sfs_bincount,i/matrix_binomial_sfs_bincount] | m, r}
                Prob <- pbinom(r2, size = r, prob = m / (matrix_binomial_sample_size * matrix_binomial_ploidy)) -
                    pbinom(r1 - 1, size = r, prob = m / (matrix_binomial_sample_size * matrix_binomial_ploidy))
                submatrix_binomial_PDF[m, i] <- Prob
            }
        }
        return(submatrix_binomial_PDF)
    }
    #---Prepare the binomial PDF table
    cat(bold(blue("Prepare binomial table...\n")))
    if (compute_parallel == FALSE) {
        pb <- txtProgressBar(
            min = 1,
            max = max_total_read,
            style = 3,
            width = 50,
            char = "+"
        )
        matrix_binomial_PDF <- array(0, dim = c(max_total_read, matrix_binomial_sample_size, matrix_binomial_sfs_bincount))
        for (r in 1:max_total_read) {
            setTxtProgressBar(pb, r)
            matrix_binomial_PDF[r, , ] <- func_submatrix_binomial_PDF(
                r = r,
                matrix_binomial_sample_size = matrix_binomial_sample_size,
                matrix_binomial_sfs_bincount = matrix_binomial_sfs_bincount,
                matrix_binomial_ploidy = matrix_binomial_ploidy,
                min_variant_read = min_variant_read,
                min_total_read = min_total_read
            )
        }
        cat("\n")
    } else {
        library(parallel)
        library(pbapply)
        #   Start parallel cluster
        numCores <- ifelse(is.null(n_cores), detectCores(), n_cores)
        cl <- makePSOCKcluster(numCores - 1)
        #   Prepare input parameters
        clusterExport(cl, varlist = c(
            "matrix_binomial_sample_size",
            "matrix_binomial_sfs_bincount",
            "matrix_binomial_ploidy",
            "min_variant_read",
            "min_total_read",
            "max_total_read"
        ), envir = environment())
        #   Compute each sub-array
        output <- pblapply(cl = cl, X = 1:max_total_read, FUN = function(r) {
            return(func_submatrix_binomial_PDF(
                r = r,
                matrix_binomial_sample_size = matrix_binomial_sample_size,
                matrix_binomial_sfs_bincount = matrix_binomial_sfs_bincount,
                matrix_binomial_ploidy = matrix_binomial_ploidy,
                min_variant_read = min_variant_read,
                min_total_read = min_total_read
            ))
        })
        stopCluster(cl)
        matrix_binomial_PDF <- array(0, dim = c(max_total_read, matrix_binomial_sample_size, matrix_binomial_sfs_bincount))
        for (r in 1:max_total_read) {
            matrix_binomial_PDF[r, , ] <- output[[r]]
        }
    }
    #---Create DECODE_binomial_matrix object
    output <- list()
    output$matrix_binomial_PDF <- matrix_binomial_PDF
    output$sample_size <- matrix_binomial_sample_size
    output$sfs_bincount <- matrix_binomial_sfs_bincount
    output$ploidy <- matrix_binomial_ploidy
    output$min_variant_read <- min_variant_read
    output$min_total_read <- min_total_read
    output$max_total_read <- max_total_read
    class(output) <- "DECODE_binomial_matrix"
    filename <- paste0(
        folder, "/DECODE_binomial_matrix_",
        matrix_binomial_sample_size, "_",
        matrix_binomial_sfs_bincount, "_",
        matrix_binomial_ploidy, "_",
        min_variant_read, "_",
        min_total_read, "_",
        max_total_read,
        ".rds"
    )
    saveRDS(output, file = filename)
    return(output)
}

build_convolution_matrix <- function(binomial_matrix,
                                     sfs_bincount,
                                     coverage_distribution,
                                     coverage_variables,
                                     sample_coverage) {
    library(progress)
    matrix_binomial_PDF <- binomial_matrix$matrix_binomial_PDF
    N_end <- binomial_matrix$sample_size
    SFS_totalsteps_base <- binomial_matrix$sfs_bincount
    min_total_read <- binomial_matrix$min_total_read
    max_total_read <- binomial_matrix$max_total_read
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
    #---Build convolution matrix to transform Griffiths-Tavare SFS to expected SFS
    vec_SFS_freq <- seq(0, 1, length.out = sfs_bincount + 1)
    mat_convolution <- matrix(0, nrow = N_end, ncol = sfs_bincount)
    pb <- txtProgressBar(
        min = 0,
        max = sfs_bincount,
        style = 3,
        width = 50,
        char = "+"
    )
    start_time <- Sys.time()
    for (i in 1:sfs_bincount) {
        setTxtProgressBar(pb, i)
        elapsed_time <- Sys.time() - start_time
        estimated_total_time <- (elapsed_time / i) * sfs_bincount
        remaining_time <- estimated_total_time - elapsed_time
        cat(sprintf(" elapsed=%02ds remaining=%02ds\r", as.integer(elapsed_time), as.integer(remaining_time)))
        j_lower <- round(SFS_totalsteps_base * vec_SFS_freq[i]) + 1 # x_1*r
        j_upper <- round(SFS_totalsteps_base * vec_SFS_freq[i + 1]) # x_2*r
        for (m in 1:N_end) {
            B_values <- sapply(sample_coverage_distribution$total_readcount, function(r) {
                sum(matrix_binomial_PDF[r, m, j_lower:j_upper])
            })
            mat_convolution[m, i] <- sum(sample_coverage_distribution$pdf * B_values)
        }
    }
    cat("\n")
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

prep_distribution_patient <- function(mutations_total_read, min_total_read, max_total_read) {
    sample_coverage <- data.frame(
        total_readcount = min_total_read:max_total_read,
        pdf = 0
    )
    #---Compute coverage distribution
    for (i in 1:nrow(sample_coverage)) {
        sample_coverage$pdf[i] <- sum(mutations_total_read == sample_coverage$total_readcount[i])
    }
    if (sum(sample_coverage$pdf) == 0) stop("No mutation with total reads within [min_total_read, max_total_read]")
    #---Normalize distribution
    sample_coverage$pdf <- sample_coverage$pdf / sum(sample_coverage$pdf)
    sample_coverage <- sample_coverage[sample_coverage$pdf != 0, ]
    #---Return results
    return(sample_coverage)
}
