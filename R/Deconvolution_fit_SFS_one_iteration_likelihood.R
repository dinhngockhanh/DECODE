fit_SFS_likelihood <- function(mutation_altcounts,
                               mutation_refcounts,
                               mutation_identity = NULL,
                               list_frequencies,
                               matrix_binomial_PDF,
                               matrix_binomial_sample_size,
                               matrix_binomial_sfs_stepcount,
                               matrix_binomial_ploidy,
                               sample_size,
                               SFS_totalsteps,
                               r_min,
                               r_max,
                               option_dist_coverage,
                               dist_coverage_var_1,
                               max_trials = 10000,
                               compute_parallel = TRUE,
                               n_cores = NULL) {
    mutation_totcounts <- mutation_refcounts + mutation_altcounts
    #---Prepare the total readcount distribution
    cat("Prepare the total readcount distribution...\n")
    TCGA_coverage_PDF <- prep_distribution_patient(mutation_totcounts)
    #---Prepare the real SFS
    cat("Prepare the real SFS...\n")
    no_mutations_total <- length(mutation_refcounts)
    vec_freq <- seq(1, SFS_totalsteps) / SFS_totalsteps
    vec_SFS_real <- rep(0, SFS_totalsteps)
    mutation_count <- 0
    for (j in 1:no_mutations_total) {
        no_variant <- mutation_altcounts[j]
        no_total <- mutation_refcounts[j] + mutation_altcounts[j]
        if (no_variant >= min_variant_read && no_total >= min_total_read) {
            mutation_count <- mutation_count + 1
            VAF <- no_variant / no_total
            pos <- which(vec_freq >= VAF)[1]
            vec_SFS_real[pos] <- vec_SFS_real[pos] + 1
        }
    }
    #---Prepare the SFS library
    cat("Prepare the SFS library...\n")
    library_SFS_component <- list()
    library_SFS_component$SFS_exact <- list()
    library_SFS_component$SFS_expected <- list()
    library_SFS_component$SFS_expected_normalized <- list()
    vec_para <- c(1)
    output <- compute_SFS(vec_para)
    library_SFS_component$SFS_exact[[1]] <- output$SFS_exact
    library_SFS_component$SFS_expected[[1]] <- output$SFS_expected
    library_SFS_component$SFS_expected_normalized[[1]] <- output$SFS_expected / sum(output$SFS_expected)
    if (compute_parallel == FALSE) {
        #---------------------------Build SFS library in sequential mode
        pb <- txtProgressBar(
            min = 0,
            max = length(list_frequencies),
            style = 3,
            width = 50,
            char = "="
        )
        for (i in seq_along(list_frequencies)) {
            setTxtProgressBar(pb, i)
            p <- list_frequencies[i]
            vec_para <- c(0, p, 1)
            output <- compute_SFS(vec_para)
            library_SFS_component$SFS_exact[[i + 1]] <- output$SFS_exact
            library_SFS_component$SFS_expected[[i + 1]] <- output$SFS_expected
            library_SFS_component$SFS_expected_normalized[[i + 1]] <- output$SFS_expected / sum(output$SFS_expected)
        }
    } else {
        #-----------------------------Build SFS library in parallel mode
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
        compute_SFS <<- compute_SFS
        SFS_Griffiths_Tavare <<- SFS_Griffiths_Tavare
        pdf_coverage <<- pdf_coverage
        N_end <<- matrix_binomial_sample_size
        SFS_totalsteps <<- SFS_totalsteps
        SFS_totalsteps_base <<- matrix_binomial_sfs_stepcount
        r_min <<- r_min
        r_max <<- r_max
        option_dist_coverage <<- option_dist_coverage
        dist_coverage_var_1 <<- dist_coverage_var_1
        matrix_binomial_PDF <<- matrix_binomial_PDF
        clusterExport(cl, varlist = c(
            "compute_SFS",
            "SFS_Griffiths_Tavare",
            "pdf_coverage",
            "N_end",
            "SFS_totalsteps",
            "SFS_totalsteps_base",
            "r_min",
            "r_max",
            "option_dist_coverage",
            "dist_coverage_var_1",
            "matrix_binomial_PDF"
        ))
        #   Build SFS library in parallel mode
        output <- pblapply(cl = cl, X = 1:length(list_frequencies), FUN = function(i) {
            p <- list_frequencies[i]
            vec_para <- c(0, p, 1)
            return(compute_SFS(vec_para))
        })
        for (i in seq_along(list_frequencies)) {
            library_SFS_component$SFS_exact[[i + 1]] <- output[[i]]$SFS_exact
            library_SFS_component$SFS_expected[[i + 1]] <- output[[i]]$SFS_expected
            library_SFS_component$SFS_expected_normalized[[i + 1]] <- output[[i]]$SFS_expected / sum(output[[i]]$SFS_expected)
        }
    }
    #---SFS deconvolution
    compute_AIC <- function(log_L, num_params) {
        return(2 * num_params - 2 * log_L)
    }
    N_humps <- 0
    AIC_best_final <- Inf
    N_fitting_rounds <- 200
    while (TRUE) {
        num_parameters <- 1 + 2 * N_humps
        AIC_best_current <- Inf
        vec_para_best_current <- c()
        #   Find best parameter set, given the number of humps
        if (N_humps == 0) {
            #   If 0 humps: the SFS consists of only neutral component
            vec_para_best_current <- 1
            log_L <- error_SFS_one_iteration_likelihood(vec_para_best_current, c(), list_frequencies, library_SFS_component, vec_SFS_real)
            AIC_current <- compute_AIC(log_L, num_parameters)
            AIC_best_current <- AIC_current
        } else {
            #   Find number of hump frequency combinations there are
            N_trials_true <- choose(length(list_frequencies), N_humps)
            N_trials <- min(N_trials_true, max_trials)
            #   Find best parameter set
            if (N_trials_true <= N_trials) {
                #   If there are not too many choices: find the best fit among all combinations
                tmp <- combn(list_frequencies, N_humps)
                for (i in 1:ncol(tmp)) {
                    vec_p <- tmp[, i]
                    results <- fit_SFS_given_humpcount_likelihood(
                        vec_SFS_real = vec_SFS_real,
                        vec_p = vec_p,
                        N_humps = N_humps,
                        list_frequencies = list_frequencies,
                        library_SFS_component = library_SFS_component
                    )
                    log_L <- results$log_L
                    vec_para <- results$parameters
                    AIC_current <- compute_AIC(log_L, num_parameters)
                    if (AIC_current < AIC_best_current) {
                        AIC_best_current <- AIC_current
                        vec_para_best_current <- vec_para
                    }
                }
            } else {
                #   If there are too many choices: try lots of random combinations and find the best fit
                for (i in 1:N_trials) {
                    results <- fit_SFS_given_humpcount_likelihood(
                        vec_SFS_real = vec_SFS_real,
                        N_humps = N_humps,
                        list_frequencies = list_frequencies,
                        library_SFS_component = library_SFS_component
                    )
                    log_L <- results$log_L
                    vec_para <- results$parameters
                    AIC_current <- compute_AIC(log_L, num_parameters)
                    if (AIC_current < AIC_best_current) {
                        AIC_best_current <- AIC_current
                        vec_para_best_current <- vec_para
                    }
                }
            }
        }
        #   Report the best fit for the current hump count
        report <- paste0(N_humps, " humps: AIC = ", round(AIC_best_current, 3), "; neutral: pi = ", round(vec_para_best_current[1], 3))
        if (length(vec_para_best_current) > 1) {
            for (i in 1:N_humps) {
                report <- paste0(report, "; pi = ", round(vec_para_best_current[2 * i + 1], 3), " at freq = ", round(vec_para_best_current[2 * i], 3))
            }
        }
        report <- paste0(report, "\n")
        cat(report)
        #   Check if the increased hump count leads to lower AIC...
        if (AIC_best_current < AIC_best_final) {
            #   ... if yes, then update the best fit and continue with 1 more hump
            AIC_best_final <- AIC_best_current
            vec_para_best_final <- vec_para_best_current
            N_humps <- N_humps + 1
        } else {
            #   ... if no, then stop
            break
        }
    }
    #---Report the best fit
    N_humps_best_final <- (length(vec_para_best_final) - 1) / 2
    report <- paste0("\n\n\nBEST FIT:\n", N_humps_best_final, " humps: AIC = ", round(AIC_best_final, 3), "; neutral: pi = ", round(vec_para_best_final[1], 3))
    if (length(vec_para_best_final) > 1) {
        for (i in 1:N_humps_best_final) {
            report <- paste0(report, "; pi = ", round(vec_para_best_final[2 * i + 1], 3), " at freq = ", round(vec_para_best_final[2 * i], 3))
        }
    }
    report <- paste0(report, "\n\n\n\n")
    cat(report)
    #---Plot the SFS deconvolution
    filename_plot <- paste0(R_workplace, "/", folder_workplace, "SFS_", n_simulation, ".png")
    png(filename_plot)
    bar_centers <- barplot(height = vec_SFS_real, names.arg = vec_freq, col = "blue", main = "SFS Fitting Results")
    vec_A_and_K <- vec_para_best_final[seq(1, length(vec_para_best_final), by = 2)]
    if (length(vec_para_best_final) == 1) {
        vec_p <- c()
    } else {
        vec_p <- vec_para_best_final[seq(2, length(vec_para_best_final) - 1, by = 2)]
    }
    vec_SFS_model <- compute_SFS_one_iteration(vec_A_and_K, vec_p, list_frequencies, library_SFS_component)
    vec_SFS_model <- vec_SFS_model / sum(vec_SFS_model) * sum(vec_SFS_real)
    lines(bar_centers, vec_SFS_model, col = "red", lwd = 2)
    dev.off()
    #---Translation to parameters of cancer evolution in the sample
    deconvolution <- data.frame()
    deconvolution[1, "Total_N"] <- sum(vec_SFS_real)
    deconvolution[1, "Tail"] <- TRUE
    deconvolution[1, "Tail_power"] <- 2
    deconvolution[1, "Tail_mutcount_observed"] <-
        vec_A_and_K[1] * sum(vec_SFS_real)
    deconvolution[1, "Tail_mutcount_predicted"] <-
        vec_A_and_K[1] * sum(vec_SFS_real) *
            sum(library_SFS_component$SFS_exact[[1]]) / sum(library_SFS_component$SFS_expected[[1]]) *
            sample_size / matrix_binomial_sample_size
    deconvolution[1, "Cluster_count"] <- N_humps_best_final
    for (k in 1:N_humps_best_final) {
        deconvolution[1, paste0("Cluster_frequency_", k)] <- vec_p[k] / matrix_binomial_ploidy
        deconvolution[1, paste0("Cluster_mutcount_observed_", k)] <-
            vec_A_and_K[k + 1] * sum(vec_SFS_real)
        deconvolution[1, paste0("Cluster_mutcount_predicted_", k)] <-
            vec_A_and_K[k + 1] * sum(vec_SFS_real) *
                sum(library_SFS_component$SFS_exact[[which(list_frequencies == vec_p[k]) + 1]]) /
                sum(library_SFS_component$SFS_expected[[which(list_frequencies == vec_p[k]) + 1]])
    }
    #---Return the SFS deconvolution results
    output <- list()
    output$vec_para_best_final <- vec_para_best_final
    output$AIC_best_final <- AIC_best_final
    output$vec_SFS_real <- vec_SFS_real
    output$vec_SFS_model <- vec_SFS_model
    output$deconvolution <- deconvolution
    return(output)
}

fit_SFS_given_humpcount_likelihood <- function(vec_SFS_real, vec_p = NULL, N_humps, list_frequencies, library_SFS_component) {
    if (is.null(vec_p)) {
        # 	Choose the hump locations at random
        vec_p <- sort(sample(list_frequencies, N_humps, replace = FALSE), decreasing = TRUE)
    } else {
        vec_p <- sort(vec_p, decreasing = TRUE)
    }
    # 	Function for parameter transformation
    parameter_transform <- function(parameters) {
        vec_A_and_K <- exp(parameters)
        vec_A_and_K <- vec_A_and_K / sum(vec_A_and_K)
        return(vec_A_and_K)
    }
    # 	Function for optimization
    func_fit <- function(parameters) {
        vec_A_and_K <- parameter_transform(parameters)
        loglikelihood <- error_SFS_one_iteration_likelihood(vec_A_and_K, vec_p, list_frequencies, library_SFS_component, vec_SFS_real)
        return(loglikelihood)
    }
    # 	Initial values for parameters
    parameters_initial <- rep(0, N_humps + 1)
    # 	Optimization using optim function
    optim_results <- optim(
        par = parameters_initial,
        fn = func_fit,
        method = "Nelder-Mead",
        control = list(fnscale = -1)
    )
    # 	Extract the results
    parameters <- optim_results$par
    vec_A_and_K <- parameter_transform(parameters)
    log_L <- optim_results$value
    # 	Prepare the parameters to be returned
    vec_para <- numeric(2 * N_humps + 1)
    vec_para[1] <- vec_A_and_K[1]
    for (i in 1:N_humps) {
        vec_para[2 * i] <- vec_p[i]
        vec_para[2 * i + 1] <- vec_A_and_K[i + 1]
    }
    output <- list()
    output$log_L <- log_L
    output$parameters <- vec_para
    return(output)
}

error_SFS_one_iteration_likelihood <- function(vec_A_and_K, vec_p, list_frequencies, library_SFS_component, vec_SFS_real) {
    #----------------Compute the SFS probability distribution from model
    vec_SFS_model <- compute_SFS_one_iteration(vec_A_and_K, vec_p, list_frequencies, library_SFS_component)
    # print(sum(vec_SFS_model))
    vec_SFS_model_normalized <- vec_SFS_model / sum(vec_SFS_model)
    # vec_SFS_model_normalized <- pmax(vec_SFS_model_normalized, 10^-10)
    #-----------------------------Compute the log-likelihood(data|model)
    loglikelihood <- sum(log10(vec_SFS_model_normalized) * vec_SFS_real)
    return(loglikelihood)
}
