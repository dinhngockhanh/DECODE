fit_SFS_likelihood <- function(mutation_table,
                               criterion = "BIC",
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
                               n_cores = NULL,
                               data_marker_colors) {
    mutation_refcounts <- mutation_table$Ref_count
    mutation_altcounts <- mutation_table$Alt_count
    mutation_totcounts <- mutation_refcounts + mutation_altcounts
    if ("Marker" %in% colnames(mutation_table)) {
        mutation_markers <- mutation_table$Marker
    } else {
        mutation_markers <- c()
    }
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
    list_frequencies_tmp <- sort(unique(c(0, list_frequencies)))
    N_end <<- matrix_binomial_sample_size
    SFS_totalsteps_base <<- matrix_binomial_sfs_stepcount
    matrix_binomial_PDF <<- matrix_binomial_PDF
    library_SFS_component <- list()
    library_SFS_component$SFS_exact <- list()
    library_SFS_component$SFS_expected <- list()
    library_SFS_component$SFS_expected_normalized <- list()
    if (compute_parallel == FALSE) {
        #---------------------------Build SFS library in sequential mode
        pb <- txtProgressBar(
            min = 0,
            max = length(list_frequencies),
            style = 3,
            width = 50,
            char = "="
        )
        for (i in seq_along(list_frequencies_tmp)) {
            setTxtProgressBar(pb, i)
            p <- list_frequencies_tmp[i]
            if (p == 0) {
                vec_para <- c(1)
            } else {
                vec_para <- c(0, p, 1)
            }
            output <- compute_SFS(vec_para)
            library_SFS_component$SFS_exact[[i]] <- output$SFS_exact
            library_SFS_component$SFS_expected[[i]] <- output$SFS_expected
            library_SFS_component$SFS_expected_normalized[[i]] <- output$SFS_expected / sum(output$SFS_expected)
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
        SFS_totalsteps <<- SFS_totalsteps
        r_min <<- r_min
        r_max <<- r_max
        option_dist_coverage <<- option_dist_coverage
        dist_coverage_var_1 <<- dist_coverage_var_1
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
        output <- pblapply(cl = cl, X = 1:length(list_frequencies_tmp), FUN = function(i) {
            p <- list_frequencies_tmp[i]
            if (p == 0) {
                vec_para <- c(1)
            } else {
                vec_para <- c(0, p, 1)
            }
            return(compute_SFS(vec_para))
        })
        for (i in seq_along(list_frequencies_tmp)) {
            library_SFS_component$SFS_exact[[i]] <- output[[i]]$SFS_exact
            library_SFS_component$SFS_expected[[i]] <- output[[i]]$SFS_expected
            library_SFS_component$SFS_expected_normalized[[i]] <- output[[i]]$SFS_expected / sum(output[[i]]$SFS_expected)
        }
        stopCluster(cl)
    }
    #---SFS deconvolution
    compute_AIC <- function(log_L, num_params) {
        return(2 * num_params - 2 * log_L)
    }
    compute_BIC <- function(log_L, num_params, num_samples) {
        return(num_params * log(num_samples) - 2 * log_L)
    }
    N_humps <- 0
    AIC_best_final <- Inf
    BIC_best_final <- Inf
    N_fitting_rounds <- 200
    while (TRUE) {
        num_parameters <- 1 + 2 * N_humps
        AIC_best_current <- Inf
        BIC_best_current <- Inf
        vec_para_AIC_best_current <- c()
        vec_para_BIC_best_current <- c()
        #   Find best parameter set, given the number of humps
        if (N_humps == 0) {
            #   If 0 humps: the SFS consists of only neutral component
            vec_para_AIC_best_current <- 1
            vec_para_BIC_best_current <- 1
            log_L <- error_SFS_one_iteration_likelihood(vec_para_AIC_best_current, c(), list_frequencies, library_SFS_component, vec_SFS_real)
            AIC_current <- compute_AIC(log_L, num_parameters)
            AIC_best_current <- AIC_current
            BIC_current <- compute_BIC(log_L, num_parameters, mutation_count)
            BIC_best_current <- BIC_current
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
                    BIC_current <- compute_BIC(log_L, num_parameters, mutation_count)
                    if (AIC_current < AIC_best_current) {
                        AIC_best_current <- AIC_current
                        vec_para_AIC_best_current <- vec_para
                    }
                    if (BIC_current < BIC_best_current) {
                        BIC_best_current <- BIC_current
                        vec_para_BIC_best_current <- vec_para
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
                        vec_para_AIC_best_current <- vec_para
                    }
                    if (BIC_current < BIC_best_current) {
                        BIC_best_current <- BIC_current
                        vec_para_BIC_best_current <- vec_para
                    }
                }
            }
        }
        #   Report the best fit for the current hump count
        if (criterion == "AIC") {
            report <- paste0(N_humps, " humps: AIC = ", round(AIC_best_current, 3), "; neutral: pi = ", round(vec_para_AIC_best_current[1], 3))
            if (length(vec_para_AIC_best_current) > 1) {
                for (i in 1:N_humps) {
                    report <- paste0(report, "; pi = ", round(vec_para_AIC_best_current[2 * i + 1], 3), " at freq = ", round(vec_para_AIC_best_current[2 * i], 3))
                }
            }
        } else if (criterion == "BIC") {
            report <- paste0(N_humps, " humps: BIC = ", round(BIC_best_current, 3), "; neutral: pi = ", round(vec_para_BIC_best_current[1], 3))
            if (length(vec_para_BIC_best_current) > 1) {
                for (i in 1:N_humps) {
                    report <- paste0(report, "; pi = ", round(vec_para_BIC_best_current[2 * i + 1], 3), " at freq = ", round(vec_para_BIC_best_current[2 * i], 3))
                }
            }
        }
        report <- paste0(report, "\n")
        cat(report)
        #   Check if the increased hump count leads to lower criterion score...
        if (criterion == "AIC") {
            if (AIC_best_current < AIC_best_final) {
                #   ... if yes, then update the best fit and continue with 1 more hump
                AIC_best_final <- AIC_best_current
                vec_para_best_final <- vec_para_AIC_best_current
                N_humps <- N_humps + 1
            } else {
                #   ... if no, then stop
                break
            }
        } else if (criterion == "BIC") {
            if (BIC_best_current < BIC_best_final) {
                #   ... if yes, then update the best fit and continue with 1 more hump
                BIC_best_final <- BIC_best_current
                vec_para_best_final <- vec_para_BIC_best_current
                N_humps <- N_humps + 1
            } else {
                #   ... if no, then stop
                break
            }
        }
    }
    #---Report the best fit
    N_humps_best_final <- (length(vec_para_best_final) - 1) / 2
    if (criterion == "AIC") {
        report <- paste0("\n\n\nBEST FIT:\n", N_humps_best_final, " humps: AIC = ", round(AIC_best_final, 3), "; neutral: pi = ", round(vec_para_best_final[1], 3))
    } else if (criterion == "BIC") {
        report <- paste0("\n\n\nBEST FIT:\n", N_humps_best_final, " humps: BIC = ", round(BIC_best_final, 3), "; neutral: pi = ", round(vec_para_best_final[1], 3))
    }
    if (length(vec_para_best_final) > 1) {
        for (i in 1:N_humps_best_final) {
            report <- paste0(report, "; pi = ", round(vec_para_best_final[2 * i + 1], 3), " at freq = ", round(vec_para_best_final[2 * i], 3))
        }
    }
    report <- paste0(report, "\n\n\n\n")
    cat(report)
    #---Plot the SFS deconvolution
    color_scheme <- c(
        data_marker_colors,
        "Fit: Neutral" = "gray",
        "Fit: Cluster 1" = "red",
        "Fit: Cluster 2" = "blue",
        "Fit: Cluster 3" = "green",
        "Fit: Cluster 4" = "purple"
    )
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
        for (j in 1:length(mutation_altcounts)) {
            no_variant <- mutation_altcounts[j]
            no_total <- mutation_refcounts[j] + mutation_altcounts[j]
            if (no_variant >= min_variant_read && no_total >= min_total_read) {
                VAF <- no_variant / no_total
                pos <- which(vec_freq >= VAF)[1]
                df_data$count[which(df_data$ID == pos & df_data$fill == mutation_markers[j])] <- df_data$count[which(df_data$ID == pos & df_data$fill == mutation_markers[j])] + 1
            }
        }
        df_data$fill <- paste0("Data: ", gsub("_", " ", df_data$fill))
    }
    #   Prepare the deconvolution inference for plotting
    N_humps <- (length(vec_para_best_final) - 1) / 2
    vec_A_and_K <- vec_para_best_final[seq(1, length(vec_para_best_final), by = 2)]
    if (N_humps == 0) {
        vec_p <- c()
    } else {
        vec_p <- vec_para_best_final[seq(2, length(vec_para_best_final) - 1, by = 2)]
        sorted_indices <- order(vec_p, decreasing = TRUE)
        vec_p <- vec_p[sorted_indices]
        vec_A_and_K <- c(vec_A_and_K[1], vec_A_and_K[sorted_indices + 1])
    }
    df_fit <- data.frame()
    SFS_neutral <- mutation_count * vec_A_and_K[1] * library_SFS_component$SFS_expected_normalized[[1]]
    df_fit <- rbind(df_fit, data.frame(frequency = vec_freq, count = SFS_neutral, fill = "Fit: Neutral"))
    if (N_humps > 0) {
        for (i in 1:N_humps) {
            SFS_hump <- mutation_count * vec_A_and_K[i + 1] * library_SFS_component$SFS_expected_normalized[[which(list_frequencies == vec_p[i]) + 1]]
            df_fit <- rbind(df_fit, data.frame(frequency = vec_freq, count = SFS_hump, fill = paste0("Fit: Cluster ", i)))
        }
    }
    filename_plot <- paste0(R_workplace, "/", folder_workplace, "DECONVOLUTION_", n_simulation, ".png")
    png(filename_plot, res = 150, width = 30, height = 15, units = "in")
    p <- ggplot() +
        geom_bar(data = df_data, aes(x = frequency, y = count, fill = fill), stat = "identity") +
        geom_area(data = df_fit, aes(x = frequency, y = count, fill = fill), position = "stack", alpha = 0.5) +
        scale_fill_manual(values = color_scheme, name = "") +
        guides(fill = guide_legend(nrow = 1, keywidth = 2, keyheight = 1)) +
        xlab("Variant Allele Frequency") +
        ylab("Mutation count") +
        theme(
            text = element_text(size = 40),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5)
        )
    print(p)
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
    output$SFS_data <- df_data
    output$SFS_fitted <- df_fit
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
