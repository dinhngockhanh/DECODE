DECODE_Yining <- function(sample_id = "",
                          mutation_table,
                          criterion = "BIC",
                          criterion_ratio = 0.999, # <<<<<<<<<<<<<<<<<<<<<<<<<<
                          neutral_power_min = 0.5,
                          neutral_power_max = 5,
                          cluster_frequency_min = 0.01,
                          cluster_frequency_max = 1,
                          max_total_read = NULL,
                          sample_size = 1000,
                          matrix_binomial_sample_size = 1000, # <<<<<<<<<<<<<<<
                          matrix_binomial_ploidy = 1, # <<<<<<<<<<<<<<<<<<<<<<<
                          sfs_bincount = 100, # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                          inference_retained_freq = 75,
                          validation_mutation_count = 5000,
                          validation_N_trials = 1000,
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
    suppressPackageStartupMessages(library(crayon))
    mutation_table$Tot_count <- mutation_table$Ref_count + mutation_table$Alt_count
    mutation_table$VAF <- mutation_table$Alt_count / mutation_table$Tot_count
    if (is.null(max_total_read)) max_total_read <- max(mutation_table$Tot_count)
    # ================================================ MAKE EXAMPLE DATA
    #---Choose mutation thresholds, get resulting SFS from data
    SFS_data_frequencies <- seq(1, sfs_bincount) / sfs_bincount
    threshold_results <- choose_mutation_thresholds(
        mutation_table = mutation_table,
        max_total_read = max_total_read,
        SFS_data_frequencies = SFS_data_frequencies,
        inference_retained_freq = inference_retained_freq,
        validation_mutation_count = validation_mutation_count,
        validation_N_trials = validation_N_trials
    )
    #---Prepare the SFS convolution matrix
    SFS_convolution_matrix <- build_convolution_matrix(
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
    #---Example unobserved true SFS
    true_SFS <- rep(0, matrix_binomial_sample_size)
    true_alpha <- 2.5
    true_A <- 10000000
    true_Ks <- c(1500, 4000)   #c(5000)
    true_ps <- c(0.2, 0.4) # c(0.2)
    true_SFS <- true_SFS + true_A / (1:matrix_binomial_sample_size)^true_alpha
    for (cluster in 1:length(true_Ks)) {
        true_SFS <- true_SFS + true_Ks[cluster] * dbinom(1:matrix_binomial_sample_size, matrix_binomial_sample_size, true_ps[cluster])
    }
    #---Example observed SFS
    observed_SFS <- rep(0, sfs_bincount)
    for (k in 1:sfs_bincount) {
        observed_SFS[k] <- sum(true_SFS * SFS_convolution_matrix$convolution_matrix[, k])
    }
    png(filename = "observed_SFS_plot.png")
    plot(1:sfs_bincount, observed_SFS)
    dev.off()
    # ============================================ TEST THE SFS FUNCTION
    # alpha <- 2.5 # ~ Uniform(neutral_power_min, neutral_power_max)
    # ps <- c(0.5) # ~ Uniform(cluster_frequency_min, cluster_frequency_max) then sort in decreasing order
    # pis <- c(0.3, 0.7) # ~ Uniform(0, 1) then renormalize
    # expected_SFS <- compute_SFS_Yining(
    #     pis = pis,
    #     alpha = alpha,
    #     ps = ps,
    #     sfs_bincount = sfs_bincount,
    #     matrix_binomial_sample_size = matrix_binomial_sample_size,
    #     SFS_convolution_matrix = SFS_convolution_matrix
    # )

    # ============================================================== MH-MCMC
    N_humps <- 2
    with_tail <- TRUE
    n_samples <- 1000
    N_trials <- 10
    burnin_porp <- 0.2
    progress_bar <- TRUE
    initial_std <- c((neutral_power_max - neutral_power_min) / 500, rep((cluster_frequency_max - cluster_frequency_min) / 500, N_humps), rep(0.01, N_humps + with_tail-1))

    # # dyanmic std
    # target_acceptance_rate <- 0.25
    # adjust_std <- function(std, acceptance_rate, target_acceptance_rate) {
    #     rate_diff <- abs(acceptance_rate - target_acceptance_rate)
    #     if (acceptance_rate > target_acceptance_rate) {
    #         return(std * (1 + rate_diff))
    #     } else {
    #         return(std * (1 - rate_diff))
    #     }
    # }
    
    # define LogLikehood function
    # params = c(alpha, p1, ..., pH, pi0, pi, ..., pi_{H-1})
    burnin <- burnin_porp *n_samples
    log_likelihood_target <- function(params) {
        alpha <- params[1]
        ps <- params[2:(1 + N_humps)]
        pis <- params[(2 + N_humps):length(params)]
        pis <- c(pis, 1-sum(pis))
        
        
        # chech if params are within the domain
        if (alpha < neutral_power_min || alpha > neutral_power_max ||
            any(ps < cluster_frequency_min) || any(ps > cluster_frequency_max) ||
            any(pis < 0) || any(pis > 1)) {
            return(-Inf)
        }
        
        expected_SFS <- compute_SFS_Yining(
                pis = pis,
                alpha = alpha,
                ps = ps,
                sfs_bincount = sfs_bincount,
                matrix_binomial_sample_size = matrix_binomial_sample_size,
                SFS_convolution_matrix = SFS_convolution_matrix
            )
        
        # check if expected_SFS < 0
        if (any(expected_SFS <= 0)) {
            warning("expected_SFS contains non-positive values")
            return(-Inf)
        }
        
        logLikelihood <- sum(observed_SFS * log(expected_SFS))
        
        # check if logLikelihood is NA or NaN
        if (is.na(logLikelihood) || is.nan(logLikelihood)) {
            warning("logLikelihood is NA or NaN")
            cat("params: ", params, "\n")
            return(-Inf)
        }
        
        return(logLikelihood)
    }  


    func_MCMC <- function() {
        library(mcmc)
        library(coda)

        # Initialization of parameters
        alpha_init <- runif(1, neutral_power_min, neutral_power_max)
        ps_init <- sort(runif(N_humps, cluster_frequency_min, cluster_frequency_max), decreasing = FALSE)
        pis_init <- runif(N_humps + with_tail, 0, 1)
        pis_init <- pis_init / sum(pis_init) # Normalize to sum to 1
        params_init <- c(alpha_init, ps_init, pis_init[-length(pis_init)]) # Exclude last value

    # MH-MCMC settings
        acceptance_rate <- 0
        ess <- 0
        count <- 0
        max_iterations <- 50
        MH_samples_all  <- matrix(NA, nrow = 0, ncol = length(params_init))  # Preallocate sample matrix

        while ((acceptance_rate < 0.15 || acceptance_rate > 0.65 ) && count < max_iterations) {
            # || ess < 10
            cat("Re-sampling:", count, "\n")
            count <- count + 1

            MH_result <- metrop(log_likelihood_target, params_init, nbatch = n_samples, scale = initial_std)
            MH_samples_all <- rbind(MH_samples_all, MH_result$batch)
            # Update acceptance rate
            acceptance_rate <- MH_result$accept
            cat("Current Acceptance Rate:", acceptance_rate, "\n")
        
            # Calculate ESS after burn-in
            MH_samples_post_burnin <- MH_samples_all[(nrow(MH_samples_all) - burnin + 1):nrow(MH_samples_all), ]
            mcmc_samples_convergence_check <- mcmc(MH_samples_post_burnin)
            all_ess <- effectiveSize(mcmc_samples_convergence_check)
            ess <- min(all_ess) # Get minimum ESS
            cat("Current ESS:", all_ess, "\n")
            }

        if (count == max_iterations) {
            warning("Maximum iterations reached without convergence.")
        }

        # Final samples post burn-in
        MH_samples <- MH_samples_all[(burnin + 1):nrow(MH_samples_all), ]

        # Adjust pis to ensure sum(pis) = 1
        pi_column <- MH_samples[, (2 + N_humps):ncol(MH_samples)]
        if (is.null(dim(pi_column))) {
            pi_column <- matrix(pi_column, ncol = 1)
        }
        MH_samples <- cbind(MH_samples, 1 - apply(pi_column, 1, sum))

        # order p
        p_col <- 2:(1 + N_humps)
        MH_samples[, p_col] <- t(apply(MH_samples[, p_col], 1, sort))
        print(head(MH_samples))

        # Create a data frame to store the final samples
        parameters <- as.data.frame(MH_samples)
        colnames(parameters) <- colnames # Assumes `colnames` is defined outside

        # Output results as a list
        output <- list(
            logLikelihood = log_likelihood_target(MH_samples),  # Calculate log likelihood
            parameters = parameters,
            acceptance_rate = acceptance_rate,
            MH_samples_all = MH_samples_all,
            ess = all_ess  # Include the final ESS value
        )

        return(output)
    }



    # func_MCMC <- function() {
    #     library(mcmc)
    #     library(coda)

    #     # initialization
    #     alpha_init <- runif(1, neutral_power_min, neutral_power_max)
    #     ps_init <- sort(runif(N_humps, cluster_frequency_min, cluster_frequency_max), decreasing = FALSE)
    #     pis_init <- runif(N_humps + with_tail, 0, 1)
    #     pis_init <- pis_init / sum(pis_init)
    #     params_init <- c(alpha_init, ps_init, pis_init[-length(pis_init)])

    #     #  MH-MCMC
    #     acceptance_rate <- 0
    #     ess <- 0
    #     count <- 0
    #     MH_samples_all  <- matrix(NA, nrow = 0, ncol = 1 + N_humps + with_tail)
    #     while (((acceptance_rate < 0.15) || (acceptance_rate > 0.65) || (ess < 100))) {
    #         cat("Re-sampling:", count, "\n")
    #             count <- count + 1           
    #         if (count > 50) {
    #             break
    #         }
    #         MH_result <- metrop(log_likelihood_target, params_init, nbatch = n_samples, scale = initial_std)
    #         MH_samples_all <- rbind(MH_samples_all, MH_result$batch)
    #         acceptance_rate <- MH_result$accept
    #     }
    #     MH_samples <- MH_samples_all[(nrow(MH_samples_all) - burnin + 1):nrow(MH_samples_all), ]

    #     # update and monitor the acceptance rate
    #     acceptance_rate <- MH_result$accept
    #     print(acceptance_rate)

    #     # check convergence： ess > 100
    #     mcmc_samples_convergence_check <- mcmc(MH_samples)
    #     ess <- effectiveSize(mcmc_samples_convergence_check)
    #     print(ess)
    #     # is_converged <- all(ess > 100)
    #     # is_converged <- TRUE

    #     # save results if converged
    #     # if (is_converged) {
    #     pi_column <- MH_samples[,(2 + N_humps):ncol(MH_samples)]
    #     if (is.null(dim(pi_column))) {
    #         pi_column <- matrix(pi_column, ncol = 1)
    #     }
    #     MH_samples <- cbind(MH_samples, 1-apply(pi_column, 1, sum))

    #     # store the result
    #     parameters <- data.frame(MH_samples, nrow = 1)
    #     colnames(parameters) <- colnames
    #     output <- list()
    #     output$logLikelihood <- log_likelihood_target(MH_samples)
    #     output$parameters <- parameters
    #     output$acceptance_rate <- acceptance_rate
    #     output$MH_samples_all <- MH_samples_all
    #     return(output)
    #     # }    
    # }


    # # ============================================================== ABC
    # progress_bar <- TRUE
    # compute_parallel <- TRUE
    # n_cores <- NULL
    # N_trials <- 100000
    # ABC_percentile <- 0.001
    # N_humps <- 1
    # with_tail <- TRUE



    # func_ABC <- function() {
    #     alpha_new <- runif(1, neutral_power_min, neutral_power_max)
    #     ps_new <- sort(runif(N_humps, cluster_frequency_min, cluster_frequency_max), decreasing = FALSE)
    #     pis_new <- runif(N_humps + with_tail, 0, 1)
    #     pis_new <- pis_new / sum(pis_new)
    #     #   compute SFS...
    #     expected_SFS_iteration <- compute_SFS_Yining(
    #         pis = pis_new,
    #         alpha = alpha_new,
    #         ps = ps_new,
    #         sfs_bincount = sfs_bincount,
    #         matrix_binomial_sample_size = matrix_binomial_sample_size,
    #         SFS_convolution_matrix = SFS_convolution_matrix
    #     )
    #     #   compute logLikelihood...
    #     logLikelihood <- sum(observed_SFS * log(expected_SFS_iteration))
    #     parameters <- data.frame(matrix(c(alpha_new, ps_new, pis_new), nrow = 1))
    #     colnames(parameters) <- colnames
    #     output <- list()
    #     output$logLikelihood <- logLikelihood
    #     output$parameters <- parameters
    #     return(output)
    # }


    colnames <- ifelse(with_tail, c("Tail_power"), c())
    colnames <- c(colnames, paste0("Cluster_frequency_", 1:N_humps))
    colnames <- c(colnames, paste0("Pi_", 0:(N_humps + with_tail - 1)))
    # if (compute_parallel == FALSE) {
    all_logLikelihood <- c()
    all_parameters <- c()
    all_acceptance_rate <- c()
    if (progress_bar) {
        pb <- txtProgressBar(
            min = 0,
            max = N_trials,
            style = 3,
            width = 50,
            char = "+"
        )
    }

    valid <- 0
    while (valid <= N_trials) {
        print(valid)
        if (progress_bar) setTxtProgressBar(pb, valid)
        result <- func_MCMC()
        if (result$acceptance_rate > 0.15 && result$acceptance_rate < 0.65) {
            all_logLikelihood <- c(all_logLikelihood, result$logLikelihood)
            all_acceptance_rate <- c(all_acceptance_rate, result$acceptance_rate)
            if (valid == 0) {
                all_parameters <- result$parameters
            }
            valid <- valid + 1
        }
    }
    
    if (progress_bar) cat("\n")
    # } else {
    #     suppressPackageStartupMessages(library(parallel))
    #     suppressPackageStartupMessages(library(pbapply))
    #     #   Start parallel cluster
    #     numCores <- ifelse(is.null(n_cores), detectCores(), n_cores)
    #     cl <- makePSOCKcluster(numCores - 1)
    #     #   Prepare input parameters
    #     clusterExport(cl, varlist = c(
    #         "compute_SFS_Yining", "build_SFS_library_Griffiths_Tavare"
    #     ), envir = environment())
    #     #   ...
    #     if (progress_bar) {
    #         output <- pblapply(cl = cl, X = 1:N_trials, FUN = function(i) func_MCMC())
    #     } else {
    #         output <- parLapply(cl = cl, X = 1:N_trials, fun = function(i) func_MCMC())
    #     }
    #     stopCluster(cl)
    #     #   Extract the results
    #     all_logLikelihood <- sapply(output, function(x) x$logLikelihood)
    #     all_parameters <- do.call(rbind, lapply(output, function(x) x$parameters))
    #     all_acceptance_rate <- sapply(output, function(x) x$acceptance_rate)
    # }
    # #   keep top percentile...
    # selected_indices <- which(all_logLikelihood >= quantile(all_logLikelihood, probs = 1 - ABC_percentile))
    # posterior_logLikelihoods <- all_logLikelihood[selected_indices]
    # posterior_parameters <- all_parameters[selected_indices, ]


    # moinitor the acceptance rate
    print(all_acceptance_rate)
    png(filename = "acceptance_rate.png")
    plot(all_acceptance_rate, type = "l", xlab = "Iteration", ylab = "Acceptance Rate")
    abline(h = mean(all_acceptance_rate), col = "red", lwd = 2, lty = 2)
    text(x = length(all_acceptance_rate), y = mean(all_acceptance_rate), labels = round(mean(all_acceptance_rate), 2), pos = 3, col = "red")
    dev.off()

    # plot for alpha
    alphas <- all_parameters$Tail_power
    alpha_post_mean <- mean(alphas)
    png(filename = "alpha_posterior.png")
    hist(alphas, main = "Histogram of Selected Alphas", xlab = "Alpha")
    abline(v = alpha_post_mean, col = "red", lwd = 2, lty = 2)
    text(x = alpha_post_mean, y = 0.9 * max(hist(alphas, plot = FALSE)$counts), labels = round(alpha_post_mean, 2), pos = 3, col = "red")
    dev.off()
    # plot for each p
    for (i in 1:N_humps) {
        ps_matrix <- all_parameters[[paste0("Cluster_frequency_", i)]]
        ps_post_mean <- mean(ps_matrix)
        png(filename = paste0("p_posterior_", i, ".png"))
        hist(ps_matrix, main = paste("Histogram of Selected Ps[", i, "]", sep = ""), xlab = paste("Ps[", i, "]", sep = ""))
        abline(v = ps_post_mean, col = "red", lwd = 2, lty = 2)
        text(x = ps_post_mean, y = 0.9 *max(hist(ps_matrix, plot = FALSE)$counts), labels = round(ps_post_mean, 2), pos = 3, col = "red")
        dev.off()
    }
    # plot for each pi
    for (i in 0:N_humps) {
        pis_matrix <- all_parameters[[paste0("Pi_", i)]]
        pi_post_mean <- mean(pis_matrix)
        png(filename = paste0("pi_posterior_", i, ".png"))
        hist(pis_matrix, main = paste("Histogram of Selected Pis[", i, "]", sep = ""), xlab = paste("Pis[", i, "]", sep = ""))
        abline(v = pi_post_mean, col = "red", lwd = 2, lty = 2)
        text(x = pi_post_mean, y = 0.9 * max(hist(pis_matrix, plot = FALSE)$counts), labels = round(pi_post_mean, 2), pos = 3, col = "red")
        dev.off()
    }

    # # print the posterior mean of each parameter
    # mean_alpha <- mean(sapply(posterior_parameters, function(x) x$alpha))
    # mean_ps <- apply(ps_matrix, 2, mean)
    # mean_pis <- apply(pis_matrix, 2, mean)
    # print posterior mean
    # print(paste("Mean Alpha:", mean_alpha))
    # print(paste("Mean Ps:", paste(mean_ps, collapse = ", ")))
    # print(paste("Mean Pis:", paste(mean_pis, collapse = ", ")))
    final_SFS <- compute_SFS_Yining(
        pis = pi_post_mean,
        alpha = alpha_post_mean,
        ps = ps_post_mean,
        sfs_bincount = sfs_bincount,
        matrix_binomial_sample_size = matrix_binomial_sample_size,
        SFS_convolution_matrix = SFS_convolution_matrix,
        plot_expected_SFS = TRUE
    )
}

compute_SFS_Yining <- function(pis, alpha, ps, sfs_bincount, matrix_binomial_sample_size, SFS_convolution_matrix, zero_cutoff = 1e-50, plot_expected_SFS=FALSE) {
    SFS_data_frequencies <- seq(1, sfs_bincount) / sfs_bincount

    matrix_GT_SFS <- c()
    vec_para <- c(1, alpha)
    matrix_GT_SFS <- rbind(matrix_GT_SFS, build_SFS_library_Griffiths_Tavare(vec_para = vec_para, N_end = matrix_binomial_sample_size))
    for (i in 1:length(ps)) {
        vec_para <- c(0, 0, 1, ps[i])
        matrix_GT_SFS <- rbind(matrix_GT_SFS, build_SFS_library_Griffiths_Tavare(vec_para = vec_para, N_end = matrix_binomial_sample_size))
    }

    matrix_expected_SFS <- c()
    for (i in 1:nrow(matrix_GT_SFS)) {
        vec_SFS_GT <- matrix_GT_SFS[i, ]
        vec_SFS_expected <- rep(0, sfs_bincount)
        for (k in 1:sfs_bincount) {
            vec_SFS_expected[k] <- sum(vec_SFS_GT * SFS_convolution_matrix$convolution_matrix[, k])
        }
        matrix_expected_SFS <- rbind(matrix_expected_SFS, vec_SFS_expected)
    }

    expected_SFS <- rep(0, sfs_bincount)
    for (i in 1:nrow(matrix_expected_SFS)) {
        expected_SFS <- expected_SFS + pis[i] * matrix_expected_SFS[i, ] / sum(matrix_expected_SFS[i, ])
    }

    expected_SFS <- expected_SFS + zero_cutoff
    expected_SFS <- expected_SFS / sum(expected_SFS)
     
    if (plot_expected_SFS) {
        png(filename = "expected_SFS_plot.png")
        plot(1:sfs_bincount, expected_SFS)
        dev.off()
    }

    return(expected_SFS)
}
