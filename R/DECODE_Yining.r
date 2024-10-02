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
    true_Ks <- c(5000)
    true_ps <- c(0.2)
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
    # ============================================================== ABC
    progress_bar <- TRUE
    compute_parallel <- TRUE
    n_cores <- NULL
    N_trials <- 100000
    ABC_percentile <- 0.001
    N_humps <- 1
    with_tail <- TRUE



    func_ABC <- function() {
        alpha_new <- runif(1, neutral_power_min, neutral_power_max)
        ps_new <- sort(runif(N_humps, cluster_frequency_min, cluster_frequency_max), decreasing = FALSE)
        pis_new <- runif(N_humps + with_tail, 0, 1)
        pis_new <- pis_new / sum(pis_new)
        #   compute SFS...
        expected_SFS_iteration <- compute_SFS_Yining(
            pis = pis_new,
            alpha = alpha_new,
            ps = ps_new,
            sfs_bincount = sfs_bincount,
            matrix_binomial_sample_size = matrix_binomial_sample_size,
            SFS_convolution_matrix = SFS_convolution_matrix
        )
        #   compute logLikelihood...
        logLikelihood <- sum(observed_SFS * log(expected_SFS_iteration))
        parameters <- data.frame(matrix(c(alpha_new, ps_new, pis_new), nrow = 1))
        colnames(parameters) <- colnames
        output <- list()
        output$logLikelihood <- logLikelihood
        output$parameters <- parameters
        return(output)
    }


    colnames <- ifelse(with_tail, c("Tail_power"), c())
    colnames <- c(colnames, paste0("Cluster_frequency_", 1:N_humps))
    colnames <- c(colnames, paste0("Pi_", 0:(N_humps + with_tail - 1)))
    if (compute_parallel == FALSE) {
        all_logLikelihood <- c()
        all_parameters <- c()
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
            result <- func_ABC()
            all_logLikelihood <- c(all_logLikelihood, result$logLikelihood)
            if (i == 1) {
                all_parameters <- result$parameters
            } else {
                all_parameters <- rbind(all_parameters, result$parameters)
            }
        }
        if (progress_bar) cat("\n")
    } else {
        suppressPackageStartupMessages(library(parallel))
        suppressPackageStartupMessages(library(pbapply))
        #   Start parallel cluster
        numCores <- ifelse(is.null(n_cores), detectCores(), n_cores)
        cl <- makePSOCKcluster(numCores - 1)
        #   Prepare input parameters
        clusterExport(cl, varlist = c(
            "compute_SFS_Yining", "build_SFS_library_Griffiths_Tavare"
        ), envir = environment())
        #   ...
        if (progress_bar) {
            output <- pblapply(cl = cl, X = 1:N_trials, FUN = function(i) func_ABC())
        } else {
            output <- parLapply(cl = cl, X = 1:N_trials, fun = function(i) func_ABC())
        }
        stopCluster(cl)
        #   Extract the results
        all_logLikelihood <- sapply(output, function(x) x$logLikelihood)
        all_parameters <- do.call(rbind, lapply(output, function(x) x$parameters))
    }
    #   keep top percentile...
    selected_indices <- which(all_logLikelihood >= quantile(all_logLikelihood, probs = 1 - ABC_percentile))
    posterior_logLikelihoods <- all_logLikelihood[selected_indices]
    posterior_parameters <- all_parameters[selected_indices, ]



    # plot for alpha
    alphas <- posterior_parameters$Tail_power
    png(filename = "alpha_posterior.png")
    hist(alphas, main = "Histogram of Selected Alphas", xlab = "Alpha")
    dev.off()
    # plot for each p
    for (i in 1:N_humps) {
        ps_matrix <- posterior_parameters[[paste0("Cluster_frequency_", i)]]
        png(filename = paste0("p_posterior_", i, ".png"))
        hist(ps_matrix, main = paste("Histogram of Selected Ps[", i, "]", sep = ""), xlab = paste("Ps[", i, "]", sep = ""))
        dev.off()
    }
    # plot for each pi
    for (i in 0:N_humps) {
        pis_matrix <- posterior_parameters[[paste0("Pi_", i)]]
        png(filename = paste0("pi_posterior_", i, ".png"))
        hist(pis_matrix, main = paste("Histogram of Selected Pis[", i, "]", sep = ""), xlab = paste("Pis[", i, "]", sep = ""))
        dev.off()
    }
    # # print the posterior mean of each parameter
    # mean_alpha <- mean(sapply(posterior_parameters, function(x) x$alpha))
    # mean_ps <- apply(ps_matrix, 2, mean)
    # mean_pis <- apply(pis_matrix, 2, mean)
    # # print posterior mean
    # print(paste("Mean Alpha:", mean_alpha))
    # print(paste("Mean Ps:", paste(mean_ps, collapse = ", ")))
    # print(paste("Mean Pis:", paste(mean_pis, collapse = ", ")))
    # final_SFS <- compute_SFS_Yining(
    #     pis = mean_pis,
    #     alpha = mean_alpha,
    #     ps = mean_ps,
    #     sfs_bincount = sfs_bincount,
    #     matrix_binomial_sample_size = matrix_binomial_sample_size,
    #     SFS_convolution_matrix = SFS_convolution_matrix
    # )
}
compute_SFS_Yining <- function(pis, alpha, ps, sfs_bincount, matrix_binomial_sample_size, SFS_convolution_matrix, zero_cutoff = 1e-50) {
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

    # png(filename = "expected_SFS_plot.png")
    # plot(1:sfs_bincount, expected_SFS)
    # dev.off()

    return(expected_SFS)
}
