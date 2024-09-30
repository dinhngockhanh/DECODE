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
    true_A <- 2000000
    true_Ks <- c(5000)
    true_ps <- c(0.5)
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
    #---Test the SFS function...
    pis <- c(0.7, 0.3)
    alpha <- 3.5
    ps <- c(0.4)
    compute_SFS_Yining(
        pis = pis,
        alpha = alpha,
        ps = ps,
        sfs_bincount = sfs_bincount,
        matrix_binomial_sample_size = matrix_binomial_sample_size,
        SFS_convolution_matrix = SFS_convolution_matrix
    )
}
compute_SFS_Yining <- function(pis, alpha, ps, sfs_bincount, matrix_binomial_sample_size, SFS_convolution_matrix) {
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

    png(filename = "expected_SFS_plot.png")
    plot(1:sfs_bincount, expected_SFS)
    dev.off()

    return(expected_SFS)
}
