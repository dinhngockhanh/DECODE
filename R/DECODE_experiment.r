DECODE_plot_experiment <- function(DECODE_result,
                                   fit = "best",
                                   mode = "inference",
                                   data_marker_colors = NULL) {
    library(ggplot2)
    if (is.null(data_marker_colors)) data_marker_colors <- c("Data" = "black")

    vec_freq <- DECODE_result$SFS_frequencies
    SFS_totalsteps <- length(vec_freq)
    mutation_table <- DECODE_result$mutational_table
    mutation_refcounts <- mutation_table$Ref_count
    mutation_altcounts <- mutation_table$Alt_count
    mutation_totcounts <- mutation_refcounts + mutation_altcounts
    tmp <- parameter_conversion_experiment(
        result = DECODE_result$final_fit,
        output_parameters_df = FALSE
    )
    if (mode == "inference") {
        vec_SFS_real <- DECODE_result$SFS_for_fitting
        min_variant_read <- DECODE_result$min_variant_read_inference
        min_total_read <- DECODE_result$min_total_read_inference
        if (fit == "best") {
            vec_para_best_final <- DECODE_result$final_fit$best_fit$parameters_inference
            tail_status_final <- DECODE_result$final_fit$best_fit$tail_status
            component_distributions_best_final <- DECODE_result$final_fit$best_fit$component_distributions_inference
        }
        vec_A <- tmp$inference$vec_A
        vec_K <- tmp$inference$vec_K
        vec_p <- tmp$inference$vec_p
        N_humps <- tmp$inference$N_humps
        tail_status <- tmp$inference$tail_status
    } else if (mode == "validation") {
        vec_SFS_real <- DECODE_result$SFS_for_validating
        min_variant_read <- DECODE_result$min_variant_read_validation
        min_total_read <- DECODE_result$min_total_read_validation
        if (fit == "best") {
            vec_para_best_final <- DECODE_result$final_fit$best_fit$parameters_validation
            tail_status_final <- DECODE_result$final_fit$best_fit$tail_status
            component_distributions_best_final <- DECODE_result$final_fit$best_fit$component_distributions_validation
        }
        vec_A <- tmp$validation$vec_A
        vec_K <- tmp$validation$vec_K
        vec_p <- tmp$validation$vec_p
        N_humps <- tmp$validation$N_humps
        tail_status <- tmp$validation$tail_status
    }
    if ("Marker" %in% colnames(mutation_table)) {
        mutation_markers <- mutation_table$Marker
    } else {
        mutation_markers <- c()
    }
    #---Plot the SFS deconvolution
    color_scheme <- c(
        data_marker_colors,
        "Neutral tail" = "#999999",
        "Cluster 1" = "#D55E00",
        "Cluster 2" = "#009E73",
        "Cluster 3" = "#0072B2",
        "Cluster 4" = "#E69F00",
        "Cluster 5" = "#56B4E9",
        "Cluster 6" = "#CC79A7",
        "Cluster 7" = "#F0E442"
    )
    mutation_count <- sum(vec_SFS_real)
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
        df_data$fill <- gsub("_", " ", df_data$fill)
    }
    #   Prepare the deconvolution inference for plotting
    df_fit <- data.frame()
    df_fit_order <- c()
    if (tail_status_final) {
        SFS_neutral <- compute_SFS(
            A = vec_A[1],
            vec_K = c(),
            component_distributions = component_distributions_best_final
        )
        SFS_neutral <- vec_A[1] * mutation_count * SFS_neutral / sum(SFS_neutral)
        df_fit <- rbind(df_fit, data.frame(frequency = vec_freq, count = SFS_neutral, fill = "Neutral tail"))
        df_fit_order <- c(df_fit_order, "Neutral tail")
    }
    if (N_humps > 0) {
        for (i in 1:N_humps) {
            vec_K_tmp <- rep(0, N_humps)
            vec_K_tmp[i] <- vec_K[i]
            SFS_hump <- compute_SFS(
                A = 0,
                vec_K = vec_K_tmp,
                component_distributions = component_distributions_best_final
            )
            SFS_hump <- vec_K[i] * mutation_count * SFS_hump / sum(SFS_hump)
            df_fit <- rbind(df_fit, data.frame(frequency = vec_freq, count = SFS_hump, fill = paste0("Cluster ", i)))
            df_fit_order <- c(df_fit_order, paste0("Cluster ", i))
        }
    }
    df_fit$fill <- factor(df_fit$fill, levels = rev(df_fit_order))
    df_data$fill <- factor(df_data$fill, levels = rev(names(data_marker_colors)))
    p <- ggplot() +
        geom_area(data = df_fit, aes(x = frequency, y = count, fill = fill), position = "stack", alpha = 0.8) +
        geom_bar(data = df_data, aes(x = frequency, y = count, fill = fill), stat = "identity", width = 0.5 / SFS_totalsteps) +
        scale_fill_manual(values = color_scheme, name = "") +
        guides(fill = guide_legend(nrow = 1, keywidth = 2, keyheight = 1)) +
        xlab("Variant Allele Frequency") +
        ylab("Mutation count") +
        theme(
            text = element_text(size = 30),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5)
        )
    return(p)
}

DECODE_plot_readcounts <- function(DECODE_result,
                                   freq_cutoff = 10,
                                   compute_parallel = TRUE,
                                   n_cores = NULL) {
    library(ggplot2)
    library(shadowtext)
    library(reshape2)
    mutation_table <- DECODE_result$mutational_table
    mutation_table$Tot_count <- mutation_table$Ref_count + mutation_table$Alt_count
    min_variant_read_inference <- max(DECODE_result$min_variant_read_inference, min(mutation_table$Alt_count))
    min_total_read_inference <- max(DECODE_result$min_total_read_inference, min(mutation_table$Tot_count))
    min_variant_read_validation <- max(DECODE_result$min_variant_read_validation, min(mutation_table$Alt_count))
    min_total_read_validation <- max(DECODE_result$min_total_read_validation, min(mutation_table$Tot_count))
    #---Find joint distribution of variant and total readcounts
    vec_min_variant_read <- min(mutation_table$Alt_count):max(mutation_table$Alt_count)
    vec_min_total_read <- min(mutation_table$Tot_count):max(mutation_table$Tot_count)
    func_dist_given_min_variant_read <- function(min_variant_read, mutation_table) {
        df <- data.frame()
        for (min_total_read in vec_min_total_read) {
            df <- rbind(
                df,
                data.frame(
                    min_total_read = min_total_read,
                    min_variant_read = min_variant_read,
                    freq = 100 * sum(mutation_table$Alt_count >= min_variant_read & mutation_table$Tot_count >= min_total_read) / nrow(mutation_table)
                )
            )
        }
        return(df)
    }
    if (compute_parallel == FALSE) {
        df_dist <- data.frame()
        pb <- txtProgressBar(
            min = 1,
            max = length(vec_min_variant_read),
            style = 3,
            width = 50,
            char = "+"
        )
        for (min_variant_read in vec_min_variant_read) {
            df_dist <- rbind(
                df_dist,
                func_dist_given_min_variant_read(min_variant_read, mutation_table)
            )
            setTxtProgressBar(pb, min_variant_read)
        }
    } else {
        library(parallel)
        library(pbapply)
        #   Start parallel cluster
        numCores <- ifelse(is.null(n_cores), detectCores(), n_cores)
        cl <- makePSOCKcluster(numCores - 1)
        #   Prepare input parameters
        clusterExport(cl, varlist = c("mutation_table"))
        #   Compute each sub-dataframe
        output <- pblapply(cl = cl, X = vec_min_variant_read, FUN = function(min_variant_read) {
            return(func_dist_given_min_variant_read(min_variant_read, mutation_table))
        })
        stopCluster(cl)
        df_dist <- do.call(rbind, output)
    }
    #---Reduce the distribution to region satisfying the frequency cutoff
    max_total_read <- vec_min_total_read[which(df_dist$freq[which(df_dist$min_variant_read == vec_min_variant_read[1])] < freq_cutoff)[1]]
    max_variant_read <- vec_min_variant_read[which(df_dist$freq[which(df_dist$min_total_read == vec_min_total_read[1])] < freq_cutoff)[1]]
    df_dist <- df_dist[which(df_dist$min_total_read <= max_total_read & df_dist$min_variant_read <= max_variant_read), ]
    #---Plot the readcount distribution
    freq_inference <- round(df_dist$freq[df_dist$min_total_read == min_total_read_inference & df_dist$min_variant_read == min_variant_read_inference], 2)
    freq_validation <- round(df_dist$freq[df_dist$min_total_read == min_total_read_validation & df_dist$min_variant_read == min_variant_read_validation], 2)
    text_inference <- paste0("Inference (", freq_inference, "% mutations retained)")
    text_validation <- paste0("Validation (", freq_validation, "% mutations retained)")
    color_inference <- "#DF536B"
    color_validation <- "#2297E6"
    p <- ggplot(df_dist, aes(x = min_total_read, y = min_variant_read, fill = freq)) +
        geom_tile() +
        geom_rect(
            aes(
                xmin = min_total_read_inference - 0.5, xmax = min_total_read_inference + 0.5,
                ymin = min_variant_read_inference - 0.5, ymax = min_variant_read_inference + 0.5
            ),
            fill = NA, color = "white", size = 1
        ) +
        geom_shadowtext(
            aes(
                x = min_total_read_inference + 2,
                y = min_variant_read_inference,
                label = text_inference
            ),
            angle = 45,
            hjust = 0,
            vjust = 0,
            size = 8,
            color = color_inference,
            bg.color = "white",
            fontface = "bold"
        ) +
        geom_rect(
            aes(
                xmin = min_total_read_validation - 0.5, xmax = min_total_read_validation + 0.5,
                ymin = min_variant_read_validation - 0.5, ymax = min_variant_read_validation + 0.5
            ),
            fill = NA, color = "white", size = 1
        ) +
        geom_shadowtext(
            aes(
                x = min_total_read_validation + 2,
                y = min_variant_read_validation,
                label = text_validation
            ),
            angle = 45,
            hjust = 0,
            vjust = 0,
            size = 8,
            color = color_validation,
            bg.color = "white",
            fontface = "bold"
        ) +
        scale_fill_gradientn(
            colors = c("#0072B2", "#56B4E9", "#009E73", "#E69F00", "#D55E00"),
            name = "% retained mutations"
        ) +
        theme_minimal() +
        labs(title = "", x = "Minimum total readcount", y = "Minimum alternative readcount") +
        theme(
            text = element_text(size = 30),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5),
            legend.key.width = unit(3.5, "cm"),
        )
    return(p)
}

DECODE_experiment <- function(sample_id = "",
                              mutation_table,
                              criterion = "BIC",
                              criterion_ratio = 1,
                              neutral_power_min = 0.5,
                              neutral_power_max = 5,
                              cluster_frequency_min = 0.01,
                              cluster_frequency_max = 1,
                              max_total_read = NULL,
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
    ####################################################################
    ####################################################################
    ####################################################################
    libPaths_binomial_table <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/DECODE_binomial_matrices"
    matrix_binomial_sample_size <- 1000
    matrix_binomial_sfs_bincount <- sfs_bincount
    matrix_binomial_ploidy <- 2
    min_variant_read_inference <- 4
    min_total_read_inference <- 50
    min_variant_read_validation <- 7
    min_total_read_validation <- 0
    if (is.null(max_total_read)) max_total_read <- max(mutation_totcounts)
    ####################################################################
    ####################################################################
    ####################################################################
    #---Prepare the real SFS
    vec_freq <- seq(1, sfs_bincount) / sfs_bincount
    vec_SFS_real_inference <- rep(0, sfs_bincount)
    for (j in 1:sfs_bincount) {
        vec_SFS_real_inference[j] <- length(which(
            mutation_altcounts >= min_variant_read_inference &
                mutation_totcounts >= min_total_read_inference &
                mutation_vaf >= ifelse(j == 1, 0, vec_freq[j - 1]) &
                mutation_vaf < vec_freq[j]
        ))
    }
    vec_SFS_real_validation <- rep(0, sfs_bincount)
    for (j in 1:sfs_bincount) {
        vec_SFS_real_validation[j] <- length(which(
            mutation_altcounts >= min_variant_read_validation &
                mutation_totcounts >= min_total_read_validation &
                mutation_vaf >= ifelse(j == 1, 0, vec_freq[j - 1]) &
                mutation_vaf < vec_freq[j]
        ))
    }
    #---Prepare the total readcount distribution
    sample_coverage_inference <- prep_distribution_patient(
        mutations_total_read = mutation_totcounts,
        min_total_read = min_total_read_inference,
        max_total_read = max_total_read
    )
    sample_coverage_validation <- prep_distribution_patient(
        mutations_total_read = mutation_totcounts,
        min_total_read = min_total_read_validation,
        max_total_read = max_total_read
    )
    #---Get DECODE binomial table
    binomial_matrix_inference <- get_binomial_matrix(
        folder = libPaths_binomial_table,
        matrix_binomial_sample_size = matrix_binomial_sample_size,
        matrix_binomial_sfs_bincount = matrix_binomial_sfs_bincount,
        matrix_binomial_ploidy = matrix_binomial_ploidy,
        min_variant_read = min_variant_read_inference,
        min_total_read = min_total_read_inference,
        max_total_read = max_total_read,
        compute_parallel = compute_parallel,
        n_cores = n_cores
    )
    binomial_matrix_validation <- get_binomial_matrix(
        folder = libPaths_binomial_table,
        matrix_binomial_sample_size = matrix_binomial_sample_size,
        matrix_binomial_sfs_bincount = matrix_binomial_sfs_bincount,
        matrix_binomial_ploidy = matrix_binomial_ploidy,
        min_variant_read = min_variant_read_validation,
        min_total_read = min_total_read_validation,
        max_total_read = max_total_read,
        compute_parallel = compute_parallel,
        n_cores = n_cores
    )
    #---Prepare the SFS convolution matrix
    cat(bold(blue("Prepare the SFS convolution matrix...\n")))
    SFS_convolution_inference <- build_convolution_matrix(
        binomial_matrix = binomial_matrix_inference,
        sfs_bincount = sfs_bincount,
        coverage_distribution = coverage_distribution,
        coverage_variables = coverage_variables,
        sample_coverage = sample_coverage_inference
    )
    SFS_convolution_validation <- build_convolution_matrix(
        binomial_matrix = binomial_matrix_validation,
        sfs_bincount = sfs_bincount,
        coverage_distribution = coverage_distribution,
        coverage_variables = coverage_variables,
        sample_coverage = sample_coverage_validation
    )
    #---DECODE
    DECODE_result <- list()
    DECODE_result$sample_id <- sample_id
    DECODE_result$sfs_bincount <- sfs_bincount
    DECODE_result$mutational_table <- mutation_table
    DECODE_result$min_variant_read_inference <- min_variant_read_inference
    DECODE_result$min_total_read_inference <- min_total_read_inference
    DECODE_result$min_variant_read_validation <- min_variant_read_validation
    DECODE_result$min_total_read_validation <- min_total_read_validation
    DECODE_result$max_total_read <- max_total_read
    DECODE_result$SFS_frequencies <- vec_freq
    DECODE_result$SFS_for_fitting <- vec_SFS_real_inference
    DECODE_result$SFS_for_validating <- vec_SFS_real_validation
    if (is.na(neutral_tail)) {
        # result_with_tail <- DECODE_given_tail_status(
        #     vec_SFS_real = vec_SFS_real,
        #     sfs_bincount = sfs_bincount,
        #     with_tail = TRUE,
        #     criterion = criterion,
        #     criterion_ratio = criterion_ratio,
        #     min_N_humps = min_N_humps,
        #     max_N_humps = max_N_humps,
        #     N_trials = N_trials,
        #     SFS_convolution = SFS_convolution,
        #     neutral_power_min = neutral_power_min,
        #     neutral_power_max = neutral_power_max,
        #     cluster_frequency_min = cluster_frequency_min,
        #     cluster_frequency_max = cluster_frequency_max,
        #     pi_cutoff = pi_cutoff,
        #     zero_cutoff = zero_cutoff,
        #     compute_parallel = compute_parallel,
        #     n_cores = n_cores
        # )
        # DECODE_result$fits_with_tail <- result_with_tail
        # result_without_tail <- DECODE_given_tail_status(
        #     vec_SFS_real = vec_SFS_real,
        #     sfs_bincount = sfs_bincount,
        #     with_tail = FALSE,
        #     criterion = criterion,
        #     criterion_ratio = criterion_ratio,
        #     min_N_humps = min_N_humps,
        #     max_N_humps = max_N_humps,
        #     N_trials = N_trials,
        #     SFS_convolution = SFS_convolution,
        #     neutral_power_min = neutral_power_min,
        #     neutral_power_max = neutral_power_max,
        #     cluster_frequency_min = cluster_frequency_min,
        #     cluster_frequency_max = cluster_frequency_max,
        #     pi_cutoff = pi_cutoff,
        #     zero_cutoff = zero_cutoff,
        #     compute_parallel = compute_parallel,
        #     n_cores = n_cores
        # )
        # DECODE_result$fits_without_tail <- result_without_tail
        # if (result_with_tail$best_fit$selected_criterion_value < result_without_tail$best_fit$selected_criterion_value) {
        #     final_result <- result_with_tail
        # } else {
        #     final_result <- result_without_tail
        # }
    } else if (neutral_tail == TRUE) {
        final_result <- DECODE_given_tail_status_experiment(
            vec_SFS_real_inference = vec_SFS_real_inference,
            vec_SFS_real_validation = vec_SFS_real_validation,
            sfs_bincount = sfs_bincount,
            with_tail = TRUE,
            criterion = criterion,
            criterion_ratio = criterion_ratio,
            min_N_humps = min_N_humps,
            max_N_humps = max_N_humps,
            N_trials = N_trials,
            SFS_convolution_inference = SFS_convolution_inference,
            SFS_convolution_validation = SFS_convolution_validation,
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
        # final_result <- DECODE_given_tail_status(
        #     vec_SFS_real = vec_SFS_real,
        #     sfs_bincount = sfs_bincount,
        #     with_tail = FALSE,
        #     criterion = criterion,
        #     criterion_ratio = criterion_ratio,
        #     min_N_humps = min_N_humps,
        #     max_N_humps = max_N_humps,
        #     N_trials = N_trials,
        #     SFS_convolution = SFS_convolution,
        #     neutral_power_min = neutral_power_min,
        #     neutral_power_max = neutral_power_max,
        #     cluster_frequency_min = cluster_frequency_min,
        #     cluster_frequency_max = cluster_frequency_max,
        #     pi_cutoff = pi_cutoff,
        #     zero_cutoff = zero_cutoff,
        #     compute_parallel = compute_parallel,
        #     n_cores = n_cores
        # )
        # DECODE_result$fits_without_tail <- final_result
    }
    DECODE_result$final_fit <- final_result
    #---Report the best fit
    tail_status_final_result <- final_result$best_fit$tail_status
    parameters_inference_final_result <- final_result$best_fit$parameters_inference
    criterion_final_result <- final_result$best_fit$selected_criterion_value
    if (tail_status_final_result) {
        N_humps_final_result <- length(parameters_inference_final_result) / 2 - 1
        report <- bold(underline(red(paste0("Best fit = neutral tail + ", N_humps_final_result, " clusters:\n"))))
        report <- paste0(report, red("Score            : "), yellow(paste0(criterion, " = ", round(criterion_final_result, 3))), "\n")
        report <- paste0(report, red("Neutral component: "), yellow(paste0("pi = ", round(parameters_inference_final_result[1], 3))), red(", "), yellow(paste0("power = ", round(parameters_inference_final_result[2], 3))), "\n")
        ii <- 0
    } else {
        N_humps_final_result <- length(parameters_inference_final_result) / 2
        report <- bold(underline(red(paste0("Best fit = no neutral tail + ", N_humps_final_result, " clusters:\n"))))
        report <- paste0(report, red("Score            : "), yellow(paste0(criterion, " = ", round(criterion_final_result, 3))), "\n")
        ii <- -1
    }
    if (N_humps_final_result > 0) {
        for (i in 1:N_humps_final_result) {
            report <- paste0(report, red(paste0("Cluster ", i, "        : ")), yellow(paste0("pi = ", round(parameters_inference_final_result[2 * (i + ii) + 1], 3))), red(", "), yellow(paste0("f = ", round(parameters_inference_final_result[2 * (i + ii) + 2], 3))), "\n")
        }
    }
    cat(report)
    #---Translation to parameters of cancer evolution in the sample
    tmp <- parameter_conversion_experiment(
        result = final_result,
        mutation_count_for_fitting = sum(vec_SFS_real_inference),
        sample_size = sample_size,
        matrix_binomial_sample_size = matrix_binomial_sample_size,
        matrix_binomial_ploidy = matrix_binomial_ploidy
    )
    DECODE_result$final_fit$parameters_df <- tmp$parameters_df
    #---Return the SFS deconvolution results
    return(DECODE_result)
}

DECODE_given_tail_status_experiment <- function(vec_SFS_real_inference,
                                                vec_SFS_real_validation,
                                                criterion,
                                                criterion_ratio,
                                                min_N_humps,
                                                max_N_humps,
                                                with_tail,
                                                N_trials,
                                                sfs_bincount,
                                                SFS_convolution_inference,
                                                SFS_convolution_validation,
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
        fit_results <- DECODE_given_tail_status_and_Ncluster_experiment(
            vec_SFS_real_inference = vec_SFS_real_inference,
            vec_SFS_real_validation = vec_SFS_real_validation,
            N_humps = N_humps,
            with_tail = with_tail,
            N_trials = N_trials,
            sfs_bincount = sfs_bincount,
            SFS_convolution_inference = SFS_convolution_inference,
            SFS_convolution_validation = SFS_convolution_validation,
            neutral_power_min = neutral_power_min,
            neutral_power_max = neutral_power_max,
            cluster_frequency_min = cluster_frequency_min,
            cluster_frequency_max = cluster_frequency_max,
            zero_cutoff = zero_cutoff,
            compute_parallel = compute_parallel,
            n_cores = n_cores
        )
        all_fits[[paste0(N_humps, "_clusters")]] <- fit_results
        parameters_inference_best_current <- fit_results$best$parameters_inference
        parameters_validation_best_current <- fit_results$best$parameters_validation
        component_distributions_inference_best_current <- fit_results$best$component_distributions_inference
        component_distributions_validation_best_current <- fit_results$best$component_distributions_validation
        criterion_all_best_current <- fit_results$best$criteria
        criterion_best_current <- criterion_all_best_current[[criterion]]
        #   Report the best fit for the current hump count
        cluster_pis <- Inf
        report <- paste0(blue("Score            : "), cyan(paste0(criterion, " = ", round(criterion_best_current, 3))), "\n")
        if (with_tail) {
            N_humps <- length(parameters_inference_best_current) / 2 - 1
            report <- paste0(report, blue("Neutral component: "), cyan(paste0("pi = ", round(parameters_inference_best_current[1], 3))), blue(", "), cyan(paste0("power = ", round(parameters_inference_best_current[2], 3))), "\n")
            ii <- 0
        } else {
            N_humps <- length(parameters_inference_best_current) / 2
            ii <- -1
        }
        if (N_humps > 0) {
            for (i in 1:N_humps) {
                report <- paste0(report, blue(paste0("Cluster ", i, "        : ")), cyan(paste0("pi = ", round(parameters_inference_best_current[2 * (i + ii) + 1], 3))), blue(", "), cyan(paste0("f = ", round(parameters_inference_best_current[2 * (i + ii) + 2], 3))), "\n")
                cluster_pis <- c(cluster_pis, parameters_inference_best_current[2 * (i + ii) + 1])
            }
        }
        cat(report)
        #   Check if the increased hump count leads to lower criterion score without tiny selective components...
        if ((N_humps == min_N_humps) | ((criterion_best_current < criterion_ratio * criterion_best_final) & (min(cluster_pis) >= pi_cutoff))) {
            #   ... if yes, then update the best fit and continue with 1 more hump
            fit_results_best_final <- fit_results
            criterion_best_final <- criterion_best_current
            criterion_all_final <- criterion_all_best_current
            parameters_inference_best_final <- parameters_inference_best_current
            parameters_validation_best_final <- parameters_validation_best_current
            component_distributions_inference_best_final <- component_distributions_inference_best_current
            component_distributions_validation_best_final <- component_distributions_validation_best_current
            N_humps <- N_humps + 1
            #   ... except if exceeding maximum number of clusters
            if (N_humps > max_N_humps) break
        } else {
            #   ... if no, then stop
            break
        }
    }
    #---Check if the neutral tail component is too tiny
    if (with_tail == TRUE & parameters_inference_best_final[1] < pi_cutoff) {
        with_tail <- FALSE
        parameters_inference_best_final[seq(3, length(parameters_inference_best_final), by = 2)] <- parameters_inference_best_final[seq(3, length(parameters_inference_best_final), by = 2)] / sum(parameters_inference_best_final[seq(3, length(parameters_inference_best_final), by = 2)])
        parameters_inference_best_final <- parameters_inference_best_final[-c(1, 2)]
        parameters_validation_best_final[seq(3, length(parameters_validation_best_final), by = 2)] <- parameters_validation_best_final[seq(3, length(parameters_validation_best_final), by = 2)] / sum(parameters_validation_best_final[seq(3, length(parameters_validation_best_final), by = 2)])
        parameters_validation_best_final <- parameters_validation_best_final[-c(1, 2)]
        component_distributions_inference_best_final$SFS_exact[1, ] <- rep(0, length(component_distributions_inference_best_final$SFS_exact[1, ]))
        component_distributions_inference_best_final$SFS_expected[1, ] <- rep(0, length(component_distributions_inference_best_final$SFS_expected[1, ]))
        component_distributions_inference_best_final$SFS_expected_normalized[1, ] <- rep(0, length(component_distributions_inference_best_final$SFS_expected_normalized[1, ]))
        component_distributions_validation_best_final$SFS_exact[1, ] <- rep(0, length(component_distributions_validation_best_final$SFS_exact[1, ]))
        component_distributions_validation_best_final$SFS_expected[1, ] <- rep(0, length(component_distributions_validation_best_final$SFS_expected[1, ]))
        component_distributions_validation_best_final$SFS_expected_normalized[1, ] <- rep(0, length(component_distributions_validation_best_final$SFS_expected_normalized[1, ]))
    }
    #---Report the best fit
    result <- list()
    result$all_fits <- all_fits
    result$best_fit <- list()
    result$best_fit$parameters_inference <- parameters_inference_best_final
    result$best_fit$parameters_validation <- parameters_validation_best_final
    result$best_fit$component_distributions_inference <- component_distributions_inference_best_final
    result$best_fit$component_distributions_validation <- component_distributions_validation_best_final
    result$best_fit$selected_criterion <- criterion
    result$best_fit$all_criteria <- criterion_all_final
    result$best_fit$selected_criterion_value <- criterion_best_final
    result$best_fit$tail_status <- with_tail
    return(result)
}

DECODE_given_tail_status_and_Ncluster_experiment <- function(vec_SFS_real_inference,
                                                             vec_SFS_real_validation,
                                                             N_humps,
                                                             with_tail,
                                                             N_trials,
                                                             sfs_bincount,
                                                             SFS_convolution_inference,
                                                             SFS_convolution_validation,
                                                             neutral_power_min,
                                                             neutral_power_max,
                                                             cluster_frequency_min,
                                                             cluster_frequency_max,
                                                             zero_cutoff,
                                                             compute_criteria = TRUE,
                                                             compute_parallel,
                                                             n_cores,
                                                             progress_bar = TRUE) {
    N_end <- SFS_convolution_inference$N_end
    SFS_convolution_matrix_inference <- SFS_convolution_inference$convolution_matrix
    SFS_convolution_matrix_validation <- SFS_convolution_validation$convolution_matrix
    #---Function to perform one trial to find A & K's
    func_one_trial <- function(with_tail,
                               compute_criteria,
                               N_humps,
                               neutral_power_min,
                               neutral_power_max,
                               cluster_frequency_min,
                               cluster_frequency_max,
                               vec_SFS_real_inference,
                               vec_SFS_real_validation,
                               sfs_bincount,
                               N_end,
                               SFS_convolution_matrix_inference,
                               SFS_convolution_matrix_validation,
                               zero_cutoff) {
        #   Sample neutral component power and cluster frequencies
        neutral_power <- ifelse(with_tail, runif(1, neutral_power_min, neutral_power_max), NA)
        cluster_frequencies <- sort(runif(N_humps, cluster_frequency_min, cluster_frequency_max), decreasing = TRUE)
        #   Build the SFS component library
        component_distributions_inference <- build_SFS_library(
            neutral_power = neutral_power,
            cluster_frequencies = cluster_frequencies,
            sfs_bincount = sfs_bincount,
            SFS_convolution_matrix = SFS_convolution_matrix_inference,
            N_end = N_end
        )
        component_distributions_validation <- build_SFS_library(
            neutral_power = neutral_power,
            cluster_frequencies = cluster_frequencies,
            sfs_bincount = sfs_bincount,
            SFS_convolution_matrix = SFS_convolution_matrix_validation,
            N_end = N_end
        )
        #   Find component proportions
        results <- DECODE_for_pis(
            vec_SFS_real = vec_SFS_real_inference,
            neutral_power = neutral_power,
            cluster_frequencies = cluster_frequencies,
            component_distributions = component_distributions_inference,
            zero_cutoff = zero_cutoff
        )
        if (compute_criteria) {
            mutation_count <- sum(vec_SFS_real_inference)
            num_parameters <- ifelse(with_tail, 2 * N_humps + 1, 2 * N_humps - 1)
            if (neutral_power_max > neutral_power_min) num_parameters <- num_parameters + 1
            criteria <- cluster_count_criteria(
                num_parameters = num_parameters,
                log_L = results$log_L,
                num_samples = mutation_count,
                vec_SFS_real = vec_SFS_real_inference,
                parameters = results$parameters,
                component_distributions = component_distributions_inference,
                zero_cutoff = zero_cutoff
            )
        }
        output <- list()
        output$parameters_inference <- results$parameters
        output$parameters_validation <- results$parameters
        if (with_tail) {
            for (i in 1:nrow(component_distributions_validation$SFS_expected)) {
                output$parameters_validation[2 * i - 1] <- output$parameters_validation[2 * i - 1] *
                    sum(component_distributions_validation$SFS_expected[i, ]) /
                    sum(component_distributions_inference$SFS_expected[i, ])
            }
        } else {
            for (i in 2:nrow(component_distributions_validation$SFS_expected)) {
                output$parameters_validation[2 * i - 3] <- output$parameters_validation[2 * i - 3] *
                    sum(component_distributions_validation$SFS_expected[i, ]) /
                    sum(component_distributions_inference$SFS_expected[i, ])
            }
        }
        output$parameters_validation[seq(1, length(output$parameters_validation), by = 2)] <-
            output$parameters_validation[seq(1, length(output$parameters_validation), by = 2)] /
                sum(output$parameters_validation[seq(1, length(output$parameters_validation), by = 2)])
        output$logLikelihood <- results$log_L
        output$component_distributions_inference <- component_distributions_inference
        output$component_distributions_validation <- component_distributions_validation
        if (compute_criteria) output$criteria <- criteria
        return(output)
    }
    #---Find best variable parameters (A & K's) for each fixed parameter set from many trials
    if (compute_parallel == FALSE) {
        all_logLikelihood <- c()
        all_parameters_inference <- c()
        all_parameters_validation <- c()
        all_component_distributions_inference <- list()
        all_component_distributions_validation <- list()
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
                vec_SFS_real_inference = vec_SFS_real_inference,
                vec_SFS_real_validation = vec_SFS_real_validation,
                sfs_bincount = sfs_bincount,
                N_end = N_end,
                SFS_convolution_matrix_inference = SFS_convolution_matrix_inference,
                SFS_convolution_matrix_validation = SFS_convolution_matrix_validation,
                zero_cutoff = zero_cutoff
            )
            all_parameters_inference <- rbind(all_para, trial_result$parameters_inference)
            all_parameters_validation <- rbind(all_para, trial_result$parameters_validation)
            all_component_distributions_inference[[i]] <- trial_result$component_distributions_inference
            all_component_distributions_validation[[i]] <- trial_result$component_distributions_validation
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
            "with_tail", "N_humps", "sfs_bincount", "N_end", "zero_cutoff", "compute_criteria",
            "vec_SFS_real_inference", "vec_SFS_real_validation",
            "SFS_convolution_matrix_inference", "SFS_convolution_matrix_validation",
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
                    vec_SFS_real_inference = vec_SFS_real_inference,
                    vec_SFS_real_validation = vec_SFS_real_validation,
                    sfs_bincount = sfs_bincount,
                    N_end = N_end,
                    SFS_convolution_matrix_inference = SFS_convolution_matrix_inference,
                    SFS_convolution_matrix_validation = SFS_convolution_matrix_validation,
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
                    vec_SFS_real_inference = vec_SFS_real_inference,
                    vec_SFS_real_validation = vec_SFS_real_validation,
                    sfs_bincount = sfs_bincount,
                    N_end = N_end,
                    SFS_convolution_matrix_inference = SFS_convolution_matrix_inference,
                    SFS_convolution_matrix_validation = SFS_convolution_matrix_validation,
                    zero_cutoff = zero_cutoff
                )
            })
        }
        stopCluster(cl)
        #   Extract the results
        all_parameters_inference <- do.call(rbind, lapply(output, function(x) x$parameters_inference))
        all_parameters_validation <- do.call(rbind, lapply(output, function(x) x$parameters_validation))
        all_component_distributions_inference <- lapply(output, function(x) x$component_distributions_inference)
        all_component_distributions_validation <- lapply(output, function(x) x$component_distributions_validation)
        all_logLikelihood <- sapply(output, function(x) x$logLikelihood)
        all_criteria <- do.call(rbind, lapply(output, function(x) x$criteria))
    }
    #---Find the best fit
    best_index <- which.max(all_logLikelihood)
    fit_results <- list()
    fit_results$all <- list()
    fit_results$all$parameters_inference <- all_parameters_inference
    fit_results$all$parameters_validation <- all_parameters_validation
    fit_results$all$logLikelihood <- all_logLikelihood
    fit_results$all$criteria <- all_criteria
    fit_results$best <- list()
    fit_results$best$parameters_inference <- all_parameters_inference[best_index, ]
    fit_results$best$parameters_validation <- all_parameters_validation[best_index, ]
    fit_results$best$logLikelihood <- all_logLikelihood[best_index]
    fit_results$best$criteria <- all_criteria[best_index, ]
    fit_results$best$component_distributions_inference <- all_component_distributions_inference[[best_index]]
    fit_results$best$component_distributions_validation <- all_component_distributions_validation[[best_index]]
    return(fit_results)
}

parameter_conversion_experiment <- function(result,
                                            output_parameters_df = TRUE,
                                            mutation_count_for_fitting,
                                            sample_size,
                                            matrix_binomial_sample_size,
                                            matrix_binomial_ploidy) {
    tail_status <- result$best_fit$tail_status
    parameters_inference <- result$best_fit$parameters_inference
    parameters_validation <- result$best_fit$parameters_validation
    component_distributions_inference <- result$best_fit$component_distributions_inference
    component_distributions_validation <- result$best_fit$component_distributions_validation
    func_get_parameters <- function(parameters, tail_status) {
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
        return(output)
    }
    output <- list()
    output$inference <- func_get_parameters(parameters_inference, tail_status)
    output$validation <- func_get_parameters(parameters_validation, tail_status)
    if (output_parameters_df) {
        parameters_df <- data.frame()
        parameters_df[1, "Mutation_count_for_fitting"] <- mutation_count_for_fitting
        parameters_df[1, "Tail"] <- tail_status
        if (tail_status) {
            parameters_df[1, "Tail_power"] <- output$inference$vec_A[2]
            parameters_df[1, "Tail_mutcount_observed"] <-
                output$inference$vec_A[1] * mutation_count_for_fitting
            parameters_df[1, "Tail_mutcount_predicted"] <-
                output$inference$vec_A[1] * mutation_count_for_fitting *
                    sum(component_distributions_inference$SFS_exact[1, ]) /
                    sum(component_distributions_inference$SFS_expected[1, ]) *
                    sample_size / matrix_binomial_sample_size
        } else {
            parameters_df[1, "Tail_power"] <- NA
            parameters_df[1, "Tail_mutcount_observed"] <- NA
            parameters_df[1, "Tail_mutcount_predicted"] <- NA
        }
        parameters_df[1, "Cluster_count"] <- output$inference$N_humps
        if (output$inference$N_humps > 0) {
            for (k in 1:output$inference$N_humps) {
                parameters_df[1, paste0("Cluster_frequency_", k)] <- output$inference$vec_p[k] / matrix_binomial_ploidy
                parameters_df[1, paste0("Cluster_mutcount_observed_", k)] <-
                    output$inference$vec_K[k] * mutation_count_for_fitting
                parameters_df[1, paste0("Cluster_mutcount_predicted_", k)] <-
                    output$inference$vec_K[k] * mutation_count_for_fitting /
                        sum(component_distributions_inference$SFS_expected[k + 1, ]) /
                        sum(component_distributions_inference$SFS_exact[k + 1, ])
            }
        }
        output$parameters_df <- parameters_df
    }
    return(output)
}
