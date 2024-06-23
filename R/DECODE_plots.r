DECODE_plot <- function(DECODE_result,
                        fit = "best",
                        data_marker_colors = NULL) {
    if (is.null(data_marker_colors)) data_marker_colors <- c("Data" = "black")

    vec_freq <- DECODE_result$SFS_frequencies
    vec_SFS_real <- DECODE_result$SFS_for_fitting
    mutation_table <- DECODE_result$mutational_table
    mutation_refcounts <- mutation_table$Ref_count
    mutation_altcounts <- mutation_table$Alt_count
    mutation_totcounts <- mutation_refcounts + mutation_altcounts

    if (fit == "best") {
        vec_para_best_final <- DECODE_result$best_result$best_fit$parameters
        tail_status_final <- DECODE_result$best_result$best_fit$tail_status
        component_distributions_best_final <- DECODE_result$best_result$best_fit$component_distributions
    }

    if ("Marker" %in% colnames(mutation_table)) {
        mutation_markers <- mutation_table$Marker
    } else {
        mutation_markers <- c()
    }

    tmp <- parameter_conversion(
        result = DECODE_result$best_result,
        output_parameters_df = FALSE
    )
    vec_A <- tmp$vec_A
    vec_K <- tmp$vec_K
    vec_p <- tmp$vec_p
    N_humps <- tmp$N_humps
    #---Plot the SFS deconvolution
    color_scheme <- c(
        data_marker_colors,
        "Neutral tail" = "black",
        "Cluster 1" = "firebrick2",
        "Cluster 2" = "springgreen3",
        "Cluster 3" = "royalblue2",
        "Cluster 4" = "darkturquoise",
        "Cluster 5" = "darkorange",
        "Cluster 6" = "magenta3",
        "Cluster 7" = "salmon4"
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
        geom_area(data = df_fit, aes(x = frequency, y = count, fill = fill), position = "stack", alpha = 0.5) +
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

analysis_synthetic_test <- function(groundtruth_df,
                                    mobster_df = NULL,
                                    decode_df = NULL,
                                    text_notation = FALSE,
                                    cluster_count = NA,
                                    tail = NA,
                                    folder_workplace = "") {
    library(ggplot2)
    tail_correct <- tail
    cluster_count_correct <- cluster_count
    is_mobster <- !is.null(mobster_df)
    is_decode <- !is.null(decode_df)
    binary_color_scheme <- c(
        "TRUE" = "#D55E00",
        "FALSE" = "#0072B2"
    )
    cluster_count_color_scheme <- c(
        "1" = "#E69F00",
        "2" = "#56B4E9",
        "3" = "#009E73",
        "4" = "#F0E442",
        "5" = "#CC79A7"
    )
    method_color_scheme <- c(
        "MOBSTER" = "darkorange2",
        "DECODE" = "magenta4"
    )
    shape_scheme <- c(
        "Cluster_1" = 15,
        "Cluster_2" = 16,
        "Cluster_3" = 17,
        "Cluster_4" = 18
    )
    shape_labels <- c(
        "Cluster_1" = "1",
        "Cluster_2" = "2",
        "Cluster_3" = "3",
        "Cluster_4" = "4"
    )
    #---Rearrange the clusters in each group in descending frequencies
    for (i in 1:(dim(groundtruth_df)[1])) {
        p_row <- groundtruth_df[i, grep("^p_", names(groundtruth_df), value = TRUE)]
        p_rowsorted_indices <- order(as.vector(unlist(p_row)), decreasing = TRUE)
        for (index_new in 1:length(p_rowsorted_indices)) {
            index_old <- p_rowsorted_indices[index_new]
            groundtruth_df[i, paste0("ordered_K_", index_new)] <- groundtruth_df[i, paste0("K_expected_", index_old)]
            groundtruth_df[i, paste0("ordered_p_", index_new)] <- groundtruth_df[i, paste0("p_", index_old)]
        }
    }
    if (is_mobster) {
        p_cols <- grep("^Cluster_frequency_", names(mobster_df), value = TRUE)
        mobster_df[p_cols][is.na(mobster_df[p_cols])] <- 0
        for (i in 1:(dim(mobster_df)[1])) {
            p_row <- mobster_df[i, p_cols]
            p_rowsorted_indices <- order(as.vector(unlist(p_row)), decreasing = TRUE)
            for (index_new in 1:length(p_rowsorted_indices)) {
                index_old <- p_rowsorted_indices[index_new]
                mobster_df[i, paste0("ordered_K_", index_new)] <- mobster_df[i, paste0("Cluster_mutcount_observed_", index_old)]
                mobster_df[i, paste0("ordered_p_", index_new)] <- mobster_df[i, paste0("Cluster_frequency_", index_old)]
            }
        }
    }
    if (is_decode) {
        p_cols <- grep("^Cluster_frequency_", names(decode_df), value = TRUE)
        decode_df[p_cols][is.na(decode_df[p_cols])] <- 0
        for (i in 1:(dim(decode_df)[1])) {
            p_row <- decode_df[i, p_cols]
            p_rowsorted_indices <- order(as.vector(unlist(p_row)), decreasing = TRUE)
            for (index_new in 1:length(p_rowsorted_indices)) {
                index_old <- p_rowsorted_indices[index_new]
                decode_df[i, paste0("ordered_K_", index_new)] <- decode_df[i, paste0("Cluster_mutcount_predicted_", index_old)]
                decode_df[i, paste0("ordered_p_", index_new)] <- decode_df[i, paste0("Cluster_frequency_", index_old)]
            }
        }
    }
    #---Find distributions of characteristics in each method with respect to sample purity and coverage
    characteristic_df <- groundtruth_df
    characteristic_df$mobster_tail <- mobster_df$Tail
    characteristic_df$mobster_cluster_count <- as.factor(mobster_df$Cluster_count)
    characteristic_df$decode_cluster_count <- as.factor(decode_df$Cluster_count)
    #---Find distributions of neutral tail power in each method
    if (is_mobster) {
        mobster_df$alpha <- mobster_df$Tail_power
    }
    #---Find distributions of neutral tail detection in each method
    df_tail_detection <- data.frame()
    if (is_mobster) {
        df_tail_detection <- rbind(
            df_tail_detection,
            data.frame(
                tail = c(TRUE, FALSE),
                frequency = c(length(which(mobster_df$Tail == TRUE)), length(which(mobster_df$Tail == FALSE))),
                method = "MOBSTER"
            )
        )
    }
    if (is_decode) {
        df_tail_detection <- rbind(
            df_tail_detection,
            data.frame(
                tail = c(TRUE, FALSE),
                frequency = c(length(which(decode_df$Tail == TRUE)), length(which(decode_df$Tail == FALSE))),
                method = "DECODE"
            )
        )
    }
    #---Find distributions of cluster counts in each method (all; or conditioned on correct cluster count & detected tail)
    df_cluster_count <- data.frame()
    cluster_count_min <- Inf
    cluster_count_max <- 0
    if (is_mobster) {
        cluster_count_min <- min(cluster_count_min, min(mobster_df$Cluster_count))
        cluster_count_max <- max(cluster_count_max, max(mobster_df$Cluster_count))
    }
    if (is_decode) {
        cluster_count_min <- min(cluster_count_min, min(decode_df$Cluster_count))
        cluster_count_max <- max(cluster_count_max, max(decode_df$Cluster_count))
    }
    if (is_mobster) {
        # for (cluster_count in unique(mobster_df$Cluster_count)) {
        for (cluster_count in cluster_count_min:cluster_count_max) {
            df_cluster_count <- rbind(
                df_cluster_count,
                data.frame(
                    cluster_count = cluster_count,
                    frequency = length(which(mobster_df$Cluster_count == cluster_count)),
                    condition = "all",
                    method = "MOBSTER"
                )
            )
        }
        if (!is.na(tail_correct)) {
            # for (cluster_count in unique(mobster_df$Cluster_count[which(mobster_df$Tail == tail_correct)])) {
            for (cluster_count in cluster_count_min:cluster_count_max) {
                df_cluster_count <- rbind(
                    df_cluster_count,
                    data.frame(
                        cluster_count = cluster_count,
                        frequency = length(which(mobster_df$Cluster_count == cluster_count & mobster_df$Tail == tail_correct)),
                        condition = "correct_tail",
                        method = "MOBSTER"
                    )
                )
            }
        }
    }
    if (is_decode) {
        # for (cluster_count in unique(decode_df$Cluster_count)) {
        for (cluster_count in cluster_count_min:cluster_count_max) {
            df_cluster_count <- rbind(
                df_cluster_count,
                data.frame(
                    cluster_count = cluster_count,
                    frequency = length(which(decode_df$Cluster_count == cluster_count)),
                    condition = "all",
                    method = "DECODE"
                )
            )
        }
        if (!is.na(tail_correct)) {
            # for (cluster_count in unique(decode_df$Cluster_count[which(decode_df$Tail == tail_correct)])) {
            for (cluster_count in cluster_count_min:cluster_count_max) {
                df_cluster_count <- rbind(
                    df_cluster_count,
                    data.frame(
                        cluster_count = cluster_count,
                        frequency = length(which(decode_df$Cluster_count == cluster_count & decode_df$Tail == tail_correct)),
                        condition = "correct_tail",
                        method = "DECODE"
                    )
                )
            }
        }
    }
    #---Find distributions of cluster parameters in each method (all results; or conditioned on correct cluster count & detected tail)
    cluster_parameter_df <- data.frame()
    if (is_mobster) {
        simulations_ids <- mobster_df$Simulation
        mobster_cluster_df <- data.frame(
            Simulation = rep(simulations_ids, 2 * cluster_count_correct),
            Parameter = c(rep("p", cluster_count_correct * length(simulations_ids)), rep("K", cluster_count_correct * length(simulations_ids))),
            Cluster_ID = paste0("Cluster_", rep(rep(as.character(1:cluster_count_correct), each = length(simulations_ids)), 2)),
            Cluster = rep(rep(as.character(1:cluster_count_correct), each = length(simulations_ids)), 2),
            Value_TRUTH = NA,
            Value_INFERRED = NA,
            Tail = rep(mobster_df$Tail, 2 * cluster_count_correct),
            Cluster_count = rep(mobster_df$Cluster_count, 2 * cluster_count_correct),
            Method = "MOBSTER"
        )
        for (row in 1:nrow(mobster_cluster_df)) {
            mobster_cluster_df$Value_TRUTH[row] <- groundtruth_df[[paste0("ordered_", mobster_cluster_df$Parameter[row], "_", mobster_cluster_df$Cluster[row])]][mobster_cluster_df$Simulation[row]]
            mobster_cluster_df$Value_INFERRED[row] <- mobster_df[[paste0("ordered_", mobster_cluster_df$Parameter[row], "_", mobster_cluster_df$Cluster[row])]][mobster_cluster_df$Simulation[row]]
        }
        cluster_parameter_df <- rbind(cluster_parameter_df, mobster_cluster_df)
    }
    if (is_decode) {
        simulations_ids <- decode_df$Simulation
        decode_cluster_df <- data.frame(
            Simulation = rep(simulations_ids, 2 * cluster_count_correct),
            Parameter = c(rep("p", cluster_count_correct * length(simulations_ids)), rep("K", cluster_count_correct * length(simulations_ids))),
            Cluster_ID = paste0("Cluster_", rep(rep(as.character(1:cluster_count_correct), each = length(simulations_ids)), 2)),
            Cluster = rep(rep(as.character(1:cluster_count_correct), each = length(simulations_ids)), 2),
            Value_TRUTH = NA,
            Value_INFERRED = NA,
            Tail = rep(decode_df$Tail, 2 * cluster_count_correct),
            Cluster_count = rep(decode_df$Cluster_count, 2 * cluster_count_correct),
            Method = "DECODE"
        )
        for (row in 1:nrow(decode_cluster_df)) {
            decode_cluster_df$Value_TRUTH[row] <- groundtruth_df[[paste0("ordered_", decode_cluster_df$Parameter[row], "_", decode_cluster_df$Cluster[row])]][decode_cluster_df$Simulation[row]]
            decode_cluster_df$Value_INFERRED[row] <- decode_df[[paste0("ordered_", decode_cluster_df$Parameter[row], "_", decode_cluster_df$Cluster[row])]][decode_cluster_df$Simulation[row]]
        }
        cluster_parameter_df <- rbind(cluster_parameter_df, decode_cluster_df)
    }
    #---Find distributions of neutral mutation count in each method
    deconvolution_neutral_df <- data.frame()
    if (is_decode) {
        simulations_ids <- decode_df$Simulation[which(decode_df$Cluster_count == cluster_count_correct & decode_df$Tail == tail_correct)]
        decode_deconvolution_neutral_df <- data.frame(
            Simulation = rep(simulations_ids, 2),
            Parameter = c(rep("A_complete", length(simulations_ids)), rep("A_observed", length(simulations_ids))),
            Value_TRUTH = NA,
            Value_INFERRED = NA,
            Method = "DECODE"
        )
        for (row in 1:length(simulations_ids)) {
            decode_deconvolution_neutral_df$Value_TRUTH[row] <- groundtruth_df[["A_total"]][decode_deconvolution_neutral_df$Simulation[row]]
            decode_deconvolution_neutral_df$Value_INFERRED[row] <- decode_df[["Tail_mutcount_predicted"]][decode_deconvolution_neutral_df$Simulation[row]]
            decode_deconvolution_neutral_df$Value_TRUTH[row + length(simulations_ids)] <- groundtruth_df[["A_observed_decode"]][decode_deconvolution_neutral_df$Simulation[row]]
            decode_deconvolution_neutral_df$Value_INFERRED[row + length(simulations_ids)] <- decode_df[["Tail_mutcount_observed"]][decode_deconvolution_neutral_df$Simulation[row]]
        }
        deconvolution_neutral_df <- rbind(deconvolution_neutral_df, decode_deconvolution_neutral_df)
    }
    if (is_mobster) {
        simulations_ids <- mobster_df$Simulation[which(mobster_df$Cluster_count == cluster_count_correct & mobster_df$Tail == tail_correct)]
        mobster_deconvolution_neutral_df <- data.frame(
            Simulation = simulations_ids,
            Parameter = rep("A_observed", length(simulations_ids)),
            Value_TRUTH = NA,
            Value_INFERRED = NA,
            Method = "MOBSTER"
        )
        for (row in 1:length(simulations_ids)) {
            mobster_deconvolution_neutral_df$Value_TRUTH[row] <- groundtruth_df[["A_observed_mobster"]][mobster_deconvolution_neutral_df$Simulation[row]]
            mobster_deconvolution_neutral_df$Value_INFERRED[row] <- mobster_df[["Tail_mutcount_observed"]][mobster_deconvolution_neutral_df$Simulation[row]]
        }
        deconvolution_neutral_df <- rbind(deconvolution_neutral_df, mobster_deconvolution_neutral_df)
    }
    #---Plot distributions of MOBSTER's tail detection w.r.t. sample purity and coverage
    png(paste0(folder_workplace, "Comparison_0_sample_info_vs_mobster_tail_detection.png"), res = 150, width = 30, height = 30, units = "in")
    p <- ggplot() +
        geom_point(
            data = characteristic_df,
            aes(x = Purity, y = Coverage, fill = mobster_tail, color = mobster_tail),
            alpha = 0.5, size = 20
        ) +
        scale_fill_manual(values = binary_color_scheme, name = "MOBSTER tail") +
        scale_color_manual(values = binary_color_scheme, name = "MOBSTER tail") +
        xlab("Purity") +
        ylab("Sequencing coverage") +
        theme(
            text = element_text(size = 120),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5)
        )
    print(p)
    dev.off()
    #---Plot distributions of MOBSTER's cluster count w.r.t. sample purity and coverage
    png(paste0(folder_workplace, "Comparison_0_sample_info_vs_mobster_cluster_count.png"), res = 150, width = 30, height = 30, units = "in")
    p <- ggplot() +
        geom_point(
            data = characteristic_df,
            aes(x = Purity, y = Coverage, fill = mobster_cluster_count, color = mobster_cluster_count),
            alpha = 0.5, size = 20
        ) +
        scale_fill_manual(values = cluster_count_color_scheme, name = "MOBSTER cluster count") +
        scale_color_manual(values = cluster_count_color_scheme, name = "MOBSTER cluster count") +
        xlab("Purity") +
        ylab("Sequencing coverage") +
        theme(
            text = element_text(size = 120),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5)
        )
    print(p)
    dev.off()
    #---Plot distributions of DECODE's cluster count w.r.t. sample purity and coverage
    png(paste0(folder_workplace, "Comparison_0_sample_info_vs_decode_cluster_count.png"), res = 150, width = 30, height = 30, units = "in")
    p <- ggplot() +
        geom_point(
            data = characteristic_df,
            aes(x = Purity, y = Coverage, fill = decode_cluster_count, color = decode_cluster_count),
            alpha = 0.5, size = 20
        ) +
        scale_fill_manual(values = cluster_count_color_scheme, name = "DECODE cluster count") +
        scale_color_manual(values = cluster_count_color_scheme, name = "DECODE cluster count") +
        xlab("Purity") +
        ylab("Sequencing coverage") +
        theme(
            text = element_text(size = 120),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5)
        )
    print(p)
    dev.off()
    #---Plot distributions of tail detection
    png(paste0(folder_workplace, "Comparison_1_tail_detection.png"), res = 150, width = 30, height = 30, units = "in")
    xticks <- sort(unique(df_tail_detection$tail))
    xticks_label <- as.character(xticks)
    if (!is.na(tail_correct)) {
        if (tail_correct) {
            xticks_label[xticks_label == "TRUE"] <- "TRUE (correct)"
        } else {
            xticks_label[xticks_label == "FALSE"] <- "FALSE (correct)"
        }
    }
    p <- ggplot() +
        geom_bar(data = df_tail_detection, aes(x = tail, y = frequency, fill = method), stat = "identity", position = "dodge") +
        scale_fill_manual(values = method_color_scheme, name = "") +
        guides(fill = guide_legend(nrow = 1, keywidth = 3, keyheight = 1)) +
        xlab("Tail detection") +
        ylab("Frequency") +
        scale_x_discrete(breaks = xticks, labels = xticks_label) +
        theme(
            text = element_text(size = 120),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5)
        )
    print(p)
    dev.off()
    #---Plot distributions of cluster counts
    png(paste0(folder_workplace, "Comparison_2_cluster_count_all.png"), res = 150, width = 30, height = 30, units = "in")
    xticks <- sort(unique(df_cluster_count$cluster_count))
    xticks_label <- as.character(xticks)
    if (!is.na(cluster_count_correct)) xticks_label[xticks_label == cluster_count_correct] <- paste0(cluster_count_correct, " (correct)")
    p <- ggplot() +
        geom_bar(data = df_cluster_count[which(df_cluster_count$condition == "all"), ], aes(x = cluster_count, y = frequency, fill = method), stat = "identity", position = "dodge") +
        scale_fill_manual(values = method_color_scheme, name = "") +
        guides(fill = guide_legend(nrow = 1, keywidth = 6, keyheight = 1)) +
        xlab("Cluster count") +
        ylab("Frequency") +
        scale_x_continuous(breaks = xticks, labels = xticks_label) +
        theme(
            text = element_text(size = 120),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5)
        )
    print(p)
    dev.off()
    if (!is.na(tail_correct)) {
        png(paste0(folder_workplace, "Comparison_2_cluster_count_conditioned.png"), res = 150, width = 30, height = 30, units = "in")
        xticks <- sort(unique(df_cluster_count$cluster_count))
        xticks_label <- as.character(xticks)
        if (!is.na(cluster_count_correct)) xticks_label[xticks_label == cluster_count_correct] <- paste0(cluster_count_correct, " (correct)")
        p <- ggplot() +
            geom_bar(data = df_cluster_count[df_cluster_count$condition == "correct_tail", ], aes(x = cluster_count, y = frequency, fill = method), stat = "identity", position = "dodge") +
            scale_fill_manual(values = method_color_scheme, name = "") +
            guides(fill = guide_legend(nrow = 1, keywidth = 6, keyheight = 1)) +
            xlab("Cluster count") +
            ylab("Frequency") +
            scale_x_continuous(breaks = xticks, labels = xticks_label) +
            theme(
                text = element_text(size = 120),
                panel.background = element_rect(fill = "white", colour = "white"),
                panel.grid.major = element_line(colour = "white"),
                panel.grid.minor = element_line(colour = "white"),
                legend.position = "top",
                legend.justification = c(0, 0.5)
            )
        print(p)
        dev.off()
    }
    #---Plot distributions of clonal frequencies
    png(paste0(folder_workplace, "Comparison_3_clonal_frequency_all.png"), res = 150, width = 30, height = 30, units = "in")
    p <- ggplot() +
        geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 2) +
        scale_fill_manual(values = method_color_scheme, name = "") +
        scale_color_manual(values = method_color_scheme, name = "") +
        scale_shape_manual(values = shape_scheme, labels = shape_labels, name = "Cluster") +
        xlab("True cluster frequency") +
        ylab("Inferred cluster frequency") +
        theme(
            text = element_text(size = 120),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5),
            plot.margin = margin(0, 2, 0, 0, "cm")
        )
    common_range <- c(
        cluster_parameter_df$Value_TRUTH[cluster_parameter_df$Parameter == "p"],
        cluster_parameter_df$Value_INFERRED[cluster_parameter_df$Parameter == "p"]
    )
    p <- p +
        geom_point(
            data = cluster_parameter_df[cluster_parameter_df$Parameter == "p", ],
            aes(x = Value_TRUTH, y = Value_INFERRED, shape = Cluster_ID, fill = Method, color = Method),
            alpha = 0.5, size = 20
        ) + xlim(range(common_range)) + ylim(range(common_range))
    if (text_notation) {
        p <- p +
            geom_text(data = cluster_parameter_df, aes(x = Value_TRUTH, y = Value_INFERRED, label = Simulation, color = Method), size = 20, vjust = 0.5, hjust = 0.5)
    }
    print(p)
    dev.off()
    if (!is.na(tail_correct) | !is.na(cluster_count_correct)) {
        if (!is.na(tail_correct) & is.na(cluster_count_correct)) {
            mini_cluster_parameter_df <- cluster_parameter_df[cluster_parameter_df$Tail == tail_correct, ]
        } else if (is.na(tail_correct) & !is.na(cluster_count_correct)) {
            mini_cluster_parameter_df <- cluster_parameter_df[cluster_parameter_df$Cluster_count == cluster_count_correct, ]
        } else {
            mini_cluster_parameter_df <- cluster_parameter_df[cluster_parameter_df$Tail == tail_correct & cluster_parameter_df$Cluster_count == cluster_count_correct, ]
        }
        png(paste0(folder_workplace, "Comparison_3_clonal_frequency_conditioned.png"), res = 150, width = 30, height = 30, units = "in")
        p <- ggplot() +
            geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 2) +
            scale_fill_manual(values = method_color_scheme, name = "") +
            scale_color_manual(values = method_color_scheme, name = "") +
            scale_shape_manual(values = shape_scheme, labels = shape_labels, name = "Cluster") +
            xlab("True cluster frequency") +
            ylab("Inferred cluster frequency") +
            theme(
                text = element_text(size = 120),
                panel.background = element_rect(fill = "white", colour = "white"),
                panel.grid.major = element_line(colour = "white"),
                panel.grid.minor = element_line(colour = "white"),
                legend.position = "top",
                legend.justification = c(0, 0.5),
                plot.margin = margin(0, 2, 0, 0, "cm")
            )
        common_range <- c(
            mini_cluster_parameter_df$Value_TRUTH[mini_cluster_parameter_df$Parameter == "p"],
            mini_cluster_parameter_df$Value_INFERRED[mini_cluster_parameter_df$Parameter == "p"]
        )
        p <- p +
            geom_point(
                data = mini_cluster_parameter_df[mini_cluster_parameter_df$Parameter == "p", ],
                aes(x = Value_TRUTH, y = Value_INFERRED, shape = Cluster_ID, fill = Method, color = Method),
                alpha = 0.5, size = 20
            ) + xlim(range(common_range)) + ylim(range(common_range))
        if (text_notation) {
            p <- p +
                geom_text(data = mini_cluster_parameter_df, aes(x = Value_TRUTH, y = Value_INFERRED, label = Simulation, color = Method), size = 20, vjust = 0.5, hjust = 0.5)
        }
        print(p)
        dev.off()
    }
    #---Plot distributions of clonal mutation counts
    png(paste0(folder_workplace, "Comparison_4_clonal_observed_mutation_count_all.png"), res = 150, width = 30, height = 30, units = "in")
    p <- ggplot() +
        geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 2) +
        scale_fill_manual(values = method_color_scheme, name = "") +
        scale_color_manual(values = method_color_scheme, name = "") +
        scale_shape_manual(values = shape_scheme, labels = shape_labels, name = "Cluster") +
        xlab("True mutation count") +
        ylab("Inferred mutation count") +
        theme(
            text = element_text(size = 120),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5),
            plot.margin = margin(0, 3, 0, 0, "cm")
        )
    common_range <- c(
        cluster_parameter_df$Value_TRUTH[cluster_parameter_df$Parameter == "K"],
        cluster_parameter_df$Value_INFERRED[cluster_parameter_df$Parameter == "K"]
    )
    p <- p +
        geom_point(
            data = cluster_parameter_df[cluster_parameter_df$Parameter == "K", ],
            aes(x = Value_TRUTH, y = Value_INFERRED, shape = Cluster_ID, fill = Method, color = Method),
            alpha = 0.5, size = 20
        ) + xlim(range(common_range)) + ylim(range(common_range))
    if (text_notation) {
        p <- p +
            geom_text(data = cluster_parameter_df, aes(x = Value_TRUTH, y = Value_INFERRED, label = Simulation, color = Method), size = 20, vjust = 0.5, hjust = 0.5)
    }
    print(p)
    dev.off()
    if (!is.na(tail_correct) | !is.na(cluster_count_correct)) {
        if (!is.na(tail_correct) & is.na(cluster_count_correct)) {
            mini_cluster_parameter_df <- cluster_parameter_df[cluster_parameter_df$Tail == tail_correct, ]
        } else if (is.na(tail_correct) & !is.na(cluster_count_correct)) {
            mini_cluster_parameter_df <- cluster_parameter_df[cluster_parameter_df$Cluster_count == cluster_count_correct, ]
        } else {
            mini_cluster_parameter_df <- cluster_parameter_df[cluster_parameter_df$Tail == tail_correct & cluster_parameter_df$Cluster_count == cluster_count_correct, ]
        }
        png(paste0(folder_workplace, "Comparison_4_clonal_observed_mutation_count_conditioned.png"), res = 150, width = 30, height = 30, units = "in")
        p <- ggplot() +
            geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 2) +
            scale_fill_manual(values = method_color_scheme, name = "") +
            scale_color_manual(values = method_color_scheme, name = "") +
            scale_shape_manual(values = shape_scheme, labels = shape_labels, name = "Cluster") +
            xlab("True mutation count") +
            ylab("Inferred mutation count") +
            theme(
                text = element_text(size = 120),
                panel.background = element_rect(fill = "white", colour = "white"),
                panel.grid.major = element_line(colour = "white"),
                panel.grid.minor = element_line(colour = "white"),
                legend.position = "top",
                legend.justification = c(0, 0.5),
                plot.margin = margin(0, 3, 0, 0, "cm")
            )
        common_range <- c(
            mini_cluster_parameter_df$Value_TRUTH[mini_cluster_parameter_df$Parameter == "K"],
            mini_cluster_parameter_df$Value_INFERRED[mini_cluster_parameter_df$Parameter == "K"]
        )
        p <- p +
            geom_point(
                data = mini_cluster_parameter_df[mini_cluster_parameter_df$Parameter == "K", ],
                aes(x = Value_TRUTH, y = Value_INFERRED, shape = Cluster_ID, fill = Method, color = Method),
                alpha = 0.5, size = 20
            ) + xlim(range(common_range)) + ylim(range(common_range))
        if (text_notation) {
            p <- p +
                geom_text(data = mini_cluster_parameter_df, aes(x = Value_TRUTH, y = Value_INFERRED, label = Simulation, color = Method), size = 20, vjust = 0.5, hjust = 0.5)
        }
        print(p)
        dev.off()
    }
    #---Plot distributions of neutral tail power
    png(paste0(folder_workplace, "Comparison_5_neutral_tail_power_raw.png"), res = 150, width = 30, height = 15, units = "in")
    tail_power <- unique(groundtruth_df$alpha)
    p <- ggplot() +
        scale_fill_manual(values = method_color_scheme, name = "") +
        guides(fill = guide_legend(nrow = 1, keywidth = 3, keyheight = 1)) +
        xlab("Neutral tail power") +
        ylab("Frequency") +
        theme(
            text = element_text(size = 60),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5)
        )
    if (length(tail_power) == 1) {
        p <- p +
            geom_vline(aes(xintercept = tail_power), color = "black", linewidth = 3)
    }
    if (is_mobster) {
        p <- p +
            geom_histogram(data = mobster_df[!is.na(mobster_df$alpha), ], aes(x = alpha, fill = "MOBSTER"), alpha = 0.5)
    }
    if (is_decode) {
        p <- p +
            geom_histogram(data = decode_df[!is.na(decode_df$Tail_power), ], aes(x = Tail_power, fill = "DECODE"), alpha = 0.5)
    }
    print(p)
    dev.off()
    #---Plot distributions of neutral tail mutation count
    png(paste0(folder_workplace, "Comparison_6_neutral_tail_observed_mutation_count.png"), res = 150, width = 30, height = 30, units = "in")
    p <- ggplot() +
        geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 2) +
        scale_fill_manual(values = method_color_scheme, name = "") +
        scale_color_manual(values = method_color_scheme, name = "") +
        xlab("True neutral mutation count") +
        ylab("Inferred neutral mutation count") +
        theme(
            text = element_text(size = 120),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5)
        )
    common_range <- c(
        deconvolution_neutral_df$Value_TRUTH[deconvolution_neutral_df$Parameter == "A_observed"],
        deconvolution_neutral_df$Value_INFERRED[deconvolution_neutral_df$Parameter == "A_observed"]
    )
    p <- p +
        geom_point(
            data = deconvolution_neutral_df[deconvolution_neutral_df$Parameter == "A_observed", ],
            aes(x = Value_TRUTH, y = Value_INFERRED, fill = Method, color = Method),
            alpha = 0.5, size = 20
        ) + xlim(range(common_range)) + ylim(range(common_range))
    if (text_notation) {
        p <- p +
            geom_text(data = deconvolution_neutral_df[deconvolution_neutral_df$Parameter == "A_observed", ], aes(x = Value_TRUTH, y = Value_INFERRED, label = Simulation, color = Method), size = 20, vjust = 0.5, hjust = 0.5)
    }
    print(p)
    dev.off()
    png(paste0(folder_workplace, "Comparison_6_neutral_tail_expected_mutation_count.png"), res = 150, width = 30, height = 30, units = "in")
    p <- ggplot() +
        geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 2) +
        scale_fill_manual(values = method_color_scheme, name = "") +
        scale_color_manual(values = method_color_scheme, name = "") +
        xlab("True neutral mutation count") +
        ylab("Inferred neutral mutation count") +
        theme(
            text = element_text(size = 120),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5)
        )
    common_range <- c(
        deconvolution_neutral_df$Value_TRUTH[deconvolution_neutral_df$Parameter == "A_complete"],
        deconvolution_neutral_df$Value_INFERRED[deconvolution_neutral_df$Parameter == "A_complete"]
    )
    p <- p +
        geom_point(
            data = deconvolution_neutral_df[deconvolution_neutral_df$Parameter == "A_complete", ],
            aes(x = Value_TRUTH, y = Value_INFERRED, fill = Method, color = Method),
            alpha = 0.5, size = 20
        ) + xlim(range(common_range)) + ylim(range(common_range))
    if (text_notation) {
        p <- p +
            geom_text(data = deconvolution_neutral_df[deconvolution_neutral_df$Parameter == "A_complete", ], aes(x = Value_TRUTH, y = Value_INFERRED, label = Simulation, color = Method), size = 20, vjust = 0.5, hjust = 0.5)
    }
    print(p)
    dev.off()
}

analysis_ICGC <- function(sample_information_df,
                          mobster_df = NULL,
                          decode_df = NULL,
                          bound_truncal_frequency_vs_sample_purity = 0.05,
                          text_notation = FALSE,
                          folder_workplace = "") {
    library(ggplot2)
    is_mobster <- !is.null(mobster_df)
    is_decode <- !is.null(decode_df)
    method_color_scheme <- c(
        "MOBSTER" = "darkorange2",
        "DECODE" = "magenta4"
    )
    shape_scheme <- c(
        "Cluster_1" = 15,
        "Cluster_2" = 16,
        "Cluster_3" = 17,
        "Cluster_4" = 18
    )
    shape_labels <- c(
        "Cluster_1" = "1",
        "Cluster_2" = "2",
        "Cluster_3" = "3",
        "Cluster_4" = "4"
    )
    #---Rearrange the clusters in each group in descending frequencies
    if (is_mobster) {
        p_cols <- grep("^Cluster_frequency_", names(mobster_df), value = TRUE)
        mobster_df[p_cols][is.na(mobster_df[p_cols])] <- 0
        for (i in 1:(dim(mobster_df)[1])) {
            p_row <- mobster_df[i, p_cols]
            p_rowsorted_indices <- order(as.vector(unlist(p_row)), decreasing = TRUE)
            for (index_new in 1:length(p_rowsorted_indices)) {
                index_old <- p_rowsorted_indices[index_new]
                mobster_df[i, paste0("ordered_K_", index_new)] <- mobster_df[i, paste0("Cluster_mutcount_observed_", index_old)]
                mobster_df[i, paste0("ordered_p_", index_new)] <- mobster_df[i, paste0("Cluster_frequency_", index_old)]
            }
        }
    }
    if (is_decode) {
        p_cols <- grep("^Cluster_frequency_", names(decode_df), value = TRUE)
        decode_df[p_cols][is.na(decode_df[p_cols])] <- 0
        for (i in 1:(dim(decode_df)[1])) {
            p_row <- decode_df[i, p_cols]
            p_rowsorted_indices <- order(as.vector(unlist(p_row)), decreasing = TRUE)
            for (index_new in 1:length(p_rowsorted_indices)) {
                index_old <- p_rowsorted_indices[index_new]
                decode_df[i, paste0("ordered_K_", index_new)] <- decode_df[i, paste0("Cluster_mutcount_predicted_", index_old)]
                decode_df[i, paste0("ordered_p_", index_new)] <- decode_df[i, paste0("Cluster_frequency_", index_old)]
            }
        }
    }
    #---Find comparison of truncal cluster frequency against sample purity
    df_truncal_frequency_vs_sample_purity <- data.frame()
    if (is_mobster) {
        mobster_df_truncal_frequency_vs_sample_purity <- data.frame(
            Sample = mobster_df$Sample,
            Sample_short = sapply(strsplit(mobster_df$Sample, "-"), `[`, 1),
            Truncal_frequency = 2 * mobster_df$ordered_p_1,
            Purity = NA,
            Method = "MOBSTER"
        )
        for (i in 1:nrow(mobster_df_truncal_frequency_vs_sample_purity)) {
            mobster_df_truncal_frequency_vs_sample_purity$Purity[i] <- sample_information_df$purity[which(sample_information_df$aliquot_id == mobster_df_truncal_frequency_vs_sample_purity$Sample[i])]
        }
        df_truncal_frequency_vs_sample_purity <- rbind(
            df_truncal_frequency_vs_sample_purity,
            mobster_df_truncal_frequency_vs_sample_purity
        )
    }

    if (is_decode) {
        decode_df_truncal_frequency_vs_sample_purity <- data.frame(
            Sample = decode_df$Sample,
            Sample_short = sapply(strsplit(decode_df$Sample, "-"), `[`, 1),
            Truncal_frequency = 2 * decode_df$ordered_p_1,
            Purity = NA,
            Method = "DECODE"
        )
        for (i in 1:nrow(decode_df_truncal_frequency_vs_sample_purity)) {
            decode_df_truncal_frequency_vs_sample_purity$Purity[i] <- sample_information_df$purity[which(sample_information_df$aliquot_id == decode_df_truncal_frequency_vs_sample_purity$Sample[i])]
        }
        df_truncal_frequency_vs_sample_purity <- rbind(
            df_truncal_frequency_vs_sample_purity,
            decode_df_truncal_frequency_vs_sample_purity
        )
    }


    df_truncal_frequency_vs_sample_purity$Within_bounds <- NA
    df_truncal_frequency_vs_sample_purity$Within_bounds[
        which(
            df_truncal_frequency_vs_sample_purity$Truncal_frequency < df_truncal_frequency_vs_sample_purity$Purity + bound_truncal_frequency_vs_sample_purity &
                df_truncal_frequency_vs_sample_purity$Truncal_frequency > df_truncal_frequency_vs_sample_purity$Purity - bound_truncal_frequency_vs_sample_purity
        )
    ] <- "Correct"
    df_truncal_frequency_vs_sample_purity$Within_bounds[
        which(
            df_truncal_frequency_vs_sample_purity$Truncal_frequency >= df_truncal_frequency_vs_sample_purity$Purity + bound_truncal_frequency_vs_sample_purity
        )
    ] <- "Above"
    df_truncal_frequency_vs_sample_purity$Within_bounds[
        which(
            df_truncal_frequency_vs_sample_purity$Truncal_frequency <= df_truncal_frequency_vs_sample_purity$Purity - bound_truncal_frequency_vs_sample_purity
        )
    ] <- "Below"
    #---Plot truncal cluster frequency against sample purity
    png(paste0(folder_workplace, "ICGC_1_MOBSTER_truncal_frequency_vs_sample_purity.png"), res = 150, width = 30, height = 30, units = "in")
    p <- ggplot() +
        geom_point(
            data = df_truncal_frequency_vs_sample_purity,
            # data = df_truncal_frequency_vs_sample_purity[which(df_truncal_frequency_vs_sample_purity$Method == "MOBSTER"), ],
            aes(x = Purity, y = Truncal_frequency, fill = Method, color = Method),
            alpha = 0.5, size = 20
        ) +
        scale_fill_manual(values = method_color_scheme, name = "") +
        scale_color_manual(values = method_color_scheme, name = "") +
        xlab("Sample purity") +
        ylab("2 \u00D7 Truncal cluster frequency") +
        theme(
            text = element_text(size = 120),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5),
            plot.margin = margin(0, 2, 0, 0, "cm")
        )
    if (text_notation) {
        p <- p +
            geom_text(data = df_truncal_frequency_vs_sample_purity[which(df_truncal_frequency_vs_sample_purity$Within_bounds != "Correct"), ], aes(x = Purity, y = Truncal_frequency, label = Sample_short, color = Method), size = 20, vjust = 0, hjust = 0, angle = 45)
    }
    p <- p +
        geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 2) +
        geom_abline(intercept = bound_truncal_frequency_vs_sample_purity, slope = 1, color = "black", linewidth = 2, linetype = "dashed") +
        geom_abline(intercept = -bound_truncal_frequency_vs_sample_purity, slope = 1, color = "black", linewidth = 2, linetype = "dashed")
    print(p)
    cat(paste0("\nMOBSTER samples with truncal cluster frequency within ", 100 * bound_truncal_frequency_vs_sample_purity, "% of sample purity: ", length(which(df_truncal_frequency_vs_sample_purity$Method == "MOBSTER" & df_truncal_frequency_vs_sample_purity$Within_bounds == "Correct")), " (", 100 * round(length(which(df_truncal_frequency_vs_sample_purity$Method == "MOBSTER" & df_truncal_frequency_vs_sample_purity$Within_bounds == "Correct")) / length(which(df_truncal_frequency_vs_sample_purity$Method == "MOBSTER")), 2), "%)", "\n"))
    cat(paste0("MOBSTER samples with truncal cluster frequency > sample purity + ", 100 * bound_truncal_frequency_vs_sample_purity, "%:       ", length(which(df_truncal_frequency_vs_sample_purity$Method == "MOBSTER" & df_truncal_frequency_vs_sample_purity$Within_bounds == "Above")), " (", 100 * round(length(which(df_truncal_frequency_vs_sample_purity$Method == "MOBSTER" & df_truncal_frequency_vs_sample_purity$Within_bounds == "Above")) / length(which(df_truncal_frequency_vs_sample_purity$Method == "MOBSTER")), 2), "%)", "\n"))
    cat(paste0("MOBSTER samples with truncal cluster frequency < sample purity - ", 100 * bound_truncal_frequency_vs_sample_purity, "%:       ", length(which(df_truncal_frequency_vs_sample_purity$Method == "MOBSTER" & df_truncal_frequency_vs_sample_purity$Within_bounds == "Below")), " (", 100 * round(length(which(df_truncal_frequency_vs_sample_purity$Method == "MOBSTER" & df_truncal_frequency_vs_sample_purity$Within_bounds == "Below")) / length(which(df_truncal_frequency_vs_sample_purity$Method == "MOBSTER")), 2), "%)", "\n\n"))
    cat(paste0("\nDECODE samples with truncal cluster frequency within ", 100 * bound_truncal_frequency_vs_sample_purity, "% of sample purity: ", length(which(df_truncal_frequency_vs_sample_purity$Method == "DECODE" & df_truncal_frequency_vs_sample_purity$Within_bounds == "Correct")), " (", 100 * round(length(which(df_truncal_frequency_vs_sample_purity$Method == "DECODE" & df_truncal_frequency_vs_sample_purity$Within_bounds == "Correct")) / length(which(df_truncal_frequency_vs_sample_purity$Method == "DECODE")), 2), "%)", "\n"))
    cat(paste0("DECODE samples with truncal cluster frequency > sample purity + ", 100 * bound_truncal_frequency_vs_sample_purity, "%:       ", length(which(df_truncal_frequency_vs_sample_purity$Method == "DECODE" & df_truncal_frequency_vs_sample_purity$Within_bounds == "Above")), " (", 100 * round(length(which(df_truncal_frequency_vs_sample_purity$Method == "DECODE" & df_truncal_frequency_vs_sample_purity$Within_bounds == "Above")) / length(which(df_truncal_frequency_vs_sample_purity$Method == "DECODE")), 2), "%)", "\n"))
    cat(paste0("DECODE samples with truncal cluster frequency < sample purity - ", 100 * bound_truncal_frequency_vs_sample_purity, "%:       ", length(which(df_truncal_frequency_vs_sample_purity$Method == "DECODE" & df_truncal_frequency_vs_sample_purity$Within_bounds == "Below")), " (", 100 * round(length(which(df_truncal_frequency_vs_sample_purity$Method == "DECODE" & df_truncal_frequency_vs_sample_purity$Within_bounds == "Below")) / length(which(df_truncal_frequency_vs_sample_purity$Method == "DECODE")), 2), "%)", "\n\n"))
    dev.off()
}
