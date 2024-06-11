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
    library_SFS_component <- DECODE_result$library_SFS_component
    list_neutral_powers <- DECODE_result$list_neutral_powers
    list_frequencies <- DECODE_result$list_frequencies

    if (fit == "best") {
        vec_para_best_final <- DECODE_result$best_fit$parameters
        tail_status_final <- DECODE_result$best_fit$tail_status
    }

    if ("Marker" %in% colnames(mutation_table)) {
        mutation_markers <- mutation_table$Marker
    } else {
        mutation_markers <- c()
    }

    tmp <- parameter_conversion(
        parameters = vec_para_best_final,
        tail_status = tail_status_final,
        parameters_df = FALSE
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
            vec_A = vec_A,
            vec_K = c(),
            vec_p = c(),
            list_neutral_powers = list_neutral_powers,
            list_frequencies = list_frequencies,
            library_SFS_component = library_SFS_component
        )
        SFS_neutral <- vec_A[1] * mutation_count * SFS_neutral / sum(SFS_neutral)
        df_fit <- rbind(df_fit, data.frame(frequency = vec_freq, count = SFS_neutral, fill = "Neutral tail"))
        df_fit_order <- c(df_fit_order, "Neutral tail")
    }
    if (N_humps > 0) {
        for (i in 1:N_humps) {
            SFS_hump <- compute_SFS(
                vec_A = c(0, list_neutral_powers[1]),
                vec_K = vec_K[i],
                vec_p = vec_p[i],
                list_neutral_powers = list_neutral_powers,
                list_frequencies = list_frequencies,
                library_SFS_component = library_SFS_component
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

comparison_synthetic_test <- function(groundtruth_df,
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
    color_scheme <- c(
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
        p_rowsorted_indices <- order(as.vector(unlist(p_row)), decreasing = TRUE) # order
        cn <- 1
        for (id in p_rowsorted_indices) {
            groundtruth_df[i, paste0("ordered_K_", id)] <- groundtruth_df[i, paste0("K_expected_", cn)]
            groundtruth_df[i, paste0("ordered_p_", id)] <- groundtruth_df[i, paste0("p_", cn)]
            cn <- cn + 1
        }
    }
    if (is_mobster) {
        p_cols <- grep("^Cluster_frequency_", names(mobster_df), value = TRUE)
        mobster_df[p_cols][is.na(mobster_df[p_cols])] <- 0
        for (i in 1:(dim(mobster_df)[1])) {
            p_row <- mobster_df[i, p_cols]
            p_rowsorted_indices <- order(as.vector(unlist(p_row)), decreasing = TRUE)
            cn <- 1
            for (id in p_rowsorted_indices) {
                mobster_df[i, paste0("ordered_K_", id)] <- mobster_df[i, paste0("Cluster_mutcount_observed_", cn)]
                mobster_df[i, paste0("ordered_p_", id)] <- mobster_df[i, paste0("Cluster_frequency_", cn)]
                cn <- cn + 1
            }
        }
    }
    if (is_decode) {
        p_cols <- grep("^Cluster_frequency_", names(decode_df), value = TRUE)
        decode_df[p_cols][is.na(decode_df[p_cols])] <- 0
        for (i in 1:(dim(decode_df)[1])) {
            p_row <- decode_df[i, p_cols]
            p_rowsorted_indices <- order(as.vector(unlist(p_row)), decreasing = TRUE)
            cn <- 1
            for (id in p_rowsorted_indices) {
                decode_df[i, paste0("ordered_K_", id)] <- decode_df[i, paste0("Cluster_mutcount_predicted_", cn)]
                decode_df[i, paste0("ordered_p_", id)] <- decode_df[i, paste0("Cluster_frequency_", cn)]
                cn <- cn + 1
            }
        }
    }
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
            decode_deconvolution_neutral_df$Value_TRUTH[row] <- groundtruth_df[["A_expected"]][decode_deconvolution_neutral_df$Simulation[row]]
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
        print(mobster_deconvolution_neutral_df)
        deconvolution_neutral_df <- rbind(deconvolution_neutral_df, mobster_deconvolution_neutral_df)
    }
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
        scale_fill_manual(values = color_scheme, name = "") +
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
        scale_fill_manual(values = color_scheme, name = "") +
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
            scale_fill_manual(values = color_scheme, name = "") +
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
        scale_fill_manual(values = color_scheme, name = "") +
        scale_color_manual(values = color_scheme, name = "") +
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
            scale_fill_manual(values = color_scheme, name = "") +
            scale_color_manual(values = color_scheme, name = "") +
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
        scale_fill_manual(values = color_scheme, name = "") +
        scale_color_manual(values = color_scheme, name = "") +
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
            scale_fill_manual(values = color_scheme, name = "") +
            scale_color_manual(values = color_scheme, name = "") +
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
        scale_fill_manual(values = color_scheme, name = "") +
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
        scale_fill_manual(values = color_scheme, name = "") +
        scale_color_manual(values = color_scheme, name = "") +
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
        scale_fill_manual(values = color_scheme, name = "") +
        scale_color_manual(values = color_scheme, name = "") +
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
    print(deconvolution_neutral_df)
    print(p)
    dev.off()
}
