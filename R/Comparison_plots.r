plot_deconvolution_components <- function(groundtruth_df,
                                          mobster_df = NULL,
                                          deconvolution_df = NULL,
                                          cluster_count,
                                          folder_workplace = "") {
    library(ggplot2)
    cluster_count_correct <- cluster_count
    is_mobster <- !is.null(mobster_df)
    is_deconvolution <- !is.null(deconvolution_df)
    color_scheme <- c(
        "MOBSTER" = "darkorange2",
        "DECONVOLUTION" = "magenta4"
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
    #---Rearrange the humps in ground truth in descending frequencies
    for (i in 1:(dim(groundtruth_df)[1])) {
        p_row <- groundtruth_df[i, grep("^p_", names(groundtruth_df), value = TRUE)]
        p_rowsorted_indices <- order(as.vector(unlist(p_row)), decreasing = TRUE) # order
        cn <- 1
        for (id in p_rowsorted_indices) {
            groundtruth_df[i, paste0("ordered_K_", id)] <- groundtruth_df[i, paste0("K_", cn)]
            groundtruth_df[i, paste0("ordered_p_", id)] <- groundtruth_df[i, paste0("p_", cn)]
            cn <- cn + 1
        }
    }
    #---Rearrange the humps in MOBSTER results in descending frequencies
    if (is_mobster) {
        p_cols <- grep("^p_", names(mobster_df), value = TRUE)
        mobster_df[p_cols][is.na(mobster_df[p_cols])] <- 0
        for (i in 1:(dim(mobster_df)[1])) {
            p_row <- mobster_df[i, p_cols]
            p_rowsorted_indices <- order(as.vector(unlist(p_row)), decreasing = TRUE)
            cn <- 1
            for (id in p_rowsorted_indices) {
                mobster_df[i, paste0("ordered_K_", id)] <- mobster_df[i, paste0("cl_num_", cn)]
                mobster_df[i, paste0("ordered_p_", id)] <- mobster_df[i, paste0("p_", cn)]
                cn <- cn + 1
            }
        }
    }
    #---Rearrange the humps in DECONVOLUTION results in descending frequencies
    if (is_deconvolution) {
        p_cols <- grep("^Cluster_frequency_", names(deconvolution_df), value = TRUE)
        deconvolution_df[p_cols][is.na(deconvolution_df[p_cols])] <- 0
        for (i in 1:(dim(deconvolution_df)[1])) {
            p_row <- deconvolution_df[i, p_cols]
            p_rowsorted_indices <- order(as.vector(unlist(p_row)), decreasing = TRUE)
            cn <- 1
            for (id in p_rowsorted_indices) {
                deconvolution_df[i, paste0("ordered_K_", id)] <- deconvolution_df[i, paste0("Cluster_mutcount_predicted_", cn)]
                deconvolution_df[i, paste0("ordered_p_", id)] <- deconvolution_df[i, paste0("Cluster_frequency_", cn)]
                cn <- cn + 1
            }
        }
    }
    #---Find distributions of neutral tail power in each method
    if (is_mobster) {
        mobster_df$alpha <- mobster_df$Tail_shape + 1
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
    #---Find distributions of cluster counts in each method
    df_cluster_count <- data.frame()
    if (is_mobster) {
        for (cluster_count in unique(mobster_df$Kbeta_cluster)) {
            df_cluster_count <- rbind(
                df_cluster_count,
                data.frame(
                    cluster_count = cluster_count,
                    frequency = length(which(mobster_df$Kbeta_cluster == cluster_count)),
                    method = "MOBSTER"
                )
            )
        }
    }
    if (is_deconvolution) {
        for (cluster_count in unique(deconvolution_df$Cluster_count)) {
            df_cluster_count <- rbind(
                df_cluster_count,
                data.frame(
                    cluster_count = cluster_count,
                    frequency = length(which(deconvolution_df$Cluster_count == cluster_count)),
                    method = "DECONVOLUTION"
                )
            )
        }
    }
    #---Find DECONVOLUTION neutral mutation count for simulations with correct cluster count & detected tail
    if (is_deconvolution) {
        simulations_ids <- deconvolution_df$Simulation[which(deconvolution_df$Cluster_count == cluster_count_correct & deconvolution_df$Tail == TRUE)]
        deconvolution_neutral_df <- data.frame(
            Simulation = simulations_ids,
            Parameter = rep("A", cluster_count_correct * length(simulations_ids)),
            Value_TRUTH = NA,
            Value_DECONVOLUTION = NA
        )
        for (row in 1:nrow(deconvolution_neutral_df)) {
            deconvolution_neutral_df$Value_TRUTH[row] <- groundtruth_df[["A"]][deconvolution_neutral_df$Simulation[row]]
            deconvolution_neutral_df$Value_DECONVOLUTION[row] <- deconvolution_df[["Tail_mutcount_predicted"]][deconvolution_neutral_df$Simulation[row]]
        }
    }
    #---Find distributions of cluster parameters in each method (conditioned on correct cluster count & detected tail)
    cluster_parameter_df <- data.frame()
    if (is_mobster) {
        simulations_ids <- mobster_df$Simulation[which(mobster_df$Kbeta_cluster == cluster_count_correct & mobster_df$Tail == TRUE)]
        mobster_cluster_df <- data.frame(
            Simulation = rep(simulations_ids, 2 * cluster_count_correct),
            Parameter = c(rep("p", cluster_count_correct * length(simulations_ids)), rep("K", cluster_count_correct * length(simulations_ids))),
            Cluster_ID = paste0("Cluster_", rep(rep(as.character(1:cluster_count_correct), each = length(simulations_ids)), 2)),
            Cluster = rep(rep(as.character(1:cluster_count_correct), each = length(simulations_ids)), 2),
            Value_TRUTH = NA,
            Value_INFERRED = NA,
            Method = "MOBSTER"
        )
        for (row in 1:nrow(mobster_cluster_df)) {
            mobster_cluster_df$Value_TRUTH[row] <- groundtruth_df[[paste0("ordered_", mobster_cluster_df$Parameter[row], "_", mobster_cluster_df$Cluster[row])]][mobster_cluster_df$Simulation[row]]
            mobster_cluster_df$Value_INFERRED[row] <- mobster_df[[paste0("ordered_", mobster_cluster_df$Parameter[row], "_", mobster_cluster_df$Cluster[row])]][mobster_cluster_df$Simulation[row]]
        }
        cluster_parameter_df <- rbind(cluster_parameter_df, mobster_cluster_df)
    }
    if (is_deconvolution) {
        simulations_ids <- deconvolution_df$Simulation[which(deconvolution_df$Cluster_count == cluster_count_correct)]
        deconvolution_cluster_df <- data.frame(
            Simulation = rep(simulations_ids, 2 * cluster_count_correct),
            Parameter = c(rep("p", cluster_count_correct * length(simulations_ids)), rep("K", cluster_count_correct * length(simulations_ids))),
            Cluster_ID = paste0("Cluster_", rep(rep(as.character(1:cluster_count_correct), each = length(simulations_ids)), 2)),
            Cluster = rep(rep(as.character(1:cluster_count_correct), each = length(simulations_ids)), 2),
            Value_TRUTH = NA,
            Value_INFERRED = NA,
            Method = "DECONVOLUTION"
        )
        for (row in 1:nrow(deconvolution_cluster_df)) {
            deconvolution_cluster_df$Value_TRUTH[row] <- groundtruth_df[[paste0("ordered_", deconvolution_cluster_df$Parameter[row], "_", deconvolution_cluster_df$Cluster[row])]][deconvolution_cluster_df$Simulation[row]]
            deconvolution_cluster_df$Value_INFERRED[row] <- deconvolution_df[[paste0("ordered_", deconvolution_cluster_df$Parameter[row], "_", deconvolution_cluster_df$Cluster[row])]][deconvolution_cluster_df$Simulation[row]]
        }
        cluster_parameter_df <- rbind(cluster_parameter_df, deconvolution_cluster_df)
    }
    cluster_parameter_df <- cluster_parameter_df[sample(nrow(cluster_parameter_df)), ]
    #---Plot distributions of tail detection
    if (is_mobster) {
        png(paste0(folder_workplace, "Comparison_tail_detection.png"), res = 150, width = 30, height = 30, units = "in")
        xticks <- sort(unique(df_tail_detection$tail))
        xticks_label <- as.character(xticks)
        xticks_label[xticks_label == "TRUE"] <- "TRUE (correct)"
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
    }
    #---Plot distributions of cluster counts
    png(paste0(folder_workplace, "Comparison_cluster_count.png"), res = 150, width = 30, height = 30, units = "in")
    xticks <- sort(unique(df_cluster_count$cluster_count))
    xticks_label <- as.character(xticks)
    xticks_label[xticks_label == cluster_count_correct] <- paste0(cluster_count_correct, " (correct)")
    p <- ggplot() +
        geom_bar(data = df_cluster_count, aes(x = cluster_count, y = frequency, fill = method), stat = "identity", position = "dodge") +
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
    #---Plot distributions of neutral tail power
    png(paste0(folder_workplace, "Comparison_neutral_tail_power.png"), res = 150, width = 30, height = 15, units = "in")
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
    if (is_deconvolution) {
        p <- p +
            geom_histogram(data = deconvolution_df[!is.na(deconvolution_df$Tail_power), ], aes(x = Tail_power, fill = "DECONVOLUTION"), alpha = 0.5)
    }
    print(p)
    dev.off()
    #---Plot distributions of neutral tail mutation count
    png(paste0(folder_workplace, "Comparison_neutral_tail_mutation_count.png"), res = 150, width = 30, height = 30, units = "in")
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
    common_range <- c()
    if (is_deconvolution) {
        common_range <- c(
            common_range,
            deconvolution_neutral_df$Value_TRUTH[deconvolution_neutral_df$Parameter == "A"],
            deconvolution_neutral_df$Value_DECONVOLUTION[deconvolution_neutral_df$Parameter == "A"]
        )
        p <- p +
            geom_point(
                data = deconvolution_neutral_df[deconvolution_neutral_df$Parameter == "A", ],
                aes(x = Value_TRUTH, y = Value_DECONVOLUTION, fill = "DECONVOLUTION", color = "DECONVOLUTION"),
                alpha = 0.5, size = 20
            ) + xlim(range(common_range)) + ylim(range(common_range))
    }
    print(p)
    dev.off()
    #---Plot distributions of clonal frequencies
    png(paste0(folder_workplace, "Comparison_clonal_frequency.png"), res = 150, width = 30, height = 30, units = "in")
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
            legend.justification = c(0, 0.5)
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
    print(p)
    dev.off()
    #---Plot distributions of clonal mutation counts
    png(paste0(folder_workplace, "Comparison_clonal_mutation_count.png"), res = 150, width = 30, height = 30, units = "in")
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
            legend.justification = c(0, 0.5)
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
    print(p)
    dev.off()
}
