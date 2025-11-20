DECODE_plot_model_selection <- function(DECODE_result,
                                        filename = NULL,
                                        filetype = "png",
                                        width = 30, height = 15, units = "in", res = 150,
                                        SFS_limit = TRUE,
                                        data_marker_colors = NULL) {
    suppressPackageStartupMessages(library(grid))
    suppressPackageStartupMessages(library(gridExtra))
    fit_colors <- c(
        "none" = "#999999"
    )
    #---Plot each SFS fit
    x_min <- min(
        min(DECODE_result$SFS_frequencies[which(DECODE_result$SFS_data_inference_A > 0)]),
        min(DECODE_result$SFS_frequencies[which(DECODE_result$SFS_data_inference_B > 0)]),
        min(DECODE_result$SFS_frequencies[which(DECODE_result$SFS_data_validation > 0)])
    )
    x_max <- max(
        max(DECODE_result$SFS_frequencies[which(DECODE_result$SFS_data_inference_A > 0)]),
        max(DECODE_result$SFS_frequencies[which(DECODE_result$SFS_data_inference_B > 0)]),
        max(DECODE_result$SFS_frequencies[which(DECODE_result$SFS_data_validation > 0)])
    )
    func_one_fit <- function(p_right,
                             p_right_widths,
                             with_tail,
                             N_clusters,
                             text,
                             box_outline = FALSE,
                             text_color = NULL) {
        color_scheme <- c(
            "Tail" = "#999999",
            "Cluster 1" = "#D55E00",
            "Cluster 2" = "#0072B2",
            "Cluster 3" = "#009E73",
            "Cluster 4" = "#CC79A7",
            "Cluster 5" = "#E69F00",
            "Cluster 6" = "#56B4E9"
        )
        parameters <- DECODE_result[[paste0("fits_", ifelse(with_tail, "with", "without"), "_tail")]]$all_fits[[paste0(N_clusters, "_clusters")]]$parameters
        if (is.na(DECODE_result$neutral_tail)) {
            p_parameters_length <- max(length(DECODE_result$fits_with_tail$all_fits), length(DECODE_result$fits_without_tail$all_fits))
        } else {
            p_parameters_length <- ifelse(DECODE_result$neutral_tail, length(DECODE_result$fits_with_tail$all_fits), length(DECODE_result$fits_without_tail$all_fits))
        }
        p_parameters_list <- vector("list", p_parameters_length)
        for (i in seq_len(p_parameters_length)) {
            p_parameters_list[[i]] <- grobTree(rectGrob(gp = gpar(col = NA, fill = "white")))
        }
        if (with_tail) {
            p_parameters_list[[1]] <- grobTree(
                rectGrob(gp = gpar(col = NA, fill = "white")),
                textGrob(paste0("\u03B1 = ", format(round(mean(parameters[["Tail_power"]]), 2), nsmall = 2), " \u00B1 ", format(round(sd(parameters[["Tail_power"]]), 2), nsmall = 2)), gp = gpar(fontsize = 20, col = color_scheme["Tail"]), y = 1, just = "top")
            )
        }
        for (N_cluster in seq_len(N_clusters)) {
            p_parameters_list[[N_cluster + with_tail]] <- grobTree(
                rectGrob(gp = gpar(col = NA, fill = "white")),
                textGrob(paste0("f = ", format(round(mean(parameters[[paste0("Cluster_", N_cluster, "_freq")]]), 2), nsmall = 2), " \u00B1 ", format(round(sd(parameters[[paste0("Cluster_", N_cluster, "_freq")]]), 2), nsmall = 2)), gp = gpar(fontsize = 20, col = color_scheme[paste0("Cluster ", N_cluster)]), y = 1, just = "top")
            )
        }
        p_parameters <- arrangeGrob(
            grobs = p_parameters_list,
            ncol = 1
        )
        p_right_inference_A <- DECODE_plot_SFS(
            DECODE_result = DECODE_result,
            with_tail = with_tail,
            N_clusters = N_clusters,
            mode = "inference_A",
            DECODE_linewidth = 1,
            text_xlab = NULL,
            text_ylab = NULL,
            notation = FALSE,
            x_min = x_min,
            x_max = x_max,
            data_marker_colors = data_marker_colors
        )
        p_right_inference_B <- DECODE_plot_SFS(
            DECODE_result = DECODE_result,
            with_tail = with_tail,
            N_clusters = N_clusters,
            mode = "inference_B",
            DECODE_linewidth = 1,
            text_xlab = NULL,
            text_ylab = NULL,
            notation = FALSE,
            x_min = x_min,
            x_max = x_max,
            data_marker_colors = data_marker_colors
        )
        p_right_validation <- DECODE_plot_SFS(
            DECODE_result = DECODE_result,
            with_tail = with_tail,
            N_clusters = N_clusters,
            mode = "validation",
            DECODE_linewidth = 1,
            text_xlab = NULL,
            text_ylab = NULL,
            notation = FALSE,
            x_min = x_min,
            x_max = x_max,
            data_marker_colors = data_marker_colors
        )
        if (length(p_right) == 0) {
            p_right[[length(p_right) + 1]] <- arrangeGrob(
                grobTree(
                    rectGrob(gp = gpar(col = NA, fill = NA))
                ),
                grobTree(
                    rectGrob(gp = gpar(col = NA, fill = NA)),
                    textGrob(paste0("Inference A\n(n=", sum(DECODE_result$SFS_data_inference_A), ")"), rot = 90, gp = gpar(fontsize = 30, col = "#DF536B"))
                ),
                grobTree(
                    rectGrob(gp = gpar(col = NA, fill = NA)),
                    textGrob(paste0("Inference B\n(n=", sum(DECODE_result$SFS_data_inference_B), ")"), rot = 90, gp = gpar(fontsize = 30, col = "#DF536B"))
                ),
                grobTree(
                    rectGrob(gp = gpar(col = NA, fill = NA)),
                    textGrob(paste0("Validation\n(n=", sum(DECODE_result$SFS_data_validation), ")"), rot = 90, gp = gpar(fontsize = 30, col = "#2297E6"))
                ),
                grobTree(
                    rectGrob(gp = gpar(col = NA, fill = NA))
                ),
                ncol = 1,
                heights = c(0.5, 1, 1, 1, 0.1)
            )
            p_right_widths <- c(p_right_widths, 0.05)
        } else {
            p_right[[length(p_right) + 1]] <- rectGrob(gp = gpar(col = NA, fill = NA))
            p_right_widths <- c(p_right_widths, 0.05)
        }
        if (box_outline) {
            p_right[[length(p_right) + 1]] <- arrangeGrob(
                grobTree(
                    rectGrob(
                        gp = gpar(col = "#A55194", fill = NA, lwd = 15)
                    ),
                    arrangeGrob(
                        p_parameters,
                        p_right_inference_A,
                        p_right_inference_B,
                        p_right_validation,
                        ncol = 1,
                        heights = c(0.3, 1, 1, 1)
                    )
                ),
                grobTree(
                    rectGrob(gp = gpar(col = NA, fill = NA)),
                    textGrob(text, gp = gpar(fontsize = 30), y = 0, just = "bottom")
                ),
                ncol = 1,
                heights = c(3, 0.1)
            )
        } else {
            p_right[[length(p_right) + 1]] <-
                arrangeGrob(
                    p_parameters,
                    p_right_inference_A,
                    p_right_inference_B,
                    p_right_validation,
                    grobTree(
                        rectGrob(gp = gpar(col = NA, fill = "white")),
                        textGrob(text, gp = gpar(fontsize = 30), y = 0, just = "bottom")
                    ),
                    ncol = 1,
                    heights = c(0.3, 1, 1, 1, 0.1)
                )
        }
        p_right_widths <- c(p_right_widths, 1)
        p_right_widths[1] <- sum(p_right_widths[2:length(p_right_widths)]) / 17
        output <- list()
        output$p_right <- p_right
        output$p_right_widths <- p_right_widths
        return(output)
    }
    p_right <- list()
    p_right_widths <- c()
    if ("fits_with_tail" %in% names(DECODE_result)) {
        all_N_clusters <- sort(as.numeric(gsub("_clusters$", "", names(DECODE_result$fits_with_tail$all_fits))))
        for (N_clusters in all_N_clusters) {
            fit_ID <- paste0("T+", N_clusters)
            output <- func_one_fit(
                p_right = p_right,
                p_right_widths = p_right_widths,
                with_tail = TRUE,
                N_clusters = N_clusters,
                text = fit_ID,
                box_outline = ifelse(DECODE_result$best_with_tail & DECODE_result$best_N_clusters == N_clusters, TRUE, FALSE)
            )
            p_right <- output$p_right
            p_right_widths <- output$p_right_widths
        }
    }
    if ("fits_without_tail" %in% names(DECODE_result)) {
        all_N_clusters <- sort(as.numeric(gsub("_clusters$", "", names(DECODE_result$fits_without_tail$all_fits))))
        for (N_clusters in all_N_clusters) {
            fit_ID <- paste0("NT+", N_clusters)
            output <- func_one_fit(
                p_right = p_right,
                p_right_widths = p_right_widths,
                with_tail = FALSE,
                N_clusters = N_clusters,
                text = fit_ID,
                box_outline = ifelse(!DECODE_result$best_with_tail & DECODE_result$best_N_clusters == N_clusters, TRUE, FALSE)
            )
            p_right <- output$p_right
            p_right_widths <- output$p_right_widths
        }
    }
    p_right <- grid.arrange(grobs = p_right, nrow = 1, widths = p_right_widths)
    #---Plot mutation threshold selection & fit selection
    p_left_threshold_selection <- DECODE_plot_readcounts(
        DECODE_result = DECODE_result
    )
    p_left_criteria <- DECODE_plot_criteria(
        DECODE_result = DECODE_result,
        fit_colors = fit_colors,
        legend = FALSE
    )
    p_left <- arrangeGrob(
        p_left_threshold_selection,
        p_left_criteria,
        ncol = 1,
        heights = c(1.8, 1.2)
    )
    #---Combine the plots
    p <- grid.arrange(
        grobs = list(p_left, p_right),
        ncol = 2,
        widths = c(1, 3),
        left = textGrob("", gp = gpar(fontsize = 1)),
        right = textGrob("", gp = gpar(fontsize = 1)),
        top = textGrob("", gp = gpar(fontsize = 1)),
        bottom = textGrob("", gp = gpar(fontsize = 1))
    )
    if (!is.null(filename)) {
        if (filetype == "png") {
            png(paste0(filename, ".png"), res = res, width = width, height = height, units = units)
        } else if (filetype == "jpeg" | filetype == "jpg") {
            jpeg(paste0(filename, ".jpg"), res = res, width = width, height = height, units = units, quality = 95)
        } else if (filetype == "svg") {
            svg(paste0(filename, ".svg"), width = width, height = height)
        } else if (filetype == "tiff" | filetype == "tif") {
            tiff(paste0(filename, ".tiff"), res = res, width = width, height = height, units = units)
        } else if (filetype == "eps") {
            setEPS()
            postscript(paste0(filename, ".eps"), width = width, height = height)
        } else if (filetype == "pdf") {
            pdf(paste0(filename, ".pdf"), width = width, height = height)
        }
        capture.output(
            {
                grid.draw(p)
                dev.off()
            },
            file = NULL
        )
        return()
    } else {
        return(p)
    }
}

DECODE_plot_criteria <- function(DECODE_result,
                                 fit_colors = NULL,
                                 legend = TRUE) {
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(reshape2))
    color_scheme <- c(
        "other" = "#999999",
        "best" = "#A55194"
    )
    #---Retrieve criterion values for all DECODE fits
    criterion_values <- data.frame(
        fit = c(),
        value = c(),
        category = c()
    )
    fit_order <- c()
    if ("fits_with_tail" %in% names(DECODE_result)) {
        all_N_clusters <- sort(as.numeric(gsub("_clusters$", "", names(DECODE_result$fits_with_tail$all_fits))))
        for (N_clusters in all_N_clusters) {
            criterion_values <- rbind(
                criterion_values,
                data.frame(
                    fit = paste0("T+", N_clusters),
                    value = DECODE_result$fits_with_tail$all_fits[[paste0(N_clusters, "_clusters")]]$criterion_values$criterion_value,
                    category = ifelse(DECODE_result$best_with_tail & DECODE_result$best_N_clusters == N_clusters, "best", "other")
                )
            )
            fit_order <- c(fit_order, paste0("T+", N_clusters))
        }
    }
    if ("fits_without_tail" %in% names(DECODE_result)) {
        all_N_clusters <- sort(as.numeric(gsub("_clusters$", "", names(DECODE_result$fits_without_tail$all_fits))))
        for (N_clusters in all_N_clusters) {
            criterion_values <- rbind(
                criterion_values,
                data.frame(
                    fit = paste0("NT+", N_clusters),
                    value = DECODE_result$fits_without_tail$all_fits[[paste0(N_clusters, "_clusters")]]$criterion_values$criterion_value,
                    category = ifelse(!DECODE_result$best_with_tail & DECODE_result$best_N_clusters == N_clusters, "best", "other")
                )
            )
            fit_order <- c(fit_order, paste0("NT+", N_clusters))
        }
    }
    criterion_values$fit <- factor(criterion_values$fit, levels = fit_order)
    #---Plot the criterion values
    p <- ggplot(criterion_values, aes(x = fit, y = value, fill = category, color = category)) +
        geom_boxplot(alpha = 0.5, outlier.size = 3, outlier.alpha = 1) +
        labs(x = NULL, y = DECODE_result$criterion, fill = NULL) +
        scale_fill_manual(values = color_scheme, name = "") +
        scale_color_manual(values = color_scheme, name = "") +
        theme(
            text = element_text(size = 30),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5),
            legend.key.width = unit(1.5, "cm"),
            axis.text.y = element_text(angle = 45),
            axis.text.x = element_text(angle = 45, hjust = 1)
        )
    if (!legend) p <- p + theme(legend.position = "none")
    return(p)
}

DECODE_plot_SFS <- function(DECODE_result,
                            filename = NULL,
                            filetype = "png",
                            width = 30, height = 15, units = "in", res = 150,
                            with_tail = "best",
                            N_clusters = "best",
                            mode = "inference_A",
                            DECODE_linewidth = 5,
                            text_xlab = "Variant Allele Frequency",
                            color_xlab = "black",
                            text_ylab = "Mutation count",
                            notation = TRUE,
                            x_min = NULL,
                            x_max = NULL,
                            data_marker_colors = NULL,
                            error_bar = TRUE) {
    suppressPackageStartupMessages(library(ggplot2))
    if (is.null(data_marker_colors)) data_marker_colors <- c("Data" = "black")
    color_scheme <- c(
        data_marker_colors,
        "Neutral tail" = "#999999",
        "Cluster 1" = "#D55E00",
        "Cluster 2" = "#0072B2",
        "Cluster 3" = "#009E73",
        "Cluster 4" = "#CC79A7",
        "Cluster 5" = "#E69F00",
        "Cluster 6" = "#56B4E9"
    )
    vec_freq <- DECODE_result$SFS_frequencies
    SFS_totalsteps <- length(vec_freq)
    mutation_table <- DECODE_result$mutational_table
    max_total_read <- DECODE_result$max_total_read
    #---Retrieve requested DECODE fit
    if (with_tail == "best") with_tail <- DECODE_result$best_with_tail
    if (N_clusters == "best") N_clusters <- DECODE_result$best_N_clusters
    SFS_fit <- DECODE_result[[paste0("fits_", ifelse(with_tail, "with", "without"), "_tail")]]$all_fits[[paste0(N_clusters, "_clusters")]]
    #---Retrieve requested data
    SFS_data <- DECODE_result[[paste0("SFS_data_", mode)]]
    min_variant_read <- DECODE_result[[paste0("min_variant_read_", mode)]]
    min_total_read <- DECODE_result[[paste0("min_total_read_", mode)]]
    #   Categorize mutations in the data if marker info is present
    if ("Marker" %in% colnames(mutation_table)) {
        mutation_markers <- mutation_table$Marker
        unique_markers <- unique(mutation_markers)
        df_data <- data.frame(
            ID = rep(1:length(vec_freq), length(unique_markers)),
            frequency = rep(vec_freq, length(unique_markers)),
            count = rep(0, SFS_totalsteps * length(unique_markers)),
            fill = rep(unique_markers, each = SFS_totalsteps)
        )
        for (j in 1:length(mutation_table$Alt_count)) {
            no_variant <- mutation_table$Alt_count[j]
            no_total <- mutation_table$Ref_count[j] + mutation_table$Alt_count[j]
            if (no_variant >= min_variant_read & no_total >= min_total_read & no_total <= max_total_read) {
                VAF <- no_variant / no_total
                pos <- which(vec_freq >= VAF)[1]
                df_data$count[which(df_data$ID == pos & df_data$fill == mutation_markers[j])] <-
                    df_data$count[which(df_data$ID == pos & df_data$fill == mutation_markers[j])] + 1
            }
        }
        df_data$fill <- gsub("_", " ", df_data$fill)
    } else {
        df_data <- data.frame(frequency = vec_freq, count = SFS_data, fill = "Data")
    }
    df_data$fill <- factor(df_data$fill, levels = rev(names(data_marker_colors)))
    #---Prepare result from each component in DECODE fit
    df_fit <- data.frame()
    df_fit_order <- c()
    if (with_tail) {
        accumulated_mean <- as.numeric(colMeans(SFS_fit[[paste0("SFS_", mode, "_tail")]]))
        df_fit <- rbind(
            df_fit,
            data.frame(
                frequency = vec_freq,
                mean = as.numeric(colMeans(SFS_fit[[paste0("SFS_", mode, "_tail")]])),
                accumulated_mean = accumulated_mean,
                sd = as.numeric(apply(SFS_fit[[paste0("SFS_", mode, "_tail")]], 2, sd)),
                fill = "Neutral tail"
            )
        )
        df_fit_order <- c(df_fit_order, "Neutral tail")
    } else {
        accumulated_mean <- rep(0, length(vec_freq))
    }
    for (i in seq_len(N_clusters)) {
        accumulated_mean <- accumulated_mean + as.numeric(colMeans(SFS_fit[[paste0("SFS_", mode, "_cluster_", i)]]))
        df_fit <- rbind(
            df_fit,
            data.frame(
                frequency = vec_freq,
                mean = as.numeric(colMeans(SFS_fit[[paste0("SFS_", mode, "_cluster_", i)]])),
                accumulated_mean = accumulated_mean,
                sd = as.numeric(apply(SFS_fit[[paste0("SFS_", mode, "_cluster_", i)]], 2, sd)),
                fill = paste("Cluster", i)
            )
        )
        df_fit_order <- c(df_fit_order, paste("Cluster", i))
    }
    df_fit$fill <- factor(df_fit$fill, levels = rev(df_fit_order))
    #---Plot SFS fitting result from DECODE
    p <- ggplot() +
        geom_area(data = df_fit, aes(x = frequency, y = mean, fill = fill), position = "stack", alpha = 0.5) +
        geom_line(data = df_fit, aes(x = frequency, y = mean, group = fill, color = fill), position = "stack", size = DECODE_linewidth, show.legend = FALSE) +
        {
            if (error_bar) {
                geom_ribbon(data = df_fit, aes(x = frequency, ymin = pmax(accumulated_mean - sd, 0), ymax = accumulated_mean + sd, fill = fill))
            }
        } +
        geom_bar(
            data = df_data, aes(x = frequency, y = count, fill = fill),
            stat = "identity", width = 0.5 / SFS_totalsteps
        ) +
        scale_fill_manual(values = color_scheme, name = "") +
        scale_color_manual(values = color_scheme, name = "") +
        guides(fill = guide_legend(nrow = 1, keywidth = 2, keyheight = 1)) +
        xlab(text_xlab) +
        ylab(text_ylab) +
        labs(title = DECODE_result$sample_id) +
        theme(
            text = element_text(size = 30),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5),
            axis.title.x = element_text(color = color_xlab)
        )
    if (!notation) {
        p <- p +
            labs(title = NULL) +
            theme(
                legend.position = "none",
                axis.ticks.x = element_blank(),
                axis.text.x = element_blank()
            )
    }
    if (is.null(text_ylab)) {
        p <- p +
            theme(
                axis.ticks.y = element_blank(),
                axis.text.y = element_blank()
            )
    }
    if (!is.null(x_min) & !is.null(x_max)) {
        p <- p + scale_x_continuous(limits = c(x_min, x_max))
    }
    if (!is.null(filename)) {
        if (filetype == "png") {
            png(paste0(filename, ".png"), res = res, width = width, height = height, units = units)
        } else if (filetype == "jpeg" | filetype == "jpg") {
            jpeg(paste0(filename, ".jpg"), res = res, width = width, height = height, units = units, quality = 95)
        } else if (filetype == "svg") {
            svg(paste0(filename, ".svg"), width = width, height = height)
        } else if (filetype == "tiff" | filetype == "tif") {
            tiff(paste0(filename, ".tiff"), res = res, width = width, height = height, units = units)
        } else if (filetype == "eps") {
            setEPS()
            postscript(paste0(filename, ".eps"), width = width, height = height)
        } else if (filetype == "pdf") {
            pdf(paste0(filename, ".pdf"), width = width, height = height)
        }
        print(p)
        dev.off()
        return()
    } else {
        return(p)
    }
}

DECODE_plot_readcounts <- function(DECODE_result) {
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(shadowtext))
    suppressPackageStartupMessages(library(reshape2))
    mutation_table <- DECODE_result$mutational_table
    vec_min_variant_read <- min(mutation_table$Alt_count):max(mutation_table$Alt_count)
    vec_min_total_read <- min(mutation_table$Tot_count):max(mutation_table$Tot_count)
    readcount_distribution <- DECODE_result$readcount_distribution
    min_variant_read_inference_A <- max(DECODE_result$min_variant_read_inference_A, min(mutation_table$Alt_count))
    min_total_read_inference_A <- max(DECODE_result$min_total_read_inference_A, min(mutation_table$Tot_count))
    min_variant_read_inference_B <- max(DECODE_result$min_variant_read_inference_B, min(mutation_table$Alt_count))
    min_total_read_inference_B <- max(DECODE_result$min_total_read_inference_B, min(mutation_table$Tot_count))
    min_variant_read_validation <- max(DECODE_result$min_variant_read_validation, min(mutation_table$Alt_count))
    min_total_read_validation <- max(DECODE_result$min_total_read_validation, min(mutation_table$Tot_count))
    # #---Reduce the distribution to region satisfying the frequency cutoff
    # max_total_read <- vec_min_total_read[which(readcount_distribution$freq[which(readcount_distribution$min_variant_read == vec_min_variant_read[1])] <= read_distribution_freq_min)[1]]
    # max_variant_read <- vec_min_variant_read[which(readcount_distribution$freq[which(readcount_distribution$min_total_read == vec_min_total_read[1])] <= read_distribution_freq_min)[1]]
    # readcount_distribution <- readcount_distribution[which(readcount_distribution$min_total_read <= max_total_read & readcount_distribution$min_variant_read <= max_variant_read), ]
    #---Plot the readcount distribution
    freq_inference_A <- round(readcount_distribution$freq[readcount_distribution$min_total_read == min_total_read_inference_A & readcount_distribution$min_variant_read == min_variant_read_inference_A])
    freq_inference_B <- round(readcount_distribution$freq[readcount_distribution$min_total_read == min_total_read_inference_B & readcount_distribution$min_variant_read == min_variant_read_inference_B])
    freq_validation <- round(readcount_distribution$freq[readcount_distribution$min_total_read == min_total_read_validation & readcount_distribution$min_variant_read == min_variant_read_validation])
    text_inference_A <- paste0("Inference A (", freq_inference_A, "%)")
    text_inference_B <- paste0("Inference B (", freq_inference_B, "%)")
    text_validation <- paste0("Validation (", freq_validation, "%)")
    color_inference_A <- "#DF536B"
    color_inference_B <- "#DF536B"
    color_validation <- "#2297E6"
    p <- ggplot(readcount_distribution, aes(x = min_total_read, y = min_variant_read, fill = freq)) +
        geom_tile() +
        geom_rect(
            aes(
                xmin = min_total_read_inference_A - 0.5, xmax = min_total_read_inference_A + 0.5,
                ymin = min_variant_read_inference_A - 0.5, ymax = min_variant_read_inference_A + 0.5
            ),
            fill = NA, color = "white", linewidth = 1
        ) +
        geom_shadowtext(
            aes(
                x = min_total_read_inference_A + 2,
                y = min_variant_read_inference_A,
                label = text_inference_A
            ),
            angle = 45,
            hjust = 0,
            vjust = 0,
            size = 8,
            color = color_inference_A,
            bg.color = "white",
            fontface = "bold"
        ) +
        geom_rect(
            aes(
                xmin = min_total_read_inference_B - 0.5, xmax = min_total_read_inference_B + 0.5,
                ymin = min_variant_read_inference_B - 0.5, ymax = min_variant_read_inference_B + 0.5
            ),
            fill = NA, color = "white", linewidth = 1
        ) +
        geom_shadowtext(
            aes(
                x = min_total_read_inference_B + 2,
                y = min_variant_read_inference_B,
                label = text_inference_B
            ),
            angle = 45,
            hjust = 0,
            vjust = 0,
            size = 8,
            color = color_inference_B,
            bg.color = "white",
            fontface = "bold"
        ) +
        geom_rect(
            aes(
                xmin = min_total_read_validation - 0.5, xmax = min_total_read_validation + 0.5,
                ymin = min_variant_read_validation - 0.5, ymax = min_variant_read_validation + 0.5
            ),
            fill = NA, color = "white", linewidth = 1
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
            name = "% mutations"
        ) +
        theme_minimal() +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        labs(title = DECODE_result$sample_id, x = "Minimum total read count", y = "Minimum variant read count") +
        theme(
            text = element_text(size = 30),
            plot.title = element_text(size = ifelse(nchar(DECODE_result$sample_id) > 20, 25, 40), face = "bold"),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.key.width = unit(2, "cm"),
        )
    return(p)
}
