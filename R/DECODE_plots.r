DECODE_plot_model_selection <- function(DECODE_result,
                                        fit = "best",
                                        mode = "inference_A",
                                        SFS_limit = TRUE,
                                        data_marker_colors = NULL) {
    suppressPackageStartupMessages(library(gridExtra))
    fit_colors <- c(
        "none" = "#999999",
        "best (T)" = "#009E73",
        "best (WT)" = "#009E73",
        "best" = "#CC79A7"
    )
    #---Plot each SFS fit
    func_one_fit <- function(p_right, fit, text, text_color) {
        if (length(p_right) == 0) {
            text_ylab_inference_A <- "Inference A"
            text_ylab_inference_B <- "Inference B"
            text_ylab_validation <- "Validation"
        } else {
            text_ylab_inference_A <- NULL
            text_ylab_inference_B <- NULL
            text_ylab_validation <- NULL
            text_xlab_validation <- NULL
        }
        p_right_inference_A <- DECODE_plot_SFS(
            DECODE_result = DECODE_result,
            fit = fit,
            mode = "inference_A",
            DECODE_linewidth = 1,
            text_xlab = NULL,
            text_ylab = text_ylab_inference_A,
            notation = FALSE,
            SFS_limit = SFS_limit,
            data_marker_colors = data_marker_colors
        )
        p_right_inference_B <- DECODE_plot_SFS(
            DECODE_result = DECODE_result,
            fit = fit,
            mode = "inference_B",
            DECODE_linewidth = 1,
            text_xlab = NULL,
            text_ylab = text_ylab_inference_B,
            notation = FALSE,
            SFS_limit = SFS_limit,
            data_marker_colors = data_marker_colors
        )
        p_right_validation <- DECODE_plot_SFS(
            DECODE_result = DECODE_result,
            fit = fit,
            mode = "validation",
            DECODE_linewidth = 1,
            text_xlab = text,
            color_xlab = text_color,
            text_ylab = text_ylab_validation,
            notation = FALSE,
            SFS_limit = SFS_limit,
            data_marker_colors = data_marker_colors
        )
        p_right[[length(p_right) + 1]] <-
            arrangeGrob(
                p_right_inference_A,
                p_right_inference_B,
                p_right_validation,
                ncol = 1
            )
        return(p_right)
    }
    p_right <- list()
    if ("fits_with_tail" %in% names(DECODE_result)) {
        criteria_with_tail <- data.frame()
        for (fit in names(DECODE_result$fits_with_tail$all_fits)) {
            fit_ID <- paste0("T+", substr(fit, 1, 1))
            criteria_new <- DECODE_result$fits_with_tail$all_fits[[fit]]$best$criteria
            criteria_new$fit <- fit_ID
            criteria_new$note <- DECODE_result$fits_with_tail$all_fits[[fit]]$note
            if (criteria_new$note == "best given tail status") criteria_new$note <- "best (T)"
            criteria_with_tail <- rbind(criteria_with_tail, criteria_new)
            p_right <- func_one_fit(
                p_right = p_right,
                fit = paste0("fits_with_tail:", fit),
                text = fit_ID,
                text_color = fit_colors[criteria_new$note]
            )
        }
    }
    if ("fits_without_tail" %in% names(DECODE_result)) {
        criteria_without_tail <- data.frame()
        for (fit in names(DECODE_result$fits_without_tail$all_fits)) {
            fit_ID <- paste0("WT+", substr(fit, 1, 1))
            criteria_new <- DECODE_result$fits_without_tail$all_fits[[fit]]$best$criteria
            criteria_new$fit <- fit_ID
            criteria_new$note <- DECODE_result$fits_without_tail$all_fits[[fit]]$note
            if (criteria_new$note == "best given tail status") criteria_new$note <- "best (WT)"
            criteria_without_tail <- rbind(criteria_without_tail, criteria_new)
            p_right <- func_one_fit(
                p_right = p_right,
                fit = paste0("fits_without_tail:", fit),
                text = fit_ID,
                text_color = fit_colors[criteria_new$note]
            )
        }
    }
    p_right <- grid.arrange(grobs = p_right, nrow = 1)
    #---Plot mutation threshold selection & fit selection
    p_left_threshold_selection <- DECODE_plot_readcounts(
        DECODE_result = DECODE_result
    )
    if (("fits_with_tail" %in% names(DECODE_result)) & ("fits_without_tail" %in% names(DECODE_result))) {
        p_left_criteria_with_tail <- DECODE_plot_criteria(
            criteria = criteria_with_tail,
            criterion = DECODE_result$criterion,
            criterion_ratio = DECODE_result$criterion_ratio,
            fit_colors = fit_colors
        )
        p_left_criteria_without_tail <- DECODE_plot_criteria(
            criteria = criteria_without_tail,
            criterion = DECODE_result$criterion,
            criterion_ratio = DECODE_result$criterion_ratio,
            fit_colors = fit_colors
        )
        p_left <- arrangeGrob(
            p_left_threshold_selection,
            p_left_criteria_with_tail,
            p_left_criteria_without_tail,
            ncol = 1,
            heights = c(4, 1, 1)
        )
    } else if ("fits_with_tail" %in% names(DECODE_result)) {
        p_left_criteria_with_tail <- DECODE_plot_criteria(
            criteria = criteria_with_tail,
            criterion = DECODE_result$criterion,
            criterion_ratio = DECODE_result$criterion_ratio,
            fit_colors = fit_colors
        )
        p_left <- arrangeGrob(
            p_left_threshold_selection,
            p_left_criteria_with_tail,
            ncol = 1,
            heights = c(2, 1)
        )
    } else if ("fits_without_tail" %in% names(DECODE_result)) {
        p_left_criteria_without_tail <- DECODE_plot_criteria(
            criteria = criteria_without_tail,
            criterion = DECODE_result$criterion,
            criterion_ratio = DECODE_result$criterion_ratio,
            fit_colors = fit_colors
        )
        p_left <- arrangeGrob(
            p_left_threshold_selection,
            p_left_criteria_without_tail,
            ncol = 1,
            heights = c(2, 1)
        )
    }
    #---Combine the plots
    p <- grid.arrange(grobs = list(p_left, p_right), ncol = 2, widths = c(1, 2.5))
    return(p)
}

DECODE_plot_criteria <- function(criteria,
                                 criterion,
                                 criterion_ratio,
                                 fit_colors,
                                 legend = TRUE) {
    #---Get ratio of each fit compared to previous fit
    criteria[[paste0(criterion, "_target")]] <- NA
    for (row in 1:nrow(criteria)) {
        fit_current <- criteria[row, "fit"]
        fit_previous <- paste0(sub("\\+\\d+$", "", fit_current), "+", as.numeric(sub(".*\\+(\\d+)$", "\\1", fit_current)) - 1)
        if (fit_previous %in% criteria$fit) {
            criteria[row, paste0(criterion, "_target")] <- criteria[criteria$fit == fit_previous, criterion] * criterion_ratio
        }
    }
    #---Plot the criterion selection
    y_min <- min(c(min(criteria[[criterion]], na.rm = TRUE), min(criteria[[paste0(criterion, "_target")]], na.rm = TRUE)))
    y_max <- max(c(max(criteria[[criterion]], na.rm = TRUE), max(criteria[[paste0(criterion, "_target")]], na.rm = TRUE)))
    p <- ggplot(criteria) +
        geom_point(aes(x = fit, y = !!sym(criterion), shape = "Actual", color = note), size = 10) +
        geom_point(aes(x = fit, y = !!sym(paste0(criterion, "_target")), shape = "Target"), color = "black", size = 6, stroke = 2, na.rm = TRUE) +
        labs(x = NULL, y = NULL, shape = NULL) +
        scale_shape_manual(values = c("Actual" = 16, "Target" = 6), labels = c(criterion, paste0(criterion, " threshold"))) +
        scale_color_manual(values = fit_colors, name = "", breaks = setdiff(names(fit_colors), "none")) +
        expand_limits(y = c(y_min - 0.2 * (y_max - y_min), y_max + 0.2 * (y_max - y_min))) +
        theme(
            text = element_text(size = 30),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5),
            legend.key.width = unit(1.5, "cm")
        )
    if (!legend) p <- p + theme(legend.position = "none")
    return(p)
}

DECODE_plot_SFS <- function(DECODE_result,
                            fit = "best",
                            mode = "inference_A",
                            DECODE_linewidth = 5,
                            text_xlab = "Variant Allele Frequency",
                            color_xlab = "black",
                            text_ylab = "Mutation count",
                            notation = TRUE,
                            SFS_limit = FALSE,
                            data_marker_colors = NULL) {
    suppressPackageStartupMessages(library(ggplot2))
    if (is.null(data_marker_colors)) data_marker_colors <- c("Data" = "black")
    vec_freq <- DECODE_result$SFS_frequencies
    SFS_totalsteps <- length(vec_freq)
    mutation_table <- DECODE_result$mutational_table
    max_total_read <- DECODE_result$max_total_read
    if (mode == "inference_A") {
        vec_SFS_real <- DECODE_result$SFS_data_inference_A
        min_variant_read <- DECODE_result$min_variant_read_inference_A
        min_total_read <- DECODE_result$min_total_read_inference_A
        if (fit == "best") {
            parameters <- DECODE_result$final_fit$best_fit$parameters_inference_A
            tail_status <- DECODE_result$final_fit$best_fit$tail_status
            component_distributions <- DECODE_result$final_fit$best_fit$component_distributions_inference_A
        } else {
            detail_tail <- sub(":.*", "", fit)
            detail_Ncluster <- sub(".*:", "", fit)
            parameters <- DECODE_result[[detail_tail]]$all_fits[[detail_Ncluster]]$best$parameters_inference_A
            tail_status <- ifelse(detail_tail == "fits_with_tail", TRUE, FALSE)
            component_distributions <- DECODE_result[[detail_tail]]$all_fits[[detail_Ncluster]]$best$component_distributions_inference_A
        }
    } else if (mode == "inference_B") {
        vec_SFS_real <- DECODE_result$SFS_data_inference_B
        min_variant_read <- DECODE_result$min_variant_read_inference_B
        min_total_read <- DECODE_result$min_total_read_inference_B
        if (fit == "best") {
            parameters <- DECODE_result$final_fit$best_fit$parameters_inference_B
            tail_status <- DECODE_result$final_fit$best_fit$tail_status
            component_distributions <- DECODE_result$final_fit$best_fit$component_distributions_inference_B
        } else {
            detail_tail <- sub(":.*", "", fit)
            detail_Ncluster <- sub(".*:", "", fit)
            parameters <- DECODE_result[[detail_tail]]$all_fits[[detail_Ncluster]]$best$parameters_inference_B
            tail_status <- ifelse(detail_tail == "fits_with_tail", TRUE, FALSE)
            component_distributions <- DECODE_result[[detail_tail]]$all_fits[[detail_Ncluster]]$best$component_distributions_inference_B
        }
    } else if (mode == "validation") {
        vec_SFS_real <- DECODE_result$SFS_data_validation
        min_variant_read <- DECODE_result$min_variant_read_validation
        min_total_read <- DECODE_result$min_total_read_validation
        if (fit == "best") {
            parameters <- DECODE_result$final_fit$best_fit$parameters_validation
            tail_status <- DECODE_result$final_fit$best_fit$tail_status
            component_distributions <- DECODE_result$final_fit$best_fit$component_distributions_validation
            criteria_validation_index <- DECODE_result$final_fit$best_fit$criteria_validation_index
        } else {
            detail_tail <- sub(":.*", "", fit)
            detail_Ncluster <- sub(".*:", "", fit)
            parameters <- DECODE_result[[detail_tail]]$all_fits[[detail_Ncluster]]$best$parameters_validation
            tail_status <- ifelse(detail_tail == "fits_with_tail", TRUE, FALSE)
            component_distributions <- DECODE_result[[detail_tail]]$all_fits[[detail_Ncluster]]$best$component_distributions_validation
            criteria_validation_index <- DECODE_result[[detail_tail]]$all_fits[[detail_Ncluster]]$best$criteria_validation_index
        }
        vec_SFS_real <- vec_SFS_real[criteria_validation_index, ]
    }
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
        "Cluster 2" = "#0072B2",
        "Cluster 3" = "#E69F00",
        "Cluster 4" = "#CC79A7",
        "Cluster 5" = "#56B4E9",
        "Cluster 6" = "#009E73",
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
        for (j in 1:length(mutation_table$Alt_count)) {
            no_variant <- mutation_table$Alt_count[j]
            no_total <- mutation_table$Ref_count[j] + mutation_table$Alt_count[j]
            if (no_variant >= min_variant_read & no_total >= min_total_read & no_total <= max_total_read) {
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
    if (tail_status) {
        SFS_neutral <- compute_SFS(
            A = vec_A[1],
            vec_K = c(),
            component_distributions = component_distributions
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
                component_distributions = component_distributions
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
        geom_line(data = df_fit, aes(x = frequency, y = count, group = fill, color = fill), position = "stack", size = DECODE_linewidth, show.legend = FALSE) +
        geom_bar(data = df_data, aes(x = frequency, y = count, fill = fill), stat = "identity", width = 0.5 / SFS_totalsteps) +
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
    if (SFS_limit) {
        x_min <- min(df_data$frequency[which(df_data$count > 0)])
        x_max <- max(df_data$frequency[which(df_data$count > 0)])
        p <- p + scale_x_continuous(limits = c(x_min, x_max))
    }
    return(p)
}

DECODE_plot_readcounts <- function(DECODE_result,
                                   freq_cutoff = 10) {
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
    #---Reduce the distribution to region satisfying the frequency cutoff
    max_total_read <- vec_min_total_read[which(readcount_distribution$freq[which(readcount_distribution$min_variant_read == vec_min_variant_read[1])] < freq_cutoff)[1]]
    max_variant_read <- vec_min_variant_read[which(readcount_distribution$freq[which(readcount_distribution$min_total_read == vec_min_total_read[1])] < freq_cutoff)[1]]
    readcount_distribution <- readcount_distribution[which(readcount_distribution$min_total_read <= max_total_read & readcount_distribution$min_variant_read <= max_variant_read), ]
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
            size = 6,
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
            size = 6,
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
            size = 6,
            color = color_validation,
            bg.color = "white",
            fontface = "bold"
        ) +
        scale_fill_gradientn(
            colors = c("#0072B2", "#56B4E9", "#009E73", "#E69F00", "#D55E00"),
            name = "% mutations retained"
        ) +
        theme_minimal() +
        labs(title = DECODE_result$sample_id, x = "Minimum total read count", y = "Minimum variant read count") +
        theme(
            text = element_text(size = 30),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5),
            legend.key.width = unit(1.5, "cm"),
        )
    return(p)
}

analysis_CINner <- function(groundtruth_df,
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
    if (is_mobster) characteristic_df$mobster_cluster_count <- as.factor(mobster_df$Cluster_count)
    if (is_decode) characteristic_df$decode_cluster_count <- as.factor(decode_df$Cluster_count)
    #---Find distributions of neutral tail power in each method
    if (is_mobster) {
        mobster_df$alpha <- mobster_df$Tail_power
        mobster_mean <- median(mobster_df$alpha, na.rm = TRUE)
    }
    if (is_decode) {
        decode_df$alpha <- decode_df$Tail_power
        decode_mean <- median(decode_df$alpha, na.rm = TRUE)
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
        cluster_count_min <- min(cluster_count_min, min(mobster_df$Cluster_count, na.rm = TRUE))
        cluster_count_max <- max(cluster_count_max, max(mobster_df$Cluster_count, na.rm = TRUE))
    }
    if (is_decode) {
        cluster_count_min <- min(cluster_count_min, min(decode_df$Cluster_count, na.rm = TRUE))
        cluster_count_max <- max(cluster_count_max, max(decode_df$Cluster_count, na.rm = TRUE))
    }
    if (is_mobster) {
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
        mobster_cluster_df <- mobster_cluster_df[which(mobster_cluster_df$Simulation %in% mobster_df$Simulation[which(mobster_df$Succeed == TRUE)]), ]
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
    if (is_mobster) {
        png(paste0(folder_workplace, "Comparison_0_sample_info_vs_mobster_tail_detection.png"), res = 150, width = 30, height = 30, units = "in")
        p <- ggplot() +
            geom_point(
                data = characteristic_df[which(!is.na(characteristic_df$mobster_tail)), ],
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
                legend.justification = c(0, 0.5),
                plot.margin = margin(0, 2, 0, 0, "cm")
            )
        print(p)
        dev.off()
    }
    #---Plot distributions of MOBSTER's cluster count w.r.t. sample purity and coverage
    if (is_mobster) {
        png(paste0(folder_workplace, "Comparison_0_sample_info_vs_mobster_cluster_count.png"), res = 150, width = 30, height = 30, units = "in")
        p <- ggplot() +
            geom_point(
                data = characteristic_df[which(!is.na(characteristic_df$mobster_cluster_count)), ],
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
                legend.justification = c(0, 0.5),
                plot.margin = margin(0, 2, 0, 0, "cm")
            )
        print(p)
        dev.off()
    }
    #---Plot distributions of DECODE's cluster count w.r.t. sample purity and coverage
    if (is_decode) {
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
                legend.justification = c(0, 0.5),
                plot.margin = margin(0, 2, 0, 0, "cm")
            )
        print(p)
        dev.off()
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
            geom_text(data = cluster_parameter_df[cluster_parameter_df$Parameter == "p", ], aes(x = Value_TRUTH, y = Value_INFERRED, label = Simulation, color = Method), size = 20, vjust = 0.5, hjust = 0.5)
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
                geom_text(data = mini_cluster_parameter_df[mini_cluster_parameter_df$Parameter == "p", ], aes(x = Value_TRUTH, y = Value_INFERRED, label = Simulation, color = Method), size = 20, vjust = 0.5, hjust = 0.5)
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
            plot.margin = margin(0, 2, 0, 0, "cm")
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
                plot.margin = margin(0, 2, 0, 0, "cm")
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
            geom_histogram(data = mobster_df[!is.na(mobster_df$alpha), ], aes(x = alpha, fill = "MOBSTER"), alpha = 0.5) +
            geom_vline(xintercept = mobster_mean, color = "orange", linetype = "dashed", size = 3)
    }
    if (is_decode) {
        p <- p +
            geom_histogram(data = decode_df[!is.na(decode_df$Tail_power), ], aes(x = Tail_power, fill = "DECODE"), alpha = 0.5) +
            geom_vline(xintercept = decode_mean, color = "purple", linetype = "dashed", size = 3)
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
            legend.justification = c(0, 0.5),
            plot.margin = margin(0, 2, 0, 0, "cm")
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
            legend.justification = c(0, 0.5),
            plot.margin = margin(0, 2, 0, 0, "cm")
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
    #---Plot DECODE neutral compartment sensitivity studies
    png(paste0(folder_workplace, "Comparison_7_DECODE_tail_sensitivity.png"), res = 150, width = 30, height = 30, units = "in")
    p <- ggplot() +
        geom_point(
            data = decode_df,
            aes(x = Tail_sensitivity_Bayesian_pi0_std, y = Tail_sensitivity_Morris_pi0_mean_abs, fill = "DECODE", color = "DECODE"),
            alpha = 0.5, size = 20
        ) +
        xlab(expression(paste("Bayesian standard deviation of ", pi[0]))) +
        ylab(expression(paste("Morris mean EE for ", pi[0]))) +
        scale_fill_manual(values = method_color_scheme, name = "") +
        scale_color_manual(values = method_color_scheme, name = "") +
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
            geom_text(
                data = decode_df,
                aes(x = Tail_sensitivity_Bayesian_pi0_std, y = Tail_sensitivity_Morris_pi0_mean_abs, label = Simulation, color = "DECODE"), size = 20, vjust = 0.5, hjust = 0.5
            )
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
    library(dplyr)
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
    cancer_type_color <- c(
        "CNS-PiloAstro" = "lightblue",
        "Liver-HCC" = "deeppink",
        "Kidney-RCC" = "forestgreen",
        "Prost-AdenoCA" = "darkviolet",
        "Lymph-BNHL" = "darkorange",
        "Panc-Endocrine" = "darkturquoise",
        "Panc-AdenoCA" = "darkslategrey",
        "CNS-Medullo" = "darksalmon",
        "Breast-AdenoCA" = "goldenrod1",
        "Stomach-AdenoCA" = "darkolivegreen",
        "Skin-Melanoma" = "sienna",
        "Lymph-CLL" = "burlywood",
        "Eso-AdenoCA" = "plum",
        "Myeloid-AML" = "palevioletred",
        "Head-SCC" = "mediumaquamarine",
        "Billary-AdenoCA" = "greenyellow",
        "Ovary-AdenoCA" = "cornflowerblue",
        "Billary-AdenoCA" = "magenta",
        "Bone-Benign" = "green",
        "Bone-Epith" = "bisque4",
        "Bone-Osteosarc" = "ivory3",
        "Breast-DCIS" = "navy",
        "Breast-LobularCA" = "lightcoral",
        "Myeloid-MDS" = "yellow3",
        "Myeloid-MPN" = "steelblue4"
    )
    cluster_count_color_scheme <- c(
        "1" = "#E69F00",
        "2" = "#56B4E9",
        "3" = "#009E73",
        "4" = "#F0E442",
        "5" = "#CC79A7"
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
    ####################################################################
    ####################################################################
    ####################################################################
    ####################################################################
    ####################################################################
    ####################################################################
    ####################################################################
    #---Plot a bar graph representing Cancer Type and mutation count
    merged_data <- merge(sample_information_df, mobster_df, by.x = "aliquot_id", by.y = "Sample") %>%
        dplyr::select(aliquot_id, histology_abbreviation, Mutation_count_in_fitting)

    average_mutation_data <- merged_data %>%
        dplyr::group_by(histology_abbreviation) %>%
        dplyr::summarise(Average_Mutation_Count = mean(Mutation_count_in_fitting, na.rm = TRUE))

    # Plot ordered in descending order of mutation count
    average_mutation_data_filtered <- average_mutation_data %>%
        dplyr::arrange(desc(Average_Mutation_Count)) %>%
        dplyr::mutate(histology_abbreviation = factor(histology_abbreviation, levels = unique(histology_abbreviation)))

    ordered_histology_abbreviation <- average_mutation_data_filtered$histology_abbreviation

    # Create the plot
    p <- ggplot(average_mutation_data_filtered, aes(x = histology_abbreviation, y = Average_Mutation_Count, fill = histology_abbreviation)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_fill_manual(values = cancer_type_color, guide = FALSE) +
        labs(x = "", y = "", title = "Average mutation count") +
        scale_y_log10() +
        theme(
            text = element_text(size = 40),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(angle = 90, hjust = 0.5),
            panel.background = element_rect(fill = "white", colour = "white"),
            plot.title = element_text(hjust = 0, size = 40),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5),
            plot.margin = margin(0.5, 0.5, 0, 2, "cm")
        )
    png(paste0(folder_workplace, "ICGC_0_mutation_count.png"), res = 150, width = 30, height = 15, units = "in")
    print(p)
    dev.off()

    #---Plot a bar graph representing Cancer Type and sample count
    # Extract Data
    cancer_type_counts <- sample_information_df %>%
        group_by(histology_abbreviation) %>%
        summarise(samplecount = n(), .groups = "drop")

    # Plot ordered in descending order of mutation count
    cancer_type_counts$histology_abbreviation <- factor(cancer_type_counts$histology_abbreviation, levels = ordered_histology_abbreviation)

    # Plotting bar graph
    p <- ggplot(cancer_type_counts, aes(x = histology_abbreviation, y = samplecount, fill = histology_abbreviation)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_fill_manual(values = cancer_type_color, guide = FALSE) +
        labs(x = "", y = "", title = "Sample count") +
        theme(
            text = element_text(size = 40),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(angle = 90, hjust = 0.5),
            panel.background = element_rect(fill = "white", colour = "white"),
            plot.title = element_text(hjust = 0, size = 40),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5),
            plot.margin = margin(0.5, 0.5, 0, 2, "cm")
        )
    # Save the plot
    png(paste0(folder_workplace, "ICGC_0_sample_count.png"), res = 150, width = 30, height = 15, units = "in")
    print(p)
    dev.off()

    #---Plot a bar graph of cancer types and tail detection - MOBSTER
    # Extract and merge required data
    merged_data <- merge(sample_information_df, mobster_df, by.x = "aliquot_id", by.y = "Sample") %>%
        dplyr::select(aliquot_id, histology_abbreviation, Tail)
    # Calculate Percentages
    merged_data <- merged_data %>%
        dplyr::group_by(histology_abbreviation) %>%
        dplyr::count(Tail) %>%
        dplyr::mutate(percentage = n / sum(n)) %>%
        dplyr::ungroup()
    # Create a data frame with all combinations of histology_abbreviation and Tail
    all_combinations <- expand.grid(
        histology_abbreviation = unique(merged_data$histology_abbreviation),
        Tail = c(TRUE, FALSE)
    )
    merged_data_complete <- merge(all_combinations, merged_data, by = c("histology_abbreviation", "Tail"), all.x = TRUE)
    merged_data_complete$percentage[is.na(merged_data_complete$percentage)] <- 0

    # Plot ordered in descending order of mutation count
    merged_data_complete$histology_abbreviation <- factor(merged_data_complete$histology_abbreviation, levels = ordered_histology_abbreviation)

    # Column for coloring for No-tail
    merged_data_complete$fill_group <- ifelse(merged_data_complete$Tail == FALSE, "Non-Tail", as.character(merged_data_complete$histology_abbreviation))
    merged_data_complete$fill_group <- factor(merged_data_complete$fill_group, levels = c("Non-Tail", unique(as.character(merged_data_complete$histology_abbreviation[merged_data_complete$Tail == TRUE]))))

    # Plotting bar graph
    p <- ggplot(merged_data_complete, aes(x = histology_abbreviation, y = percentage, fill = fill_group, alpha = Tail)) +
        geom_bar(stat = "identity", position = "stack", show.legend = TRUE) +
        scale_y_continuous(labels = scales::percent_format()) +
        scale_fill_manual(values = c(cancer_type_color, "Non-Tail" = "grey"), guide = FALSE) +
        scale_alpha_manual(
            values = c("TRUE" = 1, "FALSE" = 0.3),
            guide = guide_legend(
                title = "MOBSTER - tail detection",
                override.aes = list(fill = "grey")
            )
        ) +
        geom_text(
            aes(
                label = ifelse(percentage > 0, as.character(n), ""),
                y = ifelse(Tail, 0, 1) # Adjust y based on Tail value
            ),
            position = position_dodge(width = 0.9), size = 5
        ) +
        labs(x = "", y = "") +
        theme(
            text = element_text(size = 40),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(angle = 90, hjust = 0.5),
            panel.background = element_rect(fill = "white", colour = "white"),
            plot.title = element_text(hjust = 0),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5),
            plot.margin = margin(0.5, 0.5, 0, 2, "cm")
        )
    # Save the plot
    png(paste0(folder_workplace, "ICGC_1_tail_detection_MOBSTER.png"), res = 150, width = 30, height = 15, units = "in")
    print(p)
    dev.off()

    #---Plot a bar graph of cancer types and tail detection - DECODE
    # Extract and merge required data
    merged_data <- merge(sample_information_df, decode_df, by.x = "aliquot_id", by.y = "Sample") %>%
        dplyr::select(aliquot_id, histology_abbreviation, Tail)
    # Calculate Percentages
    merged_data <- merged_data %>%
        dplyr::group_by(histology_abbreviation) %>%
        dplyr::count(Tail) %>%
        dplyr::mutate(percentage = n / sum(n)) %>%
        dplyr::ungroup()
    # Create a data frame with all combinations of histology_abbreviation and Tail
    all_combinations <- expand.grid(
        histology_abbreviation = unique(merged_data$histology_abbreviation),
        Tail = c(TRUE, FALSE)
    )
    merged_data_complete <- merge(all_combinations, merged_data, by = c("histology_abbreviation", "Tail"), all.x = TRUE)
    merged_data_complete$percentage[is.na(merged_data_complete$percentage)] <- 0

    # Plot ordered in descending order of mutation count
    merged_data_complete$histology_abbreviation <- factor(merged_data_complete$histology_abbreviation, levels = ordered_histology_abbreviation)

    # Column for coloring for No-tail
    merged_data_complete$fill_group <- ifelse(merged_data_complete$Tail == FALSE, "Non-Tail", as.character(merged_data_complete$histology_abbreviation))
    merged_data_complete$fill_group <- factor(merged_data_complete$fill_group, levels = c("Non-Tail", unique(as.character(merged_data_complete$histology_abbreviation[merged_data_complete$Tail == TRUE]))))

    # Plotting bar graph
    p <- ggplot(merged_data_complete, aes(x = histology_abbreviation, y = percentage, fill = fill_group, alpha = Tail)) +
        geom_bar(stat = "identity", position = "stack", show.legend = TRUE) +
        scale_y_continuous(labels = scales::percent_format()) +
        scale_fill_manual(values = c(cancer_type_color, "Non-Tail" = "grey"), guide = FALSE) +
        scale_alpha_manual(
            values = c("TRUE" = 1, "FALSE" = 0.3),
            guide = guide_legend(
                title = "DECODE - tail detection",
                override.aes = list(fill = "grey")
            )
        ) +
        geom_text(
            aes(
                label = ifelse(percentage > 0, as.character(n), ""),
                y = ifelse(Tail, 0, 1) # Adjust y based on Tail value
            ),
            position = position_dodge(width = 0.9), size = 5
        ) +
        labs(x = "", y = "") +
        theme(
            text = element_text(size = 40),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(angle = 90, hjust = 0.5),
            panel.background = element_rect(fill = "white", colour = "white"),
            plot.title = element_text(hjust = 0),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5),
            plot.margin = margin(0.5, 0.5, 0, 2, "cm")
        )
    # Save the plot
    png(paste0(folder_workplace, "ICGC_1_tail_detection_DECODE.png"), res = 150, width = 30, height = 15, units = "in")
    print(p)
    dev.off()

    #---Plot a bar graph representing Cancer Type and cluster count detection - MOBSTER
    merged_data <- merge(sample_information_df, mobster_df, by.x = "aliquot_id", by.y = "Sample") %>%
        dplyr::select(aliquot_id, histology_abbreviation, Cluster_count)
    aggregated_data <- merged_data %>%
        dplyr::group_by(histology_abbreviation, Cluster_count) %>%
        dplyr::summarise(Sample_Count = n(), .groups = "drop")
    all_combinations <- expand.grid(
        histology_abbreviation = unique(aggregated_data$histology_abbreviation),
        Cluster_count = 1:3
    )
    merged_data_complete <- merge(all_combinations, aggregated_data, by = c("histology_abbreviation", "Cluster_count"), all.x = TRUE)
    merged_data_complete$Sample_Count[is.na(merged_data_complete$Sample_Count)] <- 0

    # Calculate total samples per cancer type
    aggregated_data_totals <- aggregated_data %>%
        dplyr::group_by(histology_abbreviation) %>%
        dplyr::summarise(Total_Sample_Count = sum(Sample_Count), .groups = "drop")

    # Merge to get total sample count per cancer type in the complete dataset
    merged_data_complete <- merge(merged_data_complete, aggregated_data_totals, by = "histology_abbreviation")

    # Calculate percentage
    merged_data_complete$Percentage <- (merged_data_complete$Sample_Count / merged_data_complete$Total_Sample_Count)

    # Plot ordered in descending order of mutation count
    merged_data_complete$histology_abbreviation <- factor(merged_data_complete$histology_abbreviation, levels = ordered_histology_abbreviation)

    # Create the plot
    p <- ggplot(merged_data_complete, aes(x = histology_abbreviation, y = Percentage, fill = factor(Cluster_count, levels = rev(unique(Cluster_count))))) +
        geom_bar(stat = "identity", position = "stack", show.legend = TRUE) +
        scale_y_continuous(labels = scales::percent_format()) +
        scale_fill_manual(values = cluster_count_color_scheme, name = "MOBSTER - cluster count") +
        labs(x = "", y = "") +
        theme(
            text = element_text(size = 40),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(angle = 90, hjust = 0.5),
            panel.background = element_rect(fill = "white", colour = "white"),
            plot.title = element_text(hjust = 0, size = 40),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5),
            plot.margin = margin(0.5, 0.5, 0, 2, "cm")
        )
    png(paste0(folder_workplace, "ICGC_2_cluster_count_MOBSTER.png"), res = 150, width = 30, height = 15, units = "in")
    print(p)
    dev.off()

    #---Plot a bar graph representing Cancer Type and cluster count detection - DECODE
    merged_data <- merge(sample_information_df, decode_df, by.x = "aliquot_id", by.y = "Sample") %>%
        dplyr::select(aliquot_id, histology_abbreviation, Cluster_count)

    aggregated_data <- merged_data %>%
        dplyr::group_by(histology_abbreviation, Cluster_count) %>%
        dplyr::summarise(Sample_Count = n(), .groups = "drop")

    all_combinations <- expand.grid(
        histology_abbreviation = unique(aggregated_data$histology_abbreviation),
        Cluster_count = 1:3
    )
    merged_data_complete <- merge(all_combinations, aggregated_data, by = c("histology_abbreviation", "Cluster_count"), all.x = TRUE)
    merged_data_complete$Sample_Count[is.na(merged_data_complete$Sample_Count)] <- 0

    # Calculate total samples per cancer type
    aggregated_data_totals <- aggregated_data %>%
        dplyr::group_by(histology_abbreviation) %>%
        dplyr::summarise(Total_Sample_Count = sum(Sample_Count), .groups = "drop")

    # Merge to get total sample count per cancer type in the complete dataset
    merged_data_complete <- merge(merged_data_complete, aggregated_data_totals, by = "histology_abbreviation")

    # Calculate percentage
    merged_data_complete$Percentage <- (merged_data_complete$Sample_Count / merged_data_complete$Total_Sample_Count)

    # Plot ordered in descending order of mutation count
    merged_data_complete$histology_abbreviation <- factor(merged_data_complete$histology_abbreviation, levels = ordered_histology_abbreviation)

    # Create the plot
    p <- ggplot(merged_data_complete, aes(x = histology_abbreviation, y = Percentage, fill = factor(Cluster_count, levels = rev(unique(Cluster_count))))) +
        geom_bar(stat = "identity", position = "stack", show.legend = TRUE) +
        scale_y_continuous(labels = scales::percent_format()) +
        scale_fill_manual(values = cluster_count_color_scheme, name = "DECODE - cluster count") +
        labs(x = "", y = "") +
        theme(
            text = element_text(size = 40),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(angle = 90, hjust = 0.5),
            panel.background = element_rect(fill = "white", colour = "white"),
            plot.title = element_text(hjust = 0, size = 40),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5),
            plot.margin = margin(0.5, 0.5, 0, 2, "cm")
        )
    png(paste0(folder_workplace, "ICGC_2_cluster_count_DECODE.png"), res = 150, width = 30, height = 15, units = "in")
    print(p)
    dev.off()
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
    png(paste0(folder_workplace, "ICGC_3_truncal_frequency_vs_sample_purity_MOBSTER.png"), res = 150, width = 30, height = 30, units = "in")

    common_range <- range(
        c(
            df_truncal_frequency_vs_sample_purity$Purity,
            df_truncal_frequency_vs_sample_purity$Truncal_frequency
        ),
        na.rm = TRUE
    )

    # Plot ordered in descending order of mutation count

    p <- ggplot() +
        geom_point(
            data = df_truncal_frequency_vs_sample_purity,
            aes(x = Purity, y = Truncal_frequency, fill = Method, color = Method),
            alpha = 0.5, size = 20
        ) +
        xlim(c(0, 1)) +
        ylim(c(0, 1)) +
        # xlim(common_range) +
        # ylim(common_range) +
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
            plot.margin = margin(0, 0.5, 0.5, 0.5, "cm")
        )
    if (text_notation) {
        p <- p +
            geom_text(data = df_truncal_frequency_vs_sample_purity[which(df_truncal_frequency_vs_sample_purity$Within_bounds != "Correct"), ], aes(x = Purity, y = Truncal_frequency, label = Sample_short, color = Method), size = 20, vjust = 0, hjust = 0, angle = 45)
    }
    # Extracting the calculated percentages
    percent_within <- 100 * round(length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Correct")) / nrow(df_truncal_frequency_vs_sample_purity), 2)
    percent_above <- 100 * round(length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Above")) / nrow(df_truncal_frequency_vs_sample_purity), 2)
    percent_below <- 100 * round(length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Below")) / nrow(df_truncal_frequency_vs_sample_purity), 2)

    # X position and offset
    x_pos <- 0.08
    offset <- 0.02

    p <- p +
        # Annotation for Within, Above, and Below percentages & counts
        annotate("text", x = 0, y = 0, label = paste0(percent_within, "% (n=", length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Correct")), ")"), size = 20, colour = "black", angle = 45, hjust = 0) +
        annotate("text", x = 0, y = 0.1, label = paste0(percent_above, "% (n=", length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Above")), ")"), size = 20, colour = "black", angle = 45, hjust = 0) +
        annotate("text", x = 0.1, y = 0, label = paste0(percent_below, "% (n=", length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Below")), ")"), size = 20, colour = "black", angle = 45, hjust = 0)
    # annotate("text", x = x_pos - 0.03, y = x_pos + offset - 0.03, label = paste0(percent_within, "% (n=", length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Correct")), ")"), size = 20, colour = "black", angle = 45) +
    # annotate("text", x = x_pos - 0.05, y = x_pos + bound_truncal_frequency_vs_sample_purity + offset - 0.03, label = paste0(percent_above, "% (n=", length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Above")), ")"), size = 20, colour = "black", angle = 45) +
    # annotate("text", x = x_pos + 0.02, y = x_pos - bound_truncal_frequency_vs_sample_purity - offset, label = paste0(percent_below, "% (n=", length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Below")), ")"), size = 20, colour = "black", angle = 45)

    p <- p +
        geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 3, alpha = 0.5) +
        geom_abline(intercept = bound_truncal_frequency_vs_sample_purity, slope = 1, color = "black", linewidth = 3, alpha = 0.5, linetype = "dashed") +
        geom_abline(intercept = -bound_truncal_frequency_vs_sample_purity, slope = 1, color = "black", linewidth = 3, alpha = 0.5, linetype = "dashed")
    print(p)
    dev.off()
    cat(paste0("\nSamples with truncal cluster frequency within ", 100 * bound_truncal_frequency_vs_sample_purity, "% of sample purity: ", length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Correct")), " (", 100 * round(length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Correct")) / nrow(df_truncal_frequency_vs_sample_purity), 2), "%)", "\n"))
    cat(paste0("Samples with truncal cluster frequency > sample purity + ", 100 * bound_truncal_frequency_vs_sample_purity, "%:       ", length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Above")), " (", 100 * round(length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Above")) / nrow(df_truncal_frequency_vs_sample_purity), 2), "%)", "\n"))
    write.table(df_truncal_frequency_vs_sample_purity$Sample[which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Above")], file = "MOBSTER>5%.txt", quote = FALSE, row.names = FALSE)
    cat(paste0("Samples with truncal cluster frequency < sample purity - ", 100 * bound_truncal_frequency_vs_sample_purity, "%:       ", length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Below")), " (", 100 * round(length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Below")) / nrow(df_truncal_frequency_vs_sample_purity), 2), "%)", "\n\n"))
    write.table(df_truncal_frequency_vs_sample_purity$Sample[which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Below")], file = "MOBSTER<5%.txt", quote = FALSE, row.names = FALSE)
    #---Find comparison of truncal cluster frequency against sample purity - DECODE
    df_truncal_frequency_vs_sample_purity <- data.frame()
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
    png(paste0(folder_workplace, "ICGC_3_truncal_frequency_vs_sample_purity_DECODE.png"), res = 150, width = 30, height = 30, units = "in")
    common_range <- range(
        c(
            df_truncal_frequency_vs_sample_purity$Purity,
            df_truncal_frequency_vs_sample_purity$Truncal_frequency
        ),
        na.rm = TRUE
    )

    # Plot ordered in descending order of mutation count
    p <- ggplot() +
        geom_point(
            data = df_truncal_frequency_vs_sample_purity,
            aes(x = Purity, y = Truncal_frequency, fill = Method, color = Method),
            alpha = 0.5, size = 20
        ) +
        xlim(c(0, 1)) +
        ylim(c(0, 1)) +
        # xlim(common_range) +
        # ylim(common_range) +
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
            plot.margin = margin(0, 0.5, 0.5, 0.5, "cm")
        )
    if (text_notation) {
        p <- p +
            geom_text(data = df_truncal_frequency_vs_sample_purity[which(df_truncal_frequency_vs_sample_purity$Within_bounds != "Correct"), ], aes(x = Purity, y = Truncal_frequency, label = Sample_short, color = Method), size = 20, vjust = 0, hjust = 0, angle = 45)
    }

    # Extracting the calculated percentages
    percent_within <- 100 * round(length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Correct")) / nrow(df_truncal_frequency_vs_sample_purity), 2)
    percent_above <- 100 * round(length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Above")) / nrow(df_truncal_frequency_vs_sample_purity), 2)
    percent_below <- 100 * round(length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Below")) / nrow(df_truncal_frequency_vs_sample_purity), 2)
    p <- p +
        annotate("text", x = 0.0, y = 0.0, label = paste0(percent_within, "% (n=", length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Correct")), ")"), size = 20, colour = "black", angle = 45, hjust = 0) +
        annotate("text", x = 0, y = 0.1, label = paste0(percent_above, "% (n=", length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Above")), ")"), size = 20, colour = "black", angle = 45, hjust = 0) +
        annotate("text", x = 0.1, y = 0, label = paste0(percent_below, "% (n=", length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Below")), ")"), size = 20, colour = "black", angle = 45, hjust = 0)
    p <- p +
        # geom_segment(aes(x = 0.2, y = 0.2, xend = 1, yend = 1), color = "black", linewidth = 2) +
        geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 3, alpha = 0.5) +
        geom_abline(intercept = bound_truncal_frequency_vs_sample_purity, slope = 1, color = "black", linewidth = 3, alpha = 0.5, linetype = "dashed") +
        geom_abline(intercept = -bound_truncal_frequency_vs_sample_purity, slope = 1, color = "black", linewidth = 3, alpha = 0.5, linetype = "dashed")
    print(p)
    cat(paste0("\nDECODE samples with truncal cluster frequency within ", 100 * bound_truncal_frequency_vs_sample_purity, "% of sample purity: ", length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Correct")), " (", 100 * round(length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Correct")) / nrow(df_truncal_frequency_vs_sample_purity), 2), "%)", "\n"))
    cat(paste0("DECODE samples with truncal cluster frequency > sample purity + ", 100 * bound_truncal_frequency_vs_sample_purity, "%:       ", length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Above")), " (", 100 * round(length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Above")) / nrow(df_truncal_frequency_vs_sample_purity), 2), "%)", "\n"))
    write.table(df_truncal_frequency_vs_sample_purity$Sample[which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Above")], file = "DECODE>5%.txt", quote = FALSE, row.names = FALSE)
    cat(paste0("DECODE samples with truncal cluster frequency < sample purity - ", 100 * bound_truncal_frequency_vs_sample_purity, "%:       ", length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Below")), " (", 100 * round(length(which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Below")) / nrow(df_truncal_frequency_vs_sample_purity), 2), "%)", "\n\n"))
    write.table(df_truncal_frequency_vs_sample_purity$Sample[which(df_truncal_frequency_vs_sample_purity$Within_bounds == "Below")], file = "DECODE<5%.txt", quote = FALSE, row.names = FALSE)
    dev.off()

    #---Plot a Violin Plot of Cancer type and neutral tail power MOBSTER
    merged_data <- merge(sample_information_df, mobster_df, by.x = "aliquot_id", by.y = "Sample") %>%
        dplyr::select(aliquot_id, histology_abbreviation, Tail_power)

    # Plot ordered in descending order of mutation count
    merged_data$histology_abbreviation <- factor(merged_data$histology_abbreviation, levels = ordered_histology_abbreviation)

    # Create the violin plot
    p <- ggplot(merged_data, aes(x = histology_abbreviation, y = Tail_power, fill = histology_abbreviation)) +
        geom_violin(width = 0.7, scale = "width") +
        geom_boxplot(width = 0.1, fill = "white", color = "black", alpha = 0.7) +
        scale_fill_manual(values = cancer_type_color, guide = FALSE) +
        labs(x = "", y = "", title = "MOBSTER - neutral tail power") +
        theme(
            text = element_text(size = 40),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(angle = 90, hjust = 0.5),
            panel.background = element_rect(fill = "white", colour = "white"),
            plot.title = element_text(hjust = 0, size = 40),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5),
            plot.margin = margin(0.5, 0.5, 0, 2, "cm")
        )

    p <- p +
        geom_abline(intercept = 2, slope = 0, color = "black", linewidth = 2, linetype = "dashed")
    png(paste0(folder_workplace, "ICGC_4_neutral_tail_power_MOBSTER.png"), res = 150, width = 30, height = 15, units = "in")
    print(p)
    dev.off()

    #---Plot a Violin Plot of Cancer type and neutral tail power DECODE
    merged_data <- merge(sample_information_df, decode_df, by.x = "aliquot_id", by.y = "Sample") %>%
        dplyr::select(aliquot_id, histology_abbreviation, Tail_power)

    # Plot ordered in descending order of mutation count
    merged_data$histology_abbreviation <- factor(merged_data$histology_abbreviation, levels = ordered_histology_abbreviation)

    # Create the violin plot
    p <- ggplot(merged_data, aes(x = histology_abbreviation, y = Tail_power, fill = histology_abbreviation)) +
        geom_violin(width = 0.7, scale = "width") +
        geom_boxplot(width = 0.1, fill = "white", color = "black", alpha = 0.7) +
        scale_fill_manual(values = cancer_type_color, guide = FALSE) +
        labs(x = "", y = "", title = "DECODE - neutral tail power") +
        theme(
            text = element_text(size = 40),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(angle = 90, hjust = 0.5),
            panel.background = element_rect(fill = "white", colour = "white"),
            plot.title = element_text(hjust = 0, size = 40),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5),
            plot.margin = margin(0.5, 0.5, 0, 2, "cm")
        )

    p <- p +
        geom_abline(intercept = 2, slope = 0, color = "black", linewidth = 2, linetype = "dashed")
    png(paste0(folder_workplace, "ICGC_4_neutral_tail_power_DECODE.png"), res = 150, width = 30, height = 15, units = "in")
    print(p)
    dev.off()
}
analysis_DREAM <- function(sample_information_df,
                           sample_excel_df,
                           decode_df = NULL,
                           decode_BIC_df = NULL,
                           text_notation = FALSE,
                           folder_workplace = "") {
    library(ggplot2)
    library(dplyr)
    library(readxl)
    library(readr)
    library(scales)
    is_decode <- !is.null(decode_df)
    method_color_scheme <- c(
        "MOBSTER" = "darkorange2",
        "DECODE" = "magenta4",
        "DECODE-BIC" = "blue3"
    )

    #---Plot a scatter box plot of Different algorithms and Scores - 1A
    # Convert `Raw score` to numeric for DECODE scores
    decode_scores <- read_csv("Calculated_Scores_1A.csv") %>%
        dplyr::select(Score_1A) %>%
        dplyr::rename(`Raw score` = Score_1A) %>%
        dplyr::mutate(SubChallenge = "sc1A", Entry = "DECODE") %>%
        dplyr::mutate(`Raw score` = as.numeric(`Raw score`))

    # Convert `Raw score` to numeric for DECODE-BIC scores
    decode_bic_scores <- read_csv("Calculated_Scores_1A_BIC.csv") %>%
        dplyr::select(Score_1A) %>%
        dplyr::rename(`Raw score` = Score_1A) %>%
        dplyr::mutate(SubChallenge = "sc1A", Entry = "DECODE-BIC") %>%
        dplyr::mutate(`Raw score` = as.numeric(`Raw score`))

    # Prepare other data
    sample_excel_df <- read_excel("DREAM2024_supplementary.xlsx") %>%
        dplyr::filter(SubChallenge == "sc1A") %>%
        dplyr::select(Entry, `Raw score`) %>%
        dplyr::mutate(`Raw score` = as.numeric(as.character(`Raw score`)))

    # Merge all data
    combined_data <- bind_rows(sample_excel_df, decode_scores, decode_bic_scores)

    # Ensure the Entry column is a factor with the correct levels
    combined_data$Entry <- factor(combined_data$Entry,
        levels = c(
            "random-clone", "5641737", "5964755", "6052744", "6087270", "6181022",
            "6087309", "6187040", "6087362", "6185506", "6181501", "6184520",
            "6184761", "6187127", "6204327", "SOTA1", "SOTA2", "rand_q50",
            "DECODE", "DECODE-BIC"
        )
    )

    # Convert `Raw score` to numeric
    combined_data$`Raw score` <- as.numeric(as.character(combined_data$`Raw score`))

    # Check for any conversion errors
    if (any(is.na(combined_data$`Raw score`))) {
        warning("NA introduced by coercion")
    }

    # Filter out NA values in the Entry column
    combined_data <- combined_data %>% filter(!is.na(Entry))

    # Plot
    p <- ggplot(combined_data, aes(x = Entry, y = `Raw score`, fill = Entry)) +
        geom_point(position = position_jitterdodge(), aes(color = Entry)) +
        geom_boxplot(alpha = 0.5, outlier.shape = NA) +
        scale_fill_manual(values = method_color_scheme) +
        labs(title = "Scatter Box Plot of Algorithms and Scores", x = "Algorithm", y = "Score") +
        theme_minimal() +
        theme(
            text = element_text(size = 40),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(angle = 90, hjust = 0.5),
            panel.background = element_rect(fill = "white", colour = "white"),
            plot.title = element_text(hjust = 0, size = 40),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "none",
            plot.margin = margin(0.5, 0.5, 0, 2, "cm")
        ) +
        scale_y_continuous(breaks = pretty_breaks(n = 5))

    # Save the plot
    png(paste0(folder_workplace, "DREAM2024_scatter_box_plot_Score1A.png"), res = 150, width = 30, height = 15, units = "in")
    print(p)
    dev.off()

    #---Plot a scatter box plot of Different algorithms and Scores - 1B
    # Convert `Raw score` to numeric for DECODE scores
    decode_scores <- read_csv("Calculated_Scores_1B.csv") %>%
        dplyr::select(Score_1B) %>%
        dplyr::rename(`Raw score` = Score_1B) %>%
        dplyr::mutate(SubChallenge = "sc1A", Entry = "DECODE") %>%
        dplyr::mutate(`Raw score` = as.numeric(`Raw score`))

    # Convert `Raw score` to numeric for DECODE-BIC scores
    decode_bic_scores <- read_csv("Calculated_Scores_1B_BIC.csv") %>%
        dplyr::select(Score_1B) %>%
        dplyr::rename(`Raw score` = Score_1B) %>%
        dplyr::mutate(SubChallenge = "sc1A", Entry = "DECODE-BIC") %>%
        dplyr::mutate(`Raw score` = as.numeric(`Raw score`))

    # Prepare other data
    sample_excel_df <- read_excel("DREAM2024_supplementary.xlsx") %>%
        dplyr::filter(SubChallenge == "sc1B") %>%
        dplyr::select(Entry, `Raw score`) %>%
        dplyr::mutate(`Raw score` = as.numeric(as.character(`Raw score`)))

    # Merge all data
    combined_data <- bind_rows(sample_excel_df, decode_scores, decode_bic_scores)

    # Ensure the Entry column is a factor with the correct levels
    combined_data$Entry <- factor(combined_data$Entry,
        levels = c(
            "random-clone", "5641737", "5964755", "6052744", "6087270", "6181022",
            "6087309", "6187040", "6087362", "6185506", "6181501", "6184520",
            "6184761", "6187127", "6204327", "SOTA1", "SOTA2", "rand_q50",
            "DECODE", "DECODE-BIC"
        )
    )

    # Convert `Raw score` to numeric
    combined_data$`Raw score` <- as.numeric(as.character(combined_data$`Raw score`))

    # Check for any conversion errors
    if (any(is.na(combined_data$`Raw score`))) {
        warning("NA introduced by coercion")
    }

    # Filter out NA values in the Entry column
    combined_data <- combined_data %>% filter(!is.na(Entry))

    # Plot
    p <- ggplot(combined_data, aes(x = Entry, y = `Raw score`, fill = Entry)) +
        geom_point(position = position_jitterdodge(), aes(color = Entry)) +
        geom_boxplot(alpha = 0.5, outlier.shape = NA) +
        scale_fill_manual(values = method_color_scheme) +
        labs(title = "Scatter Box Plot of Algorithms and Scores", x = "Algorithm", y = "Score") +
        theme_minimal() +
        theme(
            text = element_text(size = 40),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(angle = 90, hjust = 0.5),
            panel.background = element_rect(fill = "white", colour = "white"),
            plot.title = element_text(hjust = 0, size = 40),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "none",
            plot.margin = margin(0.5, 0.5, 0, 2, "cm")
        ) +
        scale_y_continuous(breaks = pretty_breaks(n = 5))

    # Save the plot
    png(paste0(folder_workplace, "DREAM2024_scatter_box_plot_Score1B.png"), res = 150, width = 30, height = 15, units = "in")
    print(p)
    dev.off()
}
