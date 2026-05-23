# ======================================================================
# ======================================================================
# ======================================================================
# ================================================= DECODE ILLUSTRATIONS
# ======================================================================
# ======================================================================
# ======================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
R_data <- "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/data/TCGA_COAD"
R_workplace <- "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes"
R_libPaths <- ""
R_libPaths_extra <- "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/R"
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)
setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)
library(ggplot2)
# ===============================================SET UP WORKPLACE FOLDER
folder_workplace <- "R01_figs/"
if (!dir.exists(folder_workplace)) dir.create(folder_workplace)
# ===========================================GET TCGA SAMPLE INFORMATION
sample_info <- read.table(paste0(R_data, "/sample_table.txt"), header = TRUE)
sample_IDs <- sample_info$Patient
# ==================================GET TCGA CN-BASED TIMELINE INFERENCE
# sample_timelines <- read.table("/Users/kndinh/RESEARCH AND EVERYTHING/Projects/DATASETS/PCAWG/evolution_and_heterogeneity/2018-07-24-wgdMrcaTiming.txt", header = TRUE)
# ======================================================================
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
    text_inference_A <- paste0("Inference A")
    text_inference_B <- paste0("Inference B")
    text_validation <- paste0("Validation")
    color_inference_A <- "#DF536B"
    color_inference_B <- "#DF536B"
    color_validation <- "#2297E6"
    p <- ggplot(readcount_distribution, aes(x = min_total_read, y = min_variant_read, fill = freq)) +
        geom_tile() +
        geom_rect(
            aes(
                xmin = min_total_read_inference_A - 1, xmax = min_total_read_inference_A + 1,
                ymin = min_variant_read_inference_A - 1, ymax = min_variant_read_inference_A + 1
            ),
            fill = NA, color = "white", linewidth = 2
        ) +
        geom_shadowtext(
            aes(
                x = min_total_read_inference_A + 3,
                y = min_variant_read_inference_A + 1,
                label = text_inference_A
            ),
            angle = 45,
            hjust = 0,
            vjust = 0,
            size = 12,
            color = color_inference_A,
            bg.color = "white",
            fontface = "bold"
        ) +
        geom_rect(
            aes(
                xmin = min_total_read_inference_B - 1, xmax = min_total_read_inference_B + 1,
                ymin = min_variant_read_inference_B - 1, ymax = min_variant_read_inference_B + 1
            ),
            fill = NA, color = "white", linewidth = 2
        ) +
        geom_shadowtext(
            aes(
                x = min_total_read_inference_B + 3,
                y = min_variant_read_inference_B + 1,
                label = text_inference_B
            ),
            angle = 45,
            hjust = 0,
            vjust = 0,
            size = 12,
            color = color_inference_B,
            bg.color = "white",
            fontface = "bold"
        ) +
        geom_rect(
            aes(
                xmin = min_total_read_validation - 1, xmax = min_total_read_validation + 1,
                ymin = min_variant_read_validation - 1, ymax = min_variant_read_validation + 1
            ),
            fill = NA, color = "white", linewidth = 2
        ) +
        geom_shadowtext(
            aes(
                x = min_total_read_validation + 3,
                y = min_variant_read_validation + 1,
                label = text_validation
            ),
            angle = 45,
            hjust = 0,
            vjust = 0,
            size = 12,
            color = color_validation,
            bg.color = "white",
            fontface = "bold"
        ) +
        scale_fill_gradientn(
            colors = c("#0072B2", "#56B4E9", "#009E73", "#E69F00", "#D55E00"),
            name = "% mutations retained"
        ) +
        theme_minimal() +
        labs(title = NULL, x = "Minimum total read count", y = "Minimum variant read count") +
        theme(
            text = element_text(size = 40),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5),
            legend.key.width = unit(1.3, "cm"),
        )
    return(p)
}
DECODE_plot_true_SFS_data <- function(DECODE_result,
                                      with_tail = "best",
                                      N_clusters = "best",
                                      mode = "inference_A",
                                      DECODE_linewidth = 5,
                                      text_xlab = "Variant Allele Frequency",
                                      color_xlab = "black",
                                      text_ylab = "Mutation count",
                                      notation = FALSE,
                                      x_min = NULL,
                                      x_max = NULL,
                                      data_marker_colors = NULL,
                                      error_bar = FALSE) {
    suppressPackageStartupMessages(library(ggplot2))
    if (is.null(data_marker_colors)) data_marker_colors <- c("Data" = "black")
    color_scheme <- c(
        data_marker_colors,
        "Neutral_tail" = "#999999",
        "Cluster_1" = "#D55E00",
        "Cluster_2" = "#0072B2",
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
    mutation_table$VAF <- mutation_table$Alt_count / (mutation_table$Ref_count + mutation_table$Alt_count)
    max_total_read <- DECODE_result$max_total_read
    #---Retrieve requested DECODE fit
    if (with_tail == "best") with_tail <- DECODE_result$best_with_tail
    if (N_clusters == "best") N_clusters <- DECODE_result$best_N_clusters
    SFS_fit <- DECODE_result[[paste0("fits_", ifelse(with_tail, "with", "without"), "_tail")]]$all_fits[[paste0(N_clusters, "_clusters")]]
    #---Retrieve requested data
    SFS_data_frequencies <- seq(1, DECODE_result$sfs_bincount) / DECODE_result$sfs_bincount
    VAF <- mutation_table$VAF
    lower_bound <- c(0, SFS_data_frequencies[-length(SFS_data_frequencies)])
    upper_bound <- SFS_data_frequencies
    SFS_data <- numeric(length(SFS_data_frequencies))
    for (i in seq_along(SFS_data_frequencies)) SFS_data[i] <- sum(VAF > lower_bound[i] & VAF <= upper_bound[i])
    min_variant_read <- DECODE_result[[paste0("min_variant_read_", mode)]]
    min_total_read <- DECODE_result[[paste0("min_total_read_", mode)]]
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
    df_fit_wide <- df_fit %>%
        pivot_wider(names_from = fill, values_from = mean, id_cols = frequency) %>%
        mutate(
            row_sum = `Neutral tail` + `Cluster 1` + `Cluster 2`,
            `Neutral tail` = `Neutral tail` / row_sum,
            `Cluster 1` = `Cluster 1` / row_sum,
            `Cluster 2` = `Cluster 2` / row_sum
        ) %>%
        select(-row_sum)
    #---Categorize mutations in the data
    mutation_table$Marker <- sapply(mutation_table$VAF, function(vaf) {
        # Find the closest frequency in df_fit_wide
        closest_idx <- which.min(abs(df_fit_wide$frequency - vaf))
        # Get probabilities for the three categories
        probs <- c(
            df_fit_wide$`Neutral tail`[closest_idx],
            df_fit_wide$`Cluster 1`[closest_idx],
            df_fit_wide$`Cluster 2`[closest_idx]
        )
        # Sample according to probabilities
        sample(c("Neutral_tail", "Cluster_1", "Cluster_2"), size = 1, prob = probs)
    })
    mutation_markers <- mutation_table$Marker
    unique_markers <- unique(mutation_markers)
    df_data <- data.frame(
        ID = rep(1:length(vec_freq), length(unique_markers)),
        frequency = rep(vec_freq, length(unique_markers)),
        count = rep(0, SFS_totalsteps * length(unique_markers)),
        fill = rep(unique_markers, each = SFS_totalsteps)
    )
    for (j in 1:length(mutation_table$VAF)) {
        VAF <- mutation_table$VAF[j]
        pos <- which(vec_freq >= VAF)[1]
        df_data$count[which(df_data$ID == pos & df_data$fill == mutation_markers[j])] <-
            df_data$count[which(df_data$ID == pos & df_data$fill == mutation_markers[j])] + 1
    }
    df_data$fill <- factor(df_data$fill, levels = rev(c("Neutral_tail", "Cluster_1", "Cluster_2")))
    #---Plot SFS fitting result from DECODE
    p <- ggplot() +
        {
            if (error_bar) {
                geom_ribbon(data = df_fit, aes(x = frequency, ymin = pmax(accumulated_mean - sd, 0), ymax = accumulated_mean + sd, fill = fill))
            }
        } +
        geom_bar(
            data = df_data, aes(x = frequency, y = count, fill = fill),
            stat = "identity", width = 1 / SFS_totalsteps
        ) +
        scale_fill_manual(values = color_scheme, name = "") +
        scale_color_manual(values = color_scheme, name = "") +
        scale_x_continuous(limits = c(0, 0.7), expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        guides(fill = guide_legend(nrow = 1, keywidth = 2, keyheight = 1)) +
        xlab("Variant Allele Frequency") +
        ylab("Mutation count") +
        labs(title = NULL) +
        theme(
            text = element_text(size = 50),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "none",
            axis.title.x = element_text(color = color_xlab, margin = margin(t = 20)),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank()
        )
    return(p)
}
DECODE_plot_SFS_data <- function(DECODE_result,
                                 with_tail = "best",
                                 N_clusters = "best",
                                 mode = "inference_A",
                                 DECODE_linewidth = 5,
                                 text_xlab = "Variant Allele Frequency",
                                 color_xlab = "black",
                                 text_ylab = "Mutation count",
                                 notation = FALSE,
                                 x_min = NULL,
                                 x_max = NULL,
                                 data_marker_colors = NULL,
                                 error_bar = FALSE) {
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
        {
            if (error_bar) {
                geom_ribbon(data = df_fit, aes(x = frequency, ymin = pmax(accumulated_mean - sd, 0), ymax = accumulated_mean + sd, fill = fill))
            }
        } +
        geom_bar(
            data = df_data, aes(x = frequency, y = count, fill = fill),
            stat = "identity", width = 1 / SFS_totalsteps
        ) +
        scale_fill_manual(values = color_scheme, name = "") +
        scale_color_manual(values = color_scheme, name = "") +
        guides(fill = guide_legend(nrow = 1, keywidth = 2, keyheight = 1)) +
        xlab(NULL) +
        ylab(NULL) +
        labs(title = NULL) +
        theme(
            text = element_text(size = 30),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "none",
            axis.title.x = element_text(color = color_xlab),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank()
        )
    return(p)
}
DECODE_plot_SFS_fit <- function(DECODE_result,
                                fit = "best",
                                mode = "inference_A",
                                DECODE_linewidth = 5,
                                text_xlab = "Variant Allele Frequency",
                                color_xlab = "black",
                                text_ylab = "Mutation count",
                                notation = FALSE,
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
        # geom_bar(data = df_data, aes(x = frequency, y = count, fill = fill), stat = "identity", width = 1 / SFS_totalsteps) +
        scale_fill_manual(values = color_scheme, name = "") +
        scale_color_manual(values = color_scheme, name = "") +
        guides(fill = guide_legend(nrow = 1, keywidth = 2, keyheight = 1)) +
        xlab(NULL) +
        ylab(NULL) +
        labs(title = NULL) +
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
                axis.text.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.text.y = element_blank()
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
DECODE_plot_BIC <- function(BIC) {
    p <- ggplot(data.frame(x = 1, BIC = BIC), aes(x = x, y = BIC)) +
        geom_boxplot(color = "#704D9E", fill = "#704D9E", alpha = 0.5, linewidth = 5, outlier.shape = NA) +
        xlab("GIC") +
        ylab(NULL) +
        labs(title = NULL) +
        scale_y_continuous(breaks = c(0, 1), labels = c("low", "high"), limits = c(0, 1)) +
        theme(
            text = element_text(size = 200),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
        )
}
# ================================================================DECODE
library(grid)
sample_IDs <- c("TCGA-AA-3977") # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
for (sample in sample_IDs) {
    #---Input the SFS data
    filename <- paste0(R_data, "/", sample, "_1_1.txt")
    mutation_table <- read.table(filename, sep = "\t", header = TRUE)
    mutation_table$Ref_count <- mutation_table$ref_counts
    mutation_table$Alt_count <- mutation_table$alt_counts
    if (any(is.na(mutation_table$Ref_count)) | any(is.na(mutation_table$Alt_count))) mutation_table <- mutation_table[-which(is.na(mutation_table$Ref_count) | is.na(mutation_table$Alt_count)), ]
    if (any(mutation_table$Ref_count == 0) | any(mutation_table$Alt_count == 0)) mutation_table <- mutation_table[-which(mutation_table$Ref_count == 0 | mutation_table$Alt_count == 0), ]
    #----------------------------------------------------DECODE-ORIGINAL
    # #---SFS deconvolution with DECODE
    # DECODE_result <- DECODE(
    #     sample_id = sample,
    #     neutral_tail = TRUE, min_N_clusters = 1, max_N_clusters = 3, # <<<<<<<
    #     mutation_table = mutation_table
    # )
    # save(DECODE_result, file = paste0(folder_workplace, "DECODE_", sample, ".rda"))
    #---Plot DECODE deconvolution
    # load(paste0(folder_workplace, "DECODE_", sample, ".rda"))
    # png(paste0(folder_workplace, "DECODE_", sample, ".png"), res = 150, width = 30, height = 15, units = "in")
    # print(DECODE_plot_SFS(DECODE_result = DECODE_result))
    # dev.off()
    # png(paste0(folder_workplace, "DECODE_model_selection_", sample, ".png"), res = 150, width = 30, height = 15, units = "in")
    # grid.draw(
    #     DECODE_plot_model_selection(
    #         DECODE_result = DECODE_result
    #     )
    # )
    # dev.off()
    #---Plot individual items
    # png(paste0(folder_workplace, "Readcounts.png"), res = 150, width = 10, height = 10.8, units = "in")
    # print(DECODE_plot_readcounts(DECODE_result))
    # dev.off()

    # png(paste0(folder_workplace, "SFS_true_data.png"), res = 150, width = 10, height = 10, units = "in")
    # print(DECODE_plot_true_SFS_data(DECODE_result))
    # dev.off()

    # png(paste0(folder_workplace, "SFS_data_A.png"), res = 150, width = 20, height = 10, units = "in")
    # print(DECODE_plot_SFS_data(DECODE_result, mode = "inference_A"))
    # dev.off()

    # png(paste0(folder_workplace, "SFS_data_B.png"), res = 150, width = 20, height = 10, units = "in")
    # print(DECODE_plot_SFS_data(DECODE_result, mode = "inference_B"))
    # dev.off()

    # png(paste0(folder_workplace, "SFS_data_validation.png"), res = 150, width = 20, height = 10, units = "in")
    # print(DECODE_plot_SFS_data(DECODE_result, mode = "validation"))
    # dev.off()

    # png(paste0(folder_workplace, "SFS_fit_1_A.png"), res = 150, width = 20, height = 10, units = "in")
    # print(DECODE_plot_SFS_fit(DECODE_result, mode = "inference_A", fit = "fits_with_tail:1_clusters"))
    # dev.off()

    # png(paste0(folder_workplace, "SFS_fit_1_B.png"), res = 150, width = 20, height = 10, units = "in")
    # print(DECODE_plot_SFS_fit(DECODE_result, mode = "inference_B", fit = "fits_with_tail:1_clusters"))
    # dev.off()

    # png(paste0(folder_workplace, "SFS_fit_1_validation.png"), res = 150, width = 20, height = 10, units = "in")
    # print(DECODE_plot_SFS_fit(DECODE_result, mode = "validation", fit = "fits_with_tail:1_clusters"))
    # dev.off()

    # png(paste0(folder_workplace, "SFS_fit_2_A.png"), res = 150, width = 20, height = 10, units = "in")
    # print(DECODE_plot_SFS_fit(DECODE_result, mode = "inference_A", fit = "fits_with_tail:2_clusters"))
    # dev.off()

    # png(paste0(folder_workplace, "SFS_fit_2_B.png"), res = 150, width = 20, height = 10, units = "in")
    # print(DECODE_plot_SFS_fit(DECODE_result, mode = "inference_B", fit = "fits_with_tail:2_clusters"))
    # dev.off()

    # png(paste0(folder_workplace, "SFS_fit_2_validation.png"), res = 150, width = 20, height = 10, units = "in")
    # print(DECODE_plot_SFS_fit(DECODE_result, mode = "validation", fit = "fits_with_tail:2_clusters"))
    # dev.off()

    # png(paste0(folder_workplace, "SFS_fit_3_A.png"), res = 150, width = 20, height = 10, units = "in")
    # print(DECODE_plot_SFS_fit(DECODE_result, mode = "inference_A", fit = "fits_with_tail:3_clusters"))
    # dev.off()

    # png(paste0(folder_workplace, "SFS_fit_3_B.png"), res = 150, width = 20, height = 10, units = "in")
    # print(DECODE_plot_SFS_fit(DECODE_result, mode = "inference_B", fit = "fits_with_tail:3_clusters"))
    # dev.off()

    # png(paste0(folder_workplace, "SFS_fit_3_validation.png"), res = 150, width = 20, height = 10, units = "in")
    # print(DECODE_plot_SFS_fit(DECODE_result, mode = "validation", fit = "fits_with_tail:3_clusters"))
    # dev.off()

    BIC_1 <- rbeta(1000, 10, 10)
    png(paste0(folder_workplace, "SFS_fit_1_BIC.png"), res = 150, width = 10, height = 20, units = "in")
    print(DECODE_plot_BIC(BIC_1))
    dev.off()

    BIC_2 <- rbeta(1000, 10, 25)
    png(paste0(folder_workplace, "SFS_fit_2_BIC.png"), res = 150, width = 10, height = 20, units = "in")
    print(DECODE_plot_BIC(BIC_2))
    dev.off()

    BIC_3 <- rbeta(1000, 25, 10)
    png(paste0(folder_workplace, "SFS_fit_3_BIC.png"), res = 150, width = 10, height = 20, units = "in")
    print(DECODE_plot_BIC(BIC_3))
    dev.off()
}
# ======================================================================
# ======================================================================
# ======================================================================
# ================================================ BLCA-MULTIMAP RESULTS
# ======================================================================
# ======================================================================
# ======================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
R_data <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/data/BLCA_multimap"
R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes"
R_libPaths <- ""
R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/R"
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)
setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)
suppressPackageStartupMessages({
    library(tidyverse)
    library(grid)
    library(ggpubr)
})
# ===============================================SET UP WORKPLACE FOLDER
folder_workplace <- "Results_BLCA_multimap/"
folder_plots <- "R01_figs/"
# ======================================GET TCGA-COAD SAMPLE INFORMATION
DECODE_results <- read.csv(paste0(R_data, "/BLCA_multimap_sample_information.csv"), header = TRUE)
DECODE_results <- DECODE_results[which(DECODE_results$patient_id %in% c("M23")), ] # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
sample_IDs <- DECODE_results$sample_id
# ===================================================PLOT DECODE RESULTS
BLCA_multimap_mutational_data <- readRDS(paste0(R_data, "/BLCA_multimap_mutational_data.Rds"))
#---Plot spatial DECODE decompositions
color_panel <- c(
    "NA" = "white",
    "NU" = "#333333",
    "LGIN" = "#3dbce7",
    "HGIN" = "#ffac49",
    "UC" = "#d95f02",
    "Lower_grade_neighbor" = "white",
    "local" = "#004586",
    "widespread" = "#FF420E",
    "UC>=lower-grade" = "#d95f02",
    "UC<lower-grade" = "#333333",
    "HGIN>=lower-grade" = "#ffac49",
    "HGIN<lower-grade" = "#333333",
    "LGIN>=lower-grade" = "#3dbce7",
    "LGIN<lower-grade" = "#333333"
)
# for (patient_id in unique(BLCA_multimap_mutational_data$patient_id)) {
#     # for (patient_id in "M26") {
#     filename <- paste0(folder_workplace, "DECODE_results_", patient_id, ".csv")
#     if (!file.exists(filename)) next
#     DECODE_results <- read.csv(filename, header = TRUE)
#     #---Get mutational data for each patient
#     mutational_data <- BLCA_multimap_mutational_data[BLCA_multimap_mutational_data$patient_id == patient_id, ] %>%
#         mutate(sample_id_mini = sub(".*_", "", sample_id)) %>%
#         mutate(x_label = substr(sample_id_mini, 1, 1), y_label = substr(sample_id_mini, 2, nchar(sample_id_mini)))
#     mutational_data$x_label <- factor(mutational_data$x_label, levels = sort(unique(mutational_data$x_label)))
#     mutational_data$y_label <- factor(mutational_data$y_label, levels = sort(as.numeric(unique(mutational_data$y_label)), decreasing = TRUE))
#     mutational_data$category <- factor(mutational_data$category, levels = c("widespread", "local"))
#     #---Get histology information for each sample
#     sample_data <- expand.grid(x_label = unique(mutational_data$x_label), y_label = unique(mutational_data$y_label))
#     sample_data$patient_id <- patient_id
#     sample_data$sample_id <- paste0(patient_id, "_", sample_data$x_label, sample_data$y_label)
#     sample_data <- merge(sample_data, DECODE_results[, c("sample_id", "group", "Nmutations", "local_neutral_power", "local_cluster_f", "widespread_cluster_f1", "widespread_cluster_f2")], by = "sample_id", all.x = TRUE)
#     sample_data$group[which(is.na(sample_data$Nmutations))] <- "NA"
#     #---Prepare text for each sample
#     sample_data$text_1 <- paste0(
#         # 'plain("', sample_data$group, " (n=", sample_data$Nmutations, ')")'
#         'plain("', sample_data$group, ' " (n==', sample_data$Nmutations, "))"
#     )
#     sample_data$text_2 <- ifelse(
#         !is.na(sample_data$local_neutral_power),
#         paste0("plain(alpha)==", round(sample_data$local_neutral_power, 2)),
#         NA
#     )
#     sample_data$text_3 <- NA
#     # sample_data$text_3 <- ifelse(
#     #     !is.na(sample_data$widespread_cluster_f1),
#     #     ifelse(
#     #         !is.na(sample_data$widespread_cluster_f2),
#     #         paste0(
#     #             "textstyle(atop(",
#     #             "plain(f[1]==", round(pmax(sample_data$widespread_cluster_f1, sample_data$widespread_cluster_f2), 2), "), ",
#     #             "plain(f[2]==", round(pmin(sample_data$widespread_cluster_f1, sample_data$widespread_cluster_f2), 2), ")))"
#     #         ),
#     #         paste0("plain(f)==", round(sample_data$widespread_cluster_f1, 2))
#     #     ),
#     #     NA
#     # )
#     sample_data$text_line2 <- ifelse(
#         !is.na(sample_data$text_3),
#         ifelse(
#             !is.na(sample_data$text_2),
#             paste0("textstyle(atop(", sample_data$text_2, ", ", sample_data$text_3, "))"),
#             sample_data$text_3
#         ),
#         sample_data$text_2
#     )
#     sample_data$text <- ifelse(
#         !is.na(sample_data$text_line2),
#         paste0("textstyle(atop(", sample_data$text_1, ", ", sample_data$text_line2, "))"),
#         paste0("textstyle(atop(", sample_data$text_1, ", plain('')))")
#     )
#     #---Plot SFS for each patient
#     p <- ggplot(mutational_data, aes(x = VAF, fill = category)) +
#         geom_rect(
#             data = sample_data,
#             aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = group),
#             alpha = 0.2, inherit.aes = FALSE
#         ) +
#         geom_histogram(
#             binwidth = 0.01, color = NA, position = "stack"
#         ) +
#         geom_text(
#             data = sample_data,
#             aes(x = Inf, y = Inf, label = text, color = group),
#             hjust = 1, vjust = 1, size = 10, inherit.aes = FALSE, parse = TRUE
#         ) +
#         facet_grid(y_label ~ x_label, drop = TRUE) +
#         labs(x = NULL, y = NULL) +
#         scale_x_continuous(
#             breaks = c(0, 0.25, 0.5, 0.75, 1),
#             labels = c("0", "0.25", "0.5", "0.75", "1")
#         ) +
#         theme_minimal() +
#         theme(
#             plot.background = element_rect(fill = "white"),
#             panel.grid.major = element_blank(),
#             panel.grid.minor = element_blank(),
#             strip.text = element_text(size = 30)
#         ) +
#         scale_fill_manual(values = color_panel) +
#         scale_color_manual(values = color_panel) +
#         guides(fill = FALSE, color = FALSE)
#     pdf(paste0(folder_plots, patient_id, "_spatial_SFS.pdf"), width = 16, height = 12)
#     print(p)
#     dev.off()
#     png(paste0(folder_plots, patient_id, "_spatial_SFS.png"), res = 300, width = 16, height = 12, units = "in")
#     print(p)
#     dev.off()
# }
# #---Plot comparison of neutral powers based on histology
# for (patient_id in unique(BLCA_multimap_mutational_data$patient_id)) {
#     filename <- paste0(folder_workplace, "DECODE_results_", patient_id, ".csv")
#     if (!file.exists(filename)) next
#     DECODE_results <- read.csv(filename, header = TRUE)
#     DECODE_results <- DECODE_results[which(!is.na(DECODE_results$local_neutral_power)), ]
#     DECODE_results$group <- factor(DECODE_results$group, levels = names(color_panel))
#     list_compare <- list()
#     for (group in unique(DECODE_results$group)) {
#         if (group != "NU") {
#             p.val <- wilcox.test(
#                 DECODE_results$local_neutral_power[which(DECODE_results$group == "NU")],
#                 DECODE_results$local_neutral_power[which(DECODE_results$group == group)],
#                 alternative = "two.sided"
#             )$p.val
#             if (p.val < 0.1) list_compare[[length(list_compare) + 1]] <- c("NU", group)
#         }
#     }
#     p <- ggplot(data = DECODE_results, aes(x = as.factor(group), y = local_neutral_power, group = group, color = group, fill = group)) +
#         geom_boxplot(
#             alpha = 0.3, size = 2, outlier.shape = NA
#         ) +
#         geom_jitter(
#             position = position_jitter(width = 0.2),
#             size = 8
#         ) +
#         labs(x = NULL, y = NULL, title = "Neutral tail power") +
#         theme_minimal(base_size = 30) +
#         theme(
#             plot.background = element_rect(fill = "white"),
#             panel.grid.major = element_blank(),
#             panel.grid.minor = element_blank()
#         ) +
#         scale_fill_manual(values = color_panel) +
#         scale_color_manual(values = color_panel) +
#         guides(fill = FALSE, color = FALSE, alpha = FALSE)
#     if (length(list_compare) > 0) {
#         p <- p + geom_signif(
#             comparisons = list_compare,
#             map_signif_level = TRUE, textsize = 8, step_increase = 0.05
#         )
#     }
#     pdf(paste0(folder_plots, patient_id, "_neutral_power.pdf"), width = 6, height = 12)
#     print(p)
#     dev.off()
#     png(paste0(folder_plots, patient_id, "_neutral_power.png"), res = 300, width = 6, height = 12, units = "in")
#     print(p)
#     dev.off()
# }
# #---Plot comparison of neutral powers between UC and lower-grade neighbors
# for (patient_id in unique(BLCA_multimap_mutational_data$patient_id)) {
#     filename <- paste0(folder_workplace, "DECODE_results_", patient_id, ".csv")
#     if (!file.exists(filename)) next
#     DECODE_results <- read.csv(filename, header = TRUE)
#     DECODE_results <- DECODE_results[which(!is.na(DECODE_results$local_neutral_power)), ]
#     DECODE_results$x <- sapply(DECODE_results$x_label, function(x) which(LETTERS == x))
#     DECODE_results$y <- as.numeric(DECODE_results$y_label)
#     DECODE_results_UC <- DECODE_results[DECODE_results$group == "UC", ]
#     df_plot <- data.frame()
#     for (row in 1:nrow(DECODE_results_UC)) {
#         x <- DECODE_results_UC$x[row]
#         y <- DECODE_results_UC$y[row]
#         ids <- which(
#             abs(DECODE_results$x - x) <= 1 &
#                 abs(DECODE_results$y - y) <= 1 &
#                 ((DECODE_results$x != x) | (DECODE_results$y != y)) &
#                 DECODE_results$group %in% c("HGIN", "LGIN", "NU")
#         )
#         if (length(ids) == 0) next
#         df_plot <- rbind(
#             df_plot,
#             data.frame(
#                 id = paste0(patient_id, "/", DECODE_results_UC$sample_id[row], "/", DECODE_results$sample_id[ids]),
#                 UC = DECODE_results_UC$local_neutral_power[row],
#                 Lower_grade_neighbor = DECODE_results$local_neutral_power[ids]
#             )
#         )
#     }
#     if (nrow(df_plot) == 0) next
#     df_plot$comparison <- ifelse(df_plot$UC >= df_plot$Lower_grade_neighbor, "UC>=lower-grade", "UC<lower-grade")
#     p.val <- wilcox.test(df_plot$UC, df_plot$Lower_grade_neighbor, paired = TRUE, alternative = "greater")$p.val
#     df_plot <- df_plot %>%
#         pivot_longer(
#             cols = c("UC", "Lower_grade_neighbor"),
#             names_to = "histology",
#             values_to = "neutral_power"
#         )
#     df_plot$histology <- factor(df_plot$histology, levels = names(color_panel))
#     p <- ggplot(data = df_plot, aes(x = histology, y = neutral_power, group = id, color = comparison, linewidth = comparison, alpha = comparison)) +
#         geom_line() +
#         geom_point(size = 5) +
#         scale_color_manual(values = color_panel) +
#         scale_linewidth_manual(values = c("UC>=lower-grade" = 1, "UC<lower-grade" = 1)) +
#         scale_alpha_manual(values = c("UC>=lower-grade" = 1, "UC<lower-grade" = 0.6)) +
#         labs(x = NULL, y = "Neutral tail power", title = paste0("p-value = ", formatC(p.val))) +
#         theme_minimal(base_size = 30) +
#         theme(
#             plot.background = element_rect(fill = "white"),
#             panel.grid.major = element_blank(),
#             panel.grid.minor = element_blank(),
#             axis.text.x = element_text(angle = 30, hjust = 1)
#         ) +
#         scale_x_discrete(labels = c("UC" = "UC", "Lower_grade_neighbor" = "Lower-grade neighbors"), expand = c(0.1, 0.1)) +
#         guides(fill = FALSE, color = FALSE, alpha = FALSE, linewidth = FALSE)
#     pdf(paste0(folder_plots, patient_id, "_neutral_power UC v lower grade.pdf"), width = 6, height = 12)
#     print(p)
#     dev.off()
#     png(paste0(folder_plots, patient_id, "_neutral_power UC v lower grade.png"), res = 300, width = 6, height = 12, units = "in")
#     print(p)
#     dev.off()
# }
# #---Plot comparison of neutral powers between HGIN and lower-grade neighbors
# for (patient_id in unique(BLCA_multimap_mutational_data$patient_id)) {
#     filename <- paste0(folder_workplace, "DECODE_results_", patient_id, ".csv")
#     if (!file.exists(filename)) next
#     DECODE_results <- read.csv(filename, header = TRUE)
#     DECODE_results <- DECODE_results[which(!is.na(DECODE_results$local_neutral_power)), ]
#     DECODE_results$x <- sapply(DECODE_results$x_label, function(x) which(LETTERS == x))
#     DECODE_results$y <- as.numeric(DECODE_results$y_label)
#     DECODE_results_HGIN <- DECODE_results[DECODE_results$group == "HGIN", ]
#     df_plot <- data.frame()
#     for (row in 1:nrow(DECODE_results_HGIN)) {
#         x <- DECODE_results_HGIN$x[row]
#         y <- DECODE_results_HGIN$y[row]
#         ids <- which(
#             abs(DECODE_results$x - x) <= 1 &
#                 abs(DECODE_results$y - y) <= 1 &
#                 ((DECODE_results$x != x) | (DECODE_results$y != y)) &
#                 DECODE_results$group %in% c("LGIN", "NU")
#         )
#         if (length(ids) == 0) next
#         df_plot <- rbind(
#             df_plot,
#             data.frame(
#                 id = paste0(patient_id, "/", DECODE_results_HGIN$sample_id[row], "/", DECODE_results$sample_id[ids]),
#                 HGIN = DECODE_results_HGIN$local_neutral_power[row],
#                 Lower_grade_neighbor = DECODE_results$local_neutral_power[ids]
#             )
#         )
#     }
#     if (nrow(df_plot) == 0) next
#     df_plot$comparison <- ifelse(df_plot$HGIN >= df_plot$Lower_grade_neighbor, "HGIN>=lower-grade", "HGIN<lower-grade")
#     p.val <- wilcox.test(df_plot$HGIN, df_plot$Lower_grade_neighbor, paired = TRUE, alternative = "greater")$p.val
#     df_plot <- df_plot %>%
#         pivot_longer(
#             cols = c("HGIN", "Lower_grade_neighbor"),
#             names_to = "histology",
#             values_to = "neutral_power"
#         )
#     df_plot$histology <- factor(df_plot$histology, levels = names(color_panel))
#     p <- ggplot(data = df_plot, aes(x = histology, y = neutral_power, group = id, color = comparison, linewidth = comparison, alpha = comparison)) +
#         geom_line() +
#         geom_point(size = 5) +
#         scale_color_manual(values = color_panel) +
#         scale_linewidth_manual(values = c("HGIN>=lower-grade" = 1, "HGIN<lower-grade" = 1)) +
#         scale_alpha_manual(values = c("HGIN>=lower-grade" = 1, "HGIN<lower-grade" = 0.6)) +
#         labs(x = NULL, y = "Neutral tail power", title = paste0("p-value = ", formatC(p.val))) +
#         theme_minimal(base_size = 30) +
#         theme(
#             plot.background = element_rect(fill = "white"),
#             panel.grid.major = element_blank(),
#             panel.grid.minor = element_blank(),
#             axis.text.x = element_text(angle = 30, hjust = 1)
#         ) +
#         scale_x_discrete(labels = c("HGIN" = "HGIN", "Lower_grade_neighbor" = "Lower-grade neighbors"), expand = c(0.1, 0.1)) +
#         guides(fill = FALSE, color = FALSE, alpha = FALSE, linewidth = FALSE)
#     pdf(paste0(folder_plots, patient_id, "_neutral_power HGIN v lower grade.pdf"), width = 6, height = 12)
#     print(p)
#     dev.off()
#     png(paste0(folder_plots, patient_id, "_neutral_power HGIN v lower grade.png"), res = 300, width = 6, height = 12, units = "in")
#     print(p)
#     dev.off()
# }
# #---Plot comparison of neutral powers between LGIN and lower-grade neighbors
# for (patient_id in unique(BLCA_multimap_mutational_data$patient_id)) {
#     filename <- paste0(folder_workplace, "DECODE_results_", patient_id, ".csv")
#     if (!file.exists(filename)) next
#     DECODE_results <- read.csv(filename, header = TRUE)
#     DECODE_results <- DECODE_results[which(!is.na(DECODE_results$local_neutral_power)), ]
#     DECODE_results$x <- sapply(DECODE_results$x_label, function(x) which(LETTERS == x))
#     DECODE_results$y <- as.numeric(DECODE_results$y_label)
#     DECODE_results_LGIN <- DECODE_results[DECODE_results$group == "LGIN", ]
#     df_plot <- data.frame()
#     for (row in 1:nrow(DECODE_results_LGIN)) {
#         x <- DECODE_results_LGIN$x[row]
#         y <- DECODE_results_LGIN$y[row]
#         ids <- which(
#             abs(DECODE_results$x - x) <= 1 &
#                 abs(DECODE_results$y - y) <= 1 &
#                 ((DECODE_results$x != x) | (DECODE_results$y != y)) &
#                 DECODE_results$group %in% c("NU")
#         )
#         if (length(ids) == 0) next
#         df_plot <- rbind(
#             df_plot,
#             data.frame(
#                 id = paste0(patient_id, "/", DECODE_results_LGIN$sample_id[row], "/", DECODE_results$sample_id[ids]),
#                 LGIN = DECODE_results_LGIN$local_neutral_power[row],
#                 Lower_grade_neighbor = DECODE_results$local_neutral_power[ids]
#             )
#         )
#     }
#     if (nrow(df_plot) == 0) next
#     df_plot$comparison <- ifelse(df_plot$LGIN >= df_plot$Lower_grade_neighbor, "LGIN>=lower-grade", "LGIN<lower-grade")
#     p.val <- wilcox.test(df_plot$LGIN, df_plot$Lower_grade_neighbor, paired = TRUE, alternative = "greater")$p.val
#     df_plot <- df_plot %>%
#         pivot_longer(
#             cols = c("LGIN", "Lower_grade_neighbor"),
#             names_to = "histology",
#             values_to = "neutral_power"
#         )
#     df_plot$histology <- factor(df_plot$histology, levels = names(color_panel))
#     p <- ggplot(data = df_plot, aes(x = histology, y = neutral_power, group = id, color = comparison, linewidth = comparison, alpha = comparison)) +
#         geom_line() +
#         geom_point(size = 5) +
#         scale_color_manual(values = color_panel) +
#         scale_linewidth_manual(values = c("LGIN>=lower-grade" = 1, "LGIN<lower-grade" = 1)) +
#         scale_alpha_manual(values = c("LGIN>=lower-grade" = 1, "LGIN<lower-grade" = 0.6)) +
#         labs(x = NULL, y = "Neutral tail power", title = paste0("p-value = ", formatC(p.val))) +
#         theme_minimal(base_size = 30) +
#         theme(
#             plot.background = element_rect(fill = "white"),
#             panel.grid.major = element_blank(),
#             panel.grid.minor = element_blank(),
#             axis.text.x = element_text(angle = 30, hjust = 1)
#         ) +
#         scale_x_discrete(labels = c("LGIN" = "LGIN", "Lower_grade_neighbor" = "Lower-grade neighbors"), expand = c(0.1, 0.1)) +
#         guides(fill = FALSE, color = FALSE, alpha = FALSE, linewidth = FALSE)
#     pdf(paste0(folder_plots, patient_id, "_neutral_power LGIN v lower grade.pdf"), width = 6, height = 12)
#     print(p)
#     dev.off()
#     png(paste0(folder_plots, patient_id, "_neutral_power LGIN v lower grade.png"), res = 300, width = 6, height = 12, units = "in")
#     print(p)
#     dev.off()
# }
# ======================================================================
# ======================================================================
# ======================================================================
# ==================================================== ICGC-TCGA RESULTS
# ======================================================================
# ======================================================================
# ======================================================================
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
# R_data <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/data/PCAWG"
# R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/R"
# # =======================================SET UP FOLDER PATHS & LIBRARIES
# .libPaths(R_libPaths)
# setwd(R_libPaths_extra)
# files_sources <- list.files(pattern = "\\.[rR]$")
# sapply(files_sources, source)
# setwd(R_workplace)
# library(ggplot2)
# library(dplyr)
# # ===============================================SET UP WORKPLACE FOLDER
# folder_workplace <- "Results_ICGC/"
# folder_plots <- "R01_figs/"
# # # =========================================ANALYZE DECONVOLUTION RESULTS
# #---Input plot settings
# source("plot_settings.r")
# # ======================================================================
# plot_analysis <- function(results,
#                           algorithm,
#                           cohort,
#                           algorithm_colors,
#                           cluster_shapes,
#                           cluster_labels,
#                           cluster_colors,
#                           cohort_colors,
#                           folder_workplace) {
#     library(dplyr)
#     library(ggplot2)
#     #---------------Plot distribution of cancer-specific mutation counts
#     if ("Cancer_type" %in% colnames(results)) {
#         plot_df <- results %>%
#             group_by(Cancer_type) %>%
#             mutate(median_Nmut = median(Nmut, na.rm = TRUE)) %>%
#             ungroup() %>%
#             arrange(desc(median_Nmut))
#         plot_df$Cancer_type <- factor(plot_df$Cancer_type, levels = unique(plot_df$Cancer_type))
#         p <- ggplot(plot_df, aes(x = Cancer_type, y = Nmut, fill = Cancer_type)) +
#             geom_boxplot(color = "black", size = 2) +
#             scale_y_log10() +
#             scale_fill_manual(values = ICGC_cohort_colors) +
#             xlab(NULL) +
#             ylab("Mutation count") +
#             labs(NULL) +
#             theme_bw() +
#             theme(
#                 axis.text.x = element_text(angle = 45, hjust = 1, size = 50, margin = margin(t = -30)),
#                 axis.ticks.length = unit(0, "cm"),
#                 legend.position = "none",
#                 text = element_text(size = 50),
#                 panel.background = element_rect(fill = "white", colour = "white"),
#                 panel.grid.major = element_line(colour = "white"),
#                 panel.grid.minor = element_line(colour = "white"),
#                 panel.border = element_blank(),
#                 plot.margin = margin(t = 10, r = 10, b = 10, l = 80)
#             )
#         filename <- paste0(folder_workplace, cohort, "_", algorithm, "_mutation_count.png")
#         png(filename, res = 150, width = 30, height = 15, units = "in", pointsize = 12)
#         print(p)
#         dev.off()
#         cancer_type_levels <- levels(plot_df$Cancer_type)
#     }
#     #------------------Plot distribution of cancer-specific sample sizes
#     if ("Cancer_type" %in% colnames(results)) {
#         plot_df <- results %>%
#             group_by(Cancer_type) %>%
#             summarise(sample_count = n()) %>%
#             arrange(desc(sample_count))
#         plot_df$Cancer_type <- factor(plot_df$Cancer_type, levels = cancer_type_levels)
#         p <- ggplot(plot_df, aes(x = Cancer_type, y = sample_count, fill = Cancer_type)) +
#             geom_bar(stat = "identity", color = "black", size = 2) +
#             scale_fill_manual(values = ICGC_cohort_colors) +
#             xlab(NULL) +
#             ylab("Sample count") +
#             labs(NULL) +
#             theme_bw() +
#             theme(
#                 axis.text.x = element_text(angle = 45, hjust = 1, size = 50, margin = margin(t = -30)),
#                 axis.ticks.length = unit(0, "cm"),
#                 legend.position = "none",
#                 text = element_text(size = 50),
#                 panel.background = element_rect(fill = "white", colour = "white"),
#                 panel.grid.major = element_line(colour = "white"),
#                 panel.grid.minor = element_line(colour = "white"),
#                 panel.border = element_blank(),
#                 plot.margin = margin(t = 10, r = 10, b = 10, l = 80)
#             )
#         filename <- paste0(folder_workplace, cohort, "_", algorithm, "_sample_size.png")
#         png(filename, res = 150, width = 30, height = 15, units = "in", pointsize = 12)
#         print(p)
#         dev.off()
#     }
#     #----------------Plot distribution of cancer-specific tail detection
#     if ("Cancer_type" %in% colnames(results)) {
#         plot_df <- results %>%
#             group_by(Cancer_type) %>%
#             summarise(
#                 total_count = n(),
#                 tail_count = sum(Tail == TRUE)
#             ) %>%
#             mutate(proportion = 100 * tail_count / total_count)
#         plot_df$Cancer_type <- factor(plot_df$Cancer_type, levels = cancer_type_levels)
#         # p <- ggplot(plot_df, aes(x = Cancer_type, y = proportion, fill = Cancer_type)) +
#         p <- ggplot(plot_df, aes(x = Cancer_type, y = proportion)) +
#             geom_bar(stat = "identity", color = "#990F0F", fill = "#FFB2B2", size = 2) +
#             geom_text(aes(label = paste0("n=", total_count)), vjust = -0.5, size = 15) +
#             scale_y_continuous(
#                 limits = c(0, 100),
#                 breaks = seq(0, 100, by = 25),
#                 labels = function(x) paste0(x, "%")
#             ) +
#             scale_fill_manual(values = ICGC_cohort_colors) +
#             xlab(NULL) +
#             ylab("Neutral tail detection") +
#             labs(NULL) +
#             theme_bw() +
#             theme(
#                 axis.text.x = element_text(angle = 45, hjust = 1, size = 50, margin = margin(t = -30)),
#                 axis.ticks.length = unit(0, "cm"),
#                 legend.position = "none",
#                 text = element_text(size = 50),
#                 axis.title.y = element_text(size = 80),
#                 panel.background = element_rect(fill = "white", colour = "white"),
#                 panel.grid.major = element_line(colour = "white"),
#                 panel.grid.minor = element_line(colour = "white"),
#                 panel.border = element_blank(),
#                 plot.margin = margin(t = 10, r = 10, b = 10, l = 80)
#             )
#         filename <- paste0(folder_workplace, cohort, "_", algorithm, "_tail_detection.png")
#         png(filename, res = 150, width = 30, height = 15, units = "in", pointsize = 12)
#         print(p)
#         dev.off()
#     }
#     #---------------------Plot cluster count distribution by cancer type
#     if ("Cancer_type" %in% colnames(results)) {
#         cluster_distribution_df <- results %>%
#             group_by(Cancer_type, Cluster_count) %>%
#             summarise(count = n(), .groups = "drop") %>%
#             ungroup()

#         total_counts_df <- cluster_distribution_df %>%
#             group_by(Cancer_type) %>%
#             summarise(total_count = sum(count), .groups = "drop")
#         cluster_distribution_df <- merge(cluster_distribution_df, total_counts_df, by = "Cancer_type")
#         cluster_distribution_df <- cluster_distribution_df %>%
#             mutate(percentage = (count / total_count) * 100)
#         cluster_distribution_df$Cancer_type <- factor(cluster_distribution_df$Cancer_type, levels = cancer_type_levels)
#         cluster_distribution_df$Cluster_count <- factor(cluster_distribution_df$Cluster_count, levels = rev(sort(unique(cluster_distribution_df$Cluster_count))))

#         p <- ggplot(cluster_distribution_df, aes(x = Cancer_type, y = percentage, fill = as.factor(Cluster_count))) +
#             geom_bar(stat = "identity", position = "stack", width = 1, color = "black", size = 2) +
#             xlab(NULL) +
#             ylab("Cluster count") +
#             labs(fill = NULL) +
#             scale_y_continuous(labels = scales::percent_format(scale = 1)) +
#             scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
#             scale_fill_manual(values = cluster_colors) +
#             guides(fill = guide_legend(keywidth = 2.5, keyheight = 2, reverse = TRUE)) +
#             theme_bw() +
#             theme(
#                 axis.text.x = element_text(angle = 45, hjust = 1, size = 50, margin = margin(t = -30)),
#                 axis.ticks.length = unit(0, "cm"),
#                 legend.position = "top",
#                 legend.justification = "left",
#                 text = element_text(size = 50),
#                 axis.title.y = element_text(size = 80),
#                 panel.background = element_rect(fill = "white", colour = "white"),
#                 panel.grid.major = element_line(colour = "white"),
#                 panel.grid.minor = element_line(colour = "white"),
#                 panel.border = element_blank(),
#                 plot.margin = margin(t = 10, r = -5, b = 10, l = 65)
#             )
#         filename <- paste0(folder_workplace, cohort, "_", algorithm, " _cluster_count_distribution_by_cancer_type.png")
#         png(filename, res = 150, width = 30, height = 15, units = "in", pointsize = 12)
#         print(p)
#         dev.off()
#     }
#     #-------------Plot joint distribution of truncal VAF & sample purity
#     bound_purity_error <- 0.08
#     if ("Purity" %in% colnames(results)) {
#         if (algorithm == "MOBSTER") {
#             plot_df <- results %>%
#                 rowwise() %>%
#                 mutate(
#                     max_vaf = {
#                         vafs <- c_across(starts_with("Cluster_VAF_"))
#                         valid_vafs <- vafs[!is.na(vafs)]
#                         Nmuts <- c_across(starts_with("Cluster_Nmut_"))
#                         valid_Nmuts <- Nmuts[!is.na(Nmuts)]
#                         if (length(valid_vafs) == 1) {
#                             max_vaf <- valid_vafs
#                         } else {
#                             tmp <- which.max(valid_vafs)
#                             if ((valid_Nmuts[tmp] == max(valid_Nmuts)) | (valid_Nmuts[tmp] > 0.1 * sum(valid_Nmuts))) {
#                                 max_vaf <- valid_vafs[tmp]
#                             } else {
#                                 valid_vafs <- valid_vafs[-tmp]
#                                 max_vaf <- max(valid_vafs)
#                             }
#                         }
#                         max_vaf
#                     }
#                 ) %>%
#                 # mutate(
#                 #     best_cluster = which.max(c_across(starts_with("Cluster_Nmut_"))),
#                 #     max_cluster = sub("Cluster_Nmut_", "", best_cluster),
#                 #     max_vaf = cur_data()[[paste0("Cluster_VAF_", max_cluster)]]
#                 # ) %>%
#                 # mutate(max_vaf = max(c_across(starts_with("Cluster_VAF_")), na.rm = TRUE)) %>%
#                 # mutate(max_vaf = max(c_across(starts_with("Cluster_VAF_"))[c_across(starts_with("Cluster_VAF_")) < 0.5], na.rm = TRUE)) %>%
#                 ungroup() %>%
#                 mutate(max_vaf_scaled = max_vaf * 2) %>%
#                 mutate(within_bounds = ifelse(max_vaf_scaled >= Purity - bound_purity_error & max_vaf_scaled <= Purity + bound_purity_error, "Correct", "Wrong"))
#         } else if (algorithm == "DECODE") {
#             plot_df <- results %>%
#                 rowwise() %>%
#                 mutate(
#                     max_vaf = {
#                         vafs <- c_across(starts_with("Cluster_VAF_"))
#                         valid_vafs <- vafs[!is.na(vafs)]
#                         Nmuts <- c_across(starts_with("Cluster_Nmut_inference_A_"))
#                         valid_Nmuts <- Nmuts[!is.na(Nmuts)]
#                         if (length(valid_vafs) == 1) {
#                             max_vaf <- valid_vafs
#                         } else {
#                             tmp <- which.max(valid_vafs)
#                             if ((valid_Nmuts[tmp] == max(valid_Nmuts)) | (valid_Nmuts[tmp] > 0.1 * sum(valid_Nmuts))) {
#                                 max_vaf <- valid_vafs[tmp]
#                             } else {
#                                 valid_vafs <- valid_vafs[-tmp]
#                                 max_vaf <- max(valid_vafs)
#                             }
#                         }
#                         max_vaf
#                     }
#                 ) %>%
#                 # mutate(max_vaf_1 = max(c_across(starts_with("Cluster_VAF_")), na.rm = TRUE)) %>%
#                 # mutate(
#                 #     best_cluster = which.max(c_across(starts_with("Cluster_Nmut_inference_A_"))),
#                 #     max_cluster = sub("Cluster_Nmut_inference_A_", "", best_cluster),
#                 #     max_vaf_2 = cur_data()[[paste0("Cluster_VAF_", max_cluster)]]
#                 # ) %>%
#                 ungroup() %>%
#                 mutate(max_vaf_scaled = max_vaf * 2) %>%
#                 mutate(within_bounds = ifelse(max_vaf_scaled >= Purity - bound_purity_error & max_vaf_scaled <= Purity + bound_purity_error, "Correct", "Wrong"))
#         }
#         print(plot_df$max_vaf_scaled)
#         within_bounds_percentage <- 100 * sum(plot_df$within_bounds == "Correct") / nrow(plot_df)
#         out_of_bounds_percentage <- 100 * sum(plot_df$within_bounds == "Wrong") / nrow(plot_df)
#         p <- ggplot() +
#             geom_abline(intercept = bound_purity_error, slope = 1, color = "grey", linewidth = 2) +
#             geom_abline(intercept = -bound_purity_error, slope = 1, color = "grey", linewidth = 2) +
#             geom_point(data = plot_df, aes(x = Purity, y = max_vaf_scaled, color = within_bounds), size = 10, alpha = 0.5, stroke = 2) +
#             scale_color_manual(
#                 values = c("Correct" = "#D55E00", "Wrong" = "#56B4E9"),
#                 labels = c(
#                     paste0("Within bounds (", round(within_bounds_percentage, 1), "%)"),
#                     paste0("Out of bounds (", round(out_of_bounds_percentage, 1), "%)")
#                 )
#             ) +
#             scale_x_continuous(name = "Sample purity", breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
#             scale_y_continuous(name = "Predicted purity", breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
#             theme_bw() +
#             theme(
#                 legend.position = "top",
#                 legend.justification = c(0, 1),
#                 legend.title = element_blank(),
#                 plot.title = element_blank(),
#                 aspect.ratio = 1,
#                 text = element_text(size = 50),
#                 panel.background = element_rect(fill = "white", colour = "white"),
#                 panel.grid.major = element_line(colour = "white"),
#                 panel.grid.minor = element_line(colour = "white"),
#                 panel.border = element_blank()
#             )
#         filename <- paste0(folder_workplace, cohort, "_", algorithm, "_predicted_vs_true_purity.png")
#         png(filename, res = 150, width = 15, height = 15, units = "in", pointsize = 12)
#         print(p)
#         dev.off()
#     }
# }
# # ===========================================GET ICGC SAMPLE INFORMATION
# sample_info <- read.csv(paste0(R_data, "/ICGC_sample_information.csv"))
# samples_to_delete <- c()
# for (sample_index in 1:nrow(sample_info)) {
#     print("----------")
#     print(sample_index)
#     print(nrow(sample_info))
#     if (sample_info$wgd_status[sample_index] == "wgd" | sample_info$wgd_uncertain[sample_index] == TRUE) {
#         samples_to_delete <- c(samples_to_delete, sample)
#         next
#     }
#     sample <- sample_info$aliquot_id[sample_index]
#     filename <- paste0(R_data, "/", sample, "_1_1.csv")
#     mutation_table <- read.table(filename, sep = "\t", header = TRUE)
#     mutation_table$Ref_count <- mutation_table$t_ref_count
#     mutation_table$Alt_count <- mutation_table$t_alt_count
#     if (any(is.na(mutation_table$Ref_count)) | any(is.na(mutation_table$Alt_count))) mutation_table <- mutation_table[-which(is.na(mutation_table$Ref_count) | is.na(mutation_table$Alt_count)), ]
#     if (any(mutation_table$Ref_count == 0) | any(mutation_table$Alt_count == 0)) mutation_table <- mutation_table[-which(mutation_table$Ref_count == 0 | mutation_table$Alt_count == 0), ]
#     if (nrow(mutation_table) < 100) samples_to_delete <- c(samples_to_delete, sample)
# }
# if (length(samples_to_delete) > 0) sample_info <- sample_info[!sample_info$aliquot_id %in% samples_to_delete, ]
# sample_IDs <- sample_info$aliquot_id
# # ======================================EXTRACT DECONVOLUTION PARAMETERS
# mobster_df <- data.frame()
# decode_df <- data.frame()
# for (sample in sample_IDs) {
#     print("----------")
#     print(which(sample_IDs == sample))
#     print(length(sample_IDs))
#     #---Extract MOBSTER parameters
#     filename <- paste0("/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/DECODE & MOBSTER results/2024-09-19.ICGC DECODE & MOBSTER [Yining]/MOBSTER_", sample, ".rda")
#     if (file.exists(filename)) {
#         load(filename)
#         mobster_df <- MOBSTER_summary_statistics(mobster_df, MOBSTER_result)
#     }
#     # #---Extract DECODE parameters
#     filename <- paste0("/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/DECODE & MOBSTER results/2024-09-19.ICGC DECODE & MOBSTER [Yining]/DECODE_", sample, ".rda")
#     if (file.exists(filename)) {
#         load(filename)
#         decode_df <- DECODE_summary_statistics(decode_df, DECODE_result)
#     }
# }
# mobster_df <- merge(mobster_df, sample_info[, c("aliquot_id", "purity")], by.x = "Sample", by.y = "aliquot_id", all.x = TRUE)
# names(mobster_df)[names(mobster_df) == "purity"] <- "Purity"
# mobster_df <- merge(mobster_df, sample_info[, c("aliquot_id", "histology_abbreviation")], by.x = "Sample", by.y = "aliquot_id", all.x = TRUE)
# names(mobster_df)[names(mobster_df) == "histology_abbreviation"] <- "Cancer_type"
# mobster_df <- merge(mobster_df, sample_info[, c("aliquot_id", "dcc_specimen_type")], by.x = "Sample", by.y = "aliquot_id", all.x = TRUE)
# names(mobster_df)[names(mobster_df) == "dcc_specimen_type"] <- "Specimen_type"
# mobster_df <- merge(mobster_df, sample_info[, c("aliquot_id", "tumour_grade")], by.x = "Sample", by.y = "aliquot_id", all.x = TRUE)
# names(mobster_df)[names(mobster_df) == "tumour_grade"] <- "Tumor_grade"
# mobster_df <- merge(mobster_df, sample_info[, c("aliquot_id", "specimen_donor_treatment_type")], by.x = "Sample", by.y = "aliquot_id", all.x = TRUE)
# names(mobster_df)[names(mobster_df) == "specimen_donor_treatment_type"] <- "Treatment_type"
# mobster_df <- merge(mobster_df, sample_info[, c("aliquot_id", "donor_age_at_diagnosis")], by.x = "Sample", by.y = "aliquot_id", all.x = TRUE)
# names(mobster_df)[names(mobster_df) == "donor_age_at_diagnosis"] <- "Age"
# mobster_df <- merge(mobster_df, sample_info[, c("aliquot_id", "donor_survival_time")], by.x = "Sample", by.y = "aliquot_id", all.x = TRUE)
# names(mobster_df)[names(mobster_df) == "donor_survival_time"] <- "Survival_time"


# decode_df <- merge(decode_df, sample_info[, c("aliquot_id", "purity")], by.x = "Sample", by.y = "aliquot_id", all.x = TRUE)
# names(decode_df)[names(decode_df) == "purity"] <- "Purity"
# decode_df <- merge(decode_df, sample_info[, c("aliquot_id", "histology_abbreviation")], by.x = "Sample", by.y = "aliquot_id", all.x = TRUE)
# names(decode_df)[names(decode_df) == "histology_abbreviation"] <- "Cancer_type"
# decode_df <- merge(decode_df, sample_info[, c("aliquot_id", "dcc_specimen_type")], by.x = "Sample", by.y = "aliquot_id", all.x = TRUE)
# names(decode_df)[names(decode_df) == "dcc_specimen_type"] <- "Specimen_type"
# decode_df <- merge(decode_df, sample_info[, c("aliquot_id", "tumour_grade")], by.x = "Sample", by.y = "aliquot_id", all.x = TRUE)
# names(decode_df)[names(decode_df) == "tumour_grade"] <- "Tumor_grade"
# decode_df <- merge(decode_df, sample_info[, c("aliquot_id", "specimen_donor_treatment_type")], by.x = "Sample", by.y = "aliquot_id", all.x = TRUE)
# names(decode_df)[names(decode_df) == "specimen_donor_treatment_type"] <- "Treatment_type"
# decode_df <- merge(decode_df, sample_info[, c("aliquot_id", "donor_age_at_diagnosis")], by.x = "Sample", by.y = "aliquot_id", all.x = TRUE)
# names(decode_df)[names(decode_df) == "donor_age_at_diagnosis"] <- "Age"
# decode_df <- merge(decode_df, sample_info[, c("aliquot_id", "donor_survival_time")], by.x = "Sample", by.y = "aliquot_id", all.x = TRUE)
# names(decode_df)[names(decode_df) == "donor_survival_time"] <- "Survival_time"


# decode_df <- merge(decode_df, sample_timelines[, c("uuid", "MRCA.time.branching")], by.x = "Sample", by.y = "uuid", all.x = TRUE)
# names(decode_df)[names(decode_df) == "MRCA.time.branching"] <- "MRCA_age_predict_CN_branching"
# decode_df <- merge(decode_df, sample_timelines[, c("uuid", "MRCA.time.branching.10.")], by.x = "Sample", by.y = "uuid", all.x = TRUE)
# names(decode_df)[names(decode_df) == "MRCA.time.branching.10."] <- "MRCA_age_predict_CN_branching_10"
# decode_df <- merge(decode_df, sample_timelines[, c("uuid", "MRCA.time.branching.90.")], by.x = "Sample", by.y = "uuid", all.x = TRUE)
# names(decode_df)[names(decode_df) == "MRCA.time.branching.90."] <- "MRCA_age_predict_CN_branching_90"
# decode_df <- merge(decode_df, sample_timelines[, c("uuid", "MRCA.time.linear")], by.x = "Sample", by.y = "uuid", all.x = TRUE)
# names(decode_df)[names(decode_df) == "MRCA.time.linear"] <- "MRCA_age_predict_CN_linear"
# decode_df <- merge(decode_df, sample_timelines[, c("uuid", "MRCA.time.linear.10.")], by.x = "Sample", by.y = "uuid", all.x = TRUE)
# names(decode_df)[names(decode_df) == "MRCA.time.linear.10."] <- "MRCA_age_predict_CN_linear_10"
# decode_df <- merge(decode_df, sample_timelines[, c("uuid", "MRCA.time.linear.90.")], by.x = "Sample", by.y = "uuid", all.x = TRUE)
# names(decode_df)[names(decode_df) == "MRCA.time.linear.90."] <- "MRCA_age_predict_CN_linear_90"

# #---Filter out cancer types with fewer than 10 samples
# cancer_types <- unique(mobster_df$Cancer_type)
# cancer_types_nSamples <- c()
# for (cancer_type in cancer_types) cancer_types_nSamples <- c(cancer_types_nSamples, nrow(mobster_df[which(mobster_df$Cancer_type == cancer_type), ]))
# cancer_types <- cancer_types[which(cancer_types_nSamples > 10)]
# mobster_df <- mobster_df[which(mobster_df$Cancer_type %in% cancer_types), ]
# cancer_types <- unique(decode_df$Cancer_type)
# cancer_types_nSamples <- c()
# for (cancer_type in cancer_types) cancer_types_nSamples <- c(cancer_types_nSamples, nrow(decode_df[which(decode_df$Cancer_type == cancer_type), ]))
# cancer_types <- cancer_types[which(cancer_types_nSamples > 10)]
# decode_df <- decode_df[which(decode_df$Cancer_type %in% cancer_types), ]
# write.csv(mobster_df, paste0(folder_plots, "MOBSTER_ICGC.csv"), row.names = FALSE)
# write.csv(decode_df, paste0(folder_plots, "DECODE_ICGC.csv"), row.names = FALSE)
# #---Make plots for analysis of MOBSTER results
# mobster_df <- read.csv(paste0(folder_plots, "MOBSTER_ICGC.csv"))
# plot_analysis(
#     results = mobster_df,
#     algorithm = "MOBSTER",
#     cohort = "ICGC",
#     algorithm_colors = algorithm_colors,
#     cluster_shapes = cluster_shapes,
#     cluster_labels = cluster_labels,
#     cluster_colors = cluster_colors,
#     cohort_colors = ICGC_cohort_colors,
#     folder_workplace = folder_plots
# )
# #---Make plots for analysis of DECODE results
# decode_df <- read.csv(paste0(folder_plots, "DECODE_ICGC.csv"))
# plot_analysis(
#     results = decode_df,
#     algorithm = "DECODE",
#     cohort = "ICGC",
#     algorithm_colors = algorithm_colors,
#     cluster_shapes = cluster_shapes,
#     cluster_labels = cluster_labels,
#     cluster_colors = cluster_colors,
#     cohort_colors = ICGC_cohort_colors,
#     folder_workplace = folder_plots
# )
# #---Plot current sample sizes in TCGA
# library(readxl)
# tcga_data <- read_excel(paste0(folder_plots, "TCGA_current.xlsx"))
# colnames(tcga_data) <- c("Cancer_type", "Sample_count")
# tcga_data <- tcga_data %>% arrange(desc(Sample_count))
# tcga_data$Cancer_type <- gsub("TCGA-", "", tcga_data$Cancer_type)
# tcga_data$Cancer_type <- factor(tcga_data$Cancer_type, levels = tcga_data$Cancer_type)
# p <- ggplot(tcga_data, aes(x = Cancer_type, y = Sample_count, fill = Cancer_type)) +
#     geom_bar(stat = "identity", color = "#0F6B99", fill = "#B2E5FF", size = 2) +
#     scale_fill_manual(values = ICGC_cohort_colors) +
#     xlab(NULL) +
#     ylab("Sample count") +
#     labs(NULL) +
#     theme_bw() +
#     theme(
#         axis.text.x = element_text(angle = 45, hjust = 1, size = 50, margin = margin(t = -30)),
#         axis.ticks.length = unit(0, "cm"),
#         legend.position = "none",
#         text = element_text(size = 50),
#         axis.title.y = element_text(size = 80),
#         panel.background = element_rect(fill = "white", colour = "white"),
#         panel.grid.major = element_line(colour = "white"),
#         panel.grid.minor = element_line(colour = "white"),
#         panel.border = element_blank(),
#         plot.margin = margin(t = 10, r = 10, b = 10, l = 80)
#     )
# filename <- paste0(folder_plots, "TCGA_current_sample_size.png")
# png(filename, res = 150, width = 30, height = 13, units = "in", pointsize = 12)
# print(p)
# dev.off()
# #---Plot coverage distributions in ICGC samples
# sample_info <- read.csv(paste0(R_data, "/ICGC_sample_information.csv"))
# df_avg <- data.frame(coverage = c())
# p_low <- ggplot()
# p_high <- ggplot()
# for (sample_index in 1:nrow(sample_info)) {
#     print("----------")
#     print(sample_index)
#     print(nrow(sample_info))
#     filename <- paste0(R_data, "/", sample_info$aliquot_id[sample_index], "_all.csv")
#     mutation_table <- read.table(filename, sep = ",", header = TRUE)
#     total_reads <- mutation_table$t_ref_count + mutation_table$t_alt_count
#     total_reads <- total_reads[!is.na(total_reads)]
#     df_avg <- rbind(df_avg, data.frame(coverage = mean(total_reads)))
#     if ((mean(total_reads) >= 15) & (mean(total_reads) <= 35)) {
#         df_sample <- data.frame(total_reads = total_reads, avg = mean(total_reads))
#         p_low <- p_low +
#             stat_bin(
#                 data = df_sample,
#                 aes(x = total_reads, y = after_stat(count / sum(count)), color = avg),
#                 binwidth = 5,
#                 geom = "line"
#             )
#     } else if ((mean(total_reads) >= 70) & (mean(total_reads) <= 100)) {
#         df_sample <- data.frame(total_reads = total_reads, avg = mean(total_reads))
#         p_high <- p_high +
#             stat_bin(
#                 data = df_sample,
#                 aes(x = total_reads, y = after_stat(count / sum(count)), color = avg),
#                 binwidth = 5,
#                 geom = "line"
#             )
#     }
# }
# p_low <- p_low +
#     xlim(0, 125) +
#     xlab("Coverage distribution") +
#     ylab(NULL) +
#     scale_color_gradientn(
#         name = "Mean coverage", limits = c(15, 35),
#         colors = c("#33608C", "#F5B35E", "#B81840"),
#         values = c(0, 0.5, 1),
#         guide = guide_colorbar(barwidth = unit(1, "cm"), barheight = unit(10, "cm"))
#     ) +
#     theme_bw() +
#     theme(
#         axis.ticks.length = unit(0, "cm"),
#         legend.position = c(0.9, 0.9),
#         legend.justification = c("right", "top"),
#         legend.title = element_text(margin = margin(b = 20)),
#         legend.background = element_blank(),
#         text = element_text(size = 70),
#         panel.background = element_rect(fill = "white", colour = "white"),
#         panel.grid.major = element_line(colour = "white"),
#         panel.grid.minor = element_line(colour = "white"),
#         panel.border = element_blank(),
#         plot.margin = margin(t = 0, r = -30, b = 0, l = 0)
#     )
# filename <- paste0(folder_plots, "ICGC_coverage_low.png")
# png(filename, res = 150, width = 17, height = 17, units = "in", pointsize = 12)
# print(p_low)
# dev.off()
# p_high <- p_high +
#     xlim(0, 125) +
#     xlab("Coverage distribution") +
#     ylab(NULL) +
#     scale_color_gradientn(
#         name = "Mean coverage", limits = c(70, 100),
#         colors = c("#33608C", "#F5B35E", "#B81840"),
#         values = c(0, 0.5, 1),
#         guide = guide_colorbar(barwidth = unit(1, "cm"), barheight = unit(10, "cm"))
#     ) +
#     theme_bw() +
#     theme(
#         axis.ticks.length = unit(0, "cm"),
#         legend.position = c(0.5, 0.9),
#         legend.justification = c("right", "top"),
#         legend.title = element_text(margin = margin(b = 20)),
#         legend.background = element_blank(),
#         text = element_text(size = 70),
#         panel.background = element_rect(fill = "white", colour = "white"),
#         panel.grid.major = element_line(colour = "white"),
#         panel.grid.minor = element_line(colour = "white"),
#         panel.border = element_blank(),
#         plot.margin = margin(t = 0, r = -30, b = 0, l = 0)
#     )
# filename <- paste0(folder_plots, "ICGC_coverage_high.png")
# png(filename, res = 150, width = 17, height = 17, units = "in", pointsize = 12)
# print(p_high)
# dev.off()
# p_mean_coverage <- ggplot() +
#     stat_bin(data = df_avg, aes(x = coverage), breaks = 0:125, color = "#C7C7C7", fill = "#C7C7C7") +
#     # stat_bin(data = data.frame(coverage = df_avg$coverage[which(df_avg$coverage >= 15 & df_avg$coverage <= 35)]), aes(x = coverage), breaks = 0:125, color = "#B81840", fill = "#B81840") +
#     # stat_bin(data = data.frame(coverage = df_avg$coverage[which(df_avg$coverage >= 70 & df_avg$coverage <= 100)]), aes(x = coverage), breaks = 0:125, color = "#33608C", fill = "#33608C") +
#     xlim(0, 125) +
#     xlab("Mean coverage") +
#     ylab(NULL) +
#     theme_bw() +
#     theme(
#         axis.ticks.length = unit(0, "cm"),
#         text = element_text(size = 70),
#         panel.background = element_rect(fill = "white", colour = "white"),
#         panel.grid.major = element_line(colour = "white"),
#         panel.grid.minor = element_line(colour = "white"),
#         panel.border = element_blank(),
#         plot.margin = margin(t = 0, r = -50, b = 0, l = 0)
#     )
# filename <- paste0(folder_plots, "ICGC_coverage.png")
# png(filename, res = 150, width = 17, height = 17, units = "in", pointsize = 12)
# print(p_mean_coverage)
# dev.off()
# #---Predict MRCA ages for samples with tail + 1 cluster
# decode_df <- read.csv(paste0(folder_plots, "DECODE_ICGC.csv"))
# mini_decode_df <- decode_df %>%
#     filter(
#         Cancer_type == "Liver-HCC" &
#             # Cancer_type == "Lymph-BNHL" &
#             Cluster_count == 1 &
#             Tail == TRUE
#     )

# nSampledCells <- 100000

# mini_decode_df$Mutation_rate_predict <-
#     (mini_decode_df$Tail_Nmut_predict / 1000 +
#         mini_decode_df$Cluster_Nmut_predict_1 / log(nSampledCells)) /
#         (mini_decode_df$Age / log(nSampledCells))

# # mini_decode_df$Mutation_rate_predict <-
# #     (mini_decode_df$Tail_Nmut_predict * (2 * mini_decode_df$Cluster_VAF_1) / 1000 +
# #         mini_decode_df$Cluster_Nmut_predict_1 / log(nSampledCells)) /
# #         (mini_decode_df$Age / log(nSampledCells))

# # mini_decode_df$Mutation_rate_predict <-
# #     (mini_decode_df$Tail_Nmut_predict / 1000 +
# #         2 * mini_decode_df$Cluster_VAF_1 * mini_decode_df$Cluster_Nmut_predict_1 / log(nSampledCells * 2 * mini_decode_df$Cluster_VAF_1)) /
# #         (2 * mini_decode_df$Cluster_VAF_1 * mini_decode_df$Age / log(nSampledCells * 2 * mini_decode_df$Cluster_VAF_1) +
# #             (1 - 2 * mini_decode_df$Cluster_VAF_1) * mini_decode_df$Age / log(nSampledCells * (1 - 2 * mini_decode_df$Cluster_VAF_1)))

# mini_decode_df$MRCA_age_predict <-
#     mini_decode_df$Cluster_Nmut_predict_1 /
#         mini_decode_df$Mutation_rate_predict
# mini_decode_df$Growth_rate_predict <-
#     log(nSampledCells * 2 * mini_decode_df$Cluster_VAF_1) /
#         (mini_decode_df$Age - mini_decode_df$MRCA_age_predict)

# # p <- hist(mini_decode_df$Mutation_rate_predict)
# # print(p)
# # print(mean(mini_decode_df$Mutation_rate_predict))

# # p <- ggplot(
# #     mini_decode_df[which(mini_decode_df$MRCA_age_predict_CN_linear > 0), ],
# #     aes(x = (Age - MRCA_age_predict_CN_linear) / Age, y = MRCA_age_predict / Age)
# # ) +
# #     geom_point(size = 3, color = "blue", alpha = 0.7) +
# #     theme_minimal(base_size = 20)
# # print(p)
# # print(cor.test((mini_decode_df$Age - mini_decode_df$MRCA_age_predict_CN_linear) / mini_decode_df$Age, mini_decode_df$MRCA_age_predict / mini_decode_df$Age))

# # p <- ggplot(
# #     mini_decode_df[which(mini_decode_df$MRCA_age_predict_CN_linear > 0), ],
# #     aes(x = Growth_rate_predict, y = Survival_time)
# # ) +
# #     geom_point(size = 3, color = "blue", alpha = 0.7) +
# #     theme_minimal(base_size = 20)
# # print(p)
# # print(cor.test(mini_decode_df$Survival_time, mini_decode_df$Growth_rate_predict))

# # p <- ggplot(
# #     mini_decode_df[which(mini_decode_df$MRCA_age_predict_CN_linear > 0), ],
# #     aes(x = MRCA_age_predict_CN_linear, y = Survival_time)
# # ) +
# #     geom_point(size = 3, color = "blue", alpha = 0.7) +
# #     theme_minimal(base_size = 20)
# # print(p)
# # print(cor.test(mini_decode_df$Survival_time, mini_decode_df$MRCA_age_predict_CN_linear))

# p <- ggplot(
#     mini_decode_df,
#     # mini_decode_df[which(mini_decode_df$MRCA_age_predict_CN_linear > 0), ],
#     aes(x = Mutation_rate_predict, y = Survival_time)
# ) +
#     geom_point(size = 3, color = "blue", alpha = 0.7) +
#     theme_minimal(base_size = 20)
# print(p)
# print(cor.test(mini_decode_df$Survival_time, mini_decode_df$Mutation_rate_predict))

# # p <- ggplot(
# #     mini_decode_df,
# #     # mini_decode_df[which(mini_decode_df$MRCA_age_predict_CN_linear > 0), ],
# #     aes(x = Tail_power, y = Survival_time)
# # ) +
# #     geom_point(size = 3, color = "blue", alpha = 0.7) +
# #     theme_minimal(base_size = 20)
# # print(p)
# # print(cor.test(mini_decode_df$Survival_time, mini_decode_df$Tail_power))
# # ======================================================================
# # ======================================================================
# # ======================================================================
# # ============================================== CHROMOSOMAL INSTABILITY
# # ======================================================================
# # ======================================================================
# # ======================================================================
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
# R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/R"
# # =======================================SET UP FOLDER PATHS & LIBRARIES
# .libPaths(R_libPaths)
# setwd(R_libPaths_extra)
# files_sources <- list.files(pattern = "\\.[rR]$")
# sapply(files_sources, source)
# setwd(R_workplace)
# library(Pareto)
# # ===============================================SET UP WORKPLACE FOLDER
# folder_workplace <- "R01_figs/"
# # ======================================================================
# load(file = "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/data/PCAWG/ICGC_purity_coverage.rda")
# coverage <- ICGC_purity_coverage[["sample_3"]]$coverage
# nTail <- 8 * 10^4
# nCluster <- 3 * 10^5
# purity <- 0.75
# tmp <- data.frame(
#     mutation_markers = c(rep("Tail", nTail), rep("Cluster", nCluster)),
#     CCF = c(rPareto(nTail, 0.05, 0.5, truncation = 1), rep(1, nCluster)),
#     c = sample(c(2, 1, 3), nTail + nCluster, replace = TRUE, prob = c(0.3, 0.4, 0.3)),
#     m = 1
# )
# tmp$m[which(tmp$c == 3)] <- sample(c(1, 2), sum(tmp$c == 3), replace = TRUE, prob = c(0.6, 0.4))
# tmp$mutation_markers <- ifelse(
#     tmp$mutation_markers == "Cluster",
#     paste0("m=", tmp$m, ",c=", tmp$c),
#     "Tail"
# )
# tmp$VAF <- tmp$CCF * tmp$m * purity / (tmp$c * purity + 2 * (1 - purity))
# # tmp$r <- round(pmax(rnorm(nrow(tmp), 30, 10), 0))
# tmp$r <- sample(coverage$Read_count, nrow(tmp), replace = TRUE, prob = coverage$Frequency)


# tmp$alt <- rbinom(nrow(tmp), tmp$r, tmp$VAF)
# tmp$VAF_observed <- tmp$alt / tmp$r
# if (any(is.na(tmp$VAF_observed))) {
#     tmp <- tmp[!is.na(tmp$VAF_observed), ]
# }
# if (any(tmp$VAF_observed == 0)) {
#     tmp <- tmp[tmp$VAF_observed != 0, ]
# }
# vec_freq <- seq(0, 1, by = 0.02)
# SFS_totalsteps <- length(vec_freq)
# color_scheme <- c(
#     "Tail" = "#999999",
#     "m=1,c=2" = "#D55E00",
#     "m=1,c=1" = "#6B990F",
#     "m=1,c=3" = "#6551CC",
#     "m=2,c=3" = "#BFB2FF"
# )
# unique_markers <- unique(tmp$mutation_markers)
# df_data <- data.frame(
#     ID = rep(1:length(vec_freq), length(unique_markers)),
#     frequency = rep(vec_freq, length(unique_markers)),
#     count = rep(0, SFS_totalsteps * length(unique_markers)),
#     fill = rep(unique_markers, each = SFS_totalsteps)
# )
# for (j in 1:nrow(tmp)) {
#     VAF <- tmp$VAF_observed[j]
#     pos <- which(vec_freq >= VAF)[1]
#     df_data$count[which(df_data$ID == pos & df_data$fill == tmp$mutation_markers[j])] <- df_data$count[which(df_data$ID == pos & df_data$fill == tmp$mutation_markers[j])] + 1
# }
# df_data$fill <- factor(df_data$fill, levels = rev(names(color_scheme)))
# p <- ggplot() +
#     geom_bar(data = df_data, aes(x = frequency, y = count, fill = fill), stat = "identity", width = 1 / SFS_totalsteps) +
#     scale_fill_manual(values = color_scheme, name = "") +
#     scale_color_manual(values = color_scheme, name = "") +
#     guides(fill = guide_legend(nrow = 1, keywidth = 2, keyheight = 1)) +
#     xlab(NULL) +
#     ylab(NULL) +
#     labs(title = NULL) +
#     theme(
#         text = element_text(size = 80),
#         panel.background = element_rect(fill = "white", colour = "white"),
#         panel.grid.major = element_line(colour = "white"),
#         panel.grid.minor = element_line(colour = "white"),
#         legend.position = "none",
#         axis.title.x = element_text(color = "black"),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_blank()
#     )
# png(paste0(folder_workplace, "Ex_CIN.png"), res = 150, width = 15, height = 15, units = "in")
# print(p)
# dev.off()
# # ======================================================================
# # ======================================================================
# # ======================================================================
# # =============================================================== POG570
# # ======================================================================
# # ======================================================================
# # ======================================================================
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
# R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/R"
# # =======================================SET UP FOLDER PATHS & LIBRARIES
# .libPaths(R_libPaths)
# setwd(R_libPaths_extra)
# files_sources <- list.files(pattern = "\\.[rR]$")
# sapply(files_sources, source)
# setwd(R_workplace)
# library(data.table)
# # ===============================================SET UP WORKPLACE FOLDER
# folder_workplace <- "/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/DATASETS/POG570/Ploidetect/hg38/"
# folder_plots <- "R01_figs/"
# color_scheme <- c(
#     "m=1,c=1" = "#6B990F",
#     "m=1,c=3" = "#6551CC",
#     "m=2,c=3" = "#BFB2FF"
# )
# fill_colors <- c(
#     "1" = "#6B990F",
#     "2" = "#D55E00",
#     "3" = "#6551CC",
#     "4" = "#1F77B4",
#     "5" = "#EFC000",
#     "6" = "#8C564B",
#     "7" = "#17BECF",
#     "8" = "#B22D3C",
#     "1|1" = "#6B990F",
#     "2|1" = "#D55E00",
#     "2|2" = "#a04e0f",
#     "3|3" = "#6551CC",
#     "4|2" = "#569dd1",
#     "4|3" = "#1F77B4",
#     "5|5" = "#EFC000",
#     "6|3" = "#bb938b",
#     "6|4" = "#9f746c",
#     "6|5" = "#8C564B",
#     "7|7" = "#17BECF",
#     "8|4" = "#FFB2BC",
#     "8|5" = "#E67E8A",
#     "8|6" = "#CC5260",
#     "8|7" = "#B22D3C",
#     "8|8" = "#990F20",
#     "other" = "#7F7F7F"
# )
# # =====================================PLOT DISTRIBUTION OF COPY NUMBERS
# sample_ids <- list.dirs(folder_workplace, recursive = FALSE, full.names = FALSE)
# sample_ploidies <- c()
# sample_lohs <- c()
# POG_cn <- list()
# for (sample_id in sample_ids) {
#     filepath <- paste0(folder_workplace, sample_id, "/", sample_id, "_ploidetect_cna_condensed.tsv")
#     cn <- read.table(filepath, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#     if (any(is.na(cn$state))) {
#         cn <- cn[!is.na(cn$state), ]
#     }
#     cn$width <- cn$end - cn$pos + 1
#     cn$A_rounded <- 0
#     cn$B_rounded <- 0
#     cn$A_rounded[which(cn$state > 0)] <- round(cn$A[which(cn$state > 0)] / (cn$A[which(cn$state > 0)] + cn$B[which(cn$state > 0)]) * cn$state[which(cn$state > 0)])
#     cn$B_rounded[which(cn$state > 0)] <- round(cn$B[which(cn$state > 0)] / (cn$A[which(cn$state > 0)] + cn$B[which(cn$state > 0)]) * cn$state[which(cn$state > 0)])
#     cn$minor_CN <- pmin(cn$A_rounded, cn$B_rounded)
#     cn$major_CN <- pmax(cn$A_rounded, cn$B_rounded)
#     cn$sample <- sample_id
#     POG_cn[[sample_id]] <- cn
#     sample_ploidies <- c(sample_ploidies, sum(cn$state * cn$width) / sum(cn$width))
#     sample_lohs <- c(sample_lohs, sum(cn$width[which(cn$minor_CN == 0)]) / sum(cn$width))
# }
# POG_cn <- rbindlist(POG_cn)
# POG_cn$sample <- sub(".*_", "", POG_cn$sample)
# POG_samples <- data.frame(
#     Sample_ID = sub(".*_", "", sample_ids),
#     Ploidy = sample_ploidies,
#     LOH = sample_lohs
# )
# POG_samples$WGD <- ifelse(POG_samples$Ploidy < (2.9 - 2 * POG_samples$LOH), "non-WGD", "WGD")
# p <- ggplot(data = POG_samples, aes(x = LOH, y = Ploidy, color = WGD)) +
#     geom_point(size = 10, alpha = 0.5, stroke = 2) +
#     geom_abline(intercept = 2.9, slope = -2, size = 2, color = "#666666") +
#     xlim(0, 1) +
#     labs(x = "Fraction of genome with LOH", y = "Ploidy") +
#     scale_color_manual(name = NULL, values = c("WGD" = "#631879", "non-WGD" = "#008280")) +
#     theme(
#         text = element_text(size = 65),
#         panel.background = element_rect(fill = "white"),
#         axis.line = element_line(),
#         legend.text = element_text(size = 70),
#         legend.position = c(1, 1),
#         legend.justification = c(1, 1)
#     )
# png(paste0(folder_plots, "POG_wgd.png"), res = 150, width = 15, height = 15, units = "in")
# print(p)
# dev.off()
# POG_cn_nonWGD <- POG_cn[
#     which(
#         POG_cn$sample %in% POG_samples$Sample_ID[which(POG_samples$WGD == "non-WGD")] &
#             POG_cn$chr %in% 1:22
#     ),
# ]
# POG_cn_nonWGD_state <- data.frame(
#     state = c(0:max(POG_cn_nonWGD$state))
# )
# POG_cn_nonWGD_state$frequency <- 0
# for (i in 1:nrow(POG_cn_nonWGD_state)) {
#     POG_cn_nonWGD_state$frequency[i] <- sum(POG_cn_nonWGD$width[which(POG_cn_nonWGD$state == POG_cn_nonWGD_state$state[i])]) / sum(POG_cn_nonWGD$width)
# }
# POG_cn_WGD <- POG_cn[
#     which(
#         POG_cn$sample %in% POG_samples$Sample_ID[which(POG_samples$WGD == "WGD")] &
#             POG_cn$chr %in% 1:22
#     ),
# ]
# POG_cn_WGD_state <- data.frame(
#     state = c(0:max(POG_cn_WGD$state))
# )
# POG_cn_WGD_state$frequency <- 0
# for (i in 1:nrow(POG_cn_WGD_state)) {
#     POG_cn_WGD_state$frequency[i] <- sum(POG_cn_WGD$width[which(POG_cn_WGD$state == POG_cn_WGD_state$state[i])]) / sum(POG_cn_WGD$width)
# }
# cn_keep_nonwgd <- c(1, 2, 3, 4)
# POG_cn_nonWGD_state <- rbind(
#     POG_cn_nonWGD_state[POG_cn_nonWGD_state$state %in% cn_keep_nonwgd, ],
#     data.frame(state = "other", frequency = sum(POG_cn_nonWGD_state$frequency[!(POG_cn_nonWGD_state$state %in% cn_keep_nonwgd)]))
# )
# POG_cn_nonWGD_state$label <- ifelse(as.character(POG_cn_nonWGD_state$state) == "other", "other", paste0("c=", POG_cn_nonWGD_state$state))
# POG_cn_nonWGD_state$state <- factor(as.character(POG_cn_nonWGD_state$state), levels = c("other", rev(cn_keep_nonwgd)))
# p <- ggplot(POG_cn_nonWGD_state, aes(x = factor(1), y = frequency, fill = factor(state))) +
#     geom_bar(stat = "identity", width = 0.5) +
#     geom_text(aes(label = label),
#         position = position_stack(vjust = 0.5),
#         color = "white", size = 20
#     ) +
#     scale_y_continuous(labels = scales::percent_format(scale = 100)) +
#     scale_fill_manual(values = fill_colors, name = "") +
#     labs(x = NULL, y = "Total CN frequencies (non-WGD)", title = NULL) +
#     theme_minimal() +
#     theme(
#         text = element_text(size = 65),
#         panel.background = element_rect(fill = "white"),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         legend.title = element_text(size = 40),
#         legend.position = "none",
#         panel.border = element_blank(),
#         plot.margin = margin(10, 10, 10, 10)
#     )
# png(paste0(folder_plots, "POG_state_nonwgd.png"), res = 150, width = 11, height = 15, units = "in")
# print(p)
# dev.off()
# cn_keep_wgd <- c(2, 3, 4, 5, 6, 7, 8)
# POG_cn_WGD_state <- rbind(
#     POG_cn_WGD_state[POG_cn_WGD_state$state %in% cn_keep_wgd, ],
#     data.frame(state = "other", frequency = sum(POG_cn_WGD_state$frequency[!(POG_cn_WGD_state$state %in% cn_keep_wgd)]))
# )
# POG_cn_WGD_state$label <- ifelse(as.character(POG_cn_WGD_state$state) == "other", "other", paste0("c=", POG_cn_WGD_state$state))
# POG_cn_WGD_state$state <- factor(as.character(POG_cn_WGD_state$state), levels = c("other", rev(cn_keep_wgd)))
# p <- ggplot(POG_cn_WGD_state, aes(x = factor(1), y = frequency, fill = factor(state))) +
#     geom_bar(stat = "identity", width = 0.5) +
#     geom_text(aes(label = label),
#         position = position_stack(vjust = 0.5),
#         color = "white", size = 20
#     ) +
#     scale_y_continuous(labels = scales::percent_format(scale = 100)) +
#     scale_fill_manual(values = fill_colors, name = "") +
#     labs(x = NULL, y = "Total CN frequencies (WGD)", title = NULL) +
#     theme_minimal() +
#     theme(
#         text = element_text(size = 65),
#         panel.background = element_rect(fill = "white"),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         legend.title = element_text(size = 40),
#         legend.position = "none",
#         panel.border = element_blank(),
#         plot.margin = margin(10, 10, 10, 10)
#     )
# png(paste0(folder_plots, "POG_state_wgd.png"), res = 150, width = 11, height = 15, units = "in")
# print(p)
# dev.off()
# POG_cn_nonWGD_multiplicity <- c()
# for (state in cn_keep_nonwgd) {
#     tmp <- as.data.frame(POG_cn_nonWGD)[which(POG_cn_nonWGD$state == state), ]
#     major_CNs <- unique(tmp$major_CN)
#     for (major_CN in major_CNs) {
#         POG_cn_nonWGD_multiplicity <- rbind(
#             POG_cn_nonWGD_multiplicity,
#             data.frame(
#                 state = state,
#                 major_CN = major_CN,
#                 genotype = paste0(state, "|", major_CN),
#                 alpha = major_CN / state,
#                 frequency = sum(tmp$width[which(tmp$major_CN == major_CN)]) / sum(tmp$width)
#             )
#         )
#     }
# }
# POG_cn_nonWGD_multiplicity$label <- paste0("tilde(m)==", POG_cn_nonWGD_multiplicity$major_CN)
# POG_cn_WGD_multiplicity <- c()
# for (state in cn_keep_wgd) {
#     tmp <- as.data.frame(POG_cn_WGD)[which(POG_cn_WGD$state == state), ]
#     major_CNs <- unique(tmp$major_CN)
#     for (major_CN in major_CNs) {
#         POG_cn_WGD_multiplicity <- rbind(
#             POG_cn_WGD_multiplicity,
#             data.frame(
#                 state = state,
#                 major_CN = major_CN,
#                 genotype = paste0(state, "|", major_CN),
#                 alpha = major_CN / state,
#                 frequency = sum(tmp$width[which(tmp$major_CN == major_CN)]) / sum(tmp$width)
#             )
#         )
#     }
# }
# POG_cn_WGD_multiplicity$label <- paste0("tilde(m)==", POG_cn_WGD_multiplicity$major_CN)
# POG_cn_nonWGD_multiplicity$genotype <- factor(as.character(POG_cn_nonWGD_multiplicity$genotype), levels = rev(names(fill_colors)))
# p <- ggplot(POG_cn_nonWGD_multiplicity, aes(x = factor(state), y = frequency, fill = factor(genotype))) +
#     geom_bar(stat = "identity", position = "stack") +
#     geom_text(aes(label = label, color = factor(genotype)),
#         position = position_stack(vjust = 0.5),
#         size = 20, parse = TRUE
#     ) +
#     scale_y_continuous(labels = scales::percent_format(scale = 100)) +
#     scale_x_discrete(labels = function(x) paste0("c=", x)) +
#     scale_fill_manual(values = fill_colors, name = "") +
#     scale_color_manual(values = c(
#         "1|1" = "white",
#         "2|1" = "white",
#         "2|2" = "black",
#         "3|3" = "white",
#         "4|3" = "white",
#         "4|2" = "black"
#     ), name = "") +
#     labs(x = NULL, y = "Major CN frequencies (non-WGD)") +
#     theme_minimal() +
#     theme(
#         text = element_text(size = 65),
#         panel.background = element_rect(fill = "white"),
#         legend.title = element_text(size = 40),
#         legend.position = "none",
#         panel.border = element_blank(),
#         plot.margin = margin(10, 10, 10, 10)
#     )
# png(paste0(folder_plots, "POG_major_cn_nonwgd.png"), res = 150, width = 15, height = 15, units = "in")
# print(p)
# dev.off()
# POG_cn_WGD_multiplicity$genotype <- factor(as.character(POG_cn_WGD_multiplicity$genotype), levels = rev(names(fill_colors)))
# p <- ggplot(POG_cn_WGD_multiplicity, aes(x = factor(state), y = frequency, fill = factor(genotype))) +
#     geom_bar(stat = "identity", position = "stack") +
#     geom_text(aes(label = label, color = factor(genotype)),
#         position = position_stack(vjust = 0.5),
#         size = 20, parse = TRUE
#     ) +
#     scale_y_continuous(labels = scales::percent_format(scale = 100)) +
#     scale_x_discrete(labels = function(x) paste0("c=", x)) +
#     scale_fill_manual(values = fill_colors, name = "") +
#     scale_color_manual(values = c(
#         "1|1" = "white",
#         "2|1" = "white",
#         "2|2" = "black",
#         "3|3" = "white",
#         "4|3" = "white",
#         "4|2" = "black",
#         "5|5" = "white",
#         "6|3" = "white",
#         "6|4" = "white",
#         "6|5" = "black",
#         "7|7" = "white",
#         "8|4" = "white",
#         "8|5" = "white",
#         "8|6" = "white",
#         "8|7" = "white",
#         "8|8" = "white"
#     ), name = "") +
#     labs(x = NULL, y = "Major CN frequencies (WGD)") +
#     theme_minimal() +
#     theme(
#         text = element_text(size = 65),
#         panel.background = element_rect(fill = "white"),
#         legend.title = element_text(size = 40),
#         legend.position = "none",
#         panel.border = element_blank(),
#         plot.margin = margin(10, 10, 10, 10)
#     )
# png(paste0(folder_plots, "POG_major_cn_wgd.png"), res = 150, width = 22.5, height = 15, units = "in")
# print(p)
# dev.off()
# # ======================================================PLOT EXAMPLE SFS
# func_sfs_POG <- function(sample_mut, sample_cn) {
#     vec_freq <- seq(0, 1, by = 0.02)
#     SFS_totalsteps <- length(vec_freq)
#     sample_mut$genotype <- "other"
#     for (row in 1:nrow(sample_cn)) {
#         chr <- sample_cn$chr[row]
#         start <- sample_cn$pos[row]
#         end <- sample_cn$end[row]
#         genotype <- paste0(sample_cn$state[row], "|", sample_cn$major_CN[row])
#         locs <- which(sample_mut$Chromosome == chr & sample_mut$Start_Position >= start & sample_mut$Start_Position <= end)
#         sample_mut$genotype[locs] <- genotype
#     }

#     unique_markers <- unique(sample_mut$genotype)
#     df_data <- data.frame(
#         ID = rep(1:length(vec_freq), length(unique_markers)),
#         frequency = rep(vec_freq, length(unique_markers)),
#         count = rep(0, SFS_totalsteps * length(unique_markers)),
#         fill = rep(unique_markers, each = SFS_totalsteps)
#     )
#     for (j in 1:nrow(sample_mut)) {
#         VAF <- sample_mut$t_alt_count[j] / sample_mut$t_depth[j]
#         pos <- which(vec_freq >= VAF)[1]
#         df_data$count[which(df_data$ID == pos & df_data$fill == sample_mut$genotype[j])] <- df_data$count[which(df_data$ID == pos & df_data$fill == sample_mut$genotype[j])] + 1
#     }
#     for (unique_marker in unique_markers) {
#         print("---------------------------")
#         print(unique_marker)
#         print(100 * sum(df_data$count[which(df_data$fill == unique_marker)]) / sum(df_data$count))
#     }
#     # df_data$fill <- factor(df_data$fill, levels = rev(names(fill_colors)))
#     genotype_order <- c("2|1", "2|2", "1|1", "3|3", "4|3", "4|2", "5|5", "6|3", "6|4", "6|5", "7|7", "8|4", "8|5", "8|6", "8|7", "8|8", "other")
#     df_data$fill <- factor(df_data$fill, levels = rev(genotype_order))
#     p <- ggplot() +
#         geom_bar(data = df_data, aes(x = frequency, y = count, fill = fill), stat = "identity", width = 1 / SFS_totalsteps) +
#         scale_fill_manual(values = fill_colors, name = "") +
#         scale_color_manual(values = fill_colors, name = "") +
#         guides(fill = guide_legend(nrow = 1, keywidth = 2, keyheight = 1)) +
#         xlab(NULL) +
#         ylab(NULL) +
#         labs(title = NULL) +
#         theme(
#             text = element_text(size = 80),
#             panel.background = element_rect(fill = "white", colour = "white"),
#             panel.grid.major = element_line(colour = "white"),
#             panel.grid.minor = element_line(colour = "white"),
#             legend.position = "none",
#             axis.title.x = element_text(color = "black"),
#             axis.ticks.y = element_blank(),
#             axis.text.y = element_blank()
#         )
#     return(p)
# }

# sample_ids_nonWGD <- unique(POG_cn_nonWGD$sample)
# sample_ids_WGD <- unique(POG_cn_WGD$sample)

# for (sample_id in sample_ids_nonWGD) {
#     sample_cn <- as.data.frame(POG_cn)[which(POG_cn$sample == sample_id), ]
#     filename <- paste0("/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/data/POG570/", sample_id, "_all.csv")
#     sample_mut <- read.table(filename, sep = ",", header = TRUE)
#     png(paste0(folder_plots, "tmp/POG_nonWGD_SFS_", sample_id, ".png"), res = 150, width = 30, height = 15, units = "in")
#     print(func_sfs_POG(sample_mut, sample_cn))
#     dev.off()
# }

# for (sample_id in sample_ids_WGD) {
#     sample_cn <- as.data.frame(POG_cn)[which(POG_cn$sample == sample_id), ]
#     filename <- paste0("/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/data/POG570/", sample_id, "_all.csv")
#     sample_mut <- read.table(filename, sep = ",", header = TRUE)
#     png(paste0(folder_plots, "tmp/POG_WGD_SFS_", sample_id, ".png"), res = 150, width = 30, height = 15, units = "in")
#     print(func_sfs_POG(sample_mut, sample_cn))
#     dev.off()
# }
