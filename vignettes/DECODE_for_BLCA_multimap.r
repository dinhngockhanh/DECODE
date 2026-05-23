# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
R_data <- "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/data/BLCA_multimap"
R_workplace <- "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes"
R_libPaths <- ""
R_libPaths_extra <- "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/R"
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
if (!dir.exists(folder_workplace)) dir.create(folder_workplace)
# ======================================GET TCGA-COAD SAMPLE INFORMATION
sample_info <- read.csv(paste0(R_data, "/BLCA_multimap_sample_information.csv"), header = TRUE)
patient_IDs <- unique(sample_info$patient_id)
sample_IDs <- sample_info$sample_id
# =============================DECODE FOR LOCAL AND WIDESPREAD MUTATIONS
data_marker_colors <- c(
    "Local" = "#3B3B3B",
    "Widespread" = "#A73030"
)
for (sample in sample_IDs) {
    #---Input the SFS data for local mutations
    filename <- paste0(R_data, "/", sample, "_local_1_1.csv")
    if (!file.exists(filename)) {
        filename <- paste0(R_data, "/", sample, "_local_all.csv")
        if (!file.exists(filename)) next
    }
    mutation_table_local <- read.csv(filename, header = TRUE)
    mutation_table_local$Ref_count <- mutation_table_local$ref_reads
    mutation_table_local$Alt_count <- mutation_table_local$alt_reads
    # mutation_table_local$Marker <- "Local"
    if (any(is.na(mutation_table_local$Ref_count)) | any(is.na(mutation_table_local$Alt_count))) mutation_table_local <- mutation_table_local[-which(is.na(mutation_table_local$Ref_count) | is.na(mutation_table_local$Alt_count)), ]
    if (any(mutation_table_local$Ref_count == 0) | any(mutation_table_local$Alt_count == 0)) mutation_table_local <- mutation_table_local[-which(mutation_table_local$Ref_count == 0 | mutation_table_local$Alt_count == 0), ]
    # #---Input the SFS data for widespread mutations
    # # filename <- paste0(R_data, "/", sample, "_widespread_1_1.csv")
    # # if (!file.exists(filename)) {
    # filename <- paste0(R_data, "/", sample, "_widespread_all.csv")
    # if (!file.exists(filename)) next
    # # }
    # mutation_table_widespread <- read.csv(filename, header = TRUE)
    # mutation_table_widespread$Ref_count <- mutation_table_widespread$ref_reads
    # mutation_table_widespread$Alt_count <- mutation_table_widespread$alt_reads
    # mutation_table_widespread$Marker <- "Widespread"
    # if (any(is.na(mutation_table_widespread$Ref_count)) | any(is.na(mutation_table_widespread$Alt_count))) mutation_table_widespread <- mutation_table_widespread[-which(is.na(mutation_table_widespread$Ref_count) | is.na(mutation_table_widespread$Alt_count)), ]
    # if (any(mutation_table_widespread$Ref_count == 0) | any(mutation_table_widespread$Alt_count == 0)) mutation_table_widespread <- mutation_table_widespread[-which(mutation_table_widespread$Ref_count == 0 | mutation_table_widespread$Alt_count == 0), ]
    #---Combine local and widespread mutations
    mutation_table <- mutation_table_local
    # mutation_table <- rbind(mutation_table_local, mutation_table_widespread)
    if (nrow(mutation_table) < 50) next
    #---SFS deconvolution with DECODE
    DECODE_result <- DECODE(
        sample_id = sample,
        mutation_table = mutation_table,
        neutral_tail = TRUE,
        read_distribution_freq_min = 25, ###############################
        compute_parallel = FALSE #######################################
    )
    save(DECODE_result, file = paste0(folder_workplace, "DECODE_", sample, ".rda"))
    #---Plot DECODE deconvolution
    DECODE_plot_model_selection(
        DECODE_result = DECODE_result,
        filename = file.path(folder_workplace, paste0("DECODE_", sample)),
        filetype = "png"
    )
}


# ========================================ANALYSES OF CLUSTERING RESULTS
color_scheme <- c(
    "NA" = "white",
    "NU" = "#767676",
    "LGIN" = "#109618",
    "HGIN" = "#FF9900",
    "UC" = "#DC3912",
    "Lower-grade neighbors" = "white",
    "UC>=UC neighbors" = "#DC3912",
    "UC<UC neighbors" = "#767676",
    "UC>=Lower-grade neighbors" = "#DC3912",
    "UC<Lower-grade neighbors" = "#767676",
    "HGIN>=Lower-grade neighbors" = "#FF9900",
    "HGIN<Lower-grade neighbors" = "#767676",
    "LGIN>=Lower-grade neighbors" = "#109618",
    "LGIN<Lower-grade neighbors" = "#767676"
)
#---Plot spatial DECODE decompositions
all_patient_results <- c()
for (patient in patient_IDs) {
    patient_info <- sample_info[which(sample_info$patient_id == patient), ]
    #   Get mutational data and DECODE results for each sample
    patient_mutational_data <- c()
    for (i in 1:length(patient_info$sample_id)) {
        sample <- patient_info$sample_id[i]
        if (!file.exists(paste0(folder_workplace, "DECODE_", sample, ".rda"))) next
        load(paste0(folder_workplace, "DECODE_", sample, ".rda"))
        sample_mutational_data <- DECODE_result$mutational_table
        sample_mutational_data$x_label <- patient_info$x_label[i]
        sample_mutational_data$y_label <- patient_info$y_label[i]
        patient_mutational_data <- rbind(patient_mutational_data, sample_mutational_data)
        parameters <- DECODE_result[[paste0("fits_", ifelse(DECODE_result$best_with_tail, "with", "without"), "_tail")]]$all_fits[[paste0(DECODE_result$best_N_clusters, "_clusters")]]$parameters
        parameters_summary <- parameters %>%
            summarise(across(where(is.numeric), list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE))))
        for (col in names(parameters_summary)) {
            patient_info[i, col] <- parameters_summary[[col]]
        }
    }
    if (length(patient_mutational_data) == 0) next
    patient_mutational_data$x_label <- factor(patient_mutational_data$x_label, levels = sort(unique(patient_mutational_data$x_label)))
    patient_mutational_data$y_label <- factor(patient_mutational_data$y_label, levels = sort(as.numeric(unique(patient_mutational_data$y_label)), decreasing = TRUE))
    #   Get histology information for each sample
    patient_result <- expand.grid(x_label = unique(patient_mutational_data$x_label), y_label = unique(patient_mutational_data$y_label))
    patient_result$sample_id <- paste0(patient, "_", patient_result$x_label, patient_result$y_label)
    patient_result <- merge(patient_result, patient_info[, !(colnames(patient_info) %in% c("x_label", "y_label"))], by = "sample_id", all.x = TRUE)
    patient_result$group[which(is.na(patient_result$Tail_power_mean))] <- "NA"
    if (length(all_patient_results) == 0) {
        all_patient_results <- patient_result[which(!is.na(patient_result$Tail_power_mean)), ]
    } else {
        for (col in setdiff(names(patient_result), names(all_patient_results))) all_patient_results[[col]] <- NA
        for (col in setdiff(names(all_patient_results), names(patient_result))) patient_result[[col]] <- NA
        all_patient_results <- rbind(all_patient_results, patient_result[which(!is.na(patient_result$Tail_power_mean)), ])
    }
    #   Text annotation for each sample
    patient_result$text <- paste0(
        "textstyle(atop(bold(", patient_result$group,
        "), textstyle(atop(n==", patient_result$Nmutations,
        ", plain(alpha)==", round(patient_result$Tail_power_mean, 2), " %+-% ", round(patient_result$Tail_power_sd, 2),
        "))))"
    )
    patient_result$text[which(is.na(patient_result$patient_id))] <- "NA"
    #   Plot spatial SFS decomposition result
    p <- ggplot(patient_mutational_data, aes(x = VAF)) +
        geom_rect(
            data = patient_result,
            aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = group),
            alpha = 0.2, inherit.aes = FALSE
        ) +
        geom_histogram(
            binwidth = 0.01, color = "#004586", fill = "#004586", position = "stack"
        ) +
        geom_text(
            data = patient_result,
            aes(x = Inf, y = Inf, label = text, color = group),
            hjust = 1, vjust = 1, size = 10, inherit.aes = FALSE, parse = TRUE
        ) +
        facet_grid(y_label ~ x_label, drop = TRUE) +
        labs(x = NULL, y = NULL) +
        scale_x_continuous(
            breaks = c(0, 0.25, 0.5, 0.75, 1),
            labels = c("0", "0.25", "0.5", "0.75", "1")
        ) +
        theme_minimal() +
        theme(
            plot.background = element_rect(fill = "white"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.text = element_text(size = 30),
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10)
        ) +
        scale_fill_manual(values = color_scheme) +
        scale_color_manual(values = color_scheme) +
        guides(fill = FALSE, color = FALSE)
    png(paste0(folder_workplace, "Spatial SFS - ", patient, ".png"), res = 500, width = 16, height = 12, units = "in")
    print(p)
    dev.off()
}
#---Plot distribution of neutral powers based on histology for each patient
for (patient in patient_IDs) {
    patient_result <- all_patient_results[which(all_patient_results$patient_id == patient), ]
    if (nrow(patient_result) == 0) next
    groups <- unique(patient_result$group)
    comparisons <- list()
    if ("LGIN" %in% groups & "NU" %in% groups) comparisons <- append(comparisons, list(c("LGIN", "NU")))
    if ("HGIN" %in% groups & "NU" %in% groups) comparisons <- append(comparisons, list(c("HGIN", "NU")))
    if ("UC" %in% groups & "NU" %in% groups) comparisons <- append(comparisons, list(c("UC", "NU")))
    patient_result$group <- factor(patient_result$group, levels = names(color_scheme))
    p <- ggplot(data = patient_result, aes(x = as.factor(group), y = Tail_power_mean, group = group, color = group, fill = group)) +
        geom_boxplot(
            alpha = 0.3, size = 2, outlier.shape = NA
        ) +
        geom_jitter(
            position = position_jitter(width = 0.2),
            size = 8
        ) +
        labs(x = NULL, y = "Neutral tail power", title = NULL) +
        theme_minimal(base_size = 30) +
        theme(
            plot.background = element_rect(fill = "white"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        scale_fill_manual(values = color_scheme) +
        scale_color_manual(values = color_scheme) +
        guides(fill = FALSE, color = FALSE, alpha = FALSE)
    if (length(comparisons) > 0) {
        p <- p +
            stat_compare_means(
                method = "wilcox.test",
                # label = "p.signif",
                label = "p.format",
                comparisons = comparisons,
                method.args = list(alternative = "greater"),
                hide.ns = TRUE,
                bracket.size = 1, vjust = 0, step.increase = 0.1, size = 10
            )
    }
    png(paste0(folder_workplace, "Neutral power - ", patient, ".png"), res = 500, width = 6, height = 12, units = "in")
    print(p)
    dev.off()
}
#---Plot distribution of neutral powers based on histology across all patients
all_patient_results$group <- factor(all_patient_results$group, levels = names(color_scheme))
p <- ggplot(data = all_patient_results, aes(x = as.factor(group), y = Tail_power_mean, group = group, color = group, fill = group)) +
    geom_boxplot(
        alpha = 0.3, size = 2, outlier.shape = NA
    ) +
    geom_jitter(
        position = position_jitter(width = 0.3),
        size = 6
    ) +
    stat_compare_means(
        method = "wilcox.test",
        # label = "p.signif",
        label = "p.format",
        comparisons = list(c("LGIN", "NU"), c("HGIN", "NU"), c("UC", "NU")),
        method.args = list(alternative = "greater"),
        hide.ns = TRUE,
        bracket.size = 1, vjust = 0, step.increase = 0.1, size = 10
    ) +
    labs(x = NULL, y = "Neutral tail power", title = NULL) +
    theme_minimal(base_size = 30) +
    theme(
        plot.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    scale_fill_manual(values = color_scheme) +
    scale_color_manual(values = color_scheme) +
    guides(fill = FALSE, color = FALSE, alpha = FALSE)
png(paste0(folder_workplace, "Neutral power - all.png"), res = 500, width = 6, height = 12, units = "in")
print(p)
dev.off()
#---Plot comparison of neutral powers between neighbor samples across patients
all_patient_results$x <- sapply(all_patient_results$x_label, function(x) which(LETTERS == x))
all_patient_results$y <- as.numeric(all_patient_results$y_label)
plot_alpha_comparison <- function(group1,
                                  group2,
                                  text1,
                                  text2,
                                  filename) {
    df_group1 <- all_patient_results[which(all_patient_results$group %in% group1), ]
    df_plot <- c()
    for (row in 1:nrow(df_group1)) {
        x <- df_group1$x[row]
        y <- df_group1$y[row]
        patient_id <- df_group1$patient_id[row]
        ids <- which(
            abs(all_patient_results$x - x) <= 1 &
                abs(all_patient_results$y - y) <= 1 &
                ((all_patient_results$x != x) | (all_patient_results$y != y)) &
                all_patient_results$patient_id != patient_id &
                all_patient_results$group %in% group2
        )
        if (length(ids) == 0) next
        df_plot <- rbind(
            df_plot,
            data.frame(
                id = paste0(patient_id, "/", df_group1$sample_id[row], "/", all_patient_results$sample_id[ids]),
                group1 = df_group1$Tail_power_mean[row],
                group2 = all_patient_results$Tail_power_mean[ids]
            )
        )
    }
    df_plot$comparison <- ifelse(df_plot$group1 >= df_plot$group2, paste0(text1, ">=", text2), paste0(text1, "<", text2))
    df_plot <- df_plot %>%
        pivot_longer(
            cols = c("group1", "group2"),
            names_to = "histology",
            values_to = "neutral_power"
        )
    df_plot$histology[which(df_plot$histology == "group1")] <- text1
    df_plot$histology[which(df_plot$histology == "group2")] <- text2
    df_plot$histology <- factor(df_plot$histology, levels = names(color_scheme))
    p <- ggplot(data = df_plot, aes(x = histology, y = neutral_power, group = id, color = comparison)) +
        geom_line(linewidth = 1, alpha = 0.3) +
        geom_point(size = 5) +
        stat_compare_means(
            method = "wilcox.test",
            # label = "p.signif",
            label = "p.format",
            comparisons = list(c(text1, text2)),
            paired = TRUE,
            method.args = list(alternative = "greater"),
            hide.ns = TRUE,
            bracket.size = 1, vjust = 0, step.increase = 0.1, size = 10
        ) +
        scale_color_manual(values = color_scheme) +
        labs(x = NULL, y = "Neutral tail power", title = NULL) +
        theme_minimal(base_size = 30) +
        theme(
            plot.background = element_rect(fill = "white"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 15, hjust = 1)
        ) +
        scale_x_discrete(expand = c(0.1, 0.1)) +
        guides(fill = FALSE, color = FALSE, alpha = FALSE, linewidth = FALSE)
    png(filename, res = 500, width = 6, height = 12, units = "in")
    print(p)
    dev.off()
}
plot_alpha_mean_comparison <- function(group1,
                                       group2,
                                       text1,
                                       text2,
                                       filename) {
    df_group1 <- all_patient_results[which(all_patient_results$group %in% group1), ]
    df_plot <- c()
    for (row in 1:nrow(df_group1)) {
        x <- df_group1$x[row]
        y <- df_group1$y[row]
        patient_id <- df_group1$patient_id[row]
        ids <- which(
            abs(all_patient_results$x - x) <= 1 &
                abs(all_patient_results$y - y) <= 1 &
                ((all_patient_results$x != x) | (all_patient_results$y != y)) &
                all_patient_results$patient_id != patient_id &
                all_patient_results$group %in% group2
        )
        if (length(ids) == 0) next
        df_plot <- rbind(
            df_plot,
            data.frame(
                id = paste0(patient_id, "/", df_group1$sample_id[row]),
                group1 = df_group1$Tail_power_mean[row],
                group2 = mean(all_patient_results$Tail_power_mean[ids])
            )
        )
    }
    df_plot$comparison <- ifelse(df_plot$group1 >= df_plot$group2, paste0(text1, ">=", text2), paste0(text1, "<", text2))
    df_plot <- df_plot %>%
        pivot_longer(
            cols = c("group1", "group2"),
            names_to = "histology",
            values_to = "neutral_power"
        )
    df_plot$histology[which(df_plot$histology == "group1")] <- text1
    df_plot$histology[which(df_plot$histology == "group2")] <- text2
    df_plot$histology <- factor(df_plot$histology, levels = names(color_scheme))
    p <- ggplot(data = df_plot, aes(x = histology, y = neutral_power, group = id, color = comparison)) +
        geom_line(linewidth = 1, alpha = 0.3) +
        geom_point(size = 5) +
        stat_compare_means(
            method = "wilcox.test",
            label = "p.format",
            comparisons = list(c(text1, text2)),
            paired = TRUE,
            method.args = list(alternative = "greater"),
            hide.ns = TRUE,
            bracket.size = 1, vjust = 0, step.increase = 0.1, size = 10
        ) +
        scale_color_manual(values = color_scheme) +
        labs(x = NULL, y = "Neutral tail power", title = NULL) +
        theme_minimal(base_size = 30) +
        theme(
            plot.background = element_rect(fill = "white"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 15, hjust = 1)
        ) +
        scale_x_discrete(expand = c(0.1, 0.1)) +
        guides(fill = FALSE, color = FALSE, alpha = FALSE, linewidth = FALSE)
    png(filename, res = 500, width = 6, height = 12, units = "in")
    print(p)
    dev.off()
}
plot_alpha_comparison(
    group1 = c("UC"),
    group2 = c("HGIN", "LGIN", "NU"),
    text1 = "UC",
    text2 = "Lower-grade neighbors",
    filename = paste0(folder_workplace, "Neutral power comparison - UC v lower grade.png")
)
plot_alpha_comparison(
    group1 = c("HGIN"),
    group2 = c("LGIN", "NU"),
    text1 = "HGIN",
    text2 = "Lower-grade neighbors",
    filename = paste0(folder_workplace, "Neutral power comparison - HGIN v lower grade.png")
)
plot_alpha_comparison(
    group1 = c("LGIN"),
    group2 = c("NU"),
    text1 = "LGIN",
    text2 = "Lower-grade neighbors",
    filename = paste0(folder_workplace, "Neutral power comparison - LGIN v lower grade.png")
)
plot_alpha_mean_comparison(
    group1 = c("UC"),
    group2 = c("HGIN", "LGIN", "NU"),
    text1 = "UC",
    text2 = "Lower-grade neighbors",
    filename = paste0(folder_workplace, "Mean neutral power comparison - UC v lower grade.png")
)
plot_alpha_mean_comparison(
    group1 = c("HGIN"),
    group2 = c("LGIN", "NU"),
    text1 = "HGIN",
    text2 = "Lower-grade neighbors",
    filename = paste0(folder_workplace, "Mean neutral power comparison - HGIN v lower grade.png")
)
plot_alpha_mean_comparison(
    group1 = c("LGIN"),
    group2 = c("NU"),
    text1 = "LGIN",
    text2 = "Lower-grade neighbors",
    filename = paste0(folder_workplace, "Mean neutral power comparison - LGIN v lower grade.png")
)
