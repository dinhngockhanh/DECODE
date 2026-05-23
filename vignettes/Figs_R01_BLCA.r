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
BLCA_multimap_mutational_data <- readRDS(paste0(R_data, "/BLCA_multimap_mutational_data.Rds"))
BLCA_multimap_mutational_cluster <- read.table("/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/GITHUB/Multiregion-bladder/vignettes/Results for 9 maps/all_maps.mutation_clusters.tsv", header = TRUE, sep = "\t")
# =================================GET ALL VAF DATA FOR SPECIFIC SAMPLES
patient_id <- "M24"
sample_ids_mini <- c("F6", "G6", "G7", "H7", "G9", "H9")
sample_ids <- paste0(patient_id, "_", sample_ids_mini)
# ================================================================DECODE
DECODE_parameters_half <- data.frame(
    sample_id = sample_ids,
    Tail_mutation_count = NA,
    Tail_power = NA,
    Cluster_mutation_count = NA,
    Cluster_VAF = NA
)
DECODE_parameters_full <- data.frame(
    sample_id = sample_ids,
    Tail_mutation_count = NA,
    Tail_power = NA,
    Cluster_1_mutation_count = NA,
    Cluster_1_VAF = NA,
    Cluster_2_mutation_count = NA,
    Cluster_2_VAF = NA,
    Cluster_3_mutation_count = NA,
    Cluster_3_VAF = NA
)
for (i_sample in 1:nrow(DECODE_parameters_half)) {
    sample <- DECODE_parameters_half$sample_id[i_sample]
    #---DECODE for local mutations
    #   Input the SFS data
    # filename <- paste0(R_data, "/", sample, "_local_1_1.csv")
    filename <- paste0(R_data, "/", sample, "_local_all.csv")
    if (!file.exists(filename)) {
        filename <- paste0(R_data, "/", sample, "_local_all.csv")
        if (!file.exists(filename)) next
    }
    mutation_table <- read.csv(filename, header = TRUE)
    mutation_table$Ref_count <- mutation_table$ref_reads
    mutation_table$Alt_count <- mutation_table$alt_reads
    if (any(is.na(mutation_table$Ref_count)) | any(is.na(mutation_table$Alt_count))) mutation_table <- mutation_table[-which(is.na(mutation_table$Ref_count) | is.na(mutation_table$Alt_count)), ]
    if (any(mutation_table$Ref_count == 0) | any(mutation_table$Alt_count == 0)) mutation_table <- mutation_table[-which(mutation_table$Ref_count == 0 | mutation_table$Alt_count == 0), ]
    # if (nrow(mutation_table) < 50) next
    if (nrow(mutation_table) < 10) next
    #   Perform DECODE
    N_humps <- 0
    DECODE_result <- DECODE(
        sample_id = sample,
        neutral_tail = TRUE, min_N_humps = N_humps, max_N_humps = N_humps,
        min_variant_read_inference_A = 10, # <<<<<<<<<<<<<<<<<<<<<<<<<<<<
        min_variant_read_inference_B = 10, # <<<<<<<<<<<<<<<<<<<<<<<<<<<<
        min_variant_read_validation = 10, # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        inference_retained_freq = 100, # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        make_readcount_distribution = FALSE,
        mutation_table = mutation_table
    )
    save(DECODE_result, file = paste0(folder_plots, "DECODE_", sample, "_local.rda"))
    DECODE_parameters_half$Tail_power[i_sample] <- DECODE_result$final_fit$parameters_df$Tail_power
    DECODE_parameters_half$Tail_mutation_count[i_sample] <- DECODE_result$final_fit$parameters_df$Tail_mutcount_predicted
    DECODE_parameters_full$Tail_power[i_sample] <- DECODE_result$final_fit$parameters_df$Tail_power
    DECODE_parameters_full$Tail_mutation_count[i_sample] <- DECODE_result$final_fit$parameters_df$Tail_mutcount_predicted
    #   Plot DECODE results
    png(paste0(folder_plots, "DECODE_", sample, "_local.png"), res = 150, width = 30, height = 15, units = "in")
    print(DECODE_plot_SFS(DECODE_result = DECODE_result))
    dev.off()
    #---DECODE for widespread mutations
    #   Input the SFS data
    filename <- paste0(R_data, "/", sample, "_widespread_1_1.csv")
    if (!file.exists(filename)) {
        filename <- paste0(R_data, "/", sample, "_widespread_all.csv")
        if (!file.exists(filename)) next
    }
    mutation_table <- read.csv(filename, header = TRUE)
    mutation_table$Ref_count <- mutation_table$ref_reads
    mutation_table$Alt_count <- mutation_table$alt_reads
    if (any(is.na(mutation_table$Ref_count)) | any(is.na(mutation_table$Alt_count))) mutation_table <- mutation_table[-which(is.na(mutation_table$Ref_count) | is.na(mutation_table$Alt_count)), ]
    if (any(mutation_table$Ref_count == 0) | any(mutation_table$Alt_count == 0)) mutation_table <- mutation_table[-which(mutation_table$Ref_count == 0 | mutation_table$Alt_count == 0), ]
    if (nrow(mutation_table) < 50) next
    DECODE_parameters_half$Cluster_mutation_count[i_sample] <- nrow(mutation_table)
    tmp <- mutation_table$Alt_count / (mutation_table$Alt_count + mutation_table$Ref_count)
    DECODE_parameters_half$Cluster_VAF[i_sample] <- median(tmp[which(tmp > 0)])
    #   Perform DECODE
    if (sample == "M24_F6") {
        N_humps <- 2
    } else {
        N_humps <- 3
    }
    DECODE_result <- DECODE(
        sample_id = sample,
        neutral_tail = FALSE, min_N_humps = N_humps, max_N_humps = N_humps,
        sfs_bincount = 25,
        make_readcount_distribution = FALSE,
        mutation_table = mutation_table
    )
    save(DECODE_result, file = paste0(folder_plots, "DECODE_", sample, "_widespread.rda"))
    DECODE_parameters_full$Cluster_1_mutation_count[i_sample] <- DECODE_result$final_fit$parameters_df$Cluster_mutcount_predicted_1
    DECODE_parameters_full$Cluster_1_VAF[i_sample] <- DECODE_result$final_fit$parameters_df$Cluster_frequency_1
    DECODE_parameters_full$Cluster_2_mutation_count[i_sample] <- DECODE_result$final_fit$parameters_df$Cluster_mutcount_predicted_2
    DECODE_parameters_full$Cluster_2_VAF[i_sample] <- DECODE_result$final_fit$parameters_df$Cluster_frequency_2
    if (N_humps > 2) {
        DECODE_parameters_full$Cluster_3_mutation_count[i_sample] <- DECODE_result$final_fit$parameters_df$Cluster_mutcount_predicted_3
        DECODE_parameters_full$Cluster_3_VAF[i_sample] <- DECODE_result$final_fit$parameters_df$Cluster_frequency_3
    }
    #   Plot DECODE results
    png(paste0(folder_plots, "DECODE_", sample, "_widespread.png"), res = 150, width = 30, height = 15, units = "in")
    print(DECODE_plot_SFS(DECODE_result = DECODE_result))
    dev.off()
}
write.csv(DECODE_parameters_half, file = paste0(folder_plots, "DECODE_parameters_half.csv"), row.names = FALSE)
write.csv(DECODE_parameters_full, file = paste0(folder_plots, "DECODE_parameters_full.csv"), row.names = FALSE)

# # =================================GET ALL VAF DATA FOR SPECIFIC SAMPLES
# BLCA_multimap_mutational_data <- BLCA_multimap_mutational_data %>%
#     filter(
#         sample_id %in% sample_ids,
#         # category == "widespread"
#     )
# BLCA_multimap_mutational_cluster <- BLCA_multimap_mutational_cluster[which(BLCA_multimap_mutational_cluster$patient_id == patient_id), ]
# BLCA_multimap_mutational_data <- left_join(
#     BLCA_multimap_mutational_data,
#     BLCA_multimap_mutational_cluster[, c("mutation_id", "cluster")],
#     by = "mutation_id"
# )
# mutational_data <- BLCA_multimap_mutational_data[, c("sample_id", "mutation_id", "effect", "cluster", "VAF")]
# tmp <- mutational_data %>%
#     select(mutation_id, cluster) %>%
#     distinct() %>%
#     deframe()
# unique_mutation_ids <- names(tmp)
# unique_mutation_clusters <- as.character(tmp)
# tmp <- mutational_data %>%
#     select(mutation_id, effect) %>%
#     distinct() %>%
#     deframe()
# unique_mutation_effects <- as.character(tmp)
# for (sample_id in sample_ids) {
#     if (length(unique(mutational_data$mutation_id[which(mutational_data$sample_id == sample_id)])) != length(unique_mutation_ids)) {
#         missing_indices <- which(!(unique_mutation_ids %in% mutational_data$mutation_id[mutational_data$sample_id == sample_id]))
#         mutational_data <- rbind(
#             mutational_data,
#             data.frame(
#                 sample_id = sample_id,
#                 mutation_id = unique_mutation_ids[missing_indices],
#                 cluster = unique_mutation_clusters[missing_indices],
#                 effect = unique_mutation_effects[missing_indices],
#                 VAF = 0
#             )
#         )
#     }
# }
# mutational_data <- mutational_data %>%
#     mutate(sample_id = sub(".*_", "", sample_id))
# data_wide <- mutational_data[, c("sample_id", "mutation_id", "VAF")] %>%
#     pivot_wider(names_from = sample_id, values_from = VAF)
# df_long <- data_wide %>%
#     pivot_longer(cols = -mutation_id, names_to = "sample", values_to = "VAF")
# pair_data <- inner_join(df_long, df_long, by = "mutation_id", suffix = c("_x", "_y")) %>%
#     filter(sample_x != sample_y)
# pair_data <- pair_data %>%
#     mutate(category = case_when(
#         VAF_x == 0 ~ "VAF_x=0",
#         VAF_y == 0 ~ "VAF_y=0",
#         VAF_x > 0 & VAF_y > 0 ~ "VAF>0"
#     ))
# pair_data <- left_join(
#     pair_data,
#     mutational_data[, c("mutation_id", "effect", "cluster")],
#     by = "mutation_id"
# )
# color_panel <- c(
#     "VAF_x=0" = "#17BECF",
#     "VAF_y=0" = "#BCBD22",
#     "VAF>0" = "#D62728"
# )
# # ===================================PLOT JOINT SFS FOR EACH SAMPLE PAIR
# alpha_data <- pair_data %>%
#     filter(cluster %in% c("alpha"))
# beta_data <- pair_data %>%
#     filter(cluster %in% c("beta"))
# gamma_data <- pair_data %>%
#     filter(cluster %in% c("gamma"))
# silent_data <- pair_data %>%
#     filter(effect %in% c("noncoding", "silent"))
# p <- ggplot(pair_data, aes(x = VAF_x, y = VAF_y, color = category)) +
#     geom_point(size = 3, alpha = 0.1) +
#     facet_grid(sample_y ~ sample_x) +
#     xlim(c(0, 1)) +
#     ylim(c(0, 1)) +
#     labs(
#         x = NULL,
#         y = NULL
#     ) +
#     scale_color_manual(values = color_panel) +
#     theme_minimal() +
#     theme(
#         plot.background = element_rect(fill = "white"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.text = element_text(size = 30),
#         strip.text.x = element_text(color = color_panel[["VAF_y=0"]], size = 30),
#         strip.text.y = element_text(color = color_panel[["VAF_x=0"]], size = 30),
#         legend.position = "none"
#     )
# pdf(paste0(folder_plots, patient_id, "_widespread_mutations_joint_SFS.pdf"), width = 16, height = 16)
# print(p)
# dev.off()
# p <- ggplot(silent_data, aes(x = VAF_x, y = VAF_y, color = category)) +
#     geom_point(size = 3, alpha = 0.1) +
#     facet_grid(sample_y ~ sample_x) +
#     xlim(c(0, 1)) +
#     ylim(c(0, 1)) +
#     labs(
#         x = NULL,
#         y = NULL
#     ) +
#     scale_color_manual(values = color_panel) +
#     theme_minimal() +
#     theme(
#         plot.background = element_rect(fill = "white"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.text = element_text(size = 30),
#         strip.text.x = element_text(color = color_panel[["VAF_y=0"]], size = 30),
#         strip.text.y = element_text(color = color_panel[["VAF_x=0"]], size = 30),
#         legend.position = "none"
#     )
# pdf(paste0(folder_plots, patient_id, "_widespread_mutations_joint_SFS (silent).pdf"), width = 16, height = 16)
# print(p)
# dev.off()
# p <- ggplot(alpha_data, aes(x = VAF_x, y = VAF_y, color = category)) +
#     geom_point(size = 3, alpha = 0.1) +
#     facet_grid(sample_y ~ sample_x) +
#     xlim(c(0, 1)) +
#     ylim(c(0, 1)) +
#     labs(
#         x = NULL,
#         y = NULL
#     ) +
#     scale_color_manual(values = color_panel) +
#     theme_minimal() +
#     theme(
#         plot.background = element_rect(fill = "white"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.text = element_text(size = 30),
#         strip.text.x = element_text(color = color_panel[["VAF_y=0"]], size = 30),
#         strip.text.y = element_text(color = color_panel[["VAF_x=0"]], size = 30),
#         legend.position = "none"
#     )
# pdf(paste0(folder_plots, patient_id, "_widespread_mutations_joint_SFS (alpha).pdf"), width = 16, height = 16)
# print(p)
# dev.off()
# p <- ggplot(beta_data, aes(x = VAF_x, y = VAF_y, color = category)) +
#     geom_point(size = 3, alpha = 0.1) +
#     facet_grid(sample_y ~ sample_x) +
#     xlim(c(0, 1)) +
#     ylim(c(0, 1)) +
#     labs(
#         x = NULL,
#         y = NULL
#     ) +
#     scale_color_manual(values = color_panel) +
#     theme_minimal() +
#     theme(
#         plot.background = element_rect(fill = "white"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.text = element_text(size = 30),
#         strip.text.x = element_text(color = color_panel[["VAF_y=0"]], size = 30),
#         strip.text.y = element_text(color = color_panel[["VAF_x=0"]], size = 30),
#         legend.position = "none"
#     )
# pdf(paste0(folder_plots, patient_id, "_widespread_mutations_joint_SFS (beta).pdf"), width = 16, height = 16)
# print(p)
# dev.off()
# p <- ggplot(gamma_data, aes(x = VAF_x, y = VAF_y, color = category)) +
#     geom_point(size = 3, alpha = 0.1) +
#     facet_grid(sample_y ~ sample_x) +
#     xlim(c(0, 1)) +
#     ylim(c(0, 1)) +
#     labs(
#         x = NULL,
#         y = NULL
#     ) +
#     scale_color_manual(values = color_panel) +
#     theme_minimal() +
#     theme(
#         plot.background = element_rect(fill = "white"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.text = element_text(size = 30),
#         strip.text.x = element_text(color = color_panel[["VAF_y=0"]], size = 30),
#         strip.text.y = element_text(color = color_panel[["VAF_x=0"]], size = 30),
#         legend.position = "none"
#     )
# pdf(paste0(folder_plots, patient_id, "_widespread_mutations_joint_SFS (gamma).pdf"), width = 16, height = 16)
# print(p)
# dev.off()
# pdf(paste0(folder_plots, patient_id, "_coverage.pdf"), width = 16, height = 16)
# p <- hist(BLCA_multimap_mutational_data$ref_reads + BLCA_multimap_mutational_data$alt_reads, breaks = 100, col = "grey", main = "Read depth distribution", xlab = "Read depth")
# print(p)
# dev.off()
