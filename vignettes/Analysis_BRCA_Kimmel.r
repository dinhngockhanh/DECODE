library(readxl)
# # =========================================DATA PREPARATION - 2025/08/19
# brca_data <- read_excel("/Users/khanhndinh/RESEARCH AND EVERYTHING/Projects/DATASETS/BRCA [Kimmel]/BRCA_extended_ALL_TOOLS.xlsx")
# for (sampleID in unique(brca_data$SampleID)) {
#     #---Prepare data for P1 samples
#     rows <- which(
#         brca_data$SampleID == sampleID &
#             brca_data$P1_SD_ffpe_filter == "snv" &
#             brca_data$P1_major_cn == 1 &
#             brca_data$P1_minor_cn == 1
#     )
#     if (length(rows) <= 10) next
#     cat("Number of rows for sample", sampleID, ":", length(rows), "\n")
#     subset_data_P1 <- brca_data[rows, ]
#     subset_data_P1$Ref_count <- as.integer(sub(",.*", "", subset_data_P1$P1_AD))
#     subset_data_P1$Alt_count <- as.integer(sub(".*,", "", subset_data_P1$P1_AD))
#     if (any(subset_data_P1$Ref_count == 0) | any(subset_data_P1$Alt_count == 0)) subset_data_P1 <- subset_data_P1[-which(subset_data_P1$Ref_count == 0 | subset_data_P1$Alt_count == 0), ]
#     csv_filename <- file.path(
#         "/Users/khanhndinh/RESEARCH AND EVERYTHING/Projects/DATASETS/BRCA [Kimmel]/DECODE_input_250819",
#         paste0(sampleID, "_P1.csv")
#     )
#     write.csv(subset_data_P1, file = csv_filename, row.names = FALSE)
#     #---Prepare data for L1 samples
#     rows <- which(
#         brca_data$SampleID == sampleID &
#             brca_data$L1_SD_ffpe_filter == "snv" &
#             brca_data$L1_major_cn == 1 &
#             brca_data$L1_minor_cn == 1
#     )
#     if (length(rows) <= 10) next
#     cat("Number of rows for sample", sampleID, ":", length(rows), "\n")
#     subset_data_L1 <- brca_data[rows, ]
#     subset_data_L1$Ref_count <- as.integer(sub(",.*", "", subset_data_L1$L1_AD))
#     subset_data_L1$Alt_count <- as.integer(sub(".*,", "", subset_data_L1$L1_AD))
#     if (any(subset_data_L1$Ref_count == 0) | any(subset_data_L1$Alt_count == 0)) subset_data_L1 <- subset_data_L1[-which(subset_data_L1$Ref_count == 0 | subset_data_L1$Alt_count == 0), ]
#     csv_filename <- file.path(
#         "/Users/khanhndinh/RESEARCH AND EVERYTHING/Projects/DATASETS/BRCA [Kimmel]/DECODE_input_250819",
#         paste0(sampleID, "_L1.csv")
#     )
#     write.csv(subset_data_L1, file = csv_filename, row.names = FALSE)
# }
# # ============================DATA PREPARATION - 2025/08/20 (LARGE DATA)
# dir.create("/Users/khanhndinh/RESEARCH AND EVERYTHING/Projects/DATASETS/BRCA [Kimmel]/DECODE_input_250820_large", showWarnings = FALSE, recursive = TRUE)
# brca_L1_data <- read_excel("/Users/khanhndinh/RESEARCH AND EVERYTHING/Projects/DATASETS/BRCA [Kimmel]/FFPolish_enriched/large/L1_BRCA_FFPolish_enriched_large.xlsx")
# rows <- which(
#     brca_L1_data$L1_major_cn == 1 &
#         brca_L1_data$L1_minor_cn == 1
# )
# brca_L1_data <- brca_L1_data[rows, ]
# for (sampleID in unique(brca_L1_data$SampleID)) {
#     rows <- which(brca_L1_data$SampleID == sampleID)
#     subset_data_L1 <- brca_L1_data[rows, ]
#     if (nrow(subset_data_L1) <= 10) next
#     cat("Number of rows for sample", sampleID, ":", nrow(subset_data_L1), "\n")
#     subset_data_L1$Ref_count <- as.integer(sub(",.*", "", subset_data_L1$L1_AD))
#     subset_data_L1$Alt_count <- as.integer(sub(".*,", "", subset_data_L1$L1_AD))
#     if (any(subset_data_L1$Ref_count == 0) | any(subset_data_L1$Alt_count == 0)) subset_data_L1 <- subset_data_L1[-which(subset_data_L1$Ref_count == 0 | subset_data_L1$Alt_count == 0), ]
#     if (nrow(subset_data_L1) <= 10) next
#     csv_filename <- file.path(
#         "/Users/khanhndinh/RESEARCH AND EVERYTHING/Projects/DATASETS/BRCA [Kimmel]/DECODE_input_250820_large",
#         paste0(sampleID, "_L1.csv")
#     )
#     write.csv(subset_data_L1, file = csv_filename, row.names = FALSE)
# }
# brca_P1_data <- read_excel("/Users/khanhndinh/RESEARCH AND EVERYTHING/Projects/DATASETS/BRCA [Kimmel]/FFPolish_enriched/large/P1_BRCA_FFPolish_enriched_large.xlsx")
# rows <- which(
#     brca_P1_data$P1_major_cn == 1 &
#         brca_P1_data$P1_minor_cn == 1
# )
# brca_P1_data <- brca_P1_data[rows, ]
# for (sampleID in unique(brca_P1_data$SampleID)) {
#     rows <- which(brca_P1_data$SampleID == sampleID)
#     subset_data_P1 <- brca_P1_data[rows, ]
#     if (nrow(subset_data_P1) <= 10) next
#     cat("Number of rows for sample", sampleID, ":", nrow(subset_data_P1), "\n")
#     subset_data_P1$Ref_count <- as.integer(sub(",.*", "", subset_data_P1$P1_AD))
#     subset_data_P1$Alt_count <- as.integer(sub(".*,", "", subset_data_P1$P1_AD))
#     if (any(subset_data_P1$Ref_count == 0) | any(subset_data_P1$Alt_count == 0)) subset_data_P1 <- subset_data_P1[-which(subset_data_P1$Ref_count == 0 | subset_data_P1$Alt_count == 0), ]
#     if (nrow(subset_data_P1) <= 10) next
#     csv_filename <- file.path(
#         "/Users/khanhndinh/RESEARCH AND EVERYTHING/Projects/DATASETS/BRCA [Kimmel]/DECODE_input_250820_large",
#         paste0(sampleID, "_P1.csv")
#     )
#     write.csv(subset_data_P1, file = csv_filename, row.names = FALSE)
# }
# # ============================DATA PREPARATION - 2025/08/20 (SMALL DATA)
# dir.create("/Users/khanhndinh/RESEARCH AND EVERYTHING/Projects/DATASETS/BRCA [Kimmel]/DECODE_input_250820_small", showWarnings = FALSE, recursive = TRUE)
# brca_L1_data <- read_excel("/Users/khanhndinh/RESEARCH AND EVERYTHING/Projects/DATASETS/BRCA [Kimmel]/FFPolish_enriched/small/L1_BRCA_FFPolish_enriched_small.xlsx")
# rows <- which(
#     brca_L1_data$L1_major_cn == 1 &
#         brca_L1_data$L1_minor_cn == 1
# )
# brca_L1_data <- brca_L1_data[rows, ]
# for (sampleID in unique(brca_L1_data$SampleID)) {
#     rows <- which(brca_L1_data$SampleID == sampleID)
#     subset_data_L1 <- brca_L1_data[rows, ]
#     if (nrow(subset_data_L1) <= 10) next
#     cat("Number of rows for sample", sampleID, ":", nrow(subset_data_L1), "\n")
#     subset_data_L1$Ref_count <- as.integer(sub(",.*", "", subset_data_L1$L1_AD))
#     subset_data_L1$Alt_count <- as.integer(sub(".*,", "", subset_data_L1$L1_AD))
#     if (any(subset_data_L1$Ref_count == 0) | any(subset_data_L1$Alt_count == 0)) subset_data_L1 <- subset_data_L1[-which(subset_data_L1$Ref_count == 0 | subset_data_L1$Alt_count == 0), ]
#     if (nrow(subset_data_L1) <= 10) next
#     csv_filename <- file.path(
#         "/Users/khanhndinh/RESEARCH AND EVERYTHING/Projects/DATASETS/BRCA [Kimmel]/DECODE_input_250820_small",
#         paste0(sampleID, "_L1.csv")
#     )
#     write.csv(subset_data_L1, file = csv_filename, row.names = FALSE)
# }
# brca_P1_data <- read_excel("/Users/khanhndinh/RESEARCH AND EVERYTHING/Projects/DATASETS/BRCA [Kimmel]/FFPolish_enriched/small/P1_BRCA_FFPolish_enriched_small.xlsx")
# rows <- which(
#     brca_P1_data$P1_major_cn == 1 &
#         brca_P1_data$P1_minor_cn == 1
# )
# brca_P1_data <- brca_P1_data[rows, ]
# for (sampleID in unique(brca_P1_data$SampleID)) {
#     rows <- which(brca_P1_data$SampleID == sampleID)
#     subset_data_P1 <- brca_P1_data[rows, ]
#     if (nrow(subset_data_P1) <= 10) next
#     cat("Number of rows for sample", sampleID, ":", nrow(subset_data_P1), "\n")
#     subset_data_P1$Ref_count <- as.integer(sub(",.*", "", subset_data_P1$P1_AD))
#     subset_data_P1$Alt_count <- as.integer(sub(".*,", "", subset_data_P1$P1_AD))
#     if (any(subset_data_P1$Ref_count == 0) | any(subset_data_P1$Alt_count == 0)) subset_data_P1 <- subset_data_P1[-which(subset_data_P1$Ref_count == 0 | subset_data_P1$Alt_count == 0), ]
#     if (nrow(subset_data_P1) <= 10) next
#     csv_filename <- file.path(
#         "/Users/khanhndinh/RESEARCH AND EVERYTHING/Projects/DATASETS/BRCA [Kimmel]/DECODE_input_250820_small",
#         paste0(sampleID, "_P1.csv")
#     )
#     write.csv(subset_data_P1, file = csv_filename, row.names = FALSE)
# }
# # =========================================DATA PREPARATION - 2025/08/26
# brca_data <- read_excel("/Users/khanhndinh/RESEARCH AND EVERYTHING/Projects/DATASETS/BRCA [Kimmel]/BRCA_extended_ALL_TOOLS.xlsx")
# for (sampleID in unique(brca_data$SampleID)) {
#     #---Prepare data for P1 samples
#     rows <- which(
#         brca_data$SampleID == sampleID &
#             brca_data$P1_DO_FFPE_score >= 0.25 &
#             brca_data$P1_major_cn == 1 &
#             brca_data$P1_minor_cn == 1
#     )
#     if (length(rows) <= 10) next
#     cat("Number of rows for sample", sampleID, ":", length(rows), "\n")
#     subset_data_P1 <- brca_data[rows, ]
#     subset_data_P1$Ref_count <- as.integer(sub(",.*", "", subset_data_P1$P1_AD))
#     subset_data_P1$Alt_count <- as.integer(sub(".*,", "", subset_data_P1$P1_AD))
#     if (any(subset_data_P1$Ref_count == 0) | any(subset_data_P1$Alt_count == 0)) subset_data_P1 <- subset_data_P1[-which(subset_data_P1$Ref_count == 0 | subset_data_P1$Alt_count == 0), ]
#     csv_filename <- file.path(
#         "/Users/khanhndinh/RESEARCH AND EVERYTHING/Projects/DATASETS/BRCA [Kimmel]/DECODE_input_250826",
#         paste0(sampleID, "_P1.csv")
#     )
#     write.csv(subset_data_P1, file = csv_filename, row.names = FALSE)
#     #---Prepare data for L1 samples
#     rows <- which(
#         brca_data$SampleID == sampleID &
#             brca_data$L1_DO_FFPE_score >= 0.25 &
#             brca_data$L1_major_cn == 1 &
#             brca_data$L1_minor_cn == 1
#     )
#     if (length(rows) <= 10) next
#     cat("Number of rows for sample", sampleID, ":", length(rows), "\n")
#     subset_data_L1 <- brca_data[rows, ]
#     subset_data_L1$Ref_count <- as.integer(sub(",.*", "", subset_data_L1$L1_AD))
#     subset_data_L1$Alt_count <- as.integer(sub(".*,", "", subset_data_L1$L1_AD))
#     if (any(subset_data_L1$Ref_count == 0) | any(subset_data_L1$Alt_count == 0)) subset_data_L1 <- subset_data_L1[-which(subset_data_L1$Ref_count == 0 | subset_data_L1$Alt_count == 0), ]
#     csv_filename <- file.path(
#         "/Users/khanhndinh/RESEARCH AND EVERYTHING/Projects/DATASETS/BRCA [Kimmel]/DECODE_input_250826",
#         paste0(sampleID, "_L1.csv")
#     )
#     write.csv(subset_data_L1, file = csv_filename, row.names = FALSE)
# }
# # =========================================DATA PREPARATION - 2025/12/11
# brca_data <- read_excel("/Users/kndinh/DINH LAB/DATASETS/BRCA [Kimmel]/MK_LUAD_WES_Mutect2_PASS_v6_unique_extended_ASCATsc_cn.xlsx")
# print(unique(brca_data$SampleID))
# for (sampleID in unique(brca_data$SampleID)) {
#     #---Prepare data for P1 samples
#     rows <- which(
#         brca_data$SampleID == sampleID &
#             brca_data$P1_total_cn == 2
#     )
#     if (length(rows) <= 10) next
#     cat("Number of rows for sample", sampleID, ":", length(rows), "\n")
#     subset_data_P1 <- brca_data[rows, ]
#     subset_data_P1$Ref_count <- as.integer(sub(",.*", "", subset_data_P1$P1_AD))
#     subset_data_P1$Alt_count <- as.integer(sub(".*,", "", subset_data_P1$P1_AD))
#     if (any(subset_data_P1$Ref_count == 0) | any(subset_data_P1$Alt_count == 0)) subset_data_P1 <- subset_data_P1[-which(subset_data_P1$Ref_count == 0 | subset_data_P1$Alt_count == 0), ]
#     csv_filename <- file.path(
#         "/Users/kndinh/DINH LAB/DATASETS/BRCA [Kimmel]/DECODE_input_251211",
#         paste0(sampleID, "_P1.csv")
#     )
#     write.csv(subset_data_P1, file = csv_filename, row.names = FALSE)
#     if (sampleID == "Lu03") print(subset_data_P1$Ref_count)
#     #---Prepare data for P2 samples
#     rows <- which(
#         brca_data$SampleID == sampleID &
#             brca_data$P2_total_cn == 2
#     )
#     if (length(rows) <= 10) next
#     cat("Number of rows for sample", sampleID, ":", length(rows), "\n")
#     subset_data_P2 <- brca_data[rows, ]
#     subset_data_P2$Ref_count <- as.integer(sub(",.*", "", subset_data_P2$P2_AD))
#     subset_data_P2$Alt_count <- as.integer(sub(".*,", "", subset_data_P2$P2_AD))
#     if (any(subset_data_P2$Ref_count == 0) | any(subset_data_P2$Alt_count == 0)) subset_data_P2 <- subset_data_P2[-which(subset_data_P2$Ref_count == 0 | subset_data_P2$Alt_count == 0), ]
#     csv_filename <- file.path(
#         "/Users/kndinh/DINH LAB/DATASETS/BRCA [Kimmel]/DECODE_input_251211",
#         paste0(sampleID, "_P2.csv")
#     )
#     write.csv(subset_data_P2, file = csv_filename, row.names = FALSE)
#     if (sampleID == "Lu03") print(subset_data_P2$Ref_count)
#     #---Prepare data for L1 samples
#     rows <- which(
#         brca_data$SampleID == sampleID &
#             brca_data$L1_total_cn == 2
#     )
#     if (length(rows) <= 10) next
#     cat("Number of rows for sample", sampleID, ":", length(rows), "\n")
#     subset_data_L1 <- brca_data[rows, ]
#     subset_data_L1$Ref_count <- as.integer(sub(",.*", "", subset_data_L1$L1_AD))
#     subset_data_L1$Alt_count <- as.integer(sub(".*,", "", subset_data_L1$L1_AD))
#     if (any(subset_data_L1$Ref_count == 0) | any(subset_data_L1$Alt_count == 0)) subset_data_L1 <- subset_data_L1[-which(subset_data_L1$Ref_count == 0 | subset_data_L1$Alt_count == 0), ]
#     csv_filename <- file.path(
#         "/Users/kndinh/DINH LAB/DATASETS/BRCA [Kimmel]/DECODE_input_251211",
#         paste0(sampleID, "_L1.csv")
#     )
#     write.csv(subset_data_L1, file = csv_filename, row.names = FALSE)
# }
# ====================================================DECODE - AUTOMATIC
#---Set up folder paths, workplace and libraries
R_data <- "/Users/kndinh/DINH LAB/DATASETS/BRCA [Kimmel]/DECODE_input_251211"
R_workplace <- "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes"
R_libPaths <- ""
R_libPaths_extra <- "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/R"
.libPaths(R_libPaths)
setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)
folder_workplace <- "Results_BRCA_Kimmel_without_tail/"
if (!dir.exists(folder_workplace)) dir.create(folder_workplace)
library(grid)
#---DECODE
sample_IDs <- list.files(R_data, pattern = "\\.csv$", full.names = FALSE)
sample_IDs <- sub("\\.csv$", "", sample_IDs)
for (iii in 1:length(sample_IDs)) {
    sample <- sample_IDs[iii]
    if (!grepl("_L1", sample)) next ##################################################################
    #---Input the SFS data
    filename <- paste0(R_data, "/", sample, ".csv")
    if (!file.exists(filename)) next
    mutation_table <- read.table(filename, sep = ",", header = TRUE)
    #---SFS deconvolution with DECODE
    DECODE_result <- DECODE(
        sample_id = sample,
        mutation_table = mutation_table,
        neutral_tail = FALSE,
        # min_variant_read_inference_A = 1, min_variant_read_inference_B = 1, min_variant_read_validation = 1, inference_retained_freq = 100,
        criterion_tail_weight = 0.3
    )
    save(DECODE_result, file = paste0(folder_workplace, "DECODE_", sample, ".rda"))
    #---Plot DECODE deconvolution
    DECODE_plot_SFS(
        DECODE_result = DECODE_result,
        filename = file.path(folder_workplace, paste0("DECODE_", sample)),
        filetype = "png"
    )
    DECODE_plot_model_selection(
        DECODE_result = DECODE_result,
        filename = file.path(folder_workplace, paste0("DECODE_model_selection_", sample)),
        filetype = "png"
    )
}
# # =======================================================DECODE - MANUAL
# #---Set up folder paths, workplace and libraries
# R_data <- "/Users/khanhndinh/RESEARCH AND EVERYTHING/Projects/DATASETS/BRCA [Kimmel]/DECODE_input_250826"
# R_workplace <- "/Users/khanhndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/khanhndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/R/DECODE_OLD"
# .libPaths(R_libPaths)
# setwd(R_libPaths_extra)
# files_sources <- list.files(pattern = "\\.[rR]$")
# sapply(files_sources, source)
# setwd(R_workplace)
# folder_workplace <- "Results_BRCA_Kimmel/"
# if (!dir.exists(folder_workplace)) dir.create(folder_workplace)
# library(grid)
# #---DECODE
# df_N_humps <- data.frame(matrix(NA, nrow = 0, ncol = 2))
# # df_N_humps[nrow(df_N_humps) + 1, ] <- c("G2_L1", 1)
# df_N_humps[nrow(df_N_humps) + 1, ] <- c("G4_L1", 0)
# df_N_humps[nrow(df_N_humps) + 1, ] <- c("G4_P1", 0)
# # df_N_humps[nrow(df_N_humps) + 1, ] <- c("G30_L1", 0)
# df_N_humps[nrow(df_N_humps) + 1, ] <- c("G30_P1", 0)
# df_N_humps[nrow(df_N_humps) + 1, ] <- c("G31_L1", 1)
# df_N_humps[nrow(df_N_humps) + 1, ] <- c("G31_P1", 1)
# df_N_humps[nrow(df_N_humps) + 1, ] <- c("G32_L1", 0)
# df_N_humps[nrow(df_N_humps) + 1, ] <- c("G32_P1", 1)
# df_N_humps[nrow(df_N_humps) + 1, ] <- c("G33_L1", 1)
# df_N_humps[nrow(df_N_humps) + 1, ] <- c("G33_P1", 1)
# df_N_humps[nrow(df_N_humps) + 1, ] <- c("G35_L1", 0)
# df_N_humps[nrow(df_N_humps) + 1, ] <- c("G35_P1", 1)
# df_N_humps[nrow(df_N_humps) + 1, ] <- c("G36_L1", 1)
# df_N_humps[nrow(df_N_humps) + 1, ] <- c("G36_P1", 1)
# df_N_humps[nrow(df_N_humps) + 1, ] <- c("G40_L1", 1)
# df_N_humps[nrow(df_N_humps) + 1, ] <- c("G40_P1", 1)
# df_N_humps[nrow(df_N_humps) + 1, ] <- c("G41_L1", 1)
# df_N_humps[nrow(df_N_humps) + 1, ] <- c("G41_P1", 1)
# df_N_humps[nrow(df_N_humps) + 1, ] <- c("G43_P1", 1)
# # df_N_humps[nrow(df_N_humps) + 1, ] <- c("G45_L1", 0)
# df_N_humps[nrow(df_N_humps) + 1, ] <- c("G45_P1", 0)
# df_N_humps[nrow(df_N_humps) + 1, ] <- c("G46_L1", 0)
# df_N_humps[nrow(df_N_humps) + 1, ] <- c("G46_P1", 1)
# df_N_humps[nrow(df_N_humps) + 1, ] <- c("G47_L1", 1)
# df_N_humps[nrow(df_N_humps) + 1, ] <- c("G47_P1", 1)
# df_N_humps[nrow(df_N_humps) + 1, ] <- c("G48_L1", 1)
# df_N_humps[nrow(df_N_humps) + 1, ] <- c("G48_P1", 1)
# colnames(df_N_humps) <- c("sampleID", "N_humps")
# for (iii in 1:nrow(df_N_humps)) {
#     sample <- df_N_humps$sampleID[iii]
#     #---Input the SFS data
#     filename <- paste0(R_data, "/", sample, ".csv")
#     if (!file.exists(filename)) next
#     mutation_table <- read.table(filename, sep = ",", header = TRUE)
#     #---SFS deconvolution with DECODE
#     DECODE_result <- DECODE(
#         sample_id = sample,
#         mutation_table = mutation_table,
#         min_variant_read_inference_A = 1, min_variant_read_inference_B = 1, min_variant_read_validation = 1, inference_retained_freq = 100,
#         neutral_tail = TRUE, min_N_humps = as.integer(df_N_humps$N_humps[iii]), max_N_humps = as.integer(df_N_humps$N_humps[iii])
#     )
#     save(DECODE_result, file = paste0(folder_workplace, "DECODE_", sample, ".rda"))
#     #---Plot DECODE deconvolution
#     png(paste0(folder_workplace, "DECODE_", sample, ".png"), res = 150, width = 30, height = 15, units = "in")
#     print(DECODE_plot_SFS(DECODE_result = DECODE_result))
#     dev.off()
#     png(paste0(folder_workplace, "DECODE_model_selection_", sample, ".png"), res = 150, width = 30, height = 15, units = "in")
#     grid.draw(
#         DECODE_plot_model_selection(
#             DECODE_result = DECODE_result
#         )
#     )
#     dev.off()
# }
# ============================================ANALYSIS OF DECODE RESULTS


setwd("/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes/Results_BRCA_Kimmel_without_tail/")
folder_workplace <- "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes/Results_BRCA_Kimmel_without_tail/"
sample_IDs <- list.files(pattern = "\\.rda$", full.names = FALSE)
df_N_humps <- data.frame(sampleID = sample_IDs)
df_N_humps$Tail <- NA
df_N_humps$N_humps <- NA
for (iii in 1:nrow(df_N_humps)) {
    sample <- df_N_humps$sampleID[iii]
    load(sample)
    df_N_humps$Tail[iii] <- DECODE_result$best_with_tail
    df_N_humps$N_humps[iii] <- DECODE_result$best_N_clusters
    if (df_N_humps$Tail[iii] == TRUE) {
        parameters <- DECODE_result$fits_with_tail$all_fits[[paste0(DECODE_result$best_N_clusters, "_clusters")]]$parameters
    } else {
        parameters <- DECODE_result$fits_without_tail$all_fits[[paste0(DECODE_result$best_N_clusters, "_clusters")]]$parameters
    }
    parameters <- colMeans(as.data.frame(parameters))

    if (df_N_humps$Tail[iii] == TRUE) {
        df_N_humps$Tail_power[iii] <- parameters[["Tail_power"]]
        df_N_humps$Tail_mutcount[iii] <- parameters[["Tail_Nmut_exact"]]
    }
    for (cluster_idx in 1:df_N_humps$N_humps[iii]) {
        cluster_p_name <- paste0("Cluster_", cluster_idx, "_freq")
        cluster_mutcount_name <- paste0("Cluster_", cluster_idx, "_Nmut_exact")
        df_N_humps[iii, cluster_p_name] <- parameters[[cluster_p_name]]
        df_N_humps[iii, cluster_mutcount_name] <- parameters[[cluster_mutcount_name]]
    }
}
write.csv(df_N_humps, file = paste0(folder_workplace, "DECODE_results.csv"), row.names = FALSE)


# R_libPaths_extra <- "/Users/khanhndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/R/DECODE_v00"
# setwd(R_libPaths_extra)
# files_sources <- list.files(pattern = "\\.[rR]$")
# sapply(files_sources, source)
# library(grid)
# folder_workplace <- "/Users/khanhndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes/Previous results/BRCA_Kimmel_250828_manual/"
# samples <- c("G33_L1", "G33_P1", "G40_L1", "G40_P1", "G41_L1", "G41_P1")
# for (sample in samples) {
#     load(paste0(folder_workplace, "DECODE_", sample, ".rda"))
#     print(DECODE_plot_SFS(DECODE_result = DECODE_result))
#     png(paste0(folder_workplace, "new_DECODE_", sample, ".png"), res = 150, width = 15, height = 7.5, units = "in")
#     # png(paste0(folder_workplace, "new_DECODE_", sample, ".png"), res = 150, width = 30, height = 15, units = "in")
#     print(DECODE_plot_SFS(DECODE_result = DECODE_result))
#     dev.off()
# }
