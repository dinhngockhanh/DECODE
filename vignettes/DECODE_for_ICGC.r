# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
R_data <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/data/PCAWG"
R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes"
R_libPaths <- ""
R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/R"
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Ginsburg
# R_data <- "/burg/iicd/users/knd2127/DECODE_DATA/PCAWG"
# R_workplace <- "/burg/iicd/users/knd2127/DECODE_TCGA_COAD"
# R_libPaths <- "/burg/iicd/users/knd2127/rpackages"
# R_libPaths_extra <- "/burg/iicd/users/knd2127/R_DECODE"
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Yining - Laptop
# R_data <- "D:/RESEARCH/DATA/data/PCAWG"
# R_workplace <- "C:/Users/Mayin/Documents/1GRADUATE/1. Study/41. Dinh_Lab/DECODE/vignettes"
# R_libPaths <- ""
# R_libPaths_extra <- "C:/Users/Mayin/Documents/1GRADUATE/1. Study/41. Dinh_Lab/DECODE/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Yining - Ginsburg
# R_data <- "/burg/iicd/users/ym2998/data/PCAWG"
# R_workplace <- "/burg/iicd/users/ym2998/PCAWG/PCAWG[0918]_page1"
# R_libPaths <- "/burg/iicd/users/ym2998/R_Package"
# R_libPaths_extra <- "/burg/iicd/users/ym2998/DECODE/R"
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)
setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)
# ===============================================SET UP WORKPLACE FOLDER
folder_workplace <- "Results_ICGC/"
if (!dir.exists(folder_workplace)) dir.create(folder_workplace)
# ===========================================GET ICGC SAMPLE INFORMATION
sample_info <- read.csv(paste0(R_data, "/ICGC_sample_information.csv"))
samples_to_delete <- c()
for (sample_index in 1:nrow(sample_info)) {
    if (sample_info$wgd_status[sample_index] == "wgd" | sample_info$wgd_uncertain[sample_index] == TRUE) {
        samples_to_delete <- c(samples_to_delete, sample)
        next
    }
    sample <- sample_info$aliquot_id[sample_index]
    filename <- paste0(R_data, "/", sample, "_1_1.csv")
    mutation_table <- read.table(filename, sep = "\t", header = TRUE)
    mutation_table$Ref_count <- mutation_table$t_ref_count
    mutation_table$Alt_count <- mutation_table$t_alt_count
    if (any(is.na(mutation_table$Ref_count)) | any(is.na(mutation_table$Alt_count))) mutation_table <- mutation_table[-which(is.na(mutation_table$Ref_count) | is.na(mutation_table$Alt_count)), ]
    if (any(mutation_table$Ref_count == 0) | any(mutation_table$Alt_count == 0)) mutation_table <- mutation_table[-which(mutation_table$Ref_count == 0 | mutation_table$Alt_count == 0), ]
    if (nrow(mutation_table) < 100) samples_to_delete <- c(samples_to_delete, sample)
}
if (length(samples_to_delete) > 0) sample_info <- sample_info[!sample_info$aliquot_id %in% samples_to_delete, ]
sample_IDs <- sample_info$aliquot_id

page <- 1
page_size <- 150
sample_range_min <- max(1, (page - 1) * page_size + 1)
sample_range_max <- min(page * page_size, length(sample_IDs))
sample_IDs_page <- sample_IDs[sample_range_min:sample_range_max]

# ===============================================================MOBSTER
library(mobster)
for (sample in sample_IDs_page) {
    #---Input the SFS data
    filename <- paste0(R_data, "/", sample, "_1_1.csv")
    mutation_table <- read.table(filename, sep = "\t", header = TRUE)
    mutation_table$Ref_count <- mutation_table$t_ref_count
    mutation_table$Alt_count <- mutation_table$t_alt_count
    if (any(is.na(mutation_table$Ref_count)) | any(is.na(mutation_table$Alt_count))) mutation_table <- mutation_table[-which(is.na(mutation_table$Ref_count) | is.na(mutation_table$Alt_count)), ]
    if (any(mutation_table$Ref_count == 0) | any(mutation_table$Alt_count == 0)) mutation_table <- mutation_table[-which(mutation_table$Ref_count == 0 | mutation_table$Alt_count == 0), ]
    mutation_table$VAF <- mutation_table$Alt_count / (mutation_table$Alt_count + mutation_table$Ref_count)
    MOBSTER_data <- data.frame(VAF = mutation_table$VAF)
    #---SFS deconvolution with MOBSTER
    MOBSTER_result <- mobster_fit(
        MOBSTER_data,
        description = sample,
        parallel = FALSE
    )
    save(MOBSTER_result, file = paste0(folder_workplace, "MOBSTER_", sample, ".rda"))
    #---Plot MOBSTER deconvolution
    png(paste0(folder_workplace, "MOBSTER_", sample, ".png"), res = 150, width = 15, height = 7.5, units = "in")
    print(plot(MOBSTER_result$best))
    dev.off()
    png(paste0(folder_workplace, "MOBSTER_model_selection_", sample, ".png"), res = 150, width = 15, height = 7.5, units = "in")
    print(plot_model_selection(MOBSTER_result))
    dev.off()
}
# ================================================================DECODE
library(grid)
for (sample in sample_IDs_page) {
    #---Input the SFS data
    filename <- paste0(R_data, "/", sample, "_1_1.csv")
    mutation_table <- read.table(filename, sep = "\t", header = TRUE)
    mutation_table$Ref_count <- mutation_table$t_ref_count
    mutation_table$Alt_count <- mutation_table$t_alt_count
    if (any(is.na(mutation_table$Ref_count)) | any(is.na(mutation_table$Alt_count))) mutation_table <- mutation_table[-which(is.na(mutation_table$Ref_count) | is.na(mutation_table$Alt_count)), ]
    if (any(mutation_table$Ref_count == 0) | any(mutation_table$Alt_count == 0)) mutation_table <- mutation_table[-which(mutation_table$Ref_count == 0 | mutation_table$Alt_count == 0), ]
    #---SFS deconvolution with DECODE
    DECODE_result <- DECODE(
        sample_id = sample,
        mutation_table = mutation_table
    )
    save(DECODE_result, file = paste0(folder_workplace, "DECODE_", sample, ".rda"))
    #---Plot DECODE deconvolution
    png(paste0(folder_workplace, "DECODE_", sample, ".png"), res = 150, width = 30, height = 15, units = "in")
    print(DECODE_plot_SFS(DECODE_result = DECODE_result))
    dev.off()
    png(paste0(folder_workplace, "DECODE_model_selection_", sample, ".png"), res = 150, width = 30, height = 15, units = "in")
    grid.draw(
        DECODE_plot_model_selection(
            DECODE_result = DECODE_result
        )
    )
    dev.off()
}
# ======================================EXTRACT DECONVOLUTION PARAMETERS
mobster_df <- data.frame()
decode_df <- data.frame()
for (sample in sample_IDs) {
    #---Extract MOBSTER parameters
    load(paste0(folder_workplace, "MOBSTER_", sample, ".rda"))
    mobster_df <- MOBSTER_summary_statistics(mobster_df, MOBSTER_result)
    # #---Extract DECODE parameters
    load(paste0(folder_workplace, "DECODE_", sample, ".rda"))
    decode_df <- DECODE_summary_statistics(decode_df, DECODE_result)
}
mobster_df <- merge(mobster_df, sample_info[, c("aliquot_id", "purity")], by.x = "Sample", by.y = "aliquot_id", all.x = TRUE)
names(mobster_df)[names(mobster_df) == "purity"] <- "Purity"
mobster_df <- merge(mobster_df, sample_info[, c("aliquot_id", "histology_abbreviation")], by.x = "Sample", by.y = "aliquot_id", all.x = TRUE)
names(mobster_df)[names(mobster_df) == "histology_abbreviation"] <- "Cancer_type"
decode_df <- merge(decode_df, sample_info[, c("aliquot_id", "purity")], by.x = "Sample", by.y = "aliquot_id", all.x = TRUE)
names(decode_df)[names(decode_df) == "purity"] <- "Purity"
decode_df <- merge(decode_df, sample_info[, c("aliquot_id", "histology_abbreviation")], by.x = "Sample", by.y = "aliquot_id", all.x = TRUE)
names(decode_df)[names(decode_df) == "histology_abbreviation"] <- "Cancer_type"
write.csv(mobster_df, paste0(folder_workplace, "MOBSTER_ICGC.csv"), row.names = FALSE)
write.csv(decode_df, paste0(folder_workplace, "DECODE_ICGC.csv"), row.names = FALSE)
# =========================================ANALYZE DECONVOLUTION RESULTS
#---Input plot settings
source("plot_settings.r")
#---Make plots for analysis of MOBSTER results
mobster_df <- read.csv(paste0(folder_workplace, "MOBSTER_ICGC.csv"))
plot_analysis(
    results = mobster_df,
    algorithm = "MOBSTER",
    cohort = "ICGC",
    algorithm_colors = algorithm_colors,
    cluster_shapes = cluster_shapes,
    cluster_labels = cluster_labels,
    cluster_colors = cluster_colors,
    cohort_colors = ICGC_cohort_colors,
    folder_workplace = folder_workplace
)
#---Make plots for analysis of DECODE results
decode_df <- read.csv(paste0(folder_workplace, "DECODE_ICGC.csv"))
plot_analysis(
    results = decode_df,
    algorithm = "DECODE",
    cohort = "ICGC",
    algorithm_colors = algorithm_colors,
    cluster_shapes = cluster_shapes,
    cluster_labels = cluster_labels,
    cluster_colors = cluster_colors,
    cohort_colors = ICGC_cohort_colors,
    folder_workplace = folder_workplace
)
