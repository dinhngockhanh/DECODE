# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
R_data <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/data/PCAWG"
R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes"
R_libPaths <- ""
R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/R"
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Ginsburg
# R_data <- "/burg/iicd/users/knd2127/DECODE_PCAWG/PCAWG"
# R_workplace <- "/burg/iicd/users/knd2127/DECODE_PCAWG"
# R_libPaths <- "/burg/iicd/users/knd2127/rpackages"
# R_libPaths_extra <- "/burg/iicd/users/knd2127/R_DECODE"
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
sample_info <- sample_info[
    which(
        sample_info$wgd_status == "no_wgd" &
            sample_info$wgd_uncertain == FALSE
    ),
]
sample_IDs <- sample_info$aliquot_id
# ===============================================================MOBSTER
library(mobster)
for (sample in sample_IDs) {
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
        description = sample
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
for (sample in sample_IDs) {
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
