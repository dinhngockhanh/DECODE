# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
R_data <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/MK-Cod.Analysis of the SFS/R/240325 Refitting for TCGA-COAD from Stat Sci paper"
R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes"
R_libPaths <- ""
R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/R"
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)
library(data.table)
library(stringi)
library(tidyverse)
library(R.utils)
library(parallel)
library(pbapply)
library(mobster)
library(R.matlab)
library(ggplot2)
library(grid)

setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)

folder_workplace <- "TCGA-COAD/"
dir.create(folder_workplace)
# ======================================GET TCGA-COAD SAMPLE INFORMATION
sample_IDs <- read.table(paste0(R_data, "/sample_table.txt"), header = TRUE)$Patient
# ================================================================DECODE
decode_df <- data.frame()
decode_fits <- list()
# for (sample in sample_IDs) {
for (sample in c("TCGA-AA-3514")) {
    #---Input the SFS data
    filename_2 <- paste0(R_data, "/", sample, "_1_1.txt")
    mutation_table <- read.table(filename_2, sep = "\t", header = TRUE)
    mutation_table$Ref_count <- mutation_table$ref_counts
    mutation_table$Alt_count <- mutation_table$alt_counts
    #---SFS deconvolution with DECODE
    DECODE_result <- DECODE_experiment_2(
        sample_id = sample,
        mutation_table = mutation_table,
        criterion = "BIC", # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        sfs_bincount = 100,
        compute_parallel = TRUE,
        neutral_tail = TRUE # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    )
    save(DECODE_result, file = paste0(folder_workplace, "DECODE_", sample, ".rda"))
    #---Plot DECODE deconvolution
    load(paste0(folder_workplace, "DECODE_", sample, ".rda"))
    png(paste0(folder_workplace, "DECODE_", sample, ".png"), res = 150, width = 30, height = 15, units = "in")
    invisible(grid.draw(
        DECODE_plot_model_selection_experiment_2(
            DECODE_result = DECODE_result,
            data_marker_colors = c(
                "Data" = "black",
                "Foreground 0" = rgb(0.2, 0.2, 0.2),
                "Foreground 1" = rgb(0.5, 0.5, 0.5),
                "Foreground 2" = rgb(0.7, 0.7, 0.7),
                "Background 1&2" = rgb(0.9290, 0.6940, 0.1250),
                "Background 1" = rgb(0.6350, 0.0780, 0.1840),
                "Background 2" = rgb(0.4660, 0.6740, 0.1880),
                "Truncal" = rgb(0, 0.4470, 0.7410)
            )
        )
    ))
    dev.off()
}
