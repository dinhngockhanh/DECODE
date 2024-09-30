# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
R_data <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/data/PCAWG"
R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes"
R_libPaths <- ""
R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/R"
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Yining - Laptop
# R_data <- "D:/RESEARCH/DATA/data/PCAWG"
# R_workplace <- "C:/Users/Mayin/Documents/1GRADUATE/1. Study/41. Dinh_Lab/DECODE/vignettes"
# R_libPaths <- ""
# R_libPaths_extra <- "C:/Users/Mayin/Documents/1GRADUATE/1. Study/41. Dinh_Lab/DECODE/R"
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)
setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)
# ===========================================GET ICGC SAMPLE INFORMATION
sample_info <- read.csv(paste0(R_data, "/ICGC_sample_information.csv"))
sample <- sample_info$aliquot_id[1]
filename <- paste0(R_data, "/", sample, "_1_1.csv")
mutation_table <- read.table(filename, sep = "\t", header = TRUE)
mutation_table$Ref_count <- mutation_table$t_ref_count
mutation_table$Alt_count <- mutation_table$t_alt_count
if (any(is.na(mutation_table$Ref_count)) | any(is.na(mutation_table$Alt_count))) mutation_table <- mutation_table[-which(is.na(mutation_table$Ref_count) | is.na(mutation_table$Alt_count)), ]
if (any(mutation_table$Ref_count == 0) | any(mutation_table$Alt_count == 0)) mutation_table <- mutation_table[-which(mutation_table$Ref_count == 0 | mutation_table$Alt_count == 0), ]
#---SFS deconvolution with DECODE
DECODE_Yining(
    sample_id = sample,
    mutation_table = mutation_table
)
