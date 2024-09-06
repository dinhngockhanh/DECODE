# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
R_data <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/data/TCGA_COAD"
R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes"
R_libPaths <- ""
R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/R"
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)
library(grid)
library(parallel)
library(pbapply)
setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)
# ===============================================SET UP WORKPLACE FOLDER
folder_workplace <- "TCGA-COAD/"
if (!dir.exists(folder_workplace)) dir.create(folder_workplace)
# ======================================GET TCGA-COAD SAMPLE INFORMATION
sample_IDs <- read.table(paste0(R_data, "/sample_table.txt"), header = TRUE)$Patient
# ===============================================================MOBSTER
numCores <- detectCores()
cl <- makePSOCKcluster(numCores - 1)
if (is.null(R_libPaths) == FALSE) {
    R_libPaths <<- R_libPaths
    clusterExport(cl, varlist = c("R_libPaths"))
    clusterEvalQ(cl = cl, .libPaths(R_libPaths))
}
clusterExport(cl, varlist = c(
    "sample_IDs", "folder_workplace", "R_data"
), envir = environment())
clusterEvalQ(cl, library(mobster))
df_all_mobsters <- pblapply(cl = cl, X = 1:length(sample_IDs), FUN = function(n_sample) {
    sample <- sample_IDs[n_sample]
    #---Input the SFS data
    filename <- paste0(R_data, "/", sample, "_1_1.txt")
    mutation_table <- read.csv(filename, sep = "\t", header = TRUE)
    mutation_table$Ref_count <- mutation_table$ref_counts
    mutation_table$Alt_count <- mutation_table$alt_counts
    if (any(is.na(mutation_table$Ref_count)) | any(is.na(mutation_table$Alt_count))) mutation_table <- mutation_table[-which(is.na(mutation_table$Ref_count) | is.na(mutation_table$Alt_count)), ]
    if (any(mutation_table$Ref_count == 0) | any(mutation_table$Alt_count == 0)) mutation_table <- mutation_table[-which(mutation_table$Ref_count == 0 | mutation_table$Alt_count == 0), ]
    mutation_table$VAF <- mutation_table$Alt_count / (mutation_table$Alt_count + mutation_table$Ref_count)
    MOBSTER_data <- data.frame(VAF = mutation_table$VAF)
    #---SFS deconvolution with MOBSTER
    MOBSTER_result <- mobster_fit(
        MOBSTER_data,
        # maxIter = 1000,
        parallel = FALSE,
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
    return(MOBSTER_result)
})
stopCluster(cl)
# ================================================================DECODE
for (sample in sample_IDs) {
    #---Input the SFS data
    filename <- paste0(R_data, "/", sample, "_1_1.txt")
    mutation_table <- read.table(filename, sep = "\t", header = TRUE)
    mutation_table$Ref_count <- mutation_table$ref_counts
    mutation_table$Alt_count <- mutation_table$alt_counts
    if (any(is.na(mutation_table$Ref_count)) | any(is.na(mutation_table$Alt_count))) mutation_table <- mutation_table[-which(is.na(mutation_table$Ref_count) | is.na(mutation_table$Alt_count)), ]
    if (any(mutation_table$Ref_count == 0) | any(mutation_table$Alt_count == 0)) mutation_table <- mutation_table[-which(mutation_table$Ref_count == 0 | mutation_table$Alt_count == 0), ]
    #---SFS deconvolution with DECODE
    DECODE_result <- DECODE(
        sample_id = sample,
        mutation_table = mutation_table
    )
    save(DECODE_result, file = paste0(folder_workplace, "DECODE_", sample, ".rda"))
    #---Plot DECODE deconvolution
    load(paste0(folder_workplace, "DECODE_", sample, ".rda"))
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
