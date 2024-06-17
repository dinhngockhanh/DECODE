# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
R_data <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/data/PCAWG"
R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes"
R_libPaths <- ""
R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/R"
R_libPaths_binomial_table <- ""
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)
library(dplyr)
library(VariantAnnotation)
library(GenomicRanges)
library(R.utils)
library(readxl)
library(data.table)
library(stringi)
library(tidyverse)
library(R.utils)
library(parallel)
library(pbapply)
library(mobster)
library(R.matlab)
library(ggplot2)

setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)

folder_workplace <- "ICGC-MOBSTER"
# ==========================================LIMIT TO SPECIFIC CASE STUDY
# histology <- c("Breast-AdenoCA", "Ovary-AdenoCA", "Myeloid-AML", "Myeloid-MPN", "Myeloid-MDS")
sample_df <- read.csv(paste0(R_data, "/sample_information.csv"))
sample_df <- sample_df[
    which(
        # sample_df$histology_abbreviation %in% histology &
        sample_df$wgd_status == "no_wgd" &
            sample_df$wgd_uncertain == FALSE
    ),
]
write.csv(sample_df, file = paste0(R_workplace, "/sample_information.csv"), row.names = FALSE)
# ===============================================================MOBSTER

sample_df <- sample_df[1:8, ]
print(sample_df$aliquot_id)

numCores <- detectCores()
cl <- makePSOCKcluster(numCores - 1)
if (is.null(R_libPaths) == FALSE) {
    R_libPaths <<- R_libPaths
    clusterExport(cl, varlist = c("R_libPaths"))
    clusterEvalQ(cl = cl, .libPaths(R_libPaths))
}
clusterExport(cl, varlist = c(
    "folder_workplace",
    "R_data",
    "sample_df"
))
clusterEvalQ(cl, library(mobster))
df_all_mobsters <- pblapply(cl = cl, X = 1:length(sample_df$aliquot_id), FUN = function(n_sample) {
    sample <- sample_df$aliquot_id[n_sample]
    #   Import mutational data
    filename_2 <- paste0(R_data, "/", sample, "_1_1.csv")
    mutation_table <- read.csv(filename_2, sep = "\t", header = TRUE)
    mutation_table$Ref_count <- mutation_table$t_ref_count
    mutation_table$Alt_count <- mutation_table$t_alt_count
    if (any(is.na(mutation_table$Ref_count)) | any(is.na(mutation_table$Alt_count))) {
        mutation_table <- mutation_table[-which(is.na(mutation_table$Ref_count) | is.na(mutation_table$Alt_count)), ]
    }
    if (any(mutation_table$Ref_count == 0) | any(mutation_table$Alt_count == 0)) {
        mutation_table <- mutation_table[-which(mutation_table$Ref_count == 0 | mutation_table$Alt_count == 0), ]
    }
    #   Data transformation
    mutation_table$VAF <- mutation_table$Alt_count / (mutation_table$Alt_count + mutation_table$Ref_count)
    mob_data <- as.data.frame(mutation_table$VAF)
    colnames(mob_data)[1] <- "VAF"
    #   SFS deconvolution with MOBSTER
    MOBSTER_result <- mobster_fit(
        mob_data,
        maxIter = 1000,
        parallel = FALSE, # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        description = sample
    )
    #   Find best MOBSTER model
    png(paste0(folder_workplace, "/MOBSTER_", sample, "_1_1.png"), res = 150, width = 15, height = 7.5, units = "in")
    print(plot(MOBSTER_result$best))
    dev.off()
    png(paste0(folder_workplace, "/MOBSTER_model_selection_", sample, "_1_1.png"), res = 150, width = 15, height = 7.5, units = "in")
    print(plot_model_selection(MOBSTER_result))
    dev.off()
    return(MOBSTER_result)
})
stopCluster(cl)
mobster_df <- data.frame()
mobster_fits <- list()
for (n_sample in 1:length(sample_df$aliquot_id)) {
    sample <- sample_df$aliquot_id[n_sample]
    MOBSTER_result <- df_all_mobsters[[n_sample]]
    #---Save MOBSTER results
    mobster_fits[[n_simulation]] <- MOBSTER_result
    mobster_model <- MOBSTER_result$best
    mobster_df[n_simulation, "Sample"] <- sample
    mobster_df[n_simulation, "Mutation_count_in_fitting"] <- mobster_model$N
    mobster_df[n_simulation, "Tail"] <- mobster_model$fit.tail
    mobster_df[n_simulation, "Tail_power"] <- mobster_model$shape + 1
    mobster_df[n_simulation, "Tail_pareto_shape"] <- mobster_model$shape
    mobster_df[n_simulation, "Tail_pareto_scale"] <- mobster_model$scale
    mobster_df[n_simulation, "Tail_mutcount_observed"] <- mobster_model$N.k[[1]]
    mobster_df[n_simulation, "Cluster_count"] <- mobster_model$Kbeta
    for (k in 1:mobster_model$Kbeta) {
        mobster_df[n_simulation, paste0("Cluster_mutcount_observed_", k)] <- mobster_model$N.k[[k + 1]]
        mobster_df[n_simulation, paste0("Cluster_frequency_", k)] <- mobster_model$a[[k]] / (mobster_model$a[[k]] + mobster_model$b[[k]])
        mobster_df[n_simulation, paste0("Cluster_beta_a_", k)] <- mobster_model$a[[k]]
        mobster_df[n_simulation, paste0("Cluster_beta_b_", k)] <- mobster_model$b[[k]]
    }
}
write.csv(mobster_df, paste0("Parameters_ICGC_MOBSTER"), row.names = FALSE)
save(mobster_fits, file = paste0("ICGC_MOBSTER.rda"))
# ================================================================DECODE
#   ...
#   ...
#   ...
#   ...
#   ...
#   ...
#   ...
# ==============================================================ANALYSIS
sample_information_df <- read.csv(paste0(R_data, "/sample_information.csv"))
mobster_df <- read.csv("Parameters_ICGC_MOBSTER.csv")
analysis_ICGC(
    sample_information_df = sample_information_df,
    mobster_df = mobster_df,
    # decode_df = decode_df,
    text_notation = FALSE
)
