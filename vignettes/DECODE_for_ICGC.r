# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
R_data <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/data/PCAWG"
R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes"
R_libPaths <- ""
R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/R"
R_libPaths_binomial_table <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/MK-Cod.Analysis of the SFS/Core_function_for_SFS_fitting/Binomial_tables"
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
# =====================================================DECODE parameters
#---Set DECODE parameters
# 	Total number of sampled cells in binomial table construction
matrix_binomial_sample_size <- 1000
# 	Minimum and maximum number of total reads
r_min <- 0
r_max <- 500
# 	Minimum variant read count to be accepted
min_variant_read <- 5
# 	Minimum total read count to be accepted
min_total_read <- 0
# 	Number of steps to divide SFS frequencies in [0,1]
matrix_binomial_sfs_stepcount <- 100
# 	Choice of ploidy, which changes the binomial rate
matrix_binomial_ploidy <- 2
#---Options for DECODE plotting
data_marker_colors <- c(
    "Data" = "black",
    "Foreground 0" = rgb(0.2, 0.2, 0.2),
    "Foreground 1" = rgb(0.5, 0.5, 0.5),
    "Foreground 2" = rgb(0.7, 0.7, 0.7),
    "Background 1&2" = rgb(0.9290, 0.6940, 0.1250),
    "Background 1" = rgb(0.6350, 0.0780, 0.1840),
    "Background 2" = rgb(0.4660, 0.6740, 0.1880),
    "Truncal" = rgb(0, 0.4470, 0.7410)
)
# ===========================================GET ICGC SAMPLE INFORMATION
sample_df <- read.csv(paste0(R_data, "/sample_information.csv"))
sample_df <- sample_df[
    which(
        sample_df$wgd_status == "no_wgd" &
            sample_df$wgd_uncertain == FALSE
    ),
]
write.csv(sample_df, file = paste0(R_workplace, "/sample_information.csv"), row.names = FALSE)
# # ===============================================================MOBSTER
# folder_workplace <- "ICGC-MOBSTER/"
# dir.create(folder_workplace)
# numCores <- detectCores()
# cl <- makePSOCKcluster(numCores - 1)
# if (is.null(R_libPaths) == FALSE) {
#     R_libPaths <<- R_libPaths
#     clusterExport(cl, varlist = c("R_libPaths"))
#     clusterEvalQ(cl = cl, .libPaths(R_libPaths))
# }
# clusterExport(cl, varlist = c(
#     "folder_workplace",
#     "R_data",
#     "sample_df"
# ))
# clusterEvalQ(cl, library(mobster))
# df_all_mobsters <- pblapply(cl = cl, X = 1:length(sample_df$aliquot_id), FUN = function(n_sample) {
#     sample <- sample_df$aliquot_id[n_sample]
#     #   Import mutational data
#     filename_2 <- paste0(R_data, "/", sample, "_1_1.csv")
#     mutation_table <- read.csv(filename_2, sep = "\t", header = TRUE)
#     mutation_table$Ref_count <- mutation_table$t_ref_count
#     mutation_table$Alt_count <- mutation_table$t_alt_count
#     if (any(is.na(mutation_table$Ref_count)) | any(is.na(mutation_table$Alt_count))) {
#         mutation_table <- mutation_table[-which(is.na(mutation_table$Ref_count) | is.na(mutation_table$Alt_count)), ]
#     }
#     if (any(mutation_table$Ref_count == 0) | any(mutation_table$Alt_count == 0)) {
#         mutation_table <- mutation_table[-which(mutation_table$Ref_count == 0 | mutation_table$Alt_count == 0), ]
#     }
#     #   Data transformation
#     mutation_table$VAF <- mutation_table$Alt_count / (mutation_table$Alt_count + mutation_table$Ref_count)
#     mob_data <- as.data.frame(mutation_table$VAF)
#     colnames(mob_data)[1] <- "VAF"
#     #   SFS deconvolution with MOBSTER
#     MOBSTER_result <- mobster_fit(
#         mob_data,
#         maxIter = 1000, # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#         parallel = FALSE,
#         description = sample
#     )
#     #   Find best MOBSTER model
#     png(paste0(folder_workplace, "MOBSTER_", sample, "_1_1.png"), res = 150, width = 15, height = 7.5, units = "in")
#     print(plot(MOBSTER_result$best))
#     dev.off()
#     png(paste0(folder_workplace, "MOBSTER_model_selection_", sample, "_1_1.png"), res = 150, width = 15, height = 7.5, units = "in")
#     print(plot_model_selection(MOBSTER_result))
#     dev.off()
#     return(MOBSTER_result)
# })
# stopCluster(cl)
# mobster_df <- data.frame()
# mobster_fits <- list()
# for (n_sample in 1:length(sample_df$aliquot_id)) {
#     sample <- sample_df$aliquot_id[n_sample]
#     MOBSTER_result <- df_all_mobsters[[n_sample]]
#     #---Save MOBSTER results
#     mobster_fits[[n_simulation]] <- MOBSTER_result
#     mobster_model <- MOBSTER_result$best
#     mobster_df[n_simulation, "Sample"] <- sample
#     mobster_df[n_simulation, "Mutation_count_in_fitting"] <- mobster_model$N
#     mobster_df[n_simulation, "Tail"] <- mobster_model$fit.tail
#     mobster_df[n_simulation, "Tail_power"] <- mobster_model$shape + 1
#     mobster_df[n_simulation, "Tail_pareto_shape"] <- mobster_model$shape
#     mobster_df[n_simulation, "Tail_pareto_scale"] <- mobster_model$scale
#     mobster_df[n_simulation, "Tail_mutcount_observed"] <- mobster_model$N.k[[1]]
#     mobster_df[n_simulation, "Cluster_count"] <- mobster_model$Kbeta
#     for (k in 1:mobster_model$Kbeta) {
#         mobster_df[n_simulation, paste0("Cluster_mutcount_observed_", k)] <- mobster_model$N.k[[k + 1]]
#         mobster_df[n_simulation, paste0("Cluster_frequency_", k)] <- mobster_model$a[[k]] / (mobster_model$a[[k]] + mobster_model$b[[k]])
#         mobster_df[n_simulation, paste0("Cluster_beta_a_", k)] <- mobster_model$a[[k]]
#         mobster_df[n_simulation, paste0("Cluster_beta_b_", k)] <- mobster_model$b[[k]]
#     }
# }
# write.csv(mobster_df, paste0("Parameters_ICGC_MOBSTER"), row.names = FALSE)
# save(mobster_fits, file = paste0("ICGC_MOBSTER.rda"))
# # ================================================================DECODE
# folder_workplace <- "ICGC-DECODE/"
# dir.create(folder_workplace)
# #---Input binomial table
# cat("\n==========================================================================================================================\n")
# cat(paste0("LOAD THE BINOMIAL TABLE...\n"))
# filename_1 <- paste0(
#     R_libPaths_binomial_table, "/Binomial_PDF_",
#     matrix_binomial_sample_size, "_",
#     r_max, "_",
#     min_variant_read, "_",
#     min_total_read, "_",
#     matrix_binomial_sfs_stepcount, "_",
#     matrix_binomial_ploidy, ".mat"
# )
# inputBinomialMatrix <- readMat(filename_1)
# matrix_binomial_PDF <- inputBinomialMatrix$matrix.binomial.PDF
# #---Deconvolution for each SFS
# decode_df <- data.frame()
# decode_fits <- list()
# # for (n_sample in 1:length(sample_df$aliquot_id)) {
# for (n_sample in 490:491) {
#     sample <- sample_df$aliquot_id[n_sample]
#     cat("\n==========================================================================================================================\n") # nolint
#     cat(paste0("DECODE FOR SAMPLE ", sample, "...\n"))
#     #---Input the SFS data
#     filename_2 <- paste0(R_data, "/", sample, "_1_1.csv")
#     mutation_table <- read.csv(filename_2, sep = "\t", header = TRUE)
#     mutation_table$Ref_count <- mutation_table$t_ref_count
#     mutation_table$Alt_count <- mutation_table$t_alt_count
#     if (any(is.na(mutation_table$Ref_count)) | any(is.na(mutation_table$Alt_count))) {
#         mutation_table <- mutation_table[-which(is.na(mutation_table$Ref_count) | is.na(mutation_table$Alt_count)), ]
#     }
#     if (any(mutation_table$Ref_count == 0) | any(mutation_table$Alt_count == 0)) {
#         mutation_table <- mutation_table[-which(mutation_table$Ref_count == 0 | mutation_table$Alt_count == 0), ]
#     }
#     #---SFS deconvolution with DECODE
#     DECODE_result <- DECODE(
#         mutation_table = mutation_table,
#         criterion = "ICL", # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#         neutral_power_min = 0.5,
#         neutral_power_max = 5,
#         cluster_frequency_min = 0.01,
#         cluster_frequency_max = 1,
#         matrix_binomial_PDF = matrix_binomial_PDF,
#         matrix_binomial_sample_size = matrix_binomial_sample_size,
#         matrix_binomial_sfs_stepcount = matrix_binomial_sfs_stepcount,
#         matrix_binomial_ploidy = matrix_binomial_ploidy,
#         sample_size = n_sample,
#         SFS_totalsteps = 100,
#         r_min = r_min,
#         r_max = r_max,
#         coverage_distribution = "sample-specific",
#         N_trials = 10000,
#         compute_parallel_fit = TRUE,
#         neutral_tail = TRUE # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#     )
#     save(DECODE_result, file = paste0(folder_workplace, sample, "_DECODE.rda"))
#     #---Plot DECODE deconvolution
#     load(paste0(folder_workplace, sample, "_DECODE.rda"))
#     png(paste0(folder_workplace, sample, "_DECODE.png"), res = 150, width = 15, height = 7.5, units = "in")
#     print(DECODE_plot(DECODE_result = DECODE_result, data_marker_colors = data_marker_colors))
#     dev.off()
#     # png(paste0(folder_workplace, sample, "_DECODE_model_selection.png"), res = 150, width = 15, height = 7.5, units = "in")
#     # print(DECODE_plot_model_selection(DECODE_result = DECODE_result, data_marker_colors = data_marker_colors))
#     # dev.off()
#     #---Save DECODE results
#     decode_fits[[n_sample]] <- DECODE_result
#     decode_model <- DECODE_result$best_result$parameters_df
#     decode_df[n_sample, "Sample"] <- sample
#     decode_df[n_sample, "Mutation_count_in_fitting"] <- decode_model$Mutation_count_for_fitting
#     decode_df[n_sample, "Tail"] <- decode_model$Tail
#     decode_df[n_sample, "Tail_sensitivity_Bayesian_pi0_std"] <- decode_model$Tail_sensitivity_Bayesian_pi0_std
#     decode_df[n_sample, "Tail_sensitivity_Bayesian_alpha_std"] <- decode_model$Tail_sensitivity_Bayesian_alpha_std
#     decode_df[n_sample, "Tail_sensitivity_Morris_pi0_mean"] <- decode_model$Tail_sensitivity_Morris_pi0_mean
#     decode_df[n_sample, "Tail_sensitivity_Morris_pi0_mean_abs"] <- decode_model$Tail_sensitivity_Morris_pi0_mean_abs
#     decode_df[n_sample, "Tail_sensitivity_Morris_pi0_std"] <- decode_model$Tail_sensitivity_Morris_pi0_std
#     decode_df[n_sample, "Tail_sensitivity_Morris_alpha_mean"] <- decode_model$Tail_sensitivity_Morris_alpha_mean
#     decode_df[n_sample, "Tail_sensitivity_Morris_alpha_mean_abs"] <- decode_model$Tail_sensitivity_Morris_alpha_mean_abs
#     decode_df[n_sample, "Tail_sensitivity_Morris_alpha_std"] <- decode_model$Tail_sensitivity_Morris_alpha_std
#     decode_df[n_sample, "Tail_power"] <- decode_model$Tail_power
#     decode_df[n_sample, "Tail_mutcount_observed"] <- decode_model$Tail_mutcount_observed
#     decode_df[n_sample, "Tail_mutcount_predicted"] <- decode_model$Tail_mutcount_predicted
#     decode_df[n_sample, "Cluster_count"] <- decode_model$Cluster_count
#     if (decode_model$Cluster_count >= 1) {
#         for (k in 1:decode_model$Cluster_count) {
#             decode_df[n_sample, paste0("Cluster_mutcount_observed_", k)] <- decode_model[[paste0("Cluster_mutcount_observed_", k)]]
#             decode_df[n_sample, paste0("Cluster_mutcount_predicted_", k)] <- decode_model[[paste0("Cluster_mutcount_predicted_", k)]]
#             decode_df[n_sample, paste0("Cluster_frequency_", k)] <- decode_model[[paste0("Cluster_frequency_", k)]]
#         }
#     }
# }
# write.csv(decode_df, paste0("Parameters_ICGC_DECODE.csv"), row.names = FALSE)
# save(decode_fits, file = paste0("DECODE.rda"))
# ==============================================================ANALYSIS
sample_information_df <- read.csv(paste0(R_data, "/sample_information.csv"))
mobster_df <- read.csv("Parameters_ICGC_MOBSTER.csv")
decode_df <- read.csv("Parameters_ICGC_DECODE.csv")
analysis_ICGC(
    sample_information_df = sample_information_df,
    mobster_df = mobster_df,
    decode_df = decode_df,
    text_notation = FALSE
)
