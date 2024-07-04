# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
R_PCAWG <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/data/PCAWG"
R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes"
R_libPaths <- ""
R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/R"
R_libPaths_binomial_table <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/MK-Cod.Analysis of the SFS/Core_function_for_SFS_fitting/Binomial_tables"
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

setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)

folder_workplace <- "TEST/"
n_simulations <- 100 # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
n_sample <- 100000 # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
# # ==========================================MAKE CINNER LITE SIMULATIONS
# dir.create(folder_workplace)
# load(file = paste0(R_PCAWG, "/ICGC_purity_coverage.rda"))
# # #---------------------------------------------------Set model parameters
n_selective_clones <- 0
vec_hierarchy_s_mut <- c()
t_diagnosis <- 60 # [in years]
cell_lifespan <- 10 # [in days]
# numCores <- detectCores()
# cl <- makePSOCKcluster(numCores - 1)
# if (is.null(R_libPaths) == FALSE) {
#     R_libPaths <<- R_libPaths
#     clusterExport(cl, varlist = c("R_libPaths"))
#     clusterEvalQ(cl = cl, .libPaths(R_libPaths))
# }
# clusterExport(cl, varlist = c(
#     "simulator_batch",
#     "simulator_full",
#     "simulator_one_simulation",
#     "simulation_clonal_evolution",
#     "simulation_sample_phylogeny",
#     "simulation_sequening_truth",
#     "simulation_sequencing_sc",
#     "simulation_sequencing_bulk",
#     "rand_coverage",
#     "folder_workplace",
#     "n_selective_clones",
#     "vec_hierarchy_s_mut",
#     "t_diagnosis",
#     "cell_lifespan",
#     "n_sample",
#     "ICGC_purity_coverage"
# ))
# clusterEvalQ(cl, library(data.table))
# df_all_simulation_parameters <- pblapply(cl = cl, X = 1:n_simulations, FUN = function(n_simulation) {
#     ####################################################################
#     ####################################################################
#     ####################################################################
#     #---Set simulation parameters
#     t_MRCA <- runif(1, min = 1, max = 10) # [in years]
#     ratio_theta <- runif(1, min = 1, max = 1) # <<<<<<<<<<<<<<<<<<<<<<<<
#     PCAWG_index <- sample(1:ICGC_purity_coverage$N_sample, 1)
#     purity <- ICGC_purity_coverage[[paste0("sample_", PCAWG_index)]]$purity
#     coverage <- ICGC_purity_coverage[[paste0("sample_", PCAWG_index)]]$coverage
#     mean_coverage <- sum(coverage$Read_count * coverage$Frequency) / sum(coverage$Frequency)
#     ####################################################################
#     ####################################################################
#     ####################################################################
#     t_end_time <- round(t_MRCA * 365) # [in days]
#     vec_time_points_s_mut <- t_end_time * c()
#     expected_end_population <- 10^6 # [cells at t_end_time]
#     vec_expected_percent_select <- (1 / (n_selective_clones + 1)) * rep(1, length = (n_selective_clones + 1))
#     range_population <- c(0.8, 1.2) * expected_end_population
#     range_clonal_perc <- c(20, 100)
#     ploidy <- 2
#     choice_theta <- "constant"
#     theta_normal <- 6.12 # [mutations per division]
#     vec_theta_parameters <- c(ratio_theta * theta_normal)
#     vec_theta_mean <- vec_theta_parameters
#     # bulk_coverage_model <- "binomial"
#     # bulk_coverage_variables$mean_coverage <- mean_coverage
#     bulk_coverage_model <- "custom"
#     bulk_coverage_variables <- list()
#     bulk_coverage_variables$coverage <- coverage
#     bulk_min_alt_readcounts <- 3
#     #---Make one CINner simulation
#     folder_workplace_sub <- paste0(folder_workplace, "_", n_simulation, "/")
#     dir.create(folder_workplace_sub)
#     truncal_mutations <- rpois(n = 1, lambda = theta_normal * round((t_diagnosis - t_MRCA) * 365 / cell_lifespan))
#     table_parameters <- simulator_batch(
#         n_simulations = 1,
#         t_end_time = t_end_time,
#         cell_lifespan = cell_lifespan,
#         n_selective_clones = n_selective_clones,
#         vec_time_points_s_mut = vec_time_points_s_mut,
#         vec_hierarchy_s_mut = vec_hierarchy_s_mut,
#         expected_end_population = expected_end_population,
#         vec_expected_percent_select = vec_expected_percent_select,
#         n_sample = n_sample,
#         range_population = range_population,
#         range_clonal_perc = range_clonal_perc,
#         ploidy = ploidy,
#         purity = purity,
#         truncal_mutations = truncal_mutations,
#         choice_theta = choice_theta,
#         vec_theta_parameters = vec_theta_parameters,
#         vec_theta_mean = vec_theta_mean,
#         save_rda = TRUE, # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#         save_true_mutation_table = FALSE,
#         output_bulk = TRUE,
#         output_sc = FALSE,
#         compute_parallel = FALSE,
#         bulk_coverage_model = bulk_coverage_model,
#         bulk_coverage_variables = bulk_coverage_variables,
#         bulk_min_alt_readcounts = bulk_min_alt_readcounts,
#         subfolder = folder_workplace_sub,
#         file_prefix = "",
#         R_libPaths = R_libPaths
#     )
#     #---Plot simulated cell phylogeny
#     load(paste0(folder_workplace_sub, "simulation_1.rda"))
#     simulation$cell_phylogeny_hclust$height <- simulation$cell_phylogeny_hclust$height / 365
#     png(paste0(folder_workplace, "CINner_phylogeny_", n_simulation, ".png"), res = 150, width = 30, height = 15, units = "in")
#     plot(simulation$cell_phylogeny_hclust, labels = FALSE, ylab = "year")
#     dev.off()
#     #---Save simulation parameters
#     table_parameters$Batch_ID <- NULL
#     table_parameters <- cbind(data.frame(Simulation = n_simulation), table_parameters)
#     table_parameters[["Age at diagnosis (years)"]] <- t_diagnosis
#     table_parameters[["Age of MRCA (years)"]] <- t_MRCA
#     table_parameters[["Truncal mutation count"]] <- truncal_mutations
#     table_parameters[["Cell lifespan (days)"]] <- cell_lifespan
#     table_parameters[["Normal (pre-MRCA) mutation rate"]] <- theta_normal
#     table_parameters[["Tumor mutation rate ratio"]] <- ratio_theta
#     table_parameters[["Tumor mutation rate"]] <- vec_theta_parameters
#     table_parameters[["Purity"]] <- purity
#     table_parameters[["Sequencing coverage"]] <- mean_coverage
#     table_parameters[["Expected tumor cell count"]] <- expected_end_population
#     table_parameters[["Sample cell count"]] <- n_sample
#     table_parameters[["Ploidy"]] <- ploidy
#     table_parameters[["Coverage model"]] <- bulk_coverage_model
#     table_parameters[["Limit alt count"]] <- bulk_min_alt_readcounts
#     return(table_parameters)
# })
# stopCluster(cl)
# df_all_simulation_parameters <- do.call(rbind, df_all_simulation_parameters)
# write.csv(df_all_simulation_parameters, file = "Parameters_simulation.csv", row.names = FALSE)
# #--------------------------------------------------------Clean bulk data
# for (n_simulation in 1:n_simulations) {
#     filename <- paste0(folder_workplace, "_", n_simulation, "/1_mutational_data_BULK.csv")
#     mut_table <- read.csv(filename)
#     vec_delete <- which(mut_table$Alt_count == 0 | mut_table$Ref_count == 0)
#     if (length(vec_delete) > 0) mut_table <- mut_table[-vec_delete, ]
#     filename <- paste0(folder_workplace, "_", n_simulation, "/SFS_1.txt")
#     write.table(mut_table, filename, sep = " ", row.names = FALSE, col.names = FALSE)
# }
# ========================================GROUND TRUTH FOR SFS VARIABLES
df_all_simulation_parameters <- read.csv("Parameters_simulation.csv", header = TRUE)
df <- data.frame()
for (n_simulation in 1:n_simulations) {
    #   Retrieve simulation input parameters
    purity <- df_all_simulation_parameters[["Purity"]][n_simulation]
    coverage <- df_all_simulation_parameters[["Sequencing.coverage"]][n_simulation]
    t_end_time <- df_all_simulation_parameters[["Age.of.MRCA..years."]][n_simulation] * 365
    vec_theta_parameters <- df_all_simulation_parameters[["Tumor.mutation.rate"]][n_simulation]
    n_sample <- df_all_simulation_parameters[["Sample.cell.count"]][n_simulation]
    ploidy <- df_all_simulation_parameters[["Ploidy"]][n_simulation]
    truncal_mutations <- df_all_simulation_parameters[["Truncal.mutation.count"]][n_simulation]
    #   Retrieve simulation stochastic parameters
    simulation_variables <- read.csv(paste0(folder_workplace, "_", n_simulation, "/1_simulation_variables.csv"))
    Ns <- simulation_variables$Count_in_population
    ns <- simulation_variables$Count_in_sample
    MRCA_ages <- simulation_variables$MRCA_ages
    #   Retrieve observed mutational data
    observed_mutation_table <- read.csv(paste0(folder_workplace, "_", n_simulation, "/1_mutational_data_BULK.csv"))
    #   Find observed clonal sizes in sample, including subclones
    ns_combined <- ns
    if (length(vec_hierarchy_s_mut) > 0) {
        for (i in length(vec_hierarchy_s_mut):1) {
            ns_combined[vec_hierarchy_s_mut[i] + 1] <- ns_combined[vec_hierarchy_s_mut[i] + 1] + ns_combined[i + 1]
        }
    }
    #   Find total number of neutral mutations
    folder_workplace_sub <- paste0(folder_workplace, "_", n_simulation, "/")
    load(paste0(folder_workplace_sub, "simulation_1.rda"))
    A_total <- 0
    for (i in 1:length(simulation$sample_mutational_table_truth_node_markers)) {
        if (grepl("Foreground_", simulation$sample_mutational_table_truth_node_markers[[i]])) {
            A_total <- A_total + simulation$sample_mutational_table_truth_node_mutation_counts[[i]]
        }
    }
    print("-----------------------------------------------------------")
    print(n_simulation)
    print(A_total)
    #   Find observed number of neutral mutations
    A_observed_mobster <- length(which(grepl("Foreground_", observed_mutation_table$Marker)))
    A_observed_decode <- length(which(grepl("Foreground_", observed_mutation_table$Marker) &
        observed_mutation_table$Alt_count >= min_variant_read &
        (observed_mutation_table$Ref_count + observed_mutation_table$Alt_count) >= min_total_read))
    #   Find expected power of neutral mutations
    alpha <- 2
    #   Find expected binomial hump locations
    ps <- purity * ns_combined / n_sample / ploidy
    #   Find expected number of mutations in each binomial hump
    Ks_expected <- vec_theta_parameters * MRCA_ages
    if (length(vec_hierarchy_s_mut) > 0) {
        for (i in length(vec_hierarchy_s_mut):1) {
            Ks_expected[i + 1] <- Ks_expected[i + 1] - Ks_expected[vec_hierarchy_s_mut[i] + 1]
        }
    }
    Ks_expected[1] <- Ks_expected[1] + truncal_mutations
    #   Save the results
    df <- rbind(df, c(n_simulation, purity, coverage, A_total, A_observed_mobster, A_observed_decode, alpha, ps, Ks_expected))
}
names(df) <- c("Simulation", "Purity", "Coverage", "A_total", "A_observed_mobster", "A_observed_decode", "alpha", paste0("p_", 1:(n_selective_clones + 1)), paste0("K_expected_", 1:(n_selective_clones + 1)))
write.csv(df, paste0("Parameters_true.csv"), row.names = FALSE)
# # ===============================================================MOBSTER
# numCores <- detectCores()
# cl <- makePSOCKcluster(numCores - 1)
# if (is.null(R_libPaths) == FALSE) {
#     R_libPaths <<- R_libPaths
#     clusterExport(cl, varlist = c("R_libPaths"))
#     clusterEvalQ(cl = cl, .libPaths(R_libPaths))
# }
# clusterExport(cl, varlist = c(
#     "folder_workplace"
# ))
# clusterEvalQ(cl, library(mobster))
# df_all_mobsters <- pblapply(cl = cl, X = 1:n_simulations, FUN = function(n_simulation) {
#     #---Input the SFS data
#     filename <- paste0(folder_workplace, "_", n_simulation, "/SFS_1.txt")
#     txtdata <- read.table(file = filename, header = FALSE)
#     data <- transform(txtdata, VAF = txtdata[, 2] / (txtdata[, 1] + txtdata[, 2]))
#     last_col <- ncol(data)
#     mob_data <- as.data.frame(data[, last_col])
#     colnames(mob_data)[1] <- "VAF"
#     #---SFS deconvolution with MOBSTER
#     MOBSTER_result <- try(
#         {
#             mobster_fit(
#                 mob_data,
#                 tail = c(TRUE), # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#                 parallel = FALSE, # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#                 description = paste0("Simulation ", n_simulation)
#             )
#         },
#         silent = TRUE
#     )
#     save(MOBSTER_result, file = paste0(folder_workplace, "MOBSTER_", n_simulation, ".rda"))
#     #---Plot MOBSTER deconvolution
#     if (!inherits(MOBSTER_result, "try-error")) {
#         png(paste0(folder_workplace, "MOBSTER_", n_simulation, ".png"), res = 150, width = 15, height = 7.5, units = "in")
#         print(plot(MOBSTER_result$best))
#         dev.off()
#         png(paste0(folder_workplace, "MOBSTER_model_selection_", n_simulation, ".png"), res = 150, width = 15, height = 7.5, units = "in")
#         print(plot_model_selection(MOBSTER_result))
#         dev.off()
#     }
#     return(MOBSTER_result)
# })
# stopCluster(cl)
# mobster_df <- data.frame()
# mobster_fits <- list()
# for (n_simulation in 1:n_simulations) {
#     MOBSTER_result <- df_all_mobsters[[n_simulation]]
#     #---Save MOBSTER results
#     mobster_fits[[n_simulation]] <- MOBSTER_result
#     mobster_df[n_simulation, "Simulation"] <- n_simulation
#     if (inherits(MOBSTER_result, "try-error")) {
#         mobster_df[n_simulation, "Succeed"] <- FALSE
#     } else {
#         mobster_model <- MOBSTER_result$best
#         mobster_df[n_simulation, "Succeed"] <- TRUE
#         mobster_df[n_simulation, "Mutation_count_in_fitting"] <- mobster_model$N
#         mobster_df[n_simulation, "Tail"] <- mobster_model$fit.tail
#         mobster_df[n_simulation, "Tail_power"] <- mobster_model$shape + 1
#         mobster_df[n_simulation, "Tail_pareto_shape"] <- mobster_model$shape
#         mobster_df[n_simulation, "Tail_pareto_scale"] <- mobster_model$scale
#         mobster_df[n_simulation, "Tail_mutcount_observed"] <- mobster_model$N.k[[1]]
#         mobster_df[n_simulation, "Cluster_count"] <- mobster_model$Kbeta
#         for (k in 1:mobster_model$Kbeta) {
#             mobster_df[n_simulation, paste0("Cluster_mutcount_observed_", k)] <- mobster_model$N.k[[k + 1]]
#             mobster_df[n_simulation, paste0("Cluster_frequency_", k)] <- mobster_model$a[[k]] / (mobster_model$a[[k]] + mobster_model$b[[k]])
#             mobster_df[n_simulation, paste0("Cluster_beta_a_", k)] <- mobster_model$a[[k]]
#             mobster_df[n_simulation, paste0("Cluster_beta_b_", k)] <- mobster_model$b[[k]]
#         }
#     }
# }
# write.csv(mobster_df, paste0("Parameters_MOBSTER.csv"), row.names = FALSE)
# save(mobster_fits, file = paste0("MOBSTER.rda"))
# ================================================================DECODE
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
# for (n_simulation in 1:n_simulations) {
#     cat("\n==========================================================================================================================\n") # nolint
#     cat(paste0("DECODE FOR SIMULATION ", n_simulation, "...\n"))
#     #---Input the SFS data
#     filename_2 <- paste0(folder_workplace, "_", n_simulation, "/SFS_1.txt")
#     mutation_table <- read.table(filename_2, sep = " ", header = FALSE)
#     colnames(mutation_table) <- c("Ref_count", "Alt_count", "Marker")
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
#     save(DECODE_result, file = paste0(folder_workplace, "DECODE_", n_simulation, ".rda"))
#     #---Plot DECODE deconvolution
#     load(paste0(folder_workplace, "DECODE_", n_simulation, ".rda"))
#     png(paste0(folder_workplace, "DECODE_", n_simulation, ".png"), res = 150, width = 15, height = 7.5, units = "in")
#     print(DECODE_plot(DECODE_result = DECODE_result, data_marker_colors = data_marker_colors))
#     dev.off()
#     # png(paste0(folder_workplace, "DECODE_model_selection_", n_simulation, ".png"), res = 150, width = 15, height = 7.5, units = "in")
#     # print(DECODE_plot_model_selection(DECODE_result = DECODE_result, data_marker_colors = data_marker_colors))
#     # dev.off()
#     #---Save DECODE results
#     decode_fits[[n_simulation]] <- DECODE_result
#     decode_model <- DECODE_result$best_result$parameters_df
#     decode_df[n_simulation, "Simulation"] <- n_simulation
#     decode_df[n_simulation, "Mutation_count_in_fitting"] <- decode_model$Mutation_count_for_fitting
#     decode_df[n_simulation, "Tail"] <- decode_model$Tail
#     decode_df[n_simulation, "Tail_sensitivity_Bayesian_pi0_std"] <- decode_model$Tail_sensitivity_Bayesian_pi0_std
#     decode_df[n_simulation, "Tail_sensitivity_Bayesian_alpha_std"] <- decode_model$Tail_sensitivity_Bayesian_alpha_std
#     decode_df[n_simulation, "Tail_sensitivity_Morris_pi0_mean"] <- decode_model$Tail_sensitivity_Morris_pi0_mean
#     decode_df[n_simulation, "Tail_sensitivity_Morris_pi0_mean_abs"] <- decode_model$Tail_sensitivity_Morris_pi0_mean_abs
#     decode_df[n_simulation, "Tail_sensitivity_Morris_pi0_std"] <- decode_model$Tail_sensitivity_Morris_pi0_std
#     decode_df[n_simulation, "Tail_sensitivity_Morris_alpha_mean"] <- decode_model$Tail_sensitivity_Morris_alpha_mean
#     decode_df[n_simulation, "Tail_sensitivity_Morris_alpha_mean_abs"] <- decode_model$Tail_sensitivity_Morris_alpha_mean_abs
#     decode_df[n_simulation, "Tail_sensitivity_Morris_alpha_std"] <- decode_model$Tail_sensitivity_Morris_alpha_std
#     decode_df[n_simulation, "Tail_power"] <- decode_model$Tail_power
#     decode_df[n_simulation, "Tail_mutcount_observed"] <- decode_model$Tail_mutcount_observed
#     decode_df[n_simulation, "Tail_mutcount_predicted"] <- decode_model$Tail_mutcount_predicted
#     decode_df[n_simulation, "Cluster_count"] <- decode_model$Cluster_count
#     if (decode_model$Cluster_count >= 1) {
#         for (k in 1:decode_model$Cluster_count) {
#             decode_df[n_sample, paste0("Cluster_mutcount_observed_", k)] <- decode_model[[paste0("Cluster_mutcount_observed_", k)]]
#             decode_df[n_sample, paste0("Cluster_mutcount_predicted_", k)] <- decode_model[[paste0("Cluster_mutcount_predicted_", k)]]
#             decode_df[n_sample, paste0("Cluster_frequency_", k)] <- decode_model[[paste0("Cluster_frequency_", k)]]
#         }
#     }
# }
# write.csv(decode_df, paste0("Parameters_DECODE.csv"), row.names = FALSE)
# save(decode_fits, file = paste0("DECODE.rda"))
# ==============================================================ANALYSIS
groundtruth_df <- read.csv("Parameters_true.csv")
mobster_df <- read.csv("Parameters_MOBSTER.csv")
decode_df <- read.csv("Parameters_DECODE.csv")
analysis_synthetic_test(
    groundtruth_df = groundtruth_df,
    mobster_df = mobster_df,
    decode_df = decode_df,
    text_notation = TRUE,
    cluster_count = n_selective_clones + 1,
    tail = TRUE
)
