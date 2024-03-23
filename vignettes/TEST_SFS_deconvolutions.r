# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SFS_CNA_deconvolution/vignettes"
R_libPaths <- ""
R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SFS_CNA_deconvolution/R"
R_libPaths_binomial_table <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/MK-Cod.Analysis of the SFS/Core_function_for_SFS_fitting/Binomial_tables"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Xuanwen - Laptop
# R_workplace <- "C:/5398_R/vignettes" # where the vignette folder is located
# R_libPaths <- ""
# R_libPaths_extra <- "" # where the R folder is located
# R_libPaths_binomial_table <- "D:/5398_dataset" # where the binomial tables are located
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Yining - Laptop
# R_workplace <- "C:/Users/Mayin/Documents/1GRADUATE/1. Study/2. 24Spring/5398 Dinh/DATA"
# R_libPaths <- ""
# R_libPaths_extra <- "C:/Users/Mayin/Documents/1GRADUATE/1. Study/2. 24Spring/5398 Dinh/github_clone/SFS_CNA_deconvolution-1/R"
# R_libPaths_binomial_table <- ""
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Yining - Ginsburg
# R_workplace <- "/burg/iicd/users/ym2998/MOBSTER_Test"
# R_libPaths <- "/burg/iicd/users/ym2998/R_Packages"
# R_libPaths_extra <- "/burg/iicd/users/ym2998/Mob_CINner_Function"
# R_libPaths_binomial_table <- "/burg/iicd/users/ym2998/Deconvolution" # where the binomial tables are located
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
# ==========================================MAKE CINNER LITE SIMULATIONS
#---------------------------------------------------Set model parameters
n_simulations <- 1000


t_end_time <- 1000
t_tau_step <- 1
n_selective_clones <- 0 # 1
vec_time_points_s_mut <- t_end_time * c() # c(0.6)
vec_hierarchy_s_mut <- c() # c(0)
expected_end_population <- 10^6
vec_expected_percent_select <- (1 / (n_selective_clones + 1)) * rep(1, length = (n_selective_clones + 1))
n_sample <- 100000
range_population <- c(0.8, 1.2) * expected_end_population
range_clonal_perc <- c(20, 100)
# mindiff_clonal_perc <- 10
ploidy <- 2
truncal_mutations <- 500
choice_theta <- "constant"
vec_theta_parameters <- rep(0.4, length = (n_selective_clones + 1))
vec_theta_mean <- vec_theta_parameters
bulk_coverage_model <- "binomial"
bulk_coverage_variables <- c(0, 100)
bulk_min_alt_readcounts <- 4
#------------------------------------------------Create bulk simulations
dir.create(folder_workplace)
simulator_batch(
    n_simulations = n_simulations,
    t_end_time = t_end_time,
    t_tau_step = t_tau_step,
    n_selective_clones = n_selective_clones,
    vec_time_points_s_mut = vec_time_points_s_mut,
    vec_hierarchy_s_mut = vec_hierarchy_s_mut,
    expected_end_population = expected_end_population,
    vec_expected_percent_select = vec_expected_percent_select,
    n_sample = n_sample,
    range_population = range_population,
    range_clonal_perc = range_clonal_perc,
    # mindiff_clonal_perc = mindiff_clonal_perc,
    ploidy = ploidy,
    truncal_mutations = truncal_mutations,
    choice_theta = choice_theta,
    vec_theta_parameters = vec_theta_parameters,
    vec_theta_mean = vec_theta_mean,
    save_rda = TRUE,
    save_true_mutation_table = FALSE,
    output_bulk = TRUE,
    output_sc = FALSE,
    compute_parallel = TRUE,
    bulk_coverage_model = bulk_coverage_model,
    bulk_coverage_variables = bulk_coverage_variables,
    bulk_min_alt_readcounts = bulk_min_alt_readcounts,
    subfolder = folder_workplace,
    R_libPaths = R_libPaths
)
#--------------------------------------------------------Clean bulk data
for (n_simulation in 1:n_simulations) {
    filename <- paste0(folder_workplace, "ClonalTimes=", vec_time_points_s_mut, "_ClonalHierarchy=", vec_hierarchy_s_mut, "_simulated_SFS_", n_simulation, "_mutational_data_BULK.csv")
    mut_table <- read.csv(filename)
    vec_delete <- which(mut_table$Alt_count == 0 | mut_table$Ref_count == 0)
    if (length(vec_delete) > 0) mut_table <- mut_table[-vec_delete, ]
    filename <- paste0(folder_workplace, "SFS_", n_simulation, ".txt")
    write.table(mut_table, filename, sep = " ", row.names = FALSE, col.names = FALSE)
}
# ========================================GROUND TRUTH FOR SFS VARIABLES
df <- data.frame()
for (n_simulation in 1:n_simulations) {
    #   Retrieve clonal MRCA ages and sizes in population & sample
    simulation_variables <- read.csv(paste0(folder_workplace, "ClonalTimes=", vec_time_points_s_mut, "_ClonalHierarchy=", vec_hierarchy_s_mut, "_simulated_SFS_", n_simulation, "_simulation_variables.csv"))
    Ns <- simulation_variables$Count_in_population
    ns <- simulation_variables$Count_in_sample
    MRCA_ages <- simulation_variables$MRCA_ages
    #   Find clonal sizes in sample, including subclones
    ns_combined <- ns
    if (length(vec_hierarchy_s_mut) > 0) {
        for (i in length(vec_hierarchy_s_mut):1) {
            ns_combined[vec_hierarchy_s_mut[i] + 1] <- ns_combined[vec_hierarchy_s_mut[i] + 1] + ns_combined[i + 1]
        }
    }
    #   Compute actual clonal growth rates
    growth_rates <- log(Ns) / (t_end_time - MRCA_ages)
    #   Find expected number of neutral mutations
    A <- sum(vec_theta_parameters * ns / growth_rates)
    #   Find expected power of neutral mutations
    alpha <- 2
    #   Find expected binomial hump locations
    ps <- ns_combined / n_sample / ploidy
    #   Find expected number of mutations in each binomial hump
    Ks <- vec_theta_parameters * MRCA_ages
    if (length(vec_hierarchy_s_mut) > 0) {
        for (i in length(vec_hierarchy_s_mut):1) {
            Ks[i + 1] <- Ks[i + 1] - Ks[vec_hierarchy_s_mut[i] + 1]
        }
    }
    Ks[1] <- Ks[1] + truncal_mutations
    #   Save the results
    df <- rbind(df, c(n_simulation, A, alpha, ps, Ks))
}
names(df) <- c("Simulation", "A", "alpha", paste0("p_", 1:(n_selective_clones + 1)), paste0("K_", 1:(n_selective_clones + 1)))
write.csv(df, paste0(folder_workplace, "Parameters_true.csv"), row.names = FALSE)
# ===============================================================MOBSTER
mobster_df <- data.frame()
model_list <- list()
for (i in 1:n_simulations) {
    #   Import mutational data
    filename <- paste0(folder_workplace, "SFS_", i, ".txt")
    txtdata <- read.table(file = filename, header = FALSE)
    #   Data transformation
    data <- transform(txtdata, VAF = txtdata[, 2] / (txtdata[, 1] + txtdata[, 2]))
    last_col <- ncol(data)
    mob_data <- as.data.frame(data[, last_col])
    colnames(mob_data)[1] <- "VAF"
    #   SFS deconvolution with MOBSTER
    # mobster:::template_parameters_fast_setup() # show basic setup
    fit <- mobster_fit(
        mob_data,
        auto_setup = "FAST"
    )
    #   Find best MOBSTER model
    mob_model <- fit$best
    model_list[[i]] <- mob_model # save model_list
    png(paste0(folder_workplace, "MOBSTER_", i, ".png"), res = 150, width = 15, height = 7.5, units = "in")
    print(plot(fit$best))
    dev.off()
    #   Save the results
    mobster_df[i, "Simulation"] <- i # id
    mobster_df[i, "Total_N"] <- mob_model$N # total acount
    mobster_df[i, "Tail"] <- mob_model$fit.tail # bool: if tail exists
    mobster_df[i, "Tail_Num"] <- mob_model$N.k[[1]] # number of tail
    mobster_df[i, "Tail_shape"] <- mob_model$shape # shape of tail
    mobster_df[i, "Tail_scale"] <- mob_model$scale # scale of tail
    mobster_df[i, "Kbeta_cluster"] <- mob_model$Kbeta # number of clusters
    for (k in 1:mob_model$Kbeta) {
        mobster_df[i, paste0("cl_num_", k)] <- mob_model$N.k[[k + 1]] # number of Beta
        mobster_df[i, paste0("a_", k)] <- mob_model$a[[k]] # alpha of Beta
        mobster_df[i, paste0("b_", k)] <- mob_model$b[[k]] # beta of Beta
        mobster_df[i, paste0("p_", k)] <- mobster_df[i, paste0("a_", k)] / (mobster_df[i, paste0("a_", k)] + mobster_df[i, paste0("b_", k)])
    }
}
write.csv(mobster_df, paste0(folder_workplace, "Parameters_mobster.csv"), row.names = FALSE)
# =====================================================SFS DECONVOLUTION
#---Set model parameters
# 	Total number of sampled cells in binomial table construction
matrix_binomial_sample_size <- 1000
# 	Minimum and maximum number of reads
r_min <- 0
r_max <- 500
# 	Minimum variant read count to be accepted
min_variant_read <- 5
# 	Minimum total read count to be accepted
min_total_read <- 0
# 	Number of steps to divide SFS frequencies in [0,1]
SFS_totalsteps <- 100
matrix_binomial_sfs_stepcount <- 100
# 	Choice of ploidy, which changes the binomial rate
matrix_binomial_ploidy <- 2
#   Assumption of coverage distribution
option_dist_coverage <- "binomial"
dist_coverage_var_1 <- 100
#---Options for fitting
#   Candidates for neutral tail powers
list_neutral_powers <- seq(1, 3, by = 0.01)
# 	Candidates for where the hump frequencies are
N_SFS_positions <- 100
list_frequencies <- seq(from = 1 / N_SFS_positions, to = 1, by = 1 / N_SFS_positions)
#---Options for plotting
data_marker_colors <- c(
    "Data" = "black",
    "Data: Foreground 0" = rgb(0.2, 0.2, 0.2),
    "Data: Foreground 1" = rgb(0.5, 0.5, 0.5),
    "Data: Foreground 2" = rgb(0.7, 0.7, 0.7),
    "Data: Background 1&2" = rgb(0.9290, 0.6940, 0.1250),
    "Data: Background 1" = rgb(0.6350, 0.0780, 0.1840),
    "Data: Background 2" = rgb(0.4660, 0.6740, 0.1880),
    "Data: Truncal" = rgb(0.3010, 0.7450, 0.9330)
)
#---Input binomial table
cat("\n==========================================================================================================================\n")
cat(paste0("LOAD THE BINOMIAL TABLE...\n"))
filename_1 <- paste0(
    R_libPaths_binomial_table, "/Binomial_PDF_",
    matrix_binomial_sample_size, "_",
    r_max, "_",
    min_variant_read, "_",
    min_total_read, "_",
    matrix_binomial_sfs_stepcount, "_",
    matrix_binomial_ploidy, ".mat"
)
inputBinomialMatrix <- readMat(filename_1)
matrix_binomial_PDF <- inputBinomialMatrix$matrix.binomial.PDF
#---Deconvolution for each SFS
deconvolution_df <- data.frame()
for (n_simulation in 1:n_simulations) {
    cat("\n==========================================================================================================================\n")
    cat(paste0("SFS DECONVOLUTION FOR SIMULATION-", n_simulation, "...\n"))
    #---Input the SFS data
    filename_2 <- paste0(R_workplace, "/", folder_workplace, "SFS_", n_simulation, ".txt")
    mutation_table <- read.table(filename_2, sep = " ", header = FALSE)
    colnames(mutation_table) <- c("Ref_count", "Alt_count", "Marker")
    #---Perform SFS deconvolution
    results <- SFS_deconvolution(
        mutation_table = mutation_table,
        criterion = "BIC",
        list_neutral_powers = list_neutral_powers,
        list_frequencies = list_frequencies,
        matrix_binomial_PDF = matrix_binomial_PDF,
        matrix_binomial_sample_size = matrix_binomial_sample_size,
        matrix_binomial_sfs_stepcount = matrix_binomial_sfs_stepcount,
        matrix_binomial_ploidy = matrix_binomial_ploidy,
        sample_size = n_sample,
        SFS_totalsteps = SFS_totalsteps,
        r_min = r_min,
        r_max = r_max,
        option_dist_coverage = option_dist_coverage,
        dist_coverage_var_1 = dist_coverage_var_1,
        max_trials = 100000,
        compute_parallel = TRUE,
        data_marker_colors = data_marker_colors,
        plot_filename = paste0(R_workplace, "/", folder_workplace, "DECONVOLUTION_", n_simulation, ".png")
    )
    vec_para_best_final <- results$vec_para_best_final

    deconvolution <- results$deconvolution
    if (n_simulation == 1) {
        deconvolution_df <- deconvolution
    } else {
        if (deconvolution$Cluster_count > max(deconvolution_df$Cluster_count)) {
            for (k in (deconvolution_df$Cluster_count + 1):deconvolution$Cluster_count) {
                deconvolution_df[1, paste0("Cluster_frequency_", k)] <- NA
                deconvolution_df[1, paste0("Cluster_mutcount_observed_", k)] <- NA
                deconvolution_df[1, paste0("Cluster_mutcount_predicted_", k)] <- NA
            }
        } else if (deconvolution$Cluster_count < max(deconvolution_df$Cluster_count)) {
            for (k in (deconvolution$Cluster_count + 1):max(deconvolution_df$Cluster_count)) {
                deconvolution[1, paste0("Cluster_frequency_", k)] <- NA
                deconvolution[1, paste0("Cluster_mutcount_observed_", k)] <- NA
                deconvolution[1, paste0("Cluster_mutcount_predicted_", k)] <- NA
            }
        }
        deconvolution_df <- rbind(deconvolution_df, deconvolution)
    }
    #---Store the best fit
    N_humps <- (length(vec_para_best_final) - 1) / 2
    filename <- paste0(R_workplace, "/", folder_workplace, "SFS_deconvolution_parameters_", n_simulation, ".txt")
    fileID <- file(filename, "w")
    writeLines(paste(sprintf("%.3f", vec_para_best_final), collapse = "\t"), fileID)
    close(fileID)
}
deconvolution_df <- cbind(
    data.frame(Simulation = 1:n_simulations),
    deconvolution_df
)
write.csv(deconvolution_df, paste0(folder_workplace, "Parameters_deconvolution.csv"), row.names = FALSE)
# ============================================================COMPARISON
groundtruth_df <- read.csv(paste0(folder_workplace, "Parameters_true.csv"))
mobster_df <- read.csv(paste0(folder_workplace, "Parameters_mobster.csv"))
deconvolution_df <- read.csv(paste0(folder_workplace, "Parameters_deconvolution.csv"))
plot_deconvolution_components(
    groundtruth_df = groundtruth_df,
    mobster_df = mobster_df,
    deconvolution_df = deconvolution_df,
    cluster_count = n_selective_clones + 1,
    folder_workplace = folder_workplace
)
