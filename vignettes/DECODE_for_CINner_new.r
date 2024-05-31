# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
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

# ==========================================MAKE CINNER LITE SIMULATIONS
folder_workplace <- "TEST/"
dir.create(folder_workplace)
#---------------------------------------------------Set model parameters
n_simulations <- 16 # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
n_sample <- 1000 # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
n_selective_clones <- 0
vec_hierarchy_s_mut <- c()

numCores <- detectCores()
cl <- makePSOCKcluster(numCores - 1)
if (is.null(R_libPaths) == FALSE) {
    R_libPaths <<- R_libPaths
    clusterExport(cl, varlist = c("R_libPaths"))
    clusterEvalQ(cl = cl, .libPaths(R_libPaths))
}
clusterExport(cl, varlist = c(
    "simulator_batch",
    "simulator_full",
    "simulator_one_simulation",
    "simulation_clonal_evolution",
    "simulation_sample_phylogeny",
    "simulation_sequening_truth",
    "simulation_sequencing_sc",
    "simulation_sequencing_bulk",
    "rand_coverage",
    "folder_workplace",
    # "t_end_time",
    # "t_tau_step",
    "n_selective_clones",
    # "vec_time_points_s_mut",
    "vec_hierarchy_s_mut",
    # "expected_end_population",
    # "vec_expected_percent_select",
    "n_sample"
    # "range_population",
    # "range_clonal_perc",
    # "mindiff_clonal_perc",
    # "ploidy",
    # "choice_theta",
    # "vec_theta_parameters",
    # "vec_theta_mean",
    # "bulk_coverage_model",
    # "bulk_coverage_variables",
    # "bulk_min_alt_readcounts"
))
clusterEvalQ(cl, library(data.table))
df_all_simulation_parameters <- pblapply(cl = cl, X = 1:n_simulations, FUN = function(n_simulation) {
    ####################################################################
    t_diagnosis <- 60 # [in years]
    cell_lifespan <- 10 # [in days]
    t_MRCA <- runif(1, min = 1, max = 20) # [in years]
    ratio_theta <- runif(1, min = 1, max = 10) # <<<<<<<<<<<<<<<<<<<<<<<
    purity <- runif(1, min = 0.5, max = 1)
    mean_coverage <- runif(1, min = 30, max = 60)
    ####################################################################
    t_end_time <- round(t_MRCA * 365) # [in days]
    vec_time_points_s_mut <- t_end_time * c()
    expected_end_population <- 10^6 # [cells at t_end_time]
    vec_expected_percent_select <- (1 / (n_selective_clones + 1)) * rep(1, length = (n_selective_clones + 1))
    range_population <- c(0.8, 1.2) * expected_end_population
    range_clonal_perc <- c(20, 100)
    # mindiff_clonal_perc <- 10
    ploidy <- 2

    choice_theta <- "constant"
    theta_normal <- 6.12 # [mutations per division]
    vec_theta_parameters <- c(ratio_theta * theta_normal)
    vec_theta_mean <- vec_theta_parameters
    bulk_coverage_model <- "binomial"
    bulk_coverage_variables <- c(0, mean_coverage)
    bulk_min_alt_readcounts <- 4 # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    folder_workplace_sub <- paste0(folder_workplace, "_", n_simulation, "/")
    dir.create(folder_workplace_sub)
    truncal_mutations <- rpois(n = 1, lambda = theta_normal * round((t_diagnosis - t_MRCA) * 365 / cell_lifespan))
    table_parameters <- simulator_batch(
        n_simulations = 1,
        t_end_time = t_end_time,
        cell_lifespan = cell_lifespan,
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
        purity = purity,
        truncal_mutations = truncal_mutations,
        choice_theta = choice_theta,
        vec_theta_parameters = vec_theta_parameters,
        vec_theta_mean = vec_theta_mean,
        save_rda = FALSE,
        save_true_mutation_table = FALSE,
        output_bulk = TRUE,
        output_sc = FALSE,
        compute_parallel = FALSE,
        bulk_coverage_model = bulk_coverage_model,
        bulk_coverage_variables = bulk_coverage_variables,
        bulk_min_alt_readcounts = bulk_min_alt_readcounts,
        subfolder = folder_workplace_sub,
        file_prefix = "",
        R_libPaths = R_libPaths
    )
    table_parameters$Batch_ID <- NULL
    table_parameters <- cbind(data.frame(Simulation = n_simulation), table_parameters)
    table_parameters[["Age at diagnosis (years)"]] <- t_diagnosis
    table_parameters[["Age of MRCA (years)"]] <- t_MRCA
    table_parameters[["Truncal mutation count"]] <- truncal_mutations
    table_parameters[["Cell lifespan (days)"]] <- cell_lifespan
    table_parameters[["Tumor mutation rate ratio"]] <- ratio_theta
    table_parameters[["Tumor mutation rate"]] <- vec_theta_parameters
    table_parameters[["Purity"]] <- purity
    table_parameters[["Sequencing coverage"]] <- mean_coverage
    table_parameters[["Expected tumor cell count"]] <- expected_end_population
    table_parameters[["Sample cell count"]] <- n_sample
    table_parameters[["Ploidy"]] <- ploidy
    table_parameters[["Coverage model"]] <- bulk_coverage_model
    table_parameters[["Limit alt count"]] <- bulk_min_alt_readcounts
    return(table_parameters)
})
stopCluster(cl)
df_all_simulation_parameters <- do.call(rbind, df_all_simulation_parameters)
write.csv(df_all_simulation_parameters, file = "Parameters_simulation.csv", row.names = FALSE)
#--------------------------------------------------------Clean bulk data
for (n_simulation in 1:n_simulations) {
    filename <- paste0(folder_workplace, "_", n_simulation, "/1_mutational_data_BULK.csv")
    mut_table <- read.csv(filename)
    vec_delete <- which(mut_table$Alt_count == 0 | mut_table$Ref_count == 0)
    if (length(vec_delete) > 0) mut_table <- mut_table[-vec_delete, ]
    filename <- paste0(folder_workplace, "_", n_simulation, "/SFS_1.txt")
    write.table(mut_table, filename, sep = " ", row.names = FALSE, col.names = FALSE)
}
# ========================================GROUND TRUTH FOR SFS VARIABLES
df_all_simulation_parameters <- read.csv("Parameters_simulation.csv", header = TRUE)
df <- data.frame()
for (n_simulation in 1:n_simulations) {
    purity <- df_all_simulation_parameters[["Purity"]][n_simulation]
    t_end_time <- df_all_simulation_parameters[["Age.at.diagnosis..years."]][n_simulation]
    vec_theta_parameters <- df_all_simulation_parameters[["Tumor.mutation.rate"]][n_simulation]
    n_sample <- df_all_simulation_parameters[["Sample.cell.count"]][n_simulation]
    ploidy <- df_all_simulation_parameters[["Ploidy"]][n_simulation]
    truncal_mutations <- df_all_simulation_parameters[["Truncal.mutation.count"]][n_simulation]
    #   Retrieve clonal MRCA ages and sizes in population & sample
    simulation_variables <- read.csv(paste0(folder_workplace, "_", n_simulation, "/1_simulation_variables.csv"))
    Ns <- simulation_variables$Count_in_population
    ns <- simulation_variables$Count_in_sample
    MRCA_ages <- simulation_variables$MRCA_ages / 365
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
    ps <- purity * ns_combined / n_sample / ploidy
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
write.csv(df, paste0("Parameters_true.csv"), row.names = FALSE)
# ===============================================================MOBSTER
mobster_df <- data.frame()
model_list <- list()
for (n_simulation in 1:n_simulations) {
    #   Import mutational data
    filename <- paste0(folder_workplace, "_", n_simulation, "/SFS_1.txt")
    txtdata <- read.table(file = filename, header = FALSE)
    #   Data transformation
    data <- transform(txtdata, VAF = txtdata[, 2] / (txtdata[, 1] + txtdata[, 2]))
    last_col <- ncol(data)
    mob_data <- as.data.frame(data[, last_col])
    colnames(mob_data)[1] <- "VAF"
    #   SFS deconvolution with MOBSTER
    fit <- mobster_fit(
        mob_data,
        maxIter = 1000
    )
    #   Find best MOBSTER model
    mob_model <- fit$best
    model_list[[n_simulation]] <- mob_model # save model_list
    png(paste0(folder_workplace, "MOBSTER_", n_simulation, ".png"), res = 150, width = 15, height = 7.5, units = "in")
    print(plot(fit$best))
    dev.off()
    #   Save the results
    mobster_df[n_simulation, "Simulation"] <- n_simulation # id
    mobster_df[n_simulation, "Total_N"] <- mob_model$N # total acount
    mobster_df[n_simulation, "Tail"] <- mob_model$fit.tail # bool: if tail exists
    mobster_df[n_simulation, "Tail_Num"] <- mob_model$N.k[[1]] # number of tail
    mobster_df[n_simulation, "Tail_shape"] <- mob_model$shape # shape of tail
    mobster_df[n_simulation, "Tail_scale"] <- mob_model$scale # scale of tail
    mobster_df[n_simulation, "Kbeta_cluster"] <- mob_model$Kbeta # number of clusters
    for (k in 1:mob_model$Kbeta) {
        mobster_df[n_simulation, paste0("cl_num_", k)] <- mob_model$N.k[[k + 1]] # number of Beta
        mobster_df[n_simulation, paste0("a_", k)] <- mob_model$a[[k]] # alpha of Beta
        mobster_df[n_simulation, paste0("b_", k)] <- mob_model$b[[k]] # beta of Beta
        mobster_df[n_simulation, paste0("p_", k)] <- mobster_df[n_simulation, paste0("a_", k)] / (mobster_df[n_simulation, paste0("a_", k)] + mobster_df[n_simulation, paste0("b_", k)])
    }
}
write.csv(mobster_df, paste0("Parameters_MOBSTER.csv"), row.names = FALSE)
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
coverage_distribution <- "sample-specific"
#   Maximum number of trials for fitting each hump count
max_trials <- 10000
#---Options for fitting
#   Candidates for neutral tail powers
list_neutral_powers <- seq(1, 3, by = 0.01)
# 	Candidates for where the hump frequencies are
N_SFS_positions <- 100
list_frequencies <- seq(from = 1 / N_SFS_positions, to = 1, by = 1 / N_SFS_positions)
#---Options for plotting
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
    cat(paste0("DECODE FOR SIMULATION ", n_simulation, "...\n"))
    #---Input the SFS data
    filename_2 <- paste0(folder_workplace, "_", n_simulation, "/SFS_1.txt")
    mutation_table <- read.table(filename_2, sep = " ", header = FALSE)
    colnames(mutation_table) <- c("Ref_count", "Alt_count", "Marker")
    #---Perform SFS deconvolution
    results <- DECODE(
        mutation_table = mutation_table,
        criterion = "ICL", # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
        coverage_distribution = coverage_distribution,
        max_trials = max_trials,
        neutral_tail = TRUE, # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        compute_parallel_fit = FALSE, # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        data_marker_colors = data_marker_colors,
        plot_filename = paste0(folder_workplace, "DECODE_", n_simulation, ".png")
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
    filename <- paste0(paste0(folder_workplace, "_", n_simulation, "/SFS_DECODE_parameters_1.txt"))
    fileID <- file(filename, "w")
    writeLines(paste(sprintf("%.3f", vec_para_best_final), collapse = "\t"), fileID)
    close(fileID)
}
deconvolution_df <- cbind(
    data.frame(Simulation = 1:n_simulations),
    deconvolution_df
)
write.csv(deconvolution_df, paste0("Parameters_DECODE.csv"), row.names = FALSE)
