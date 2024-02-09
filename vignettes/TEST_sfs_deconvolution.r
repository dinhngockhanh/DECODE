# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SFS_CNA_deconvolution/vignettes"
R_libPaths <- ""
R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SFS_CNA_deconvolution/R"
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)
library(data.table)
library(stringi)
library(tidyverse)
library(R.utils)
library(parallel)
library(pbapply)

setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "*.r$")
sapply(files_sources, source)
setwd(R_workplace)

folder_workplace <- "TEST_SFS_DECONVOLUTION/"
#---------------------------------------------------Set model parameters
n_simulations <- 100
t_end_time <- 3000
t_tau_step <- 1
n_selective_clones <- 1
vec_time_points_s_mut <- t_end_time * c(0.6)
vec_hierarchy_s_mut <- c(0)
expected_end_population <- 10^6
vec_expected_percent_select <- (1 / (n_selective_clones + 1)) * rep(1, length = (n_selective_clones + 1))
n_sample <- 10000
range_population <- c(0.8, 1.2) * expected_end_population
range_clonal_perc <- c(30, 70)
choice_theta <- "constant"
vec_theta_parameters <- rep(0.4, length = (n_selective_clones + 1))
vec_theta_mean <- vec_theta_parameters
bulk_coverage_model <- "binomial"
bulk_coverage_variables <- c(0, 100)
bulk_min_alt_readcounts <- 0
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
    subfolder = folder_workplace
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
