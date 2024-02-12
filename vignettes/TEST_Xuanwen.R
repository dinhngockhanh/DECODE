# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Xuanwen - Laptop
R_workplace <- "C:/5398_R/vignettes" # where the vignette folder is located
R_libPaths <- ""
R_libPaths_extra <- "" # where the R folder is located
R_libPaths_binomial_table <- "D:/5398_dataset" # where the binomial tables are located
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SFS_CNA_deconvolution/vignettes"
R_libPaths <- ""
R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SFS_CNA_deconvolution/R"
R_libPaths_binomial_table <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/MK-Cod.Analysis of the SFS/Core_function_for_SFS_fitting/Binomial_tables"
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)
library(R.matlab)
library(parallel)
library(pbapply)

setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)

folder_workplace <- "TEST_SFS_DECONVOLUTION/"
# ==========================================MAKE CINNER LITE SIMULATIONS
#---------------------------------------------------Set model parameters
n_simulations <- 100



t_end_time <- 1000
t_tau_step <- 1
n_selective_clones <- 1
vec_time_points_s_mut <- t_end_time * c(0.6)
vec_hierarchy_s_mut <- c(0)
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
# =====================================================SFS DECONVOLUTION
#---------------------------------------------------Set model parameters
# Number of cells at T_N
N_end <- 1000
# Minimum and maximum number of reads
r_min <- 0
r_max <- 500
# Minimum variant read count to be accepted
min_variant_read <- 5
# Minimum total read count to be accepted
min_total_read <- 0
# Number of steps to divide SFS frequencies in [0,1]
SFS_totalsteps <- 25
SFS_totalsteps_base <- 100
# Choice of ploidy, which changes the binomial rate
option_ploidy <- 1
# Maximum number of humps
max_hump_count <- 1
option_dist_coverage <- "binomial"
dist_coverage_var_1 <- 100
#----------------------------------------------------Options for fitting
N_SFS_positions <- 20
N_fitting_rounds <- 200
threshold_stop <- 0.95
threshold_error <- 0.2
#---------------------------------------------------Input binomial table
filename_1 <- paste0(
    R_libPaths_binomial_table, "/Binomial_PDF_",
    N_end, "_",
    r_max, "_",
    min_variant_read, "_",
    min_total_read, "_",
    SFS_totalsteps_base, "_",
    option_ploidy, ".mat"
)
inputBinomialMatrix <- readMat(filename_1)
matrix_binomial_PDF <- inputBinomialMatrix$matrix.binomial.PDF
#---------------------------------------------Deconvolution for each SFS
table_parameters <- matrix(0, nrow = n_simulations, ncol = 2 * max_hump_count + 1)
for (n_simulation in 1:n_simulations) {
    #   Input the SFS data
    filename_2 <- paste0(R_workplace, "/", folder_workplace, "SFS_", n_simulation, ".txt")
    data <- read.table(filename_2, sep = " ", header = FALSE)
    vec_refcount <- as.numeric(data[, 1])
    vec_altcount <- as.numeric(data[, 2])
    vec_totcount <- vec_refcount + vec_altcount
    # 	Prepare the total readcount distribution
    L <- max(vec_totcount)
    TCGA_coverage_values <- 1:L
    TCGA_coverage_PDF <- rep(0, L)
    TCGA_coverage_PDF <- prep_distribution_patient(vec_totcount)
    # 	Prepare the SFS library
    vec_SFS_positions <- seq(from = 1 / N_SFS_positions, to = 1, by = 1 / N_SFS_positions)
    library_SFS_component <- vector("list", 2 * N_SFS_positions)
    dim(library_SFS_component) <- c(2, N_SFS_positions)
    vec_para <- c(1)
    vec_SFS <- SFS_expected(vec_para)
    library_SFS_component[[1, 1]] <- vec_SFS
    for (i in seq_along(vec_SFS_positions)) {
        p <- vec_SFS_positions[i]
        vec_para <- c(0, p, 1)
        vec_SFS <- SFS_expected(vec_para)
        library_SFS_component[[2, i]] <- vec_SFS
    }
    # 	Prepare the real SFS
    no_mutations_total <- length(vec_refcount)
    vec_freq <- seq(1, SFS_totalsteps) / SFS_totalsteps
    vec_SFS_real <- rep(0, SFS_totalsteps)
    mutation_count <- 0
    for (j in 1:no_mutations_total) {
        no_variant <- vec_altcount[j]
        no_total <- vec_refcount[j] + vec_altcount[j]
        if (no_variant >= min_variant_read && no_total >= min_total_read) {
            mutation_count <- mutation_count + 1
            VAF <- no_variant / no_total
            pos <- which(vec_freq >= VAF)[1]
            vec_SFS_real[pos] <- vec_SFS_real[pos] + 1
        }
    }
    #   SFS deconvolution
    N_humps <- -1
    ratio_error <- 0
    err_best_final <- Inf
    while (ratio_error < threshold_stop && err_best_final > threshold_error && N_humps < max_hump_count) {
        N_humps <- N_humps + 1
        err_best_current <- Inf
        vec_para_best_current <- c() # An empty numeric vector in R
        N_fitting_rounds_current <- factorial(N_humps + 1) * N_fitting_rounds
        for (i in 1:N_fitting_rounds_current) {
            results <- fit_SFS_one_iteration(vec_SFS_real, N_humps, vec_SFS_positions, library_SFS_component)
            err <- results$err
            vec_para <- results$vec_para
            if (err < err_best_current) {
                err_best_current <- err
                vec_para_best_current <- vec_para
            }
        }
        if (N_humps == 0) {
            ratio_error <- 0
            err_best_final <- err_best_current
            vec_para_best_final <- vec_para_best_current
        } else {
            ratio_error <- err_best_current / err_best_final
            if (ratio_error < threshold_stop) {
                err_best_final <- err_best_current
                vec_para_best_final <- vec_para_best_current
            }
        }
        cat(sprintf("%d humps: error = %f; error ratio = %f:      ", N_humps, err_best_current, ratio_error))
        cat(sprintf("%.3f   ", vec_para_best_current))
        cat("\n")
    }
    # Store text result
    table_parameters[n_simulation, ] <- vec_para_best_final
    filename <- paste0(R_workplace, "/", folder_workplace, "Deconvolution_parameters_", n_simulation, ".txt")
    fileID <- file(filename, "w")
    writeLines(paste(sprintf("%.3f", vec_para_best_final), collapse = "\t"), fileID)
    close(fileID)
    # Store pictures
    filename <- paste0(R_workplace, "/", folder_workplace, "Deconvolution_", n_simulation, ".png")
    png(filename)
    barplot(height = vec_SFS_real, names.arg = vec_freq, col = "blue", main = "SFS Fitting Results")
    bar_centers <- barplot(height = vec_SFS_real, names.arg = vec_freq, col = "blue", main = "SFS Fitting Results")
    vec_A_and_K <- vec_para_best_final[seq(1, length(vec_para_best_final), by = 2)]
    vec_p <- vec_para_best_final[seq(2, length(vec_para_best_final) - 1, by = 2)]
    vec_SFS_model <- compute_SFS_one_iteration(vec_A_and_K, vec_p, vec_SFS_positions, library_SFS_component)
    lines(bar_centers, vec_SFS_model, col = "red", lwd = 2)
    dev.off()
}
