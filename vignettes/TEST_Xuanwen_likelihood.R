library(R.matlab)
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

option_dist_coverage <- "binomial"
dist_coverage_var_1 <- 100

#----------------------------------------------------Options for fitting
N_SFS_positions <- 100
N_fitting_rounds <- 200
filename_1 <- paste0(
    "D:/5398_dataset/Binomial_PDF_",
    N_end, "_",
    r_max, "_",
    min_variant_read, "_",
    min_total_read, "_",
    SFS_totalsteps_base, "_",
    option_ploidy, ".mat"
)

inputBinomialMatrix <- readMat(filename_1)
matrix_binomial_PDF <- inputBinomialMatrix$matrix.binomial.PDF

n_simulations <- 100

# table_parameters <- matrix(0, nrow = n_simulations, ncol = 2 * max_hump_count + 1)
for (n_simulation in 1:n_simulations) {
    filename_2 <- paste0("C:/5398_R/vignettes/TEST_CINner_lite/TEST_CINner_lite_", n_simulation, ".txt")
    data <- read.table(filename_2, sep = " ", header = FALSE)
    vec_refcount <- as.numeric(data[, 1])
    vec_altcount <- as.numeric(data[, 2])
    vec_totcount <- vec_refcount + vec_altcount
    #
    L <- max(vec_totcount)
    TCGA_coverage_values <- 1:L
    TCGA_coverage_PDF <- rep(0, L)
    TCGA_coverage_PDF <- prep_distribution_patient(vec_totcount)

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

    # Main optimization loop
    N_humps <- -1
    AIC_old <- Inf
    AIC_best_current <- 0
    # AIC_threshold <- 0.01
    N_fitting_rounds <- 200

    # while (is.infinite(AIC_old) || (AIC_old - AIC_new >= AIC_threshold)) {
    while (AIC_best_current < AIC_old) {
        N_humps <- N_humps + 1
        num_parameters <- 1 + 2 * N_humps
        N_fitting_rounds_current <- factorial(N_humps + 1) * N_fitting_rounds
        AIC_best_current <- Inf
        for (i in 1:N_fitting_rounds_current) {
            results <- fit_SFS_one_iteration_likelihood(vec_SFS_real, N_humps, vec_SFS_positions, library_SFS_component)
            logL <- results$log_L
            AIC_current <- 2 * num_parameters - 2 * -results$log_L
            if (AIC_current < AIC_best_current) {
                AIC_best_current <- AIC_current
                vec_para_best_current <- results$parameters
            }
        }
        if (AIC_best_current < AIC_old) {
            AIC_old <- AIC_best_current
            vec_para_best_final <- vec_para_best_current
        }
        # if (is.infinite(AIC_old) || (AIC_old - AIC_new >= AIC_threshold)) {
        #     AIC_old <- AIC_new
        # } else {
        #     vec_para_best_final <- results$parameters
        #     break
        # }
    }
    vec_A_and_K <- vec_para_best_final[seq(1, length(vec_para_best_final), by = 2)]
    vec_p <- vec_para_best_final[seq(2, length(vec_para_best_final) - 1, by = 2)]

    cat(sprintf("Final number of humps: %d, AIC: %f\n", N_humps, AIC_old))
}
