# devtools::install_github("mg14/mg14")
# devtools::install_github("gerstung-lab/MutationTimeR")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
R_workplace <- "/Users/apple/Desktop/MutationTimeR/SFS_CNA_deconvolution/vignettes"
R_libPaths <- ""
R_libPaths_extra <- "/Users/apple/Desktop/MutationTimeR/SFS_CNA_deconvolution/R"
R_libPaths_CINner <- "/Users/apple/Desktop/MutationTimeR/SFS_CNA_deconvolution/R_CINner"
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)

library("MutationTimeR")

setwd(R_libPaths_CINner)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)

setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)

model_name <- "MUTATIONTIMER"
folder_workplace <- "MUTATIONTIMER"
# ==================================================SET MODEL PARAMETERS
n_simulations <- 8
#---Passenger mutation rate (per nucleotide per cell division)
passenger_mutation_rate <- 1e-9 #1e-9 # try increase this
#---Probabilities of CNA
prob_CN_missegregation <- 3e-4 # try to increase this 
prob_CN_chrom_arm_missegregation <- 0 # 3e-4 # set this to 0 for now
#---Parameters for chromosome arm selection rates
s_rate_max <- 1.2
s_rate_prob_loss <- 0.2
#---Sequencing coverage distribution
n_sample <- 1000
bulk_coverage_model <- "binomial"
bulk_coverage_variables <- c(n_sample, 100)
bulk_min_alt_readcounts <- 0



# =========================================================SET UP CINNER
cell_lifespan <- 30
T_0 <- list(0, "year")
T_end <- list(80, "year")
Table_sample <- data.frame(Sample_ID = c("SA01"), Cell_count = c(n_sample), Age_sample = c(80))
selection_model <- "chrom-arm-selection"
#---Viability thresholds
bound_maximum_CN <- 4
bound_average_ploidy <- 4.5
#---Population dynamics
vec_time <- T_0[[1]]:T_end[[1]]
L <- 10000
t_0 <- 20
k <- 0.3
vec_cell_count <- L / (1 + exp(-k * (vec_time - t_0)))
table_population_dynamics <- cbind(vec_time, vec_cell_count)
#---Initialize model variables
model_variables <- BUILD_general_variables(
    cell_lifespan = cell_lifespan,
    T_0 = T_0, T_end = T_end,
    CN_arm_level = FALSE,
    Table_sample = Table_sample,
    prob_CN_missegregation = prob_CN_missegregation,
    prob_CN_chrom_arm_missegregation = prob_CN_chrom_arm_missegregation,
    selection_model = selection_model,
    bound_maximum_CN = bound_maximum_CN,
    bound_average_ploidy = bound_average_ploidy,
    table_population_dynamics = table_population_dynamics
)
#---Set up (randomized) chromosome arm selection rates
arm_id <- c(
    paste0(model_variables$cn_info$Chromosome, "p"),
    paste0(model_variables$cn_info$Chromosome, "q")
)
arm_chromosome <- rep(model_variables$cn_info$Chromosome, 2)
arm_start <- c(
    rep(1, length(model_variables$cn_info$Chromosome)),
    model_variables$cn_info$Centromere_location + 1
)
arm_end <- c(
    model_variables$cn_info$Centromere_location,
    model_variables$cn_info$Bin_count
)
arm_s <- rep(1, length(arm_id))
set.seed(1)
for (i in 1:length(arm_s)) {
    arm_s[i] <- runif(1, 1, s_rate_max)
    if (runif(1) < s_rate_prob_loss) arm_s[i] <- 1 / arm_s[i]
    # if (grepl("q$", arm_id[i])) {
    #     arm_s[i] <- 1
    # } else if (grepl("p$", arm_id[i])) {
    #     arm_s[i] <- runif(1, 1, s_rate_max)
    #     if (runif(1) < s_rate_prob_loss) arm_s[i] <- 1 / arm_s[i]
    # }
}
set.seed(NULL)
table_arm_selection_rates <- data.frame(Arm_ID = arm_id, Chromosome = arm_chromosome, Bin_start = arm_start, Bin_end = arm_end, s_rate = arm_s)
model_variables <- BUILD_driver_library(model_variables = model_variables, table_arm_selection_rates = table_arm_selection_rates, )
#---Set up initial cell population
cell_count <- 1
CN_matrix <- BUILD_cn_normal_autosomes(model_variables$cn_info)
drivers <- list()
model_variables <- BUILD_initial_population(model_variables = model_variables, cell_count = cell_count, CN_matrix = CN_matrix, drivers = drivers)
#---Check model variables
model_variables <- CHECK_model_variables(model_variables)
# =============================================CREATE CINNER SIMULATIONS
simulator_full_program(
    model = model_variables, model_prefix = model_name,
    n_simulations = n_simulations,
    stage_final = 3,
    compute_parallel = FALSE,
    folder_workplace = folder_workplace,
    R_libPaths = R_libPaths
)
# ========================================SUPPLEMENT CINNER SIMULATIONS:
# ================================ADD PASSENGER MUTATIONS & SEGMENTED CN
# ============================================AND OUTPUT TRUE CNA TIMING
CINner_postprocessing(
    n_simulations = n_simulations,
    model_name = model_name,
    folder_workplace = folder_workplace,
    passenger_mutation_rate = passenger_mutation_rate,
    bulk_coverage_model = bulk_coverage_model,
    bulk_coverage_variables = bulk_coverage_variables,
    bulk_min_alt_readcounts = bulk_min_alt_readcounts,
    compute_parallel = FALSE
)
# =========================================================MUTATIONTIMER
for (i_simulation in 1:n_simulations) {
    cat("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
    cat(paste0("MUTATIONTIMER FOR SIMULATION-", i_simulation, "..."))
    #   Load the simulated mutation and copy number data
    cn_filename <- paste0(folder_workplace, "/", model_name, "_simulation_", i_simulation, "_cn.RData")
    mut_filename <- paste0(folder_workplace, "/", model_name, "_simulation_", i_simulation, "_mut.vcf")
    mut_vcf <- readVcf(mut_filename)
    load(cn_filename)
    #   Run MutationTimeR
    mutationtimer <- mutationTime(mut_vcf, cn_granges, n.boot = 10)
    mcols(cn_granges) <- cbind(mcols(cn_granges), mutationtimer$T)
    #   Plot MutationTimeR results
    png_filename <- paste0(folder_workplace, "/", model_name, "_simulation_", i_simulation, "_mutationtimer.png")
    png(png_filename, width = 500, height = 750)
    plotSample(mut_vcf, cn_granges)
    dev.off()
}
