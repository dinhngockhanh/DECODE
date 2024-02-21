# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
# R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SFS_CNA_deconvolution/vignettes"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SFS_CNA_deconvolution/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Yining - Laptop
R_workplace <- "C:/Users/Mayin/Documents/1GRADUATE/1. Study/2. 24Spring/5398 Dinh/DATA"
R_libPaths <- ""
R_libPaths_extra <- "C:/Users/Mayin/Documents/1GRADUATE/1. Study/2. 24Spring/5398 Dinh/github_clone/SFS_CNA_deconvolution-1/R"
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)
library(data.table)
library(stringi)
library(tidyverse)
library(R.utils)
library(parallel)
library(pbapply)
library(mobster)

setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)

folder_workplace <- "04_TEST_SFS_DECONVOLUTION/"
#---------------------------------------------------Set model parameters
n_simulations <- 8



t_end_time <- 1000
t_tau_step <- 1
n_selective_clones <- 1
vec_time_points_s_mut <- t_end_time * c(0.6)
vec_hierarchy_s_mut <- c(0)
expected_end_population <- 10^6
vec_expected_percent_select <- (1 / (n_selective_clones + 1)) * rep(1, length = (n_selective_clones + 1))
n_sample <- 10000 # CHANGE from 100000 to 10000
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
    for (i in length(vec_hierarchy_s_mut):1) {
        ns_combined[vec_hierarchy_s_mut[i] + 1] <- ns_combined[vec_hierarchy_s_mut[i] + 1] + ns_combined[i + 1]
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
    for (i in length(vec_hierarchy_s_mut):1) {
        Ks[i + 1] <- Ks[i + 1] - Ks[vec_hierarchy_s_mut[i] + 1]
    }
    Ks[1] <- Ks[1] + truncal_mutations
    #   Save the results
    df <- rbind(df, c(n_simulation, A, alpha, ps, Ks))
}
names(df) <- c("Simulation", "A", "alpha", paste0("p_", 1:(n_selective_clones + 1)), paste0("K_", 1:(n_selective_clones + 1)))
write.csv(df, paste0(folder_workplace, "Parameters_true.csv"), row.names = FALSE)

# ===============================================================MOBSTER
mob_df <- data.frame()
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
    # png(paste0(folder_workplace, "MOBSTER_", i, ".png"))
    # plot(fit$best)
    # dev.off()
    #   Save the results
    mob_df[i, "Simulation"] <- i # id
    mob_df[i, "Total_N"] <- mob_model$N # total acount
    mob_df[i, "Tail"] <- mob_model$fit.tail # bool: if tail exists
    mob_df[i, "Tail_Num"] <- mob_model$N.k[[1]] # number of tail
    mob_df[i, "Tail_shape"] <- mob_model$shape # shape of tail
    mob_df[i, "Tail_scale"] <- mob_model$scale # scale of tail
    mob_df[i, "Kbeta_cluster"] <- mob_model$Kbeta # number of clusters
    for (k in 1:mob_model$Kbeta) {
        mob_df[i, paste0("cl_num_", k)] <- mob_model$N.k[[k + 1]] # number of Beta
        mob_df[i, paste0("a_", k)] <- mob_model$a[[k]] # alpha of Beta
        mob_df[i, paste0("b_", k)] <- mob_model$b[[k]] # beta of Beta
    }
}
write.csv(mob_df, paste0(folder_workplace, "Parameters_mobster.csv"), row.names = FALSE)
# # ===============================================================COMPARISON
filename_1 <- paste0(folder_workplace, "Parameters_true.csv")
df <- read.csv(filename_1)
filename_2 <- paste0(folder_workplace, "Parameters_mobster.csv")
data_2 <- read.csv(filename_2)
mob_df <- data_2

com_df <- mob_df # initialize com_df with mob_df

#== Rearrange Clusters Part ==
## calculate p_i: can be inserted in Mobster part
for (k in 1:(((dim(com_df)[2])-7)/3)){
    com_df[paste0("p_", k)] <-com_df[paste0("a_", k)]/(com_df[paste0("a_", k)]+com_df[paste0("b_", k)])
}

# rerange cl_num, (a, b) based on p
p_cols <- grep("^p_", names(com_df), value = TRUE) # get columns represent p
com_df[p_cols][is.na(com_df[p_cols])] <- -1 # replace NA with -1 for p
for (i in 1:(dim(com_df)[1])){
    p_row <- com_df[i, p_cols]
    p_rowsorted_indices <- order(as.vector(unlist(p_row)), decreasing = TRUE) # order
    cn <- 1
    for (id in p_rowsorted_indices){
        com_df[i, paste0("K_", id)] <- com_df[i, paste0("cl_num_", cn)]
        # com_df[i, paste0("re_a_", id)] <- com_df[i, paste0("a_", cn)]
        # com_df[i, paste0("re_b_", id)] <- com_df[i, paste0("b_", cn)]
        cn <- cn + 1
  }
}

cols_to_drop <- grep("^(a_|b_|cl_num_)", names(com_df), value = TRUE)  # drop old columns
com_df <- com_df[, !(names(com_df) %in% cols_to_drop)]

#== Compare Part ==
## Number of clusters
freq_kbeta <- table(com_df$Kbeta_cluster)
barplot(freq_kbeta, main = "The Number of Clusters from MOBSTER", xlab = "Values", ylab = "Number of Clusters",  border = "black")

# compre p1
plot(com_df$p_1, df$p_1, xlab = "p_1 from MOBSTER", ylab = "p_1 from Ground Truth", main = "Comparison of p_1", pch = 16, col = "blue")
abline(a = 0, b = 1, lty = 2, main = "y = x")

# compare p_2 还要改
plot(com_df$p_2, df$p_2, xlab = "p_2 from MOBSTER", ylab = "p_2 from Ground Truth", main = "Comparison of p_2", pch = 16, col = "blue")
abline(a = 0, b = 1, lty = 2, main = "y = x")

# compare K_1
plot(com_df$K_1, df$K_1, xlab = "K_1 from MOBSTER", ylab = "K_1 from Ground Truth", main = "Comparison of K_1", pch = 16, col = "blue")
abline(a = 0, b = 1, lty = 2, main = "y = x")

# compare K_2 还要改
plot(com_df$K_2, df$K_2, xlab = "K_2 from MOBSTER", ylab = "K_2 from Ground Truth", main = "Comparison of K_2", pch = 16, col = "blue")
abline(a = 0, b = 1, lty = 1, main = "y = x")

#== Calculate Tail Part ==

(2.41*(0.007246^2.41))/((0.125/3)^3.41)
x <- seq(0.001, 10, by = 0.001)
plot(x, (1.41 * (0.007246^1.41)) / (x^2.41), type = "l", xlab = "x", ylab = "y")
y <- (2.41 * (0.007246^2.41)) / (x^3.41)
plot(x, y, type = "l", xlab = "x", ylab = "y")

x <- seq(0, 0.25, by = 0.001)
y1 <- (1.41 * (0.007246^1.41)) / (x^2.41)
y2 <- (2.41 * (0.007246^2.41)) / (x^3.41)
plot(x, y1, type = "l", ylim = c(-1, 1), xlim = c(0, 0.25), xlab = "x", ylab = "y")
lines(x, y2, col = "red")

1.41*0.007246/0.41
