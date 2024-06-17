# n_selective_clones <- 0
# for (n_simulation in 1:n_simulations) {
#     load(paste0("/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes/TEST/_", n_simulation, "/simulation_1.rda"))
#     Count_in_sample <- rep(0, n_selective_clones + 1)
#     for (clone in 0:n_selective_clones) Count_in_sample[clone + 1] <- sum(simulation$sample_genotype == clone)
#     simulation_variables <- data.frame(
#         Clone_ID = paste0("Clone_", 0:n_selective_clones),
#         MRCA_ages = simulation$MRCA_ages,
#         Count_in_population = simulation$record_vec_populations[nrow(simulation$record_vec_populations), ],
#         Count_in_sample = Count_in_sample
#     )
#     filename <- paste0("/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes/TEST/_", n_simulation, "/1_simulation_variables.csv")
#     write.csv(simulation_variables, file = filename)
# }
#-----------------------------------------------------------------------
# for (n_simulation in 1:n_simulations) {
#     load(paste0("/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes/TEST/_", n_simulation, "/simulation_1.rda"))

#     simulation$cell_phylogeny_hclust$height <- simulation$cell_phylogeny_hclust$height / 365

#     png(paste0("/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes/TEST/CINner_phylogeny_", n_simulation, ".png"), res = 150, width = 30, height = 15, units = "in")
#     p <- plot(simulation$cell_phylogeny_hclust, labels = FALSE, ylab = "year")
#     print(p)
#     dev.off()
# }
#-----------------------------------------------------------------------
# for (n_simulation in 1:n_simulations) {
#     load(paste0("/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes/TEST/_", n_simulation, "/simulation_1.rda"))
#     N_neutral <- 0
#     for (i in 1:length(simulation$sample_mutational_table_truth_node_markers)) {
#         if (grepl("Foreground_", simulation$sample_mutational_table_truth_node_markers[[i]])) {
#             N_neutral <- N_neutral + simulation$sample_mutational_table_truth_node_mutation_counts[[i]]
#         }
#     }
#     simulation$sample_mutational_table_truth_node_markers
# }
#-----------------------------------------------------------------------
# df_all <- c()
# for (batch in 1:10) {
#     # for (batch in 1:1) {
#     df <- read.csv(paste0("/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes/ICGC-MOBSTER/", batch, "_Parameters_mobster.csv"))
#     df_new <- data.frame(
#         Sample = df$Sample,
#         Mutation_count_in_fitting = df$Total_N,
#         Tail = df$Tail,
#         Tail_power = df$Tail_shape + 1,
#         Tail_pareto_shape = df$Tail_shape,
#         Tail_pareto_scale = df$Tail_scale,
#         Tail_mutcount_observed = df$Tail_Num,
#         Cluster_count = df$Kbeta_cluster
#     )
#     for (i in 1:max(df_new$Cluster_count)) {
#         df_new[[paste0("Cluster_mutcount_observed_", i)]] <- df[[paste0("cl_num_", i)]]
#         df_new[[paste0("Cluster_frequency_", i)]] <- df[[paste0("a_", i)]] / (df[[paste0("a_", i)]] + df[[paste0("b_", i)]])
#         df_new[[paste0("Cluster_beta_a_", i)]] <- df[[paste0("a_", i)]]
#         df_new[[paste0("Cluster_beta_b_", i)]] <- df[[paste0("b_", i)]]
#     }
#     if (batch > 1 & max(df_all$Cluster_count) < max(df_new$Cluster_count)) {
#         for (i in (max(df_all$Cluster_count) + 1):max(df_new$Cluster_count)) {
#             df_all[[paste0("Cluster_mutcount_observed_", i)]] <- NA
#             df_all[[paste0("Cluster_frequency_", i)]] <- NA
#             df_all[[paste0("Cluster_beta_a_", i)]] <- NA
#             df_all[[paste0("Cluster_beta_b_", i)]] <- NA
#         }
#     } else if (batch > 1 & max(df_all$Cluster_count) > max(df_new$Cluster_count)) {
#         for (i in (max(df_new$Cluster_count) + 1):max(df_all$Cluster_count)) {
#             df_new[[paste0("Cluster_mutcount_observed_", i)]] <- NA
#             df_new[[paste0("Cluster_frequency_", i)]] <- NA
#             df_new[[paste0("Cluster_beta_a_", i)]] <- NA
#             df_new[[paste0("Cluster_beta_b_", i)]] <- NA
#         }
#     }
#     df_all <- rbind(df_all, df_new)
# }
# write.csv(df_all, file = "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes/Parameters_ICGC_MOBSTER.csv", row.names = FALSE)
# length(which(df_all$Tail == TRUE)) / length(df_all$Tail)
#-----------------------------------------------------------------------
filename <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes/trial/SFS_1.txt"
txtdata <- read.table(file = filename, header = FALSE)
data <- transform(txtdata, VAF = txtdata[, 2] / (txtdata[, 1] + txtdata[, 2]))
last_col <- ncol(data)
mob_data <- as.data.frame(data[, last_col])
colnames(mob_data)[1] <- "VAF"


png(paste0("VAF_histogram.png"), res = 150, width = 15, height = 7.5, units = "in")
hist(data$VAF, main = "VAF Histogram", xlab = "VAF", ylab = "Frequency")
dev.off()


#---SFS deconvolution with MOBSTER
MOBSTER_result <- try(
    {
        mobster_fit(
            mob_data,
            tail = c(TRUE),
            parallel = FALSE,
            description = paste0("Simulation ", n_simulation)
        )
    },
    silent = TRUE
)
# MOBSTER_result <- mobster_fit(
#     mob_data,
#     tail = c(TRUE), # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#     parallel = FALSE, # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#     seed = 200,
#     description = paste0("Simulation ", n_simulation)
# )
#---Plot MOBSTER deconvolution
if (!inherits(MOBSTER_result, "try-error")) {
    png(paste0("MOBSTER_", n_simulation, ".png"), res = 150, width = 15, height = 7.5, units = "in")
    print(plot(MOBSTER_result$best))
    dev.off()
    png(paste0("MOBSTER_model_selection_", n_simulation, ".png"), res = 150, width = 15, height = 7.5, units = "in")
    print(plot_model_selection(MOBSTER_result))
    dev.off()
}
