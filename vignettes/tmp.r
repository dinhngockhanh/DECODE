n_selective_clones <- 0
for (n_simulation in 1:n_simulations) {
    load(paste0("/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes/TEST/_", n_simulation, "/simulation_1.rda"))
    Count_in_sample <- rep(0, n_selective_clones + 1)
    for (clone in 0:n_selective_clones) Count_in_sample[clone + 1] <- sum(simulation$sample_genotype == clone)
    simulation_variables <- data.frame(
        Clone_ID = paste0("Clone_", 0:n_selective_clones),
        MRCA_ages = simulation$MRCA_ages,
        Count_in_population = simulation$record_vec_populations[nrow(simulation$record_vec_populations), ],
        Count_in_sample = Count_in_sample
    )
    filename <- paste0("/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes/TEST/_", n_simulation, "/1_simulation_variables.csv")
    write.csv(simulation_variables, file = filename)
}
