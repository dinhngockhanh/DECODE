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

setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)

folder_workplace <- "03_TEST_SFS_DECONVOLUTION/"
#---------------------------------------------------Set model parameters
n_simulations <- 100



t_end_time <- 1000
t_tau_step <- 1
n_selective_clones <- 1
vec_time_points_s_mut <- t_end_time * c(0.6)
vec_hierarchy_s_mut <- c(0)
expected_end_population <- 10^6
vec_expected_percent_select <- (1 / (n_selective_clones + 1)) * rep(1, length = (n_selective_clones + 1))
n_sample <- 10000  # CHANGE from 100000 to 10000
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
# ===============================================================MOBSTER
library(mobster) # Load the mobster package
df <- data.frame() # Create an empty data frame
for (i in 1:n_simulations) {
    print(paste("Dealing with file: ", i))
    # import data
    # base_filename <- paste(folder_workplace , sep = " ", " ")"C:/Users/Mayin/Documents/1GRADUATE/1. Study/2. 24Spring/5398 Dinh/DATA/02_TEST_SFS_DECONVOLUTION/"
    filename <- paste0(folder_workplace,"SFS_", i, ".txt")
    txtdata <- read.table(file = filename, header = FALSE)

    # data transformation
    data <- transform(txtdata, VAF = txtdata[,2] / (txtdata[,1] + txtdata[,2]))
    last_col <- ncol(data)
    mob_data <- as.data.frame(data[,last_col])
    colnames(mob_data)[1] <- "VAF"

    # # data check
    # head(mob_data,5)
    # class(mob_data)
    # hist(mob_data[,1])
    # Sys.sleep(3)

    # mobster
    # mobster:::template_parameters_fast_setup() # show basic setup
    fit <-  mobster_fit(
    mob_data,    
    auto_setup = "FAST"
    )

    print("Start fitting")
    mob_model <- fit$best # best model
    print("Fitting done")

    # plot_path <- paste0("C:/Users/Mayin/Desktop/02_Mobster_Test_plot/", i, ".png")
    # png(plot_path) # Save the plot
    # plot(mob_model)
    # dev.off()


    
     print("Start saving")
    df[i,"id"] <- i # id
    df[i, "Total_N"] <- mob_model$N # total acount

    # tail part 
    df[i, "Tail"] <- mob_model$fit.tail # bool: if tail exists
    df[i, "Tail_Num"] <- mob_model$N.k[[1]] # number of tail
    df[i, "Tail_shape"] <- mob_model$shape # shape of tail
    df[i, "Tail_scale"] <- mob_model$scale # scale of tail

    # cluster part
    df[i, "Kbeta_cluster"] <- mob_model$Kbeta # number of clusters
    for (k in 1:mob_model$Kbeta){
    df[i, paste0("cl_num_", k)] <- mob_model$N.k[[k+1]] # alpha of Beta
    df[i, paste0("a_", k)] <- mob_model$a[[k]] # alpha of Beta
    df[i, paste0("b_", k)] <- mob_model$b[[k]] # beta of Beta
    }

    print("Saving done")
    print("=====================================")

    #print finished
    if (i == n_simulations){
        print("=====================================")
        print(paste("||       All done ", n_simulations, " files!      ||"))
        print("=====================================")
    }
}

# Summary of the result
print(paste("Shape of the result: ", dim(df)[1],",", dim(df)[2]))
print(head(df, 5))
csv_name <- paste0(folder_workplace, "Mobster_Test.csv")
write.csv(df, file = csv_name)



    # ...
# ...
# ...
# ...
# ...

# ===============================================================GROUND TRUTH
for (n_simulation in 1:n_simulations) {
    # ...
# ...
# ...
# ...
# ...
}
# ===============================================================COMPARISON
