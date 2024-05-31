# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
R_data <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SFS_CNA_deconvolution/data/PCAWG"
R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SFS_CNA_deconvolution/vignettes"
R_libPaths <- ""
R_libPaths_extra <- "C:/Users/Mayin/Documents/1GRADUATE/1. Study/2. 24Spring/5398 Dinh/github_clone/SFS_CNA_deconvolution-1/R"
R_libPaths_binomial_table <- ""
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)
library(dplyr)
library(VariantAnnotation)
library(GenomicRanges)
library(R.utils)
library(readxl)
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
# # ====================================EXTRACT MUTATIONAL DATA FROM PCAWG
# #--------------------------------------------------Unzip mutational data
# directory <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/DATASETS/PCAWG/consensus_snv_indel/final_consensus_snv_indel_passonly_icgc.public/indel"
# indel_files <- list.files(directory, pattern = "\\.vcf\\.gz$")
# for (file in indel_files) {
#     gunzip(paste0("/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/DATASETS/PCAWG/consensus_snv_indel/final_consensus_snv_indel_passonly_icgc.public/indel/", file), remove = FALSE, overwrite = TRUE)
# }
# directory <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/DATASETS/PCAWG/consensus_snv_indel/final_consensus_snv_indel_passonly_icgc.public/snv_mnv"
# snv_mnv_files <- list.files(directory, pattern = "\\.vcf\\.gz$")
# for (file in snv_mnv_files) {
#     gunzip(paste0("/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/DATASETS/PCAWG/consensus_snv_indel/final_consensus_snv_indel_passonly_icgc.public/snv_mnv/", file), remove = FALSE, overwrite = TRUE)
# }
# #---------------------------------------Prepare sample information table
# consensus_20170218_purity_ploidy <- read.table("/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/DATASETS/PCAWG/consensus_cnv/consensus.20170218.purity.ploidy.txt", header = TRUE, sep = "\t")
# consensus_20170218_purity_ploidy <- read.table("consensus.20170218.purity.ploidy.txt", header = TRUE, sep = "\t")
# pcawg_sample_sheet <- read.table("/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/DATASETS/PCAWG/pcawg_sample_sheet.tsv", header = TRUE, sep = "\t")
# pcawg_specimen_histology_August2016_v9 <- as.data.frame(read_excel("/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/DATASETS/PCAWG/pcawg_specimen_histology_August2016_v9.xlsx"))
# sample_df <- data.frame(
#     aliquot_id = unique(
#         c(
#             sapply(indel_files, function(x) strsplit(x, "\\.")[[1]][1]),
#             sapply(snv_mnv_files, function(x) strsplit(x, "\\.")[[1]][1])
#         )
#     )
# )
# sample_df <- merge(sample_df, pcawg_sample_sheet, by = "aliquot_id")
# sample_df <- merge(sample_df, pcawg_specimen_histology_August2016_v9[which(pcawg_specimen_histology_August2016_v9$specimen_library_strategy == "WGS"), ], by = "icgc_specimen_id")
# sample_df <- merge(sample_df, consensus_20170218_purity_ploidy, by.x = "aliquot_id", by.y = "samplename")
# write.csv(sample_df, file = "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SFS_CNA_deconvolution/data/PCAWG/sample_information.csv", row.names = FALSE)
# #--------------------------------Extract mutational data for each sample
# pb <- txtProgressBar(
#     min = 0,
#     max = length(sample_df$aliquot_id),
#     style = 3,
#     width = 50,
#     char = "+"
# )
# for (i_sample in 1:length(sample_df$aliquot_id)) {
#     sample <- sample_df$aliquot_id[i_sample]
#     setTxtProgressBar(pb, i_sample)
#     #---Get INSERTION/DELETION data
#     indel_filename <- paste0("/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/DATASETS/PCAWG/consensus_snv_indel/final_consensus_snv_indel_passonly_icgc.public/indel/", sample, ".consensus.20161006.somatic.indel.vcf")
#     if (file.exists(indel_filename)) {
#         indel_vcf <- readVcf(indel_filename)
#         if (nrow(indel_vcf) > 0) {
#             indel_df <- data.frame(
#                 chromosome = as.character(seqnames(indel_vcf)),
#                 start = start(ranges(indel_vcf)),
#                 end = end(ranges(indel_vcf)),
#                 width = width(ranges(indel_vcf)),
#                 type = "indel"
#             )
#             indel_df <- cbind(indel_df, info(indel_vcf))
#             indel_flag <- TRUE
#         } else {
#             indel_flag <- FALSE
#         }
#     } else {
#         indel_flag <- FALSE
#     }
#     #---Get SNV/MNV data
#     snv_mnv_filename <- paste0("/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/DATASETS/PCAWG/consensus_snv_indel/final_consensus_snv_indel_passonly_icgc.public/snv_mnv/", sample, ".consensus.20160830.somatic.snv_mnv.vcf")
#     if (file.exists(snv_mnv_filename)) {
#         snv_mnv_vcf <- readVcf(snv_mnv_filename)
#         if (nrow(snv_mnv_vcf) > 0) {
#             snv_mnv_df <- data.frame(
#                 chromosome = as.character(seqnames(snv_mnv_vcf)),
#                 start = start(ranges(snv_mnv_vcf)),
#                 end = end(ranges(snv_mnv_vcf)),
#                 width = width(ranges(snv_mnv_vcf)),
#                 type = "snv_mnv"
#             )
#             snv_mnv_df <- cbind(snv_mnv_df, info(snv_mnv_vcf))
#             snv_mnv_flag <- TRUE
#         } else {
#             snv_mnv_flag <- FALSE
#         }
#     } else {
#         snv_mnv_flag <- FALSE
#     }
#     #---Combine all mutations
#     if (indel_flag == FALSE) {
#         all_mutations_df <- snv_mnv_df
#     } else if (snv_mnv_flag == FALSE) {
#         all_mutations_df <- indel_df
#     } else {
#         all_mutations_df <- bind_rows(indel_df, snv_mnv_df)
#     }
#     tmp <- names(all_mutations_df)[which(sapply(all_mutations_df, class) == "list")]
#     if (length(tmp) > 0) {
#         for (element in tmp) {
#             all_mutations_df[[element]] <- sapply(all_mutations_df[[element]], paste, collapse = ",")
#         }
#     }
#     #---Get CN for each mutation
#     cns <- read.table(paste0("/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/DATASETS/PCAWG/consensus_cnv/consensus.20170119.somatic.cna.annotated/", sample, ".consensus.20170119.somatic.cna.annotated.txt"), header = TRUE, sep = "\t")
#     all_mutations_df$total_cn <- NA
#     all_mutations_df$major_cn <- NA
#     all_mutations_df$minor_cn <- NA
#     for (row in 1:nrow(cns)) {
#         if (is.na(cns$total_cn[row]) | is.na(cns$major_cn[row]) | is.na(cns$minor_cn[row])) next
#         chr <- cns$chromosome[row]
#         start <- cns$start[row]
#         end <- cns$end[row]
#         tmp <- which(all_mutations_df$chromosome == chr & all_mutations_df$start >= start & all_mutations_df$end <= end)
#         all_mutations_df$total_cn[tmp] <- cns$total_cn[row]
#         all_mutations_df$major_cn[tmp] <- cns$major_cn[row]
#         all_mutations_df$minor_cn[tmp] <- cns$minor_cn[row]
#     }
#     if (length(which(is.na(all_mutations_df$total_cn))) > 0) {
#         all_mutations_df <- all_mutations_df[-which(is.na(all_mutations_df$total_cn)), ]
#     }
#     all_mutations_df$karyotype <- paste0(pmax(all_mutations_df$major_cn, all_mutations_df$minor_cn), "_", pmin(all_mutations_df$major_cn, all_mutations_df$minor_cn))
#     #---Save all mutations
#     write.csv(all_mutations_df, file = paste0("/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SFS_CNA_deconvolution/data/PCAWG/", sample, "_all.csv"), row.names = FALSE)
#     #---Get CN-specific mutations
#     for (karyotype in unique(all_mutations_df$karyotype)) {
#         muts_karyotype <- all_mutations_df[all_mutations_df$karyotype == karyotype, ]
#         filename <- paste0("/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SFS_CNA_deconvolution/data/PCAWG/", sample, "_", karyotype, ".csv")
#         write.table(muts_karyotype, filename, sep = "\t", row.names = FALSE)
#     }
# }

# ==========================================LIMIT TO SPECIFIC CASE STUDY
histology <- c("Breast-AdenoCA", "Ovary-AdenoCA", "Myeloid-AML", "Myeloid-MPN", "Myeloid-MDS")
sample_df <- read.csv(paste0(R_data, "/sample_information.csv"))
sample_df <- sample_df[
    which(sample_df$histology_abbreviation %in% histology &
        sample_df$wgd_status == "no_wgd" &
        sample_df$wgd_uncertain == FALSE),
]
write.csv(sample_df, file = paste0(R_workplace, "/sample_information.csv"), row.names = FALSE)
# # ===============================================================MOBSTER
# mobster_df <- data.frame()
# i <- 0
# # for (sample in sample_df$aliquot_id) {
# for (sample in c("097a7d36-905b-72be-e050-11ac0d482c9a")) {
#     i <- i + 1
#     #   Import mutational data
#     filename_2 <- paste0(R_data, "/", sample, "_1_1.csv")
#     mutation_table <- read.csv(filename_2, sep = "\t", header = TRUE)
#     mutation_table$Ref_count <- mutation_table$t_ref_count
#     mutation_table$Alt_count <- mutation_table$t_alt_count
#     if (any(is.na(mutation_table$Ref_count)) | any(is.na(mutation_table$Alt_count))) {
#         mutation_table <- mutation_table[-which(is.na(mutation_table$Ref_count) | is.na(mutation_table$Alt_count)), ]
#     }
#     if (any(mutation_table$Ref_count == 0) | any(mutation_table$Alt_count == 0)) {
#         mutation_table <- mutation_table[-which(mutation_table$Ref_count == 0 | mutation_table$Alt_count == 0), ]
#     }
#     #   Data transformation
#     mutation_table$VAF <- mutation_table$Alt_count / (mutation_table$Alt_count + mutation_table$Ref_count)
#     mob_data <- as.data.frame(mutation_table$VAF)
#     colnames(mob_data)[1] <- "VAF"
#     #   SFS deconvolution with MOBSTER
#     fit <- mobster_fit(
#         mob_data,
#         maxIter = 1000
#     )
#     #   Find best MOBSTER model
#     mob_model <- fit$best
#     png(paste0(R_workplace, "/", sample, "_1_1_MOBSTER.png"), res = 150, width = 15, height = 7.5, units = "in")
#     print(plot(fit$best))
#     dev.off()
#     png(paste0(sample, "_1_1_MOBSTER_model_selection.png"), res = 150, width = 15, height = 7.5, units = "in")
#     print(plot_model_selection(fit))
#     dev.off()
#     #   Save the results
#     mobster_df[i, "Sample"] <- sample # id
#     mobster_df[i, "Total_N"] <- mob_model$N # total acount
#     mobster_df[i, "Tail"] <- mob_model$fit.tail # bool: if tail exists
#     mobster_df[i, "Tail_Num"] <- mob_model$N.k[[1]] # number of tail
#     mobster_df[i, "Tail_shape"] <- mob_model$shape # shape of tail
#     mobster_df[i, "Tail_scale"] <- mob_model$scale # scale of tail
#     mobster_df[i, "Kbeta_cluster"] <- mob_model$Kbeta # number of clusters
#     for (k in 1:mob_model$Kbeta) {
#         mobster_df[i, paste0("cl_num_", k)] <- mob_model$N.k[[k + 1]] # number of Beta
#         mobster_df[i, paste0("a_", k)] <- mob_model$a[[k]] # alpha of Beta
#         mobster_df[i, paste0("b_", k)] <- mob_model$b[[k]] # beta of Beta
#         mobster_df[i, paste0("p_", k)] <- mobster_df[i, paste0("a_", k)] / (mobster_df[i, paste0("a_", k)] + mobster_df[i, paste0("b_", k)])
#     }
#     #---Store the best fit
#     filename <- paste0(R_workplace, "/", sample, "_1_1_MOBSTER_parameters.txt")
#     write.table(mobster_df[i, ], file = filename, sep = "\t", row.names = FALSE, col.names = TRUE)
# }
# write.csv(mobster_df, "Parameters_mobster.csv", row.names = FALSE)
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
SFS_totalsteps <- 100 # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
matrix_binomial_sfs_stepcount <- 100
# 	Choice of ploidy, which changes the binomial rate
matrix_binomial_ploidy <- 2
#   Assumption of coverage distribution
coverage_distribution <- "sample-specific"
#   Maximum number of trials for fitting each hump count
max_trials <- 10000
#   Number of cells in the sample (for calibrating neutral mutation counts)
n_sample <- 1000000
#---Options for fitting
#   Candidates for neutral tail powers
list_neutral_powers <- seq(1, 3, by = 0.01)
# 	Candidates for where the hump frequencies are
N_SFS_positions <- 500
list_frequencies <- seq(from = 1 / N_SFS_positions, to = 1, by = 1 / N_SFS_positions)
# #---Input binomial table
# cat("\n==========================================================================================================================\n")
# cat(paste0("LOAD THE BINOMIAL TABLE...\n"))
# filename_1 <- paste0(
#     R_libPaths_binomial_table, "/Binomial_PDF_",
#     matrix_binomial_sample_size, "_",
#     r_max, "_",
#     min_variant_read, "_",
#     min_total_read, "_",
#     matrix_binomial_sfs_stepcount, "_",
#     matrix_binomial_ploidy, ".mat"
# )
# inputBinomialMatrix <- readMat(filename_1)
# matrix_binomial_PDF <- inputBinomialMatrix$matrix.binomial.PDF
#---Deconvolution for each SFS
for (sample in sample_df$aliquot_id) {
    # for (sample in c("2b40a733-7a63-4bb8-a953-95a4ee28f962")) {
    # for (sample in sample_df$aliquot_id[3]) {
    cat("\n==========================================================================================================================\n")
    cat(paste0("SFS DECONVOLUTION FOR SAMPLE ", sample, "...\n"))
    #---Input the SFS data
    filename_2 <- paste0(R_data, "/", sample, "_1_1.csv")
    mutation_table <- read.table(filename_2, sep = "\t", header = TRUE)
    mutation_table$Ref_count <- mutation_table$t_ref_count
    mutation_table$Alt_count <- mutation_table$t_alt_count
    if (any(is.na(mutation_table$Ref_count)) | any(is.na(mutation_table$Alt_count))) {
        mutation_table <- mutation_table[-which(is.na(mutation_table$Ref_count) | is.na(mutation_table$Alt_count)), ]
    }
    #---Perform SFS deconvolution
    results <- DECODE(
        mutation_table = mutation_table,
        # criterion = "BIC", # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        criterion = "ICL", # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # criterion = "ICL_MAP", # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
        min_N_humps = 1, # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        parameter_filename = paste0(R_workplace, "/", sample, "_1_1_DECONVOLUTION.txt"),
        plot_filename = paste0(R_workplace, "/", sample, "_1_1_DECONVOLUTION.png")
    )
    vec_para_best_final <- results$vec_para_best_final
    deconvolution <- results$deconvolution
    # #---Store the best fit
    # N_humps <- (length(vec_para_best_final) - 1) / 2
    # filename <- paste0(R_workplace, "/", sample, "_1_1_DECONVOLUTION_parameters.txt")

    # print(filename)

    # fileID <- file(filename, "w")
    # writeLines(paste(sprintf("%.3f", vec_para_best_final), collapse = "\t"), fileID)
    # close(fileID)
}
# deconvolution_df <- cbind(
#     data.frame(Sample = sample_IDs),
#     deconvolution_df
# )
# write.csv(deconvolution_df, "Parameters_deconvolution.csv", row.names = FALSE)
