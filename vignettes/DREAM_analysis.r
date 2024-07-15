# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
R_DREAM_raw_data <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/DATASETS/DREAM Challenge/truthFiles"
R_DREAM_processed_data <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/data/DREAM"
R_libPaths <- ""
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)
library(R.utils)
library(VariantAnnotation)
library(dplyr)
# ====================================EXTRACT MUTATIONAL DATA FROM DREAM
sample_IDs <- paste0(c(paste0("P", 1:25), paste0("S", 2:10), paste0("T", 0:16)), "-noXY")
########################################################################
########################################################################
########################################################################
######################################################   2 TETRABYTES!!!
########################################################################
########################################################################
########################################################################
#----------------------------------Unzip ground truth 2B and 3B datasets
# for (sample_ID in sample_IDs) {
#     if (file.exists(paste0(R_DREAM_raw_data, "/", sample_ID, "/", sample_ID, ".truth.2B.gz"))) {
#         gunzip(paste0(R_DREAM_raw_data, "/", sample_ID, "/", sample_ID, ".truth.2B.gz"), remove = TRUE, overwrite = TRUE)
#     }
#     if (file.exists(paste0(R_DREAM_raw_data, "/", sample_ID, "/", sample_ID, ".truth.3B.gz"))) {
#         gunzip(paste0(R_DREAM_raw_data, "/", sample_ID, "/", sample_ID, ".truth.3B.gz"), remove = TRUE, overwrite = TRUE)
#     }
#     if (file.exists(paste0(R_DREAM_raw_data, "/", sample_ID, "/", sample_ID, ".truth.2B"))) {
#         file.rename(paste0(R_DREAM_raw_data, "/", sample_ID, "/", sample_ID, ".truth.2B"), paste0(R_DREAM_raw_data, "/", sample_ID, "/", sample_ID, ".truth.2B.txt"))
#     }
#     if (file.exists(paste0(R_DREAM_raw_data, "/", sample_ID, "/", sample_ID, ".truth.3B"))) {
#         file.rename(paste0(R_DREAM_raw_data, "/", sample_ID, "/", sample_ID, ".truth.3B"), paste0(R_DREAM_raw_data, "/", sample_ID, "/", sample_ID, ".truth.3B.txt"))
#     }
# }
########################################################################
########################################################################
########################################################################
######################################################   2 TETRABYTES!!!
########################################################################
########################################################################
########################################################################
# #---------------------------------------Prepare sample information table
# sample_df <- data.frame()
# for (sample_ID in sample_IDs) {
#     sample_df[nrow(sample_df) + 1, "Sample"] <- sample_ID
#     sample_df[nrow(sample_df), "Purity"] <- read.table(paste0(R_DREAM_raw_data, "/", sample_ID, "/", sample_ID, ".truth.1A.txt"))$V1
#     sample_df[nrow(sample_df), "Ploidy"] <- read.table(paste0(R_DREAM_raw_data, "/", sample_ID, "/", sample_ID, ".cellularity_ploidy.txt"), header = TRUE)$ploidy
#     snv_filename <- paste0(R_DREAM_raw_data, "/", sample_ID, "/", sample_ID, ".truth.scoring_vcf.vcf")
#     snv_vcf <- readVcf(snv_filename)
#     sample_df[nrow(sample_df), "Mutcount"] <- nrow(snv_vcf)
#     sample_df[nrow(sample_df), "Subclone_count"] <- read.table(paste0(R_DREAM_raw_data, "/", sample_ID, "/", sample_ID, ".truth.1B.txt"))$V1
#     clone_phylo <- read.table(paste0(R_DREAM_raw_data, "/", sample_ID, "/", sample_ID, ".truth.3A.txt"))
#     clone_phylo_string <- ""
#     for (i in 1:nrow(clone_phylo)) {
#         clone_phylo_string <- paste0(clone_phylo_string, paste0(clone_phylo$V1[i], "from", clone_phylo$V2[i]))
#         if (i < nrow(clone_phylo)) clone_phylo_string <- paste0(clone_phylo_string, ",")
#     }
#     sample_df[nrow(sample_df), "Subclone_phylogeny"] <- clone_phylo_string
#     clone_info <- read.table(paste0(R_DREAM_raw_data, "/", sample_ID, "/", sample_ID, ".truth.1C.txt"))
#     for (k in 1:nrow(clone_info)) {
#         sample_df[nrow(sample_df), paste0("Clone_mutcount_", k)] <- clone_info[k, "V2"]
#         sample_df[nrow(sample_df), paste0("Clone_frequency_", k)] <- clone_info[k, "V3"]
#     }
# }
# write.csv(sample_df, file = paste0(R_DREAM_processed_data, "/DREAM_sample_information.csv"), row.names = FALSE)
# #--------------------------------Extract mutational data for each sample
# for (sample in sample_IDs) {
#     print(sample)
#     #---Get mutational data
#     snv_filename <- paste0(R_DREAM_raw_data, "/", sample, "/", sample, ".mutect.vcf")
#     snv_vcf <- readVcf(snv_filename)
#     snv_df <- data.frame(
#         chromosome = as.character(seqnames(snv_vcf)),
#         start = start(ranges(snv_vcf)),
#         end = end(ranges(snv_vcf)),
#         width = width(ranges(snv_vcf)),
#         REF = as.vector(rowRanges(snv_vcf)$REF),
#         ALT = as.character(unlist(rowRanges(snv_vcf)$ALT)),
#         QUAL = as.vector(rowRanges(snv_vcf)$QUAL),
#         FILTER = as.vector(rowRanges(snv_vcf)$FILTER)
#     )
#     snv_df <- cbind(snv_df, info(snv_vcf))
#     snv_df$mutation_ID <- rownames(snv_df)
#     snv_df <- snv_df[, c("mutation_ID", names(snv_df)[-which(names(snv_df) == "mutation_ID")])]
#     #   Get genotype data for normal and tumor samples
#     snv_genotype_df <- as.data.frame(geno(snv_vcf)$GT)
#     colnames(snv_genotype_df) <- c("GT_normal", "GT_tumor")
#     snv_genotype_df$mutation_ID <- rownames(snv_genotype_df)
#     snv_df <- snv_df %>%
#         left_join(snv_genotype_df, by = "mutation_ID")
#     #   Get allelic depths for normal and tumor samples
#     snv_allelic_depth <- geno(snv_vcf)$AD
#     snv_allelic_depth_normal <- snv_allelic_depth[, "normal"]
#     snv_allelic_depth_tumor <- snv_allelic_depth[, "tumor"]
#     snv_allelic_depth_df <- cbind(
#         do.call(rbind, lapply(snv_allelic_depth_normal, function(x) {
#             data.frame(ref_normal = x[1], alt_normal = x[2])
#         })),
#         do.call(rbind, lapply(snv_allelic_depth_tumor, function(x) {
#             data.frame(ref_tumor = x[1], alt_tumor = x[2])
#         }))
#     )
#     snv_allelic_depth_df <- snv_allelic_depth_df %>%
#         tibble::rownames_to_column(var = "mutation_ID")
#     snv_df <- snv_df %>%
#         left_join(snv_allelic_depth_df, by = "mutation_ID")
#     #   Get average base quality for normal and tumor samples
#     snv_bq_df <- as.data.frame(geno(snv_vcf)$BQ)
#     colnames(snv_bq_df) <- c("BQ_normal", "BQ_tumor")
#     snv_bq_df$mutation_ID <- rownames(snv_bq_df)
#     snv_df <- snv_df %>%
#         left_join(snv_bq_df, by = "mutation_ID")
#     snv_df$BQ_normal <- unlist(snv_df$BQ_normal)
#     snv_df$BQ_tumor <- unlist(snv_df$BQ_tumor)
#     #   Get approximate read depth for normal and tumor samples
#     snv_dp_df <- as.data.frame(geno(snv_vcf)$DP)
#     colnames(snv_dp_df) <- c("DP_normal", "DP_tumor")
#     snv_dp_df$mutation_ID <- rownames(snv_dp_df)
#     snv_df <- snv_df %>%
#         left_join(snv_dp_df, by = "mutation_ID")
#     #   Get allele fraction of alternate allele w.r.t. reference for normal and tumor samples
#     snv_fa_df <- as.data.frame(geno(snv_vcf)$FA)
#     colnames(snv_fa_df) <- c("FA_normal", "FA_tumor")
#     snv_fa_df$mutation_ID <- rownames(snv_fa_df)
#     snv_df <- snv_df %>%
#         left_join(snv_fa_df, by = "mutation_ID")
#     snv_df$FA_normal <- unlist(snv_df$FA_normal)
#     snv_df$FA_tumor <- unlist(snv_df$FA_tumor)
#     #   Get genotype quality for normal and tumor samples
#     snv_gq_df <- as.data.frame(geno(snv_vcf)$GQ)
#     colnames(snv_gq_df) <- c("GQ_normal", "GQ_tumor")
#     snv_gq_df$mutation_ID <- rownames(snv_gq_df)
#     snv_df <- snv_df %>%
#         left_join(snv_gq_df, by = "mutation_ID")
#     #   Get variant status relative to non-adjacent for normal and tumor samples
#     snv_ss_df <- as.data.frame(geno(snv_vcf)$SS)
#     colnames(snv_ss_df) <- c("SS_normal", "SS_tumor")
#     snv_ss_df$mutation_ID <- rownames(snv_ss_df)
#     snv_df <- snv_df %>%
#         left_join(snv_ss_df, by = "mutation_ID")
#     #---Get CN for each mutation
#     all_mutations_df <- snv_df
#     cns <- read.table(paste0(R_DREAM_raw_data, "/", sample, "/", sample, ".battenberg.txt"), header = TRUE, sep = "\t")
#     all_mutations_df$clonal_cn <- NA
#     all_mutations_df$total_cn <- NA
#     all_mutations_df$major_cn <- NA
#     all_mutations_df$minor_cn <- NA
#     for (row in 1:nrow(cns)) {
#         # if (is.na(cns$total_cn[row]) | is.na(cns$major_cn[row]) | is.na(cns$minor_cn[row])) next
#         chr <- cns$chr[row]
#         start <- cns$startpos[row]
#         end <- cns$endpos[row]
#         clonal_cn <- ifelse(cns$frac1_A[row] == 1, TRUE, FALSE)
#         tmp <- which(all_mutations_df$chromosome == chr & all_mutations_df$start >= start & all_mutations_df$end <= end)
#         all_mutations_df$clonal_cn[tmp] <- clonal_cn
#         all_mutations_df$total_cn[tmp] <- cns$nMaj1_A[row] + cns$nMin1_A[row]
#         all_mutations_df$major_cn[tmp] <- cns$nMaj1_A[row]
#         all_mutations_df$minor_cn[tmp] <- cns$nMin1_A[row]
#     }
#     if (length(which(is.na(all_mutations_df$total_cn))) > 0) {
#         all_mutations_df <- all_mutations_df[-which(is.na(all_mutations_df$total_cn)), ]
#     }
#     all_mutations_df$karyotype <- paste0(pmax(all_mutations_df$major_cn, all_mutations_df$minor_cn), "_", pmin(all_mutations_df$major_cn, all_mutations_df$minor_cn))
#     #---Save all mutations
#     write.csv(all_mutations_df, file = paste0(R_DREAM_processed_data, "/", sample, "_all.csv"), row.names = FALSE)
#     #---Get CN-specific mutations
#     for (karyotype in unique(all_mutations_df$karyotype)) {
#         muts_karyotype <- all_mutations_df[all_mutations_df$karyotype == karyotype, ]
#         filename <- paste0(R_DREAM_processed_data, "/", sample, "_", karyotype, ".csv")
#         write.csv(muts_karyotype, file = filename, row.names = FALSE)
#     }
# }


for (sample in sample_IDs) {
    if (!file.exists(paste0(R_DREAM_processed_data, "/", sample, "_1_1.csv"))) next
    data <- read.csv(paste0(R_DREAM_processed_data, "/", sample, "_1_1.csv"))
    data$VAF <- data$alt_tumor / (data$ref_tumor + data$alt_tumor)

    hist(data$VAF, main = "Histogram of VAF", xlab = "VAF", ylab = "Frequency")

    # Save the figure
    dev.copy(png, file = paste0(R_DREAM_processed_data, "/", sample, "_1_1.png"))
    dev.off()
}
