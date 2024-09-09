# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
R_ICGC_raw_data <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/DATASETS/PCAWG"
R_ICGC_processed_data <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/data/PCAWG"
R_libPaths <- ""
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)
library(R.utils)
# ====================================EXTRACT MUTATIONAL DATA FROM PCAWG
#--------------------------------------------------Unzip mutational data
directory <- paste0(R_ICGC_raw_data, "/consensus_snv_indel/final_consensus_snv_indel_passonly_icgc.public/indel")
indel_files <- list.files(directory, pattern = "\\.vcf\\.gz$")
for (file in indel_files) {
    gunzip(paste0(R_ICGC_raw_data, "/consensus_snv_indel/final_consensus_snv_indel_passonly_icgc.public/indel/", file), remove = FALSE, overwrite = TRUE)
}
directory <- paste0(R_ICGC_raw_data, "/consensus_snv_indel/final_consensus_snv_indel_passonly_icgc.public/snv_mnv")
snv_mnv_files <- list.files(directory, pattern = "\\.vcf\\.gz$")
for (file in snv_mnv_files) {
    gunzip(paste0(R_ICGC_raw_data, "/consensus_snv_indel/final_consensus_snv_indel_passonly_icgc.public/snv_mnv/", file), remove = FALSE, overwrite = TRUE)
}
#---------------------------------------Prepare sample information table
consensus_20170218_purity_ploidy <- read.table(paste0(R_ICGC_raw_data, "/consensus_cnv/consensus.20170218.purity.ploidy.txt"), header = TRUE, sep = "\t")
pcawg_sample_sheet <- read.table(paste0(R_ICGC_raw_data, "/pcawg_sample_sheet.tsv"), header = TRUE, sep = "\t")
pcawg_specimen_histology_August2016_v9 <- as.data.frame(read_excel(paste0(R_ICGC_raw_data, "/pcawg_specimen_histology_August2016_v9.xlsx")))
sample_df <- data.frame(
    aliquot_id = unique(
        c(
            sapply(indel_files, function(x) strsplit(x, "\\.")[[1]][1]),
            sapply(snv_mnv_files, function(x) strsplit(x, "\\.")[[1]][1])
        )
    )
)
sample_df <- merge(sample_df, pcawg_sample_sheet, by = "aliquot_id")
sample_df <- merge(sample_df, pcawg_specimen_histology_August2016_v9[which(pcawg_specimen_histology_August2016_v9$specimen_library_strategy == "WGS"), ], by = "icgc_specimen_id")
sample_df <- merge(sample_df, consensus_20170218_purity_ploidy, by.x = "aliquot_id", by.y = "samplename")
write.csv(sample_df, file = paste0(R_ICGC_processed_data, "/ICGC_sample_information.csv"), row.names = FALSE)
#--------------------------------Extract mutational data for each sample
pb <- txtProgressBar(
    min = 0,
    max = length(sample_df$aliquot_id),
    style = 3,
    width = 50,
    char = "+"
)
for (i_sample in 1:length(sample_df$aliquot_id)) {
    sample <- sample_df$aliquot_id[i_sample]
    setTxtProgressBar(pb, i_sample)
    #---Get INSERTION/DELETION data
    indel_filename <- paste0(R_ICGC_raw_data, "/consensus_snv_indel/final_consensus_snv_indel_passonly_icgc.public/indel/", sample, ".consensus.20161006.somatic.indel.vcf")
    if (file.exists(indel_filename)) {
        indel_vcf <- readVcf(indel_filename)
        if (nrow(indel_vcf) > 0) {
            indel_df <- data.frame(
                chromosome = as.character(seqnames(indel_vcf)),
                start = start(ranges(indel_vcf)),
                end = end(ranges(indel_vcf)),
                width = width(ranges(indel_vcf)),
                type = "indel"
            )
            indel_df <- cbind(indel_df, info(indel_vcf))
            indel_flag <- TRUE
        } else {
            indel_flag <- FALSE
        }
    } else {
        indel_flag <- FALSE
    }
    #---Get SNV/MNV data
    snv_mnv_filename <- paste0(R_ICGC_raw_data, "/consensus_snv_indel/final_consensus_snv_indel_passonly_icgc.public/snv_mnv/", sample, ".consensus.20160830.somatic.snv_mnv.vcf")
    if (file.exists(snv_mnv_filename)) {
        snv_mnv_vcf <- readVcf(snv_mnv_filename)
        if (nrow(snv_mnv_vcf) > 0) {
            snv_mnv_df <- data.frame(
                chromosome = as.character(seqnames(snv_mnv_vcf)),
                start = start(ranges(snv_mnv_vcf)),
                end = end(ranges(snv_mnv_vcf)),
                width = width(ranges(snv_mnv_vcf)),
                type = "snv_mnv"
            )
            snv_mnv_df <- cbind(snv_mnv_df, info(snv_mnv_vcf))
            snv_mnv_flag <- TRUE
        } else {
            snv_mnv_flag <- FALSE
        }
    } else {
        snv_mnv_flag <- FALSE
    }
    #---Combine all mutations
    if (indel_flag == FALSE) {
        all_mutations_df <- snv_mnv_df
    } else if (snv_mnv_flag == FALSE) {
        all_mutations_df <- indel_df
    } else {
        all_mutations_df <- bind_rows(indel_df, snv_mnv_df)
    }
    tmp <- names(all_mutations_df)[which(sapply(all_mutations_df, class) == "list")]
    if (length(tmp) > 0) {
        for (element in tmp) {
            all_mutations_df[[element]] <- sapply(all_mutations_df[[element]], paste, collapse = ",")
        }
    }
    #---Get CN for each mutation
    cns <- read.table(paste0(R_ICGC_raw_data, "/consensus_cnv/consensus.20170119.somatic.cna.annotated/", sample, ".consensus.20170119.somatic.cna.annotated.txt"), header = TRUE, sep = "\t")
    all_mutations_df$total_cn <- NA
    all_mutations_df$major_cn <- NA
    all_mutations_df$minor_cn <- NA
    for (row in 1:nrow(cns)) {
        if (is.na(cns$total_cn[row]) | is.na(cns$major_cn[row]) | is.na(cns$minor_cn[row])) next
        chr <- cns$chromosome[row]
        start <- cns$start[row]
        end <- cns$end[row]
        tmp <- which(all_mutations_df$chromosome == chr & all_mutations_df$start >= start & all_mutations_df$end <= end)
        all_mutations_df$total_cn[tmp] <- cns$total_cn[row]
        all_mutations_df$major_cn[tmp] <- cns$major_cn[row]
        all_mutations_df$minor_cn[tmp] <- cns$minor_cn[row]
    }
    if (length(which(is.na(all_mutations_df$total_cn))) > 0) {
        all_mutations_df <- all_mutations_df[-which(is.na(all_mutations_df$total_cn)), ]
    }
    all_mutations_df$karyotype <- paste0(pmax(all_mutations_df$major_cn, all_mutations_df$minor_cn), "_", pmin(all_mutations_df$major_cn, all_mutations_df$minor_cn))
    #---Save all mutations
    write.csv(all_mutations_df, file = paste0(R_ICGC_processed_data, "/", sample, "_all.csv"), row.names = FALSE)
    #---Get CN-specific mutations
    for (karyotype in unique(all_mutations_df$karyotype)) {
        muts_karyotype <- all_mutations_df[all_mutations_df$karyotype == karyotype, ]
        filename <- paste0(R_ICGC_processed_data, "/", sample, "_", karyotype, ".csv")
        write.table(muts_karyotype, filename, sep = "\t", row.names = FALSE)
    }
}
#-----------------------------------------Extract purity & coverage data
ICGC_sample_info <- read.csv(file.path(R_ICGC_processed_data, "/ICGC_sample_information.csv"))
ICGC_purity_coverage <- list()
ICGC_purity_coverage$N_sample <- nrow(ICGC_sample_info)
for (i in 1:nrow(ICGC_sample_info)) {
    print(i)
    sample_id <- ICGC_sample_info$aliquot_id[i]
    purity <- ICGC_sample_info$purity[i]
    tmp <- read.csv(file.path(R_ICGC_processed_data, paste0(sample_id, "_all.csv")))
    tot_counts <- tmp$t_alt_count + tmp$t_ref_count
    tmp <- read.csv(file.path(R_ICGC_processed_data, paste0(sample_id, "_all.csv")))
    tot_counts <- tmp$t_alt_count + tmp$t_ref_count
    tot_counts <- tot_counts[!is.na(tot_counts)]
    unique_counts <- table(tot_counts)
    coverage <- data.frame(Read_count = as.numeric(names(unique_counts)), Frequency = as.numeric(unique_counts))
    coverage <- coverage[order(coverage$Read_count), ]
    ICGC_purity_coverage[[paste0("sample_", i)]] <- list()
    ICGC_purity_coverage[[paste0("sample_", i)]]$purity <- purity
    ICGC_purity_coverage[[paste0("sample_", i)]]$coverage <- coverage
}
save(ICGC_purity_coverage, file = paste0(R_ICGC_processed_data, "/ICGC_purity_coverage.rda"))
