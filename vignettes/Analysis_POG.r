# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
R_POG570_raw_data <- "/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/DATASETS/POG570"
R_POG570_processed_data <- "/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/data/POG570"
R_libPaths <- ""
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)
library(readxl)
library(dplyr)
# ===================================EXTRACT MUTATIONAL DATA FROM POG570
#---------------------------------------------------Load mutational data
lines <- readLines(paste0(R_POG570_raw_data, "/cbioportal/POG570_mutations.txt"))
split_lines <- strsplit(lines, "\t")
snv_files <- do.call(rbind, split_lines)
snv_files <- as.data.frame(snv_files, stringsAsFactors = FALSE)
colnames(snv_files) <- split_lines[[1]]
snv_files <- snv_files[-1, ]
#----------------------------------------------------------Load CNV data
cnv_files <- read.table(paste0(R_POG570_raw_data, "/POG570_cnv_segments.txt"), header = TRUE, sep = "\t")
cnv_files$length <- cnv_files$end_position - cnv_files$start_position + 1
#---------------------------------------Prepare sample information table
Table_S1_Demographics <- read_excel(paste0(R_POG570_raw_data, "/Table_S1_Demographics.xlsx"))
Table_S2_Treatment <- lapply(excel_sheets(paste0(R_POG570_raw_data, "/Table_S2_Treatment.xlsx")), function(sheet) {
    read_excel(paste0(R_POG570_raw_data, "/Table_S2_Treatment.xlsx"), sheet = sheet)
})
names(Table_S2_Treatment) <- excel_sheets(paste0(R_POG570_raw_data, "/Table_S2_Treatment.xlsx"))
sample_df <- Table_S1_Demographics
colnames(sample_df) <- c(
    "Patient_ID", "Sample_ID_DNA", "Sample_ID_RNA", "Age", "Gender", "Tumour_type",
    "Histological_type", "Biopsy_site", "Biopsy_cohort", "Analysis_cohort",
    "Primary_site", "Metastatic_or_recurrence", "Tumour_content", "EGAD_ID"
)
sample_df$Treatment <- NA
sample_df$Average_CN <- NA
sample_df$Biopsy_date <- NA
sample_df$Overall_survival_days <- NA
sample_df$Alive_0_Death_1 <- NA
for (Patient_ID in sample_df$Patient_ID) {
    sub_df_Treatment <- Table_S2_Treatment$Treatment %>%
        filter(Patient_ID == !!Patient_ID) %>%
        select(-Patient_ID)
    if (nrow(sub_df_Treatment) == 0) {
        Treatment <- NA
    } else {
        Treatments <- c()
        for (i in 1:nrow(sub_df_Treatment)) {
            row <- sub_df_Treatment[i, ]
            drug_name <- row$Drug_name
            other_columns <- row[-which(names(row) == "Drug_name")]
            details <- paste(names(other_columns), other_columns, sep = "=", collapse = ",")
            Treatments <- c(Treatments, paste(drug_name, details, sep = ":"))
        }
        Treatment <- paste(Treatments, collapse = ";")
    }
    sample_df$Treatment[sample_df$Patient_ID == Patient_ID] <- Treatment
    sub_df_cnv_files <- cnv_files %>% filter(pog_id == !!Patient_ID)
    sample_df$Average_CN[sample_df$Patient_ID == Patient_ID] <- sum(sub_df_cnv_files$copy_number * sub_df_cnv_files$length) / sum(sub_df_cnv_files$length)
    sample_df$Biopsy_date[sample_df$Patient_ID == Patient_ID] <- Table_S2_Treatment$Biopsy %>%
        filter(Patient_ID == !!Patient_ID) %>%
        pull(Biopsy_date)
    sample_df$Overall_survival_days[sample_df$Patient_ID == Patient_ID] <- Table_S2_Treatment$Survival %>%
        filter(Patient_ID == !!Patient_ID) %>%
        pull(Overall_survival_days)
    sample_df$Alive_0_Death_1[sample_df$Patient_ID == Patient_ID] <- Table_S2_Treatment$Survival %>%
        filter(Patient_ID == !!Patient_ID) %>%
        pull(Alive_0_Death_1)
}
sample_df <- sample_df[order(sample_df$Patient_ID), ]
write.csv(sample_df, file = paste0(R_POG570_processed_data, "/POG570_sample_information.csv"), row.names = FALSE)
#--------------------------------Extract mutational data for each sample
pb <- txtProgressBar(
    min = 0,
    max = length(sample_df$Patient_ID),
    style = 3,
    width = 50,
    char = "+"
)
tmp_num <- 0
for (i_sample in 1:length(sample_df$Patient_ID)) {
    Patient_ID <- sample_df$Patient_ID[i_sample]
    setTxtProgressBar(pb, i_sample)
    #---Get MUTATION data
    mutation_df <- snv_files %>% filter(Tumor_Sample_Barcode == !!paste0(Patient_ID, "_T"))
    tmp_num <- tmp_num + nrow(mutation_df)
    print(tmp_num)
    #---Get CN for each mutation
    cns <- cnv_files %>% filter(pog_id == !!Patient_ID)
    mutation_df$total_cn <- NA
    for (row in 1:nrow(cns)) {
        chr <- cns$chromosome[row]
        start <- cns$start_position[row]
        end <- cns$end_position[row]
        tmp <- which(mutation_df$Chromosome == chr & mutation_df$Start_Position >= start & mutation_df$End_Position <= end)
        mutation_df$total_cn[tmp] <- cns$copy_number[row]
    }
    if (length(which(is.na(mutation_df$total_cn))) > 0) {
        mutation_df <- mutation_df[-which(is.na(mutation_df$total_cn)), ]
    }
    #---Save all mutations
    write.csv(mutation_df, file = paste0(R_POG570_processed_data, "/", Patient_ID, "_all.csv"), row.names = FALSE)
    #---Get CN-specific mutations
    for (total_cn in unique(mutation_df$total_cn)) {
        muts_karyotype <- mutation_df[mutation_df$total_cn == total_cn, ]
        filename <- paste0(R_POG570_processed_data, "/", Patient_ID, "_", total_cn, ".csv")
        write.table(muts_karyotype, filename, sep = "\t", row.names = FALSE)
    }
}
for (i_sample in 1:length(sample_df$Patient_ID)) {
    Patient_ID <- sample_df$Patient_ID[i_sample]
    filename <- paste0(R_POG570_processed_data, "/", Patient_ID, "_2.csv")
    table <- read.table(filename, header = TRUE, sep = "\t")
    table$Ref_count <- table$t_depth - table$t_alt_count
    table$Alt_count <- table$t_alt_count
    table$VAF <- table$Alt_count / table$t_depth
    png(filename = paste0(R_POG570_processed_data, "/", Patient_ID, "_2.png"))
    hist(table$VAF)
    dev.off()
}
