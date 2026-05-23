# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
R_BLCA_multimap_raw_data <- "/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/GITHUB/Multiregion-bladder/vignettes"
R_BLCA_multimap_processed_data <- "/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/data/BLCA_multimap"
R_libPaths <- ""
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)
suppressPackageStartupMessages({
    library(tidyverse)
    library(readxl)
})
patient_ids <- c("M17", "M19", "M20", "M21", "M22", "M23", "M24", "M25", "M26")
# ============================EXTRACT MUTATIONAL DATA FROM BLCA-MULTIMAP
#-------------------------------------------Get sample data for all maps
directory <- paste0(R_BLCA_multimap_raw_data, "/FINAL9MAPS_filtered_variants.Rds")
BLCA_multimap_mutational_data <- readRDS(directory)
#--------------------------Get mutational data for all filtered variants
directory <- paste0(R_BLCA_multimap_raw_data, "/FINAL9MAPS_samples_data.Rds")
BLCA_multimap_sample_data <- readRDS(directory)
#-----------------------------------------------Get CN data for all maps
directory <- paste0(R_BLCA_multimap_raw_data, "/MDA_BLD1_CNV-Illumina_OncoSNP_v2a.xlsx")
BLCA_multimap_CN_data <- read_excel(directory)
sample_id_1 <- sub("_.*", "", BLCA_multimap_CN_data$sample)
sample_id_1 <- gsub("Map24", "M24", sample_id_1)
sample_id_1 <- gsub("Map19", "M19", sample_id_1)
sample_id_2 <- sub("^[^_]*_([^_]*)_.*$", "\\1", BLCA_multimap_CN_data$sample)
BLCA_multimap_CN_data$sample_id <- paste0(sample_id_1, "_", sample_id_2)
#-----------------------------------------Add CN info to mutational data
BLCA_multimap_mutational_data$CN <- "NA"
for (row in 1:nrow(BLCA_multimap_CN_data)) {
    sample_id <- BLCA_multimap_CN_data[row, ]$sample_id
    chromosome <- BLCA_multimap_CN_data[row, ]$Chromosome
    start <- BLCA_multimap_CN_data[row, ]$StartPosition
    end <- BLCA_multimap_CN_data[row, ]$EndPosition
    CN_major <- BLCA_multimap_CN_data[row, ]$MajorCopyNumber
    CN_minor <- BLCA_multimap_CN_data[row, ]$MinorCopyNumber
    CN <- paste0(max(CN_major, CN_minor), "_", min(CN_major, CN_minor))
    BLCA_multimap_mutational_data$CN[
        which(
            BLCA_multimap_mutational_data$sample_id == sample_id &
                BLCA_multimap_mutational_data$chrom == paste0("chr", chromosome) &
                BLCA_multimap_mutational_data$pos >= start &
                BLCA_multimap_mutational_data$pos <= end
        )
    ] <- CN
}
#--------------------------Classify each mutation as local or widespread
BLCA_multimap_mutational_data <- BLCA_multimap_mutational_data |>
    left_join(
        BLCA_multimap_mutational_data |>
            group_by(patient_id, mutation_id) |>
            summarize(count = n()) |>
            mutate(category = if_else(count == 1, "local", "widespread")),
        by = c("patient_id", "mutation_id")
    )
saveRDS(BLCA_multimap_mutational_data, file = paste0(R_BLCA_multimap_processed_data, "/BLCA_multimap_mutational_data.Rds"))
#----------------------Plot sample-specific SFS with spatial information
histology_colors <- c(
    "NA" = "white",
    "NU" = "#B0BEC5",
    "LGIN" = "#78909C",
    "HGIN" = "#455A64",
    "UC" = "#263238",
    "local" = "#004586",
    "widespread" = "#FF420E"
)
for (patient_id in patient_ids) {
    #---Get mutational data for each patient
    mutational_data <- BLCA_multimap_mutational_data[BLCA_multimap_mutational_data$patient_id == patient_id, ] %>%
        mutate(sample_id_mini = sub(".*_", "", sample_id)) %>%
        mutate(x_label = substr(sample_id_mini, 1, 1), y_label = substr(sample_id_mini, 2, nchar(sample_id_mini)))
    mutational_data$x_label <- factor(mutational_data$x_label, levels = sort(unique(mutational_data$x_label)))
    mutational_data$y_label <- factor(mutational_data$y_label, levels = sort(as.numeric(unique(mutational_data$y_label)), decreasing = TRUE))
    mutational_data$category <- factor(mutational_data$category, levels = c("widespread", "local"))
    #---Get histology information for each sample
    sample_data <- expand.grid(x_label = unique(mutational_data$x_label), y_label = unique(mutational_data$y_label))
    sample_data$patient_id <- patient_id
    sample_data$sample_id <- paste0(patient_id, "_", sample_data$x_label, sample_data$y_label)
    sample_data <- merge(sample_data, BLCA_multimap_sample_data[, c("sample_id", "group4")], by = "sample_id", all.x = TRUE)
    names(sample_data)[names(sample_data) == "group4"] <- "group"
    sample_data[which(is.na(sample_data$group)), "group"] <- "NA"
    sample_data <- sample_data %>%
        left_join(mutational_data %>% count(sample_id), by = "sample_id") %>%
        rename(Nmutations = n)
    sample_data$group[which(is.na(sample_data$Nmutations))] <- "NA"
    if (patient_id == patient_ids[1]) {
        sample_df <- sample_data[which(!is.na(sample_data$Nmutations)), ]
    } else {
        sample_df <- rbind(sample_df, sample_data[which(!is.na(sample_data$Nmutations)), ])
    }
    #---Plot SFS for each patient
    pdf(paste0(R_BLCA_multimap_processed_data, "/SFS_", patient_id, ".pdf"), width = 16, height = 12)
    p <- ggplot(mutational_data, aes(x = VAF, fill = category)) +
        geom_rect(
            data = sample_data,
            aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = group),
            alpha = 0.5, inherit.aes = FALSE
        ) +
        geom_histogram(
            binwidth = 0.01, color = NA, position = "stack"
        ) +
        geom_text(
            data = sample_data,
            aes(x = Inf, y = Inf, label = paste0(group, " (n=", Nmutations, ")"), color = group),
            hjust = 1.1, vjust = 1.1, size = 5, inherit.aes = FALSE
        ) +
        facet_grid(y_label ~ x_label, drop = TRUE) +
        labs(x = NULL, y = NULL) +
        theme_minimal() +
        theme(
            plot.background = element_rect(fill = "white"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.text = element_text(size = 20)
        ) +
        scale_fill_manual(values = histology_colors) +
        scale_color_manual(values = histology_colors) +
        guides(fill = FALSE, color = FALSE)
    print(p)
    dev.off()
}
write.csv(sample_df, file = paste0(R_BLCA_multimap_processed_data, "/BLCA_multimap_sample_information.csv"), row.names = FALSE)
#-------------------------------------------------Output mutational data
for (sample_id in sample_df$sample_id) {
    mutational_data <- BLCA_multimap_mutational_data[BLCA_multimap_mutational_data$sample_id == sample_id, ]
    write.csv(mutational_data[which(mutational_data$category == "local"), ], file = paste0(R_BLCA_multimap_processed_data, "/", sample_id, "_local_all.csv"), row.names = FALSE)
    write.csv(mutational_data[which(mutational_data$category == "widespread"), ], file = paste0(R_BLCA_multimap_processed_data, "/", sample_id, "_widespread_all.csv"), row.names = FALSE)
    for (CN in unique(mutational_data$CN[which(mutational_data$category == "local")])) {
        write.csv(mutational_data[which(mutational_data$category == "local" & mutational_data$CN == CN), ], file = paste0(R_BLCA_multimap_processed_data, "/", sample_id, "_local_", CN, ".csv"), row.names = FALSE)
    }
    for (CN in unique(mutational_data$CN[which(mutational_data$category == "widespread")])) {
        write.csv(mutational_data[which(mutational_data$category == "widespread" & mutational_data$CN == CN), ], file = paste0(R_BLCA_multimap_processed_data, "/", sample_id, "_widespread_", CN, ".csv"), row.names = FALSE)
    }
}
