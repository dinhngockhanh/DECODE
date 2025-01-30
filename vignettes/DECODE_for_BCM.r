R_data <- "/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/DATASETS/BCM-WES+scRNA/PhylExInputFinal/"
R_workplace <- "/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/DATASETS/BCM-WES+scRNA/PhylExInputFinal/"
R_libPaths <- ""
R_libPaths_extra <- "/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/R/"
.libPaths(R_libPaths)
setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)
library(ggplot2)
library(grid)



# sample_IDs <- c(52, 57, 63, 71, 104, 115, 221)
sample_IDs <- c(115)



for (sample_ID in sample_IDs) {
    bulk <- read.table(paste0("/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/DATASETS/BCM-WES+scRNA/PhylExInputFinal/bulk_", sample_ID, ".txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    colnames(bulk) <- c("ID", "WES_alt_count", "WES_tot_count", "Major_CN", "Minor_CN")
    loci <- read.table(paste0("/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/DATASETS/BCM-WES+scRNA/PhylExInputFinal/loci_", sample_ID, ".txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    colnames(loci) <- c("ID", "Chrom", "Position", "Ref", "Alt")
    data_WES <- merge(bulk, loci, by = "ID")
    data_scRNA <- read.table(paste0("/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/DATASETS/BCM-WES+scRNA/PhylExInputFinal/sc_", sample_ID, ".txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    colnames(data_scRNA) <- c("ID", "Cell", "scRNA_ref_count", "scRNA_tot_count")

    data_all <- data_WES
    data_all$WES_ref_count <- data_all$WES_tot_count - data_all$WES_alt_count
    data_all$WES_VAF <- data_all$WES_alt_count / data_all$WES_tot_count
    data_all$scRNA_ref_count <- NA
    data_all$scRNA_tot_count <- NA
    for (row in 1:nrow(data_all)) {
        data_all$scRNA_ref_count[row] <- sum(data_scRNA$scRNA_ref_count[data_scRNA$ID == data_all$ID[row]])
        data_all$scRNA_tot_count[row] <- sum(data_scRNA$scRNA_tot_count[data_scRNA$ID == data_all$ID[row]])
    }
    data_all$scRNA_alt_count <- data_all$scRNA_tot_count - data_all$scRNA_ref_count
    data_all$scRNA_VAF <- data_all$scRNA_alt_count / data_all$scRNA_tot_count
    ggplot() +
        geom_point(data = data_all, aes(x = WES_VAF, y = scRNA_VAF), alpha = 0.5) +
        geom_smooth(method = "lm", se = FALSE, color = "blue") +
        labs(
            title = NULL,
            x = "WES VAF",
            y = "scRNA VAF"
        ) +
        theme(
            text = element_text(size = 20),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5),
            legend.key.width = unit(1.5, "cm")
        )
    ggsave(filename = paste0("WES_vs_scRNA_VAF_", sample_ID, ".png"), plot = last_plot(), path = "/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/DATASETS/BCM-WES+scRNA/PhylExInputFinal/")

    ggplot(data_all, aes(x = WES_VAF)) +
        geom_histogram(binwidth = 0.01, fill = "blue", color = "black", alpha = 0.7) +
        labs(
            title = "Histogram of WES VAF",
            x = "WES VAF",
            y = "Frequency"
        ) +
        xlim(0, 1) +
        theme(
            text = element_text(size = 30),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5),
            legend.key.width = unit(1.5, "cm")
        )
    ggsave(filename = paste0("WES_VAF_histogram_", sample_ID, ".png"), plot = last_plot(), path = "/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/DATASETS/BCM-WES+scRNA/PhylExInputFinal/")
}



for (sample_ID in sample_IDs) {
    bulk <- read.table(paste0("/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/DATASETS/BCM-WES+scRNA/PhylExInputFinal/bulk_", sample_ID, ".txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    colnames(bulk) <- c("ID", "Alt_count", "Tot_count", "Major_CN", "Minor_CN")
    loci <- read.table(paste0("/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/DATASETS/BCM-WES+scRNA/PhylExInputFinal/loci_", sample_ID, ".txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    colnames(loci) <- c("ID", "Chrom", "Position", "Ref", "Alt")
    mutation_table <- merge(bulk, loci, by = "ID")


    mutation_table <- mutation_table[which(mutation_table$Major_CN == 1 & mutation_table$Minor_CN == 1), ]


    mutation_table$Ref_count <- mutation_table$Tot_count - mutation_table$Alt_count
    if (any(is.na(mutation_table$Ref_count)) | any(is.na(mutation_table$Alt_count))) mutation_table <- mutation_table[-which(is.na(mutation_table$Ref_count) | is.na(mutation_table$Alt_count)), ]
    if (any(mutation_table$Ref_count == 0) | any(mutation_table$Alt_count == 0)) mutation_table <- mutation_table[-which(mutation_table$Ref_count == 0 | mutation_table$Alt_count == 0), ]
    #---SFS deconvolution with DECODE
    DECODE_result <- DECODE(
        sample_id = sample_ID,
        mutation_table = mutation_table,
        neutral_tail = TRUE, # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        min_N_humps = 2, max_N_humps = 2 # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    )
    save(DECODE_result, file = paste0(R_data, "DECODE_", sample_ID, ".rda"))
    #---Plot DECODE deconvolution
    png(paste0(R_data, "DECODE_", sample_ID, ".png"), res = 150, width = 30, height = 15, units = "in")
    print(DECODE_plot_SFS(DECODE_result = DECODE_result))
    dev.off()
    png(paste0(R_data, "DECODE_model_selection_", sample_ID, ".png"), res = 150, width = 30, height = 15, units = "in")
    grid.draw(
        DECODE_plot_model_selection(
            DECODE_result = DECODE_result
        )
    )
    dev.off()
}



DECODE_result$final_fit$best_fit$parameters_inference_A[1] *
    sum(DECODE_result$SFS_data_inference_A) *
    sum(DECODE_result$final_fit$best_fit$component_distributions_inference_A$SFS_exact[1, ]) /
    sum(DECODE_result$final_fit$best_fit$component_distributions_inference_A$SFS_expected[1, ])

sample_ID <- sample_IDs[1]
true_VAF <- data.frame(
    Cell = 1:length(DECODE_result$final_fit$best_fit$component_distributions_inference_A$SFS_exact[1, ]),
    VAF = DECODE_result$final_fit$best_fit$component_distributions_inference_A$SFS_exact[1, ]
)
library(openxlsx)
write.xlsx(true_VAF, file = paste0(R_data, "neutral_tail_true_VAF_", sample_ID, ".xlsx"), rowNames = FALSE)
observed_VAF <- data.frame(
    Bin = paste0(
        (0:(length(DECODE_result$final_fit$best_fit$component_distributions_inference_A$SFS_expected[1, ]) - 1)) / length(DECODE_result$final_fit$best_fit$component_distributions_inference_A$SFS_expected[1, ]), "-",
        (1:length(DECODE_result$final_fit$best_fit$component_distributions_inference_A$SFS_expected[1, ])) / length(DECODE_result$final_fit$best_fit$component_distributions_inference_A$SFS_expected[1, ])
    ),
    VAF = DECODE_result$final_fit$best_fit$component_distributions_inference_A$SFS_expected[1, ]
)
write.xlsx(observed_VAF, file = paste0(R_data, "neutral_tail_observed_VAF_", sample_ID, ".xlsx"), rowNames = FALSE)
