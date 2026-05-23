# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
R_data <- "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/data/PCAWG"
R_workplace <- "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/vignettes"
R_libPaths <- ""
R_libPaths_extra <- "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DECODE/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Ginsburg
# R_data <- "/burg/iicd/users/knd2127/DECODE_DATA/PCAWG"
# R_workplace <- "/burg/iicd/users/knd2127/DECODE_TCGA_COAD"
# R_libPaths <- "/burg/iicd/users/knd2127/rpackages"
# R_libPaths_extra <- "/burg/iicd/users/knd2127/R_DECODE"
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)
setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)
library(dplyr)
library(dplyr)
library(ggplot2)
library(shadowtext)
# ===============================================SET UP WORKPLACE FOLDER
folder_workplace <- "Results_ICGC/"
if (!dir.exists(folder_workplace)) dir.create(folder_workplace)
# ===========================================GET ICGC SAMPLE INFORMATION
sample_info <- read.csv(paste0(R_data, "/ICGC_sample_information.csv"))
# samples_to_delete <- c()
# for (sample_index in 1:nrow(sample_info)) {
#     if (sample_info$wgd_status[sample_index] == "wgd" | sample_info$wgd_uncertain[sample_index] == TRUE) {
#         samples_to_delete <- c(samples_to_delete, sample)
#         next
#     }
#     sample <- sample_info$aliquot_id[sample_index]
#     filename <- paste0(R_data, "/", sample, "_1_1.csv")
#     mutation_table <- read.table(filename, sep = "\t", header = TRUE)
#     mutation_table$Ref_count <- mutation_table$t_ref_count
#     mutation_table$Alt_count <- mutation_table$t_alt_count
#     if (any(is.na(mutation_table$Ref_count)) | any(is.na(mutation_table$Alt_count))) mutation_table <- mutation_table[-which(is.na(mutation_table$Ref_count) | is.na(mutation_table$Alt_count)), ]
#     if (any(mutation_table$Ref_count == 0) | any(mutation_table$Alt_count == 0)) mutation_table <- mutation_table[-which(mutation_table$Ref_count == 0 | mutation_table$Alt_count == 0), ]
#     if (nrow(mutation_table) < 100) samples_to_delete <- c(samples_to_delete, sample)
# }
# if (length(samples_to_delete) > 0) sample_info <- sample_info[!sample_info$aliquot_id %in% samples_to_delete, ]
sample_IDs <- sample_info$aliquot_id
# # ===============================================================MOBSTER
# library(mobster)
# for (sample in sample_IDs) {
#     #---Input the SFS data
#     filename <- paste0(R_data, "/", sample, "_1_1.csv")
#     mutation_table <- read.table(filename, sep = "\t", header = TRUE)
#     mutation_table$Ref_count <- mutation_table$t_ref_count
#     mutation_table$Alt_count <- mutation_table$t_alt_count
#     if (any(is.na(mutation_table$Ref_count)) | any(is.na(mutation_table$Alt_count))) mutation_table <- mutation_table[-which(is.na(mutation_table$Ref_count) | is.na(mutation_table$Alt_count)), ]
#     if (any(mutation_table$Ref_count == 0) | any(mutation_table$Alt_count == 0)) mutation_table <- mutation_table[-which(mutation_table$Ref_count == 0 | mutation_table$Alt_count == 0), ]
#     mutation_table$VAF <- mutation_table$Alt_count / (mutation_table$Alt_count + mutation_table$Ref_count)
#     MOBSTER_data <- data.frame(VAF = mutation_table$VAF)
#     #---SFS deconvolution with MOBSTER
#     MOBSTER_result <- mobster_fit(
#         MOBSTER_data,
#         description = sample,
#         parallel = FALSE
#     )
#     save(MOBSTER_result, file = paste0(folder_workplace, "MOBSTER_", sample, ".rda"))
#     #---Plot MOBSTER deconvolution
#     png(paste0(folder_workplace, "MOBSTER_", sample, ".png"), res = 150, width = 15, height = 7.5, units = "in")
#     print(plot(MOBSTER_result$best))
#     dev.off()
#     png(paste0(folder_workplace, "MOBSTER_model_selection_", sample, ".png"), res = 150, width = 15, height = 7.5, units = "in")
#     print(plot_model_selection(MOBSTER_result))
#     dev.off()
# }
# ================================================================DECODE
library(grid)
# sample_IDs <- c(
#     "0ab4d782-9a50-48b9-96e4-6ce42b2ea034",
#     "0be08326-c623-11e3-bf01-24c6515278c0",
#     "0bfd1043-7ec1-aaec-e050-11ac0c482f39",
#     "0bfd1043-7344-fdd0-e050-11ac0c484cab",
#     "0bfd1068-3fc5-a95b-e050-11ac0c4860c3",
#     "0c0038ff-6cc4-b0b0-e050-11ac0d483d73",
#     "1a4633c4-72a0-4e30-8c4c-345e04337627",
#     "1bd47e40-d708-4ca2-b4b3-eb8d996c916b",
#     "2b4feb84-89e4-4c38-8561-5ffab02c8132",
#     "2bd9ccca-3fae-4b66-a762-6f30d6276222",
#     "3ed783cf-2248-44a1-a2a2-d6b6519b91ef",
#     "3fb8f017-576f-4901-b8bf-3a58e5d43de3",
#     "4b91ece6-c9b2-4889-b18c-c63eb58eb061",
#     "4c4aa1b1-fda3-4c5b-b588-68aa727500ad",
#     "05c487aa-72d8-42e6-aa2b-b9b5ce273f5c",
#     "7a95af21-ca7c-4596-9c83-66d11ca0c417",
#     "7ae9b843-488f-459c-8c0d-c81dcae57f99",
#     "7ae3671f-bf98-4693-8f35-3b762c9121d4",
#     "8bbe4006-be0a-4cd5-91f6-529100d4f06e",
#     "8c443bd0-7987-44de-b312-5b859e6d13a9"
# )
# sample_IDs <- c(
#     "0bfd1043-8181-e3e4-e050-11ac0c4860c5",
#     "4255582e-c622-11e3-bf01-24c6515278c0",
#     "05616329-e7ba-4efd-87b1-d79cd0f7af3d",
#     "6867811f-ac89-47da-b5dc-1270033c36e7",
#     "8282283d-247a-431d-9421-0fcc52f0a897",
#     "9078333d-73d3-496a-9fc3-a94353b7e107",
#     "9650640f-154d-4696-aa96-3611c6fcee7b",
#     "45614404-2149-4468-848c-0796e3757d62",
#     "47050918-c623-11e3-bf01-24c6515278c0",
#     "51800588-c622-11e3-bf01-24c6515278c0",
#     "61973578-4c0d-4a3f-b9c4-f96ceab24629",
#     "68956108-2606-4696-b038-462b6c432398",
#     "a3edc9cc-f54a-4459-a5d0-097879c811e5",
#     "a4ca18dc-c622-11e3-bf01-24c6515278c0",
#     "a9a240f3-d237-4bb8-b968-e4a3cc7c2633",
#     "a47c2012-c13d-48ac-88b6-e09bfd50122b",
#     "a335b03d-41ac-4d41-a2a9-3134b5b0a0a7",
#     "a3210fd0-344c-468e-8ff2-2d0869a2fb75",
#     "a6045753-60bb-4e65-bc89-1ef0b47aab35",
#     "aa4a868a-df23-4eef-a618-e945aa2ce98a"
# )
for (sample in sample_IDs) {
    # if (file.exists(paste0(folder_workplace, "DECODE_model_selection_", sample, ".png"))) next
    if (file.exists(paste0(folder_workplace, "DECODE_", sample, ".rda"))) next
    tryCatch(
        {
            #---Input the SFS data
            filename <- paste0(R_data, "/", sample, "_1_1.csv")
            mutation_table <- read.table(filename, sep = "\t", header = TRUE)
            mutation_table$Ref_count <- mutation_table$t_ref_count
            mutation_table$Alt_count <- mutation_table$t_alt_count
            if (any(is.na(mutation_table$Ref_count)) | any(is.na(mutation_table$Alt_count))) mutation_table <- mutation_table[-which(is.na(mutation_table$Ref_count) | is.na(mutation_table$Alt_count)), ]
            if (any(mutation_table$Ref_count == 0) | any(mutation_table$Alt_count == 0)) mutation_table <- mutation_table[-which(mutation_table$Ref_count == 0 | mutation_table$Alt_count == 0), ]
            #---SFS deconvolution with DECODE
            DECODE_result <- DECODE(
                sample_id = sample,
                mutation_table = mutation_table,
                n_SMCRF_particles = rep(1000, 5),
                criterion_tail_weight = 0.3
            )
            save(DECODE_result, file = paste0(folder_workplace, "DECODE_", sample, ".rda"))
            #---Plot DECODE deconvolution
            # png(paste0(folder_workplace, "DECODE_", sample, ".png"), res = 150, width = 30, height = 15, units = "in")
            # print(DECODE_plot_SFS(DECODE_result = DECODE_result))
            # dev.off()
            png(paste0(folder_workplace, "DECODE_model_selection_", sample, ".png"), res = 150, width = 30, height = 15, units = "in")
            grid.draw(
                DECODE_plot_model_selection(
                    DECODE_result = DECODE_result
                )
            )
            dev.off()
        },
        error = function(e) {
            cat("\n\n\n")
            cat("Failed for sample:     ", sample, "\n")
            cat("Error message:         ", e$message, "\n")
        }
    )
}
# # =============================================EXTRACT DECODE PARAMETERS
# decode_df <- data.frame(
#     Sample = c(),
#     Nmut = c(),
#     Tail = c(),
#     Nclusters = c()
# )
# for (sample in sample_IDs) {
#     cat(which(sample_IDs == sample), "\n")
#     if (!file.exists(paste0(folder_workplace, "DECODE_", sample, ".rda"))) next
#     load(paste0(folder_workplace, "DECODE_", sample, ".rda"))
#     Nmut <- nrow(DECODE_result$mutational_table)
#     with_tail <- DECODE_result$best_with_tail
#     Nclusters <- DECODE_result$best_N_clusters
#     parameters <- DECODE_result[[paste0("fits_", ifelse(with_tail, "with", "without"), "_tail")]]$all_fits[[paste0(Nclusters, "_clusters")]]$parameters
#     parameters_summary <- parameters %>%
#         summarise(across(where(is.numeric), list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE))))

#     decode_df[nrow(decode_df) + 1, ] <- NA
#     decode_df[nrow(decode_df), "Sample"] <- sample
#     decode_df[nrow(decode_df), "Nmut"] <- Nmut
#     decode_df[nrow(decode_df), "Tail"] <- with_tail
#     decode_df[nrow(decode_df), "Nclusters"] <- Nclusters
#     for (col in names(parameters_summary)) {
#         decode_df[nrow(decode_df), col] <- parameters_summary[[col]]
#     }
# }
# write.csv(decode_df, file = paste0(folder_workplace, "DECODE_ICGC.csv"), row.names = FALSE)
# # ============================================EXTRACT MOBSTER PARAMETERS
# mobster_df <- data.frame(
#     Sample = c(),
#     Nmut = c(),
#     Tail = c(),
#     Nclusters = c()
# )
# for (sample in sample_IDs) {
#     cat(which(sample_IDs == sample), "\n")
#     if (!file.exists(paste0(folder_workplace, "MOBSTER_", sample, ".rda"))) next
#     load(paste0(folder_workplace, "MOBSTER_", sample, ".rda"))
#     Nmut <- MOBSTER_result$best$N
#     with_tail <- MOBSTER_result$best$fit.tail
#     Nclusters <- MOBSTER_result$best$Kbeta
#     components_Nmut <- MOBSTER_result$best$N.k
#     components_Beta_a <- MOBSTER_result$best$a
#     components_Beta_b <- MOBSTER_result$best$b
#     components_VAF <- components_Beta_a / (components_Beta_a + components_Beta_b)
#     cluster_order <- names(components_VAF)[order(as.numeric(components_VAF), decreasing = TRUE)]
#     mobster_df[nrow(mobster_df) + 1, ] <- NA
#     mobster_df[nrow(mobster_df), "Sample"] <- sample
#     mobster_df[nrow(mobster_df), "Nmut"] <- Nmut
#     mobster_df[nrow(mobster_df), "Tail"] <- with_tail
#     mobster_df[nrow(mobster_df), "Nclusters"] <- Nclusters
#     if (with_tail) {
#         mobster_df[nrow(mobster_df), "Tail_power"] <- MOBSTER_result$best$shape + 1
#         mobster_df[nrow(mobster_df), "Tail_Pareto_shape"] <- MOBSTER_result$best$shape
#         mobster_df[nrow(mobster_df), "Tail_Pareto_scale"] <- MOBSTER_result$best$scale
#         mobster_df[nrow(mobster_df), "Tail_Nmut"] <- MOBSTER_result$best$N.k[[1]]
#     }
#     for (i in seq_len(Nclusters)) {
#         mobster_df[nrow(mobster_df), paste0("Cluster_", i, "_freq")] <- components_VAF[[cluster_order[i]]]
#         mobster_df[nrow(mobster_df), paste0("Cluster_", i, "_Beta_a")] <- components_Beta_a[[cluster_order[i]]]
#         mobster_df[nrow(mobster_df), paste0("Cluster_", i, "_Beta_b")] <- components_Beta_b[[cluster_order[i]]]
#         mobster_df[nrow(mobster_df), paste0("Cluster_", i, "_Nmut")] <- components_Nmut[[cluster_order[i]]]
#     }
# }
# write.csv(mobster_df, file = paste0(folder_workplace, "MOBSTER_ICGC.csv"), row.names = FALSE)
# # ========================================ANALYSES OF CLUSTERING RESULTS
# folder_workplace <- "/Users/kndinh/DINH LAB/RESULTS/2025-09-09.Yanjie Chen.DECODE for ICGC/rda_files-tail=NA/"
# decode_df <- read.csv(paste0(folder_workplace, "DECODE_ICGC.csv"))
# decode_df <- merge(decode_df, sample_info, by.x = "Sample", by.y = "aliquot_id", all.x = TRUE)
# mobster_df <- read.csv(paste0(folder_workplace, "MOBSTER_ICGC.csv"))
# mobster_df <- merge(mobster_df, sample_info, by.x = "Sample", by.y = "aliquot_id", all.x = TRUE)
# #---Predicted versus observed sample purity
# plot_predicted_vs_observed_purity <- function(df,
#                                               filename,
#                                               error_threshold = 0.05,
#                                               color_condition = "logNmut") {
#     if (color_condition == "logNmut") df$logNmut <- log10(df$Nmut + 1)
#     freq_correct <- 100 * sum(abs(df$purity - df$Predicted_purity) <= error_threshold) / nrow(df)
#     freq_over <- 100 * sum(df$Predicted_purity - df$purity > error_threshold) / nrow(df)
#     freq_under <- 100 * sum(df$purity - df$Predicted_purity > error_threshold) / nrow(df)
#     p <- ggplot() +
#         geom_point(data = df, aes(x = purity, y = Predicted_purity, color = !!sym(color_condition)), size = 10, alpha = 0.5, stroke = 2) +
#         geom_abline(intercept = error_threshold, slope = 1, color = "white", linewidth = 4) +
#         geom_abline(intercept = error_threshold, slope = 1, color = "#999999", linewidth = 2) +
#         geom_abline(intercept = -error_threshold, slope = 1, color = "white", linewidth = 4) +
#         geom_abline(intercept = -error_threshold, slope = 1, color = "#999999", linewidth = 2) +
#         geom_shadowtext(
#             aes(
#                 x = 0.05,
#                 y = 0.05,
#                 label = paste0(round(freq_correct), "%")
#             ),
#             angle = 45,
#             hjust = 0,
#             vjust = 0.5,
#             size = 20,
#             color = "#999999",
#             bg.color = "white",
#             fontface = "bold"
#         ) +
#         geom_shadowtext(
#             aes(
#                 x = 0,
#                 y = error_threshold + 0.05,
#                 label = paste0(round(freq_over), "%")
#             ),
#             angle = 45,
#             hjust = 0,
#             vjust = 0.5,
#             size = 20,
#             color = "#999999",
#             bg.color = "white",
#             fontface = "bold"
#         ) +
#         geom_shadowtext(
#             aes(
#                 x = error_threshold + 0.05,
#                 y = 0,
#                 label = paste0(round(freq_under), "%")
#             ),
#             angle = 45,
#             hjust = 0,
#             vjust = 0.5,
#             size = 20,
#             color = "#999999",
#             bg.color = "white",
#             fontface = "bold"
#         ) +
#         scale_x_continuous(name = "Sample purity", breaks = seq(0, 1, by = 0.2), limits = c(-0.05, 1.12), expand = c(0, 0)) +
#         scale_y_continuous(name = "Predicted purity", breaks = seq(0, 1, by = 0.2), limits = c(-0.05, 1.12), expand = c(0, 0)) +
#         theme_bw() +
#         theme(
#             legend.position = "top",
#             legend.justification = c(0, 1),
#             plot.title = element_blank(),
#             aspect.ratio = 1,
#             text = element_text(size = 50),
#             panel.background = element_rect(fill = "white", colour = "white"),
#             panel.grid.major = element_line(colour = "white"),
#             panel.grid.minor = element_line(colour = "white"),
#             panel.border = element_blank()
#         )
#     if (color_condition == "logNmut") {
#         p <- p +
#             scale_color_gradientn(
#                 colors = c("#1E1E1E", "#3B99B1", "#EAC728", "#F5191C"),
#                 values = scales::rescale(c(2, 3, 4, 5)),
#                 labels = scales::label_math(expr = 10^.x),
#                 guide = guide_colorbar(barwidth = 40, barheight = 1)
#             ) +
#             labs(color = "Mutation count")
#     } else if (color_condition %in% c("Nclusters", "Cluster_count")) {
#         p <- p +
#             scale_color_manual(
#                 values = c(
#                     "1" = "#D55E00",
#                     "2" = "#0072B2",
#                     "3" = "#009E73",
#                     "4" = "#CC79A7",
#                     "5" = "#E69F00",
#                     "6" = "#56B4E9"
#                 ),
#                 breaks = c(1, 2, 3, 4, 5, 6)
#             ) +
#             labs(color = "Cluster count")
#     }
#     png(filename, res = 500, width = 15, height = 15.7, units = "in", pointsize = 12)
#     print(p)
#     dev.off()
# }
# decode_df$Predicted_purity <- 2 * decode_df$Cluster_1_freq_mean
# decode_df$Nclusters <- as.factor(decode_df$Nclusters)
# mobster_df$Predicted_purity <- 2 * mobster_df$Cluster_1_freq
# mobster_df$Cluster_count <- as.factor(mobster_df$Nclusters)
# plot_predicted_vs_observed_purity(
#     df = decode_df,
#     color_condition = "logNmut",
#     filename = paste0(folder_workplace, "ICGC_DECODE_purity_by_Nmut.png")
# )
# plot_predicted_vs_observed_purity(
#     df = mobster_df,
#     color_condition = "logNmut",
#     filename = paste0(folder_workplace, "ICGC_MOBSTER_purity_by_Nmut.png")
# )
# plot_predicted_vs_observed_purity(
#     df = decode_df,
#     color_condition = "Nclusters",
#     filename = paste0(folder_workplace, "ICGC_DECODE_purity_by_Nclusters.png")
# )
# plot_predicted_vs_observed_purity(
#     df = mobster_df,
#     color_condition = "Cluster_count",
#     filename = paste0(folder_workplace, "ICGC_MOBSTER_purity_by_Nclusters.png")
# )
# #---Save results with clinical data
# write.csv(mobster_df, paste0(folder_workplace, "MOBSTER_ICGC_with_clinical_data.csv"), row.names = FALSE)
# write.csv(decode_df, paste0(folder_workplace, "DECODE_ICGC_with_clinical_data.csv"), row.names = FALSE)
