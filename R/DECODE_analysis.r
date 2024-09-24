MOBSTER_summary_statistics <- function(mobster_df, MOBSTER_result) {
    mobster_df[nrow(mobster_df) + 1, "Sample"] <- MOBSTER_result$best$description
    mobster_df[nrow(mobster_df), "Nmut"] <- MOBSTER_result$best$N
    mobster_df[nrow(mobster_df), "Tail"] <- MOBSTER_result$best$fit.tail
    mobster_df[nrow(mobster_df), "Tail_Nmut"] <- MOBSTER_result$best$N.k[[1]]
    mobster_df[nrow(mobster_df), "Tail_power"] <- MOBSTER_result$best$shape + 1
    mobster_df[nrow(mobster_df), "Tail_Pareto_shape"] <- MOBSTER_result$best$shape
    mobster_df[nrow(mobster_df), "Tail_Pareto_scale"] <- MOBSTER_result$best$scale
    mobster_df[nrow(mobster_df), "Cluster_count"] <- MOBSTER_result$best$Kbeta
    for (k in 1:MOBSTER_result$best$Kbeta) {
        mobster_df[nrow(mobster_df), paste0("Cluster_Nmut_", k)] <- MOBSTER_result$best$N.k[[k + 1]]
        mobster_df[nrow(mobster_df), paste0("Cluster_VAF_", k)] <- MOBSTER_result$best$a[[k]] / (MOBSTER_result$best$a[[k]] + MOBSTER_result$best$b[[k]])
        mobster_df[nrow(mobster_df), paste0("Cluster_Beta_a_", k)] <- MOBSTER_result$best$a[[k]]
        mobster_df[nrow(mobster_df), paste0("Cluster_Beta_b_", k)] <- MOBSTER_result$best$b[[k]]
    }
    return(mobster_df)
}

DECODE_summary_statistics <- function(decode_df, DECODE_result) {
    parameters_inference_A <- DECODE_result$final_fit$best_fit$parameters_inference_A
    parameters_inference_B <- DECODE_result$final_fit$best_fit$parameters_inference_B
    if (DECODE_result$final_fit$best_fit$tail_status) {
        Tail_Nmut_inference_A <- parameters_inference_A[1] * sum(DECODE_result$SFS_data_inference_A)
        Tail_Nmut_inference_B <- parameters_inference_B[1] * sum(DECODE_result$SFS_data_inference_B)
        Tail_power <- parameters_inference_A[2]
        Cluster_Nmut_inference_A <- parameters_inference_A[seq(3, length(parameters_inference_A), 2)] * sum(DECODE_result$SFS_data_inference_A)
        Cluster_Nmut_inference_B <- parameters_inference_B[seq(3, length(parameters_inference_B), 2)] * sum(DECODE_result$SFS_data_inference_B)
        Cluster_VAF <- parameters_inference_A[seq(4, length(parameters_inference_A), 2)]
    } else {
        Tail_Nmut_inference_A <- 0
        Tail_Nmut_inference_B <- 0
        Tail_power <- NA
        Cluster_Nmut_inference_A <- parameters_inference_A[seq(1, length(parameters_inference_A), 2)] * sum(DECODE_result$SFS_data_inference_A)
        Cluster_Nmut_inference_B <- parameters_inference_B[seq(1, length(parameters_inference_B), 2)] * sum(DECODE_result$SFS_data_inference_B)
        Cluster_VAF <- parameters_inference_A[seq(2, length(parameters_inference_A), 2)]
    }
    decode_df[nrow(decode_df) + 1, "Sample"] <- sample
    decode_df[nrow(decode_df), "Nmut"] <- nrow(DECODE_result$mutational_table)
    decode_df[nrow(decode_df), "Nmut_inference_A"] <- sum(DECODE_result$SFS_data_inference_A)
    decode_df[nrow(decode_df), "Nmut_inference_B"] <- sum(DECODE_result$SFS_data_inference_B)
    decode_df[nrow(decode_df), "Tail"] <- DECODE_result$final_fit$best_fit$tail_status
    decode_df[nrow(decode_df), "Tail_Nmut_inference_A"] <- Tail_Nmut_inference_A
    decode_df[nrow(decode_df), "Tail_Nmut_inference_B"] <- Tail_Nmut_inference_B
    decode_df[nrow(decode_df), "Tail_power"] <- Tail_power
    decode_df[nrow(decode_df), "Cluster_count"] <- length(Cluster_VAF)
    for (k in 1:(length(Cluster_VAF))) {
        decode_df[nrow(decode_df), paste0("Cluster_Nmut_inference_A_", k)] <- Cluster_Nmut_inference_A[k]
        decode_df[nrow(decode_df), paste0("Cluster_Nmut_inference_B_", k)] <- Cluster_Nmut_inference_B[k]
        decode_df[nrow(decode_df), paste0("Cluster_VAF_", k)] <- Cluster_VAF[k]
    }
    return(decode_df)
}
