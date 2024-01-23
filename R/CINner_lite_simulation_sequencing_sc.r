simulation_sequencing_sc <- function(simulation = list(),
                                     sc_rate_false_positive = 0,
                                     sc_rate_false_negative = 0,
                                     sc_rate_unknown = 0) {
    sample_mutational_table_full <- simulation$sample_mutational_table_truth
    sample_mutational_table_id <- sample_mutational_table_full[, 1:3]
    sample_mutational_table_sc <- sample_mutational_table_full[, 4:ncol(sample_mutational_table_full)]
    sample_mutational_table_sc_mat <- data.matrix(sample_mutational_table_sc)
    #-------------------------Simulate False Positive and False Negative
    #   Find all locations of true negatives
    mat_true_neg <- which(sample_mutational_table_sc_mat == 0, arr.ind = TRUE)
    row_true_neg <- mat_true_neg[, 1]
    col_true_neg <- mat_true_neg[, 2]
    #   Find all locations of true positives
    mat_true_pos <- which(sample_mutational_table_sc_mat == 1, arr.ind = TRUE)
    row_true_pos <- mat_true_pos[, 1]
    col_true_pos <- mat_true_pos[, 2]
    #   Simmulate the false positives
    n_false_positive <- rbinom(
        n = 1,
        size = length(row_true_neg),
        prob = sc_rate_false_positive
    )
    loc_false_pos <- sample(
        x = length(row_true_neg),
        size = n_false_positive,
        replace = FALSE
    )
    row_false_pos <- row_true_neg[loc_false_pos]
    col_false_pos <- col_true_neg[loc_false_pos]
    ind_false_pos <- row_false_pos +
        nrow(sample_mutational_table_sc_mat) * (col_false_pos - 1)
    sample_mutational_table_sc_mat[ind_false_pos] <- 1
    #   Simulate the false negatives
    n_false_negative <- rbinom(
        n = 1,
        size = length(row_true_pos),
        prob = sc_rate_false_negative
    )
    loc_false_neg <- sample(
        x = length(row_true_pos),
        size = n_false_negative,
        replace = FALSE
    )
    row_false_neg <- row_true_pos[loc_false_neg]
    col_false_neg <- col_true_pos[loc_false_neg]
    ind_false_neg <- row_false_neg +
        nrow(sample_mutational_table_sc_mat) * (col_false_neg - 1)
    sample_mutational_table_sc_mat[ind_false_neg] <- 0
    #---------------------------------------------------Simulate Unknown
    mat_all <- which(sample_mutational_table_sc_mat >= 0, arr.ind = TRUE)
    row_all <- mat_all[, 1]
    col_all <- mat_all[, 2]
    #   Simulate the unknowns
    n_unknown <- rbinom(
        n = 1,
        size = length(row_all),
        prob = sc_rate_unknown
    )
    loc_unknown <- sample(
        x = length(row_all),
        size = n_unknown,
        replace = FALSE
    )
    row_unknown <- row_all[loc_unknown]
    col_unknown <- col_all[loc_unknown]
    ind_unknown <- row_unknown +
        nrow(sample_mutational_table_sc_mat) * (col_unknown - 1)
    sample_mutational_table_sc_mat[ind_unknown] <- NA
    #---------------------------------------------Prepare output package
    sample_mutational_table_sc_tmp <- data.frame(sample_mutational_table_sc_mat)
    colnames(sample_mutational_table_sc_tmp) <- colnames(sample_mutational_table_sc)
    sample_mutational_table_sc <- cbind(sample_mutational_table_id, sample_mutational_table_sc_tmp)
    simulation$sample_mutational_table_sc <- sample_mutational_table_sc
    return(simulation)
}
