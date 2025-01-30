smcrf_multi_param <- function(statistics_target = NULL,
                              model,
                              bounds = NULL,
                              parameters_initial = NULL,
                              nParticles,
                              parallel,
                              save_model = TRUE,
                              save_rds = FALSE,
                              filename_rds = "ABCSMCDRF.rds",
                              splitting.rule = "CART",
                              smcrf_multi_param_results = NULL,
                              ...) {
    library(drf)
    library(matrixStats)
    library(Hmisc)
    #---Obtain information from previous SMC-DRF
    if (!is.null(smcrf_multi_param_results)) {
        SMCDRF <- smcrf_multi_param_results
        old_nIterations <- SMCDRF[["nIterations"]]
        begin_iteration <- old_nIterations + 1
        nIterations <- length(nParticles) + old_nIterations
        nParticles <- c(SMCDRF[["nParticles"]], nParticles)
        statistics_target <- SMCDRF[["statistics_target"]]
        parameters_ids <- colnames(SMCDRF[["Iteration_1"]]$parameters)
    } else {
        #---Initialize information for SMC-DRF
        SMCDRF <- list()
        parameters_ids <- colnames(parameters_initial)
        nIterations <- length(nParticles)
        begin_iteration <- 1
    }
    SMCDRF[["method"]] <- "smcrf-multi-param"
    SMCDRF[["nIterations"]] <- nIterations
    SMCDRF[["nParticles"]] <- nParticles
    SMCDRF[["statistics_target"]] <- statistics_target
    SMCDRF[["parameters_labels"]] <- data.frame(parameter = parameters_ids)
    SMCDRF[["statistics_labels"]] <- data.frame(ID = colnames(statistics_target))
    for (iteration in begin_iteration:nIterations) {
        cat(paste0(blue("ABC-SMC-DRF iteration", iteration, "* sampling parameters...\n")))
        #---Sample prior parameters for this round of iteration...
        if (!is.null(smcrf_multi_param_results)) {
            if (iteration == (old_nIterations + 1)) {
                SMCDRF[[paste0("Iteration_", (old_nIterations + 1))]] <- c()
                DRF_weights <- SMCDRF[[paste0("Iteration_", old_nIterations)]]$weights
                parameters <- SMCDRF[[paste0("Iteration_", old_nIterations)]]$parameters
            }
        }
        if (iteration == 1) {
            #   ... For iteration 1: sample from initial parameters
            parameters <- data.frame(parameters_initial[1:nParticles[iteration], ])
            colnames(parameters) <- parameters_ids
            parameters_unperturbed <- parameters
        } else {
            #   ... For later iterations:
            ifelse(iteration == (nIterations + 1), nrow <- nParticles[nIterations], nrow <- nParticles[iteration])
            Beaumont_variances <- data.frame(matrix(0, nrow = 1, ncol = 0))
            for (parameter_id in parameters_ids) {
                Beaumont_variances[[parameter_id]] <- 2 * as.numeric(var(data.frame(parameters[sample(nrow(parameters), size = 10000, prob = DRF_weights[, 1], replace = T), parameter_id])))
            }
            parameters_unperturbed <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids)))
            colnames(parameters_unperturbed) <- parameters_ids
            parameters_next <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids)))
            colnames(parameters_next) <- parameters_ids
            invalid_indices <- 1:nrow
            while (length(invalid_indices) > 0) {
                #   Sample parameters from previous posterior distribution
                parameter_replace <- data.frame(parameters[sample(nrow(parameters), size = length(invalid_indices), prob = DRF_weights[, 1], replace = T), ])
                colnames(parameter_replace) <- parameters_ids
                if (length(invalid_indices) == nrow) parameters_unperturbed <- parameter_replace
                #   Perturb parameters
                parameter_replace <- perturb(
                    parameters = parameter_replace,
                    Beaumont_variances = Beaumont_variances,
                    parameters_ids = parameters_ids
                )
                #   Check if parameters are within bounds, otherwise redo the failed parameters
                parameters_next[invalid_indices, ] <- parameter_replace
                if (is.null(bounds)) {
                    invalid_indices <- which(apply(parameters_next, 1, function(x) any(is.na(x))))
                } else {
                    invalid_indices <- c()
                    for (parameter_id in parameters_ids) {
                        invalid_indices <- union(invalid_indices, which(
                            is.na(parameters_next[, parameter_id]) |
                                parameters_next[, parameter_id] < bounds$min[which(bounds$parameter == parameter_id)] |
                                parameters_next[, parameter_id] > bounds$max[which(bounds$parameter == parameter_id)]
                        ))
                    }
                }
            }
            parameters <- parameters_next
        }
        #---Simulate statistics
        cat(paste0(blue("ABC-SMC-DRF iteration", iteration, "* computing SFS...\n")))
        reference <- model(parameters = parameters)
        statistics <- data.frame(reference[, colnames(reference)[!colnames(reference) %in% parameters_ids]])
        parameters <- data.frame(reference[, colnames(reference)[colnames(reference) %in% parameters_ids]])
        colnames(parameters) <- parameters_ids
        colnames(statistics) <- colnames(reference)[!colnames(reference) %in% parameters_ids]
        #---Run DRF for all parameter
        cat(paste0(blue("ABC-SMC-DRF iteration", iteration, "* predicting posterior distributions...\n")))
        Xdrf <- statistics
        Ydrf <- reference[, parameters_ids]
        drfmodel <- drf(Xdrf, Ydrf, splitting.rule = splitting.rule, ...)
        def_pred <- predict(drfmodel, statistics_target)
        DRF_weights <- as.vector(get_sample_weights(drfmodel, statistics_target))
        DRF_weights <- data.frame(matrix(rep(DRF_weights, length(parameters_ids)), ncol = length(parameters_ids)))
        colnames(DRF_weights) <- parameters_ids
        ################################################################
        ################################################################
        ################################################################
        ################################################################
        ################################################################
        cat(red("+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"))
        cat(red(paste0("+++++++++++++++++++++++++++++++++++++++++++ Iteration ", iteration, "\n")))
        cat(red("+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"))
        reference_tmp <- reference[sample(x = 1:nrow(reference), prob = as.numeric(DRF_weights[, "neutral_power"]), size = 10000, replace = TRUE), ]
        cat(paste0(red("Neutral component power:        "), yellow(paste0(format(round(mean(reference_tmp[, "neutral_power"]), 3), nsmall = 3), " \u00B1 ", format(round(sd(reference_tmp[, "neutral_power"]), 3), nsmall = 3))), "\n"))
        cat("\n")
        N_humps <- length(grep("cluster_frequency_", colnames(reference)))
        for (i in 1:N_humps) {
            cat(paste0(red(paste0("Cluster ", i, " frequency:            ")), yellow(paste0(format(round(mean(reference_tmp[, paste0("cluster_frequency_", i)]), 3), nsmall = 3), " \u00B1 ", format(round(sd(reference_tmp[, paste0("cluster_frequency_", i)]), 3), nsmall = 3))), "\n"))
        }
        cat("\n")
        cat(paste0(red("Neutral component proportion:   "), yellow(paste0(format(round(mean(reference_tmp[, "cluster_proportion_inference_A_1"]), 3), nsmall = 3), " \u00B1 ", format(round(sd(reference_tmp[, "cluster_proportion_inference_A_1"]), 3), nsmall = 3))), "\n"))
        for (i in 1:N_humps) {
            cat(paste0(red(paste0("Cluster ", i, " proportion:           ")), yellow(paste0(format(round(mean(reference_tmp[, paste0("cluster_proportion_inference_A_", i + 1)]), 3), nsmall = 3), " \u00B1 ", format(round(sd(reference_tmp[, paste0("cluster_proportion_inference_A_", i + 1)]), 3), nsmall = 3))), "\n"))
        }
        ################################################################
        ################################################################
        ################################################################
        ################################################################
        ################################################################
        # ---Save SMC-DRF results from this iteration
        SMCDRF_iteration <- list()
        SMCDRF_iteration$reference <- reference
        SMCDRF_iteration$parameters <- parameters
        SMCDRF_iteration$parameters_unperturbed <- parameters_unperturbed
        SMCDRF_iteration$statistics <- statistics
        SMCDRF_iteration$weights <- DRF_weights
        if (save_model == TRUE) {
            SMCDRF_iteration$rf_model <- drfmodel
            SMCDRF_iteration$rf_predict <- def_pred
        }
        SMCDRF[[paste0("Iteration_", iteration)]] <- SMCDRF_iteration
        if (save_rds == TRUE) {
            saveRDS(SMCDRF, file = filename_rds)
        }
    }
    return(SMCDRF)
}
smcrf_single_param <- function(statistics_target = NULL,
                               model,
                               bounds = NULL,
                               parameters_initial = NULL,
                               nParticles,
                               parallel,
                               save_model = TRUE,
                               save_rds = FALSE,
                               filename_rds = "ABCSMCRF.rds",
                               smcrf_single_param_results = NULL,
                               ...) {
    library(abcrf)
    library(Hmisc)
    # nSimulations <<- 0
    # cat("\n\n==============================================================================================================================================\n")
    #---Obtain information from previous SMC-RF
    if (!is.null(smcrf_single_param_results)) {
        SMCRF <- smcrf_single_param_results
        old_nIterations <- SMCRF[["nIterations"]]
        begin_iteration <- old_nIterations + 1
        nIterations <- length(nParticles) + old_nIterations
        nParticles <- c(SMCRF[["nParticles"]], nParticles)
        statistics_target <- SMCRF[["statistics_target"]]
        parameters_ids <- colnames(SMCRF[["Iteration_1"]]$parameters)
    } else {
        #---Initialize information for SMC-RF
        SMCRF <- list()
        parameters_ids <- colnames(parameters_initial)
        nIterations <- length(nParticles)
        begin_iteration <- 1
    }
    SMCRF[["method"]] <- "smcrf-single-param"
    SMCRF[["nIterations"]] <- nIterations
    SMCRF[["nParticles"]] <- nParticles
    SMCRF[["statistics_target"]] <- statistics_target
    SMCRF[["parameters_labels"]] <- data.frame(parameter = parameters_ids)
    SMCRF[["statistics_labels"]] <- data.frame(ID = colnames(statistics_target))
    for (iteration in begin_iteration:nIterations) {
        cat(paste0(blue("ABC-SMC-RF iteration", iteration, "* sampling parameters...\n")))
        #---Sample prior parameters for this round of iteration...
        if (!is.null(smcrf_single_param_results)) {
            if (iteration == (old_nIterations + 1)) {
                SMCRF[[paste0("Iteration_", (old_nIterations + 1))]] <- c()
                ABCRF_weights <- SMCRF[[paste0("Iteration_", old_nIterations)]]$weights
                parameters <- SMCRF[[paste0("Iteration_", old_nIterations)]]$parameters
            }
        }
        if (iteration == 1) {
            #   ... For iteration 1: sample from initial parameters
            parameters <- data.frame(parameters_initial[1:nParticles[iteration], ])
            colnames(parameters) <- parameters_ids
            parameters_unperturbed <- parameters
        } else {
            #   ... For later iterations:
            ifelse(iteration == (nIterations + 1), nrow <- nParticles[nIterations], nrow <- nParticles[iteration])
            Beaumont_variances <- data.frame(matrix(0, nrow = 1, ncol = 0))
            for (parameter_id in parameters_ids) {
                Beaumont_variances[[parameter_id]] <- 2 * as.numeric(var(data.frame(parameters[sample(nrow(parameters), size = 10000, prob = ABCRF_weights[, parameter_id], replace = T), parameter_id])))
            }
            parameters_unperturbed <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids)))
            colnames(parameters_unperturbed) <- parameters_ids
            parameters_next <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids)))
            colnames(parameters_next) <- parameters_ids
            invalid_indices <- 1:nrow
            while (length(invalid_indices) > 0) {
                #   Sample parameters from previous posterior distribution
                parameter_replace <- data.frame(matrix(NA, nrow = length(invalid_indices), ncol = 0))
                for (parameter_id in parameters_ids) {
                    parameter_replace[[parameter_id]] <- as.numeric(sample(parameters[, parameter_id], size = length(invalid_indices), prob = ABCRF_weights[, parameter_id], replace = TRUE))
                }
                if (length(invalid_indices) == nrow) parameters_unperturbed <- parameter_replace
                #   Perturb parameters
                parameter_replace <- perturb(
                    parameters = parameter_replace,
                    Beaumont_variances = Beaumont_variances,
                    parameters_ids = parameters_ids
                )
                #   Check if parameters are within bounds, otherwise redo the failed parameters
                parameters_next[invalid_indices, ] <- parameter_replace
                if (is.null(bounds)) {
                    invalid_indices <- which(apply(parameters_next, 1, function(x) any(is.na(x))))
                } else {
                    invalid_indices <- c()
                    for (parameter_id in parameters_ids) {
                        invalid_indices <- union(invalid_indices, which(
                            is.na(parameters_next[, parameter_id]) |
                                parameters_next[, parameter_id] < bounds$min[which(bounds$parameter == parameter_id)] |
                                parameters_next[, parameter_id] > bounds$max[which(bounds$parameter == parameter_id)]
                        ))
                    }
                }
            }
            parameters <- parameters_next
        }
        #---Simulate statistics
        cat(paste0(blue("ABC-SMC-RF iteration", iteration, "* computing SFS...\n")))
        reference <- model(parameters = parameters)
        statistics <- data.frame(reference[, colnames(reference)[!colnames(reference) %in% parameters_ids]])
        parameters <- data.frame(reference[, colnames(reference)[colnames(reference) %in% parameters_ids]])
        colnames(parameters) <- parameters_ids
        colnames(statistics) <- colnames(reference)[!colnames(reference) %in% parameters_ids]
        #---Run ABCRF for each parameter
        cat(paste0(blue("ABC-SMC-RF iteration", iteration, "* predicting posterior distributions...\n")))
        ABCRF_weights <- data.frame(matrix(NA, nrow = nParticles[iteration], ncol = 0))
        RFmodels <- list()
        posterior_gamma_RFs <- list()
        for (parameter_id in parameters_ids) {
            mini_reference <- reference[, c(parameter_id, colnames(reference)[!colnames(reference) %in% parameters_ids])]
            colnames(mini_reference)[1] <- "para"
            f <- as.formula("para ~.")
            RFmodel <- regAbcrf(
                formula = f,
                data = mini_reference,
                paral = parallel,
                ...
            )
            posterior_gamma_RF <- predict(
                object = RFmodel,
                obs = statistics_target,
                training = mini_reference,
                paral = parallel,
                rf.weights = T
            )
            ABCRF_weights[, parameter_id] <- posterior_gamma_RF$weights
            if (save_model == TRUE) {
                RFmodels[[parameter_id]] <- RFmodel
                posterior_gamma_RFs[[parameter_id]] <- posterior_gamma_RF
            }
        }
        ################################################################
        ################################################################
        ################################################################
        ################################################################
        ################################################################
        cat(red("+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"))
        cat(red(paste0("+++++++++++++++++++++++++++++++++++++++++++ Iteration ", iteration, "\n")))
        cat(red("+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"))
        neutral_power_sample <- sample(x = as.numeric(reference[, "neutral_power"]), prob = as.numeric(ABCRF_weights[, "neutral_power"]), size = 10000, replace = TRUE)
        cat(paste0(red("Neutral component power:        "), yellow(paste0(format(round(mean(neutral_power_sample), 3), nsmall = 3), " \u00B1 ", format(round(sd(neutral_power_sample), 3), nsmall = 3))), "\n"))
        cat("\n")
        N_humps <- length(grep("cluster_frequency_", colnames(reference)))
        for (i in 1:N_humps) {
            cluster_frequency_sample <- sample(x = as.numeric(reference[, paste0("cluster_frequency_", i)]), prob = as.numeric(ABCRF_weights[, paste0("cluster_frequency_", i)]), size = 10000, replace = TRUE)
            cat(paste0(red(paste0("Cluster ", i, " frequency:            ")), yellow(paste0(format(round(mean(cluster_frequency_sample), 3), nsmall = 3), " \u00B1 ", format(round(sd(cluster_frequency_sample), 3), nsmall = 3))), "\n"))
        }
        cat("\n")
        neutral_proportion_sample <- sample(x = as.numeric(reference[, "cluster_proportion_inference_A_1"]), prob = as.numeric(ABCRF_weights[, "cluster_proportion_inference_A_1"]), size = 10000, replace = TRUE)
        cat(paste0(red("Neutral component proportion:   "), yellow(paste0(format(round(mean(neutral_proportion_sample), 3), nsmall = 3), " \u00B1 ", format(round(sd(neutral_proportion_sample), 3), nsmall = 3))), "\n"))
        for (i in 1:N_humps) {
            cluster_proportion_sample <- sample(x = as.numeric(reference[, paste0("cluster_proportion_inference_A_", i + 1)]), prob = as.numeric(ABCRF_weights[, paste0("cluster_proportion_inference_A_", i + 1)]), size = 10000, replace = TRUE)
            cat(paste0(red(paste0("Cluster ", i, " proportion:           ")), yellow(paste0(format(round(mean(cluster_proportion_sample), 3), nsmall = 3), " \u00B1 ", format(round(sd(cluster_proportion_sample), 3), nsmall = 3))), "\n"))
        }
        ################################################################
        ################################################################
        ################################################################
        ################################################################
        ################################################################
        #---Save SMC-RF results from this iteration
        SMCRF_iteration <- list()
        SMCRF_iteration$reference <- reference
        SMCRF_iteration$parameters <- parameters
        SMCRF_iteration$parameters_unperturbed <- parameters_unperturbed
        SMCRF_iteration$statistics <- statistics
        SMCRF_iteration$weights <- ABCRF_weights
        if (save_model == TRUE) {
            SMCRF_iteration$rf_model <- RFmodels
            SMCRF_iteration$rf_predict <- posterior_gamma_RFs
        }
        SMCRF[[paste0("Iteration_", iteration)]] <- SMCRF_iteration
        if (save_rds == TRUE) {
            saveRDS(SMCRF, file = filename_rds)
        }
    }
    return(SMCRF)
}
perturb <- function(parameters, Beaumont_variances, parameters_ids) {
    library(gtools)
    # for (parameter_id in parameters_ids) {
    #     parameters[[parameter_id]] <- rnorm(
    #         n = nrow(parameters),
    #         mean = parameters[[parameter_id]],
    #         sd = sqrt(Beaumont_variances[[parameter_id]])
    #     )
    # }
    #---Perturb parameters using Beaumont method for non-pi columns
    non_tmp <- grep("cluster_proportion_inference_A_", colnames(parameters), invert = TRUE)
    non_pi_ids <- parameters_ids[non_tmp]
    for (parameter_id in non_pi_ids) {
        parameters[[parameter_id]] <- rnorm(
            n = nrow(parameters),
            mean = parameters[[parameter_id]],
            sd = sqrt(Beaumont_variances[[parameter_id]])
        )
    }
    #---Perturb parameters using Beaumont method for pi columns
    tmp <- grep("cluster_proportion_inference_A_", colnames(parameters))
    #   Find the cluster proportion with the highest Beaumont variance
    scale_id <- colnames(Beaumont_variances)[tmp][which.max(Beaumont_variances[tmp])]
    pi_ids <- parameters_ids[tmp]
    if (length(tmp) > 0) {
        parameters[pi_ids] <- t(apply(parameters[pi_ids], 1, function(alpha) {
            #   Find Dirichlet scale so that the cluster proportion matches Beaumont variance
            scale <- alpha[scale_id] * (1 - alpha[scale_id]) / Beaumont_variances[scale_id] - 1
            alpha <- as.numeric(scale) * alpha
            new_alpha <- rdirichlet(1, as.numeric(alpha))
            return(new_alpha)
            # alpha <- alpha / sum(alpha) + Beaumont_variances[pi_ids]
            # return(rdirichlet(1, as.numeric(alpha)))
        }))
    } else {
        message("No matching columns found for 'cluster_proportion_inference_A_'")
    }
    # #---Normalize cluster proportions
    # # tmp <- grep("cluster_proportion_inference_A_", colnames(parameters))
    # # parameters[, tmp] <- t(apply(parameters[, tmp], 1, function(x) x / sum(x)))
    # if (length(tmp) > 0) {
    #     parameters[, tmp] <- t(apply(parameters[, tmp], 1, function(x) x - 1 / length(tmp)))
    # } else {
    #     message("No matching columns found for 'cluster_proportion_inference_A_'")
    # }
    #---Order clusters by decreasing frequency
    nClusters <- max(as.numeric(sub("cluster_frequency_", "", grep("cluster_frequency_", colnames(parameters), value = TRUE))))
    if (nClusters > 1) {
        cluster_frequency_cols <- grep("cluster_frequency_", colnames(parameters), value = TRUE)
        cluster_proportion_cols <- grep("cluster_proportion_inference_A_", colnames(parameters), value = TRUE)
        sort_columns <- function(cols, prefix) {
            numeric_part <- as.numeric(sub(prefix, "", cols))
            sorted_indices <- order(numeric_part)
            return(cols[sorted_indices])
        }
        cluster_frequency_cols <- sort_columns(cluster_frequency_cols, "cluster_frequency_")
        cluster_proportion_cols <- sort_columns(cluster_proportion_cols, "cluster_proportion_inference_A_")
        if ("neutral_power" %in% colnames(parameters)) cluster_proportion_cols <- cluster_proportion_cols[cluster_proportion_cols != "cluster_proportion_inference_A_1"]
        reorder_clusters <- function(row) {
            cluster_frequencies <- row[cluster_frequency_cols]
            cluster_proportions <- row[cluster_proportion_cols]
            order_decreasing <- order(cluster_frequencies, decreasing = TRUE)
            row[cluster_frequency_cols] <- cluster_frequencies[order_decreasing]
            row[cluster_proportion_cols] <- cluster_proportions[order_decreasing]
            return(row)
        }
        parameters <- t(apply(parameters, 1, reorder_clusters))
        parameters <- as.data.frame(parameters)
    }
    return(parameters)
}

perturb_KHANH <- function(parameters, Beaumont_variances, parameters_ids) {
    #---Perturb parameters using Beaumont method
    for (parameter_id in parameters_ids) {
        parameters[[parameter_id]] <- rnorm(
            n = nrow(parameters),
            mean = parameters[[parameter_id]],
            sd = sqrt(Beaumont_variances[[parameter_id]])
        )
    }
    #---Normalize cluster proportions
    tmp <- grep("cluster_proportion_inference_A_", colnames(parameters))
    parameters[, tmp] <- t(apply(parameters[, tmp], 1, function(x) x / sum(x)))
    #---Order clusters by decreasing frequency
    nClusters <- max(as.numeric(sub("cluster_frequency_", "", grep("cluster_frequency_", colnames(parameters), value = TRUE))))
    if (nClusters > 1) {
        cluster_frequency_cols <- grep("cluster_frequency_", colnames(parameters), value = TRUE)
        cluster_proportion_cols <- grep("cluster_proportion_inference_A_", colnames(parameters), value = TRUE)
        sort_columns <- function(cols, prefix) {
            numeric_part <- as.numeric(sub(prefix, "", cols))
            sorted_indices <- order(numeric_part)
            return(cols[sorted_indices])
        }
        cluster_frequency_cols <- sort_columns(cluster_frequency_cols, "cluster_frequency_")
        cluster_proportion_cols <- sort_columns(cluster_proportion_cols, "cluster_proportion_inference_A_")
        if ("neutral_power" %in% colnames(parameters)) cluster_proportion_cols <- cluster_proportion_cols[cluster_proportion_cols != "cluster_proportion_inference_A_1"]
        reorder_clusters <- function(row) {
            cluster_frequencies <- row[cluster_frequency_cols]
            cluster_proportions <- row[cluster_proportion_cols]
            order_decreasing <- order(cluster_frequencies, decreasing = TRUE)
            row[cluster_frequency_cols] <- cluster_frequencies[order_decreasing]
            row[cluster_proportion_cols] <- cluster_proportions[order_decreasing]
            return(row)
        }
        parameters <- t(apply(parameters, 1, reorder_clusters))
        parameters <- as.data.frame(parameters)
    }
    return(parameters)
}
