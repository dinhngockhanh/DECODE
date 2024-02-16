fit_SFS_one_iteration_likelihood <- function(vec_SFS_real, N_humps, vec_SFS_positions, library_SFS_component) {
    # 	Choose the hump locations at random
    vec_p <- sort(sample(vec_SFS_positions, N_humps))
    # 	Define the fit function
    func_fit <- function(vec_A_and_K) {
        error_SFS_one_iteration_likelihood(vec_A_and_K, vec_p, vec_SFS_positions, library_SFS_component, vec_SFS_real)
    }
    # 	Initial values for A and K's
    vec_A_and_K_initial <- c(1, rep(100, N_humps))
    # 	Optimization using optim function
    optim_results <- optim(
        par = vec_A_and_K_initial, fn = func_fit,
        method = "Nelder-Mead",
        control = list(maxit = 100000),
    )
    # 	Extract the results
    vec_A_and_K <- optim_results$par
    log_L <- optim_results$value
    # 	Prepare the parameters to be returned
    vec_para <- numeric(2 * N_humps + 1)
    vec_para[1] <- vec_A_and_K[1]
    for (i in 1:N_humps) {
        vec_para[2 * i] <- vec_p[i]
        vec_para[2 * i + 1] <- vec_A_and_K[i + 1]
    }
    output <- list()
    output$log_L <- log_L
    output$parameters <- vec_para
    return(output)
}
