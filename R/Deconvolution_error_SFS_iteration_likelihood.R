error_SFS_one_iteration_likelihood <- function(vec_A_and_K, vec_p, vec_SFS_positions, library_SFS_component, vec_SFS_real) {
    #----------------------------------Check point to make sure A, K > 0
    if (min(vec_A_and_K) < 0) {
        return(Inf)
    }
    #----------------Compute the SFS probability distribution from model
    vec_SFS_model <- compute_SFS_one_iteration(vec_A_and_K, vec_p, vec_SFS_positions, library_SFS_component)
    vec_SFS_model_normalized <- vec_SFS_model / sum(vec_SFS_model)
    vec_SFS_model_normalized <- pmax(vec_SFS_model_normalized, 10^-10)
    #-----------------------------Compute the log-likelihood(data|model)
    loglikelihood <- sum(log10(vec_SFS_model_normalized) * vec_SFS_real)
    output <- -loglikelihood
    return(output)
}
