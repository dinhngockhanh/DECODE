error_SFS_one_iteration <- function(vec_A_and_K, vec_p, vec_SFS_positions, library_SFS_component, vec_SFS_real) {
  #--------------------------------------Check point to make sure A, K > 0
  if (min(vec_A_and_K) < 0) {
    return(Inf)
  }
  
  #-------------------------Compute the expected SFS with input parameters
  vec_SFS_model <- compute_SFS_one_iteration(vec_A_and_K, vec_p, vec_SFS_positions, library_SFS_component)
  
  #-----------------Compute the 1-norm error between expected and real SFS
  output <- sum(abs(vec_SFS_model - vec_SFS_real)) / sum(abs(vec_SFS_model))
  
  return(output)
}
