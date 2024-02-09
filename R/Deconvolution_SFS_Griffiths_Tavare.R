SFS_Griffiths_Tavare <- function(vec_para) {
  #-----------------------------------------------------Get the parameters
  no_hump <- (length(vec_para) - 1) / 2
  para_A <- vec_para[1]
  
  para_K <- if (no_hump > 0) numeric(no_hump) else NULL
  para_P <- if (no_hump > 0) numeric(no_hump) else NULL
  
  for (i in 1:no_hump) {
    para_P[i] <- vec_para[2 * i]
    para_K[i] <- vec_para[2 * i + 1]
  }
  
  #---------------------------------------Compute the Griffiths-Tavare SFS
  vec_SFS_GT <- numeric(N_end)
  
  for (m in 2:N_end) {
    vec_SFS_GT[m] <- para_A * N_end / (m * (m - 1))
    if (no_hump > 0) {
      for (i in 1:no_hump) {
        K <- para_K[i]
        P <- para_P[i]
        #dbinom change the binopdf in matlab to calculate the binomial prob.
        vec_SFS_GT[m] <- vec_SFS_GT[m] + K * dbinom(m, N_end, P)
      }
    }
  }
  
  return(vec_SFS_GT)
}
