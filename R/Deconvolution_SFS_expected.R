SFS_expected <- function(vec_para, package_input) {
  vec_SFS_GT <- SFS_Griffiths_Tavare(vec_para)
  
  vec_SFS_freq <- seq(0, 1, length.out = SFS_totalsteps + 1)
  vec_SFS_expected <- numeric(SFS_totalsteps)
  
  for (i in 1:SFS_totalsteps) {
    j_lower <- round(SFS_totalsteps_base * vec_SFS_freq[i]) + 1
    j_upper <- round(SFS_totalsteps_base * vec_SFS_freq[i + 1])
    omega <- 0
    
    for (r in max(r_min, 1):r_max) {
      PDF_coverage <- pdf_coverage(r) 
      
      if (PDF_coverage > 0) {
        Sum <- 0
        
        for (m in 1:N_end) {
          q_m <- vec_SFS_GT[m]
          
          if (q_m > 0) {
            PDF_success <- 0
            
            for (j in j_lower:j_upper) {
              PDF_success <- PDF_success + matrix_binomial_PDF[r, m, j]
            }
            
            Sum <- Sum + q_m * PDF_success
          }
        }
        
        omega <- omega + PDF_coverage * Sum
      }
    }
    
    vec_SFS_expected[i] <- omega
  }
  
  vec_SFS_expected
}
