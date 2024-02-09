pdf_coverage <- function(r) {
  
  # Compute the probability of a read number
  # based on choice of sampling coverage distribution
  if (option_dist_coverage == 'uniform') {
    if (r < r_min || r > r_max || r_min > r_max) {
      phi_r <- 0
    } else if (r_min == r_max) {
      if (dist_coverage_var_1 <= r_min && dist_coverage_var_2 >= r_min) {
        phi_r <- 1
      } else {
        phi_r <- 0
      }
    } else {
      phi_r <- 1 / (dist_coverage_var_2 - dist_coverage_var_1)
    }
  } else if (option_dist_coverage == 'binomial') {
    D <- dist_coverage_var_1
    if (r < r_min || r > r_max || r_min > r_max) {
      phi_r <- 0
    } else if (D > 0) {
      phi_r <- dbinom(r, size = N_end, prob = D/N_end)
    } else {
      phi_r <- -1
    }
  } else if (option_dist_coverage == 'TCGA') {
    pos <- which(TCGA_coverage_values == r)
    if (length(pos) == 0 || r < r_min || r > r_max || r_min > r_max) {
      phi_r <- 0
    } else {
      phi_r <- TCGA_coverage_PDF[pos]
    }
  }
  
  return(phi_r)
}

