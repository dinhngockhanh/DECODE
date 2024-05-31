pdf_coverage <- function(r, sample_coverage) {
    # Compute the probability of a read number
    # based on choice of sampling coverage distribution
    if (coverage_distribution == "uniform") {
        if (r < r_min || r > r_max || r_min > r_max) {
            phi_r <- 0
        } else if (r_min == r_max) {
            if (coverage_variables <= r_min && dist_coverage_var_2 >= r_min) {
                phi_r <- 1
            } else {
                phi_r <- 0
            }
        } else {
            phi_r <- 1 / (dist_coverage_var_2 - coverage_variables)
        }
    } else if (coverage_distribution == "binomial") {
        D <- coverage_variables
        if (r < r_min || r > r_max || r_min > r_max) {
            phi_r <- 0
        } else if (D > 0) {
            phi_r <- dbinom(r, size = N_end, prob = D / N_end)
        } else {
            phi_r <- -1
        }
    } else if (coverage_distribution == "sample-specific") {
        phi_r <- sample_coverage$pdf[sample_coverage$total_readcount == r]
    }
    return(phi_r)
}
