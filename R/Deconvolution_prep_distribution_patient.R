prep_distribution_patient <- function(vec_totcount) {
    L <- max(vec_totcount)
    TCGA_coverage_values <- 1:L
    TCGA_coverage_PDF <- rep(0, L)
    #   Compute coverage distribution
    for (i in 1:length(vec_totcount)) {
        pos <- vec_totcount[i]
        TCGA_coverage_PDF[pos] <- TCGA_coverage_PDF[pos] + 1
    }
    TCGA_coverage_PDF <- TCGA_coverage_PDF / sum(TCGA_coverage_PDF)
    #   If the coverage is over the range
    if (sum(TCGA_coverage_PDF[1:min(r_max, length(TCGA_coverage_PDF))]) <= 0) {
        TCGA_coverage_PDF[r_max] <- 1
    }
    #   Return results
    return(TCGA_coverage_PDF)
}
