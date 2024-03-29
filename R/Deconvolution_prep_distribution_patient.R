prep_distribution_patient <- function(vec_totcount) {
    L <- max(vec_totcount)
    sample_coverage <- data.frame(
        total_readcount = 1:max(L, r_max),
        pdf = 0
    )
    #   Compute coverage distribution
    for (i in 1:length(vec_totcount)) {
        pos <- which(sample_coverage$total_readcount == vec_totcount[i])
        sample_coverage$pdf[pos] <- sample_coverage$pdf[pos] + 1
    }
    #   Delete coverage outside minimum and maximum number of reads
    locs <- which(sample_coverage$total_readcount < r_min | sample_coverage$total_readcount > r_max)
    if (length(locs) > 0) {
        sample_coverage <- sample_coverage[-locs, ]
    }
    if (sum(sample_coverage$pdf) == 0) {
        stop("No reads remain after imposing range of [r_min, r_max]")
    }
    #   Normalize distribution
    sample_coverage$pdf <- sample_coverage$pdf / sum(sample_coverage$pdf)
    #   Return results
    return(sample_coverage)
}
