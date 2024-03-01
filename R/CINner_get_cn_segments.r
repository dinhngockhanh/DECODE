CINner_get_cn_segments <- function(simulation) {
    library(dplyr)
    library(data.table)
    options(scipen = 999)
    #-----------------------------------------Input the clonal evolution
    genotype_list_ploidy_chrom <- simulation$clonal_evolution$genotype_list_ploidy_chrom
    genotype_list_ploidy_block <- simulation$clonal_evolution$genotype_list_ploidy_block
    genotype_list_ploidy_allele <- simulation$clonal_evolution$genotype_list_ploidy_allele
    #---------------------------------------------------Input the sample
    sample_genotype_unique <- simulation$sample$sample_genotype_unique
    sample_genotype_unique_profile <- simulation$sample$sample_genotype_unique_profile
    all_sample_genotype <- simulation$sample$all_sample_genotype
    N_sample <- length(simulation$sample$sample_cell_ID)
    N_clones <- length(sample_genotype_unique)
    #--------------------------------------Find cell count of each clone
    sample_genotype_unique_count <- sapply(sample_genotype_unique, function(x) sum(all_sample_genotype == x))
    #------------------------------Find all clonal and subclonal CN bins
    segmented_copynumber <- list()
    pb <- txtProgressBar(
        min = 1, max = nrow(sample_genotype_unique_profile[[1]]),
        style = 3, width = 50, char = "="
    )
    for (row in 1:nrow(sample_genotype_unique_profile[[1]])) {
        setTxtProgressBar(pb, row)
        CN_block <- sample_genotype_unique_profile[[1]][row, ]
        CN_block_Cell_count <- sample_genotype_unique_count[1]
        if (N_clones > 1) {
            for (clone in 2:N_clones) {
                CN_segment <- sample_genotype_unique_profile[[clone]][row, ]
                #   Compare the CN segment from this clone against others
                CN_block_str <- apply(CN_block[, c("Min", "Maj")], 1, paste, collapse = ",")
                CN_segment_str <- paste(CN_segment[, c("Min", "Maj")], collapse = ",")
                identical_row_index <- which(CN_block_str == CN_segment_str)
                #   If same as another clone, add the cell count.
                #   Otherwise, record the new CN segment
                if (length(identical_row_index) > 0) {
                    CN_block_Cell_count[identical_row_index] <- CN_block_Cell_count[identical_row_index] + sample_genotype_unique_count[clone]
                } else {
                    CN_block <- rbind(CN_block, CN_segment)
                    CN_block_Cell_count <- c(CN_block_Cell_count, sample_genotype_unique_count[clone])
                }
            }
        }
        CN_block$Cell_count <- CN_block_Cell_count
        segmented_copynumber[[row]] <- CN_block
    }
    cat("\n")
    segmented_copynumber <- rbindlist(lapply(segmented_copynumber, as.data.table))
    #----------------------------------------Segmentation of the CN bins
    #   Segment the CN bins
    segmentation <- function(df) {
        group <- numeric(nrow(df))
        counter <- 1
        group[1] <- counter
        for (i in 2:nrow(df)) {
            if (df$chr[i] == df$chr[i - 1] & df$start[i] == (df$end[i - 1] + 1) &
                df$Min[i] == df$Min[i - 1] & df$Maj[i] == df$Maj[i - 1] &
                df$Cell_count[i] == df$Cell_count[i - 1]) {
                group[i] <- counter
            } else {
                counter <- counter + 1
                group[i] <- counter
            }
        }
        return(group)
    }
    segmented_copynumber <- segmented_copynumber %>% arrange(desc(Cell_count))
    segmented_copynumber$group <- segmentation(segmented_copynumber)
    segmented_copynumber <- segmented_copynumber %>%
        group_by(chr, Min, Maj, Cell_count, group) %>%
        summarize(start = min(start), end = max(end), .groups = "drop") %>%
        select(-group)
    #   Reorder the CN segments based on chromosomes
    segmented_copynumber <- as.data.frame(segmented_copynumber)
    chr_order <- c(as.character(1:22), "X", "Y")
    segmented_copynumber$chr <- factor(segmented_copynumber$chr, levels = chr_order)
    segmented_copynumber <- segmented_copynumber %>%
        arrange(chr, start, Cell_count) %>%
        select(chr, start, end, Maj, Min, Cell_count)
    segmented_copynumber$Clonal_frequency <- segmented_copynumber$Cell_count / N_sample
    #---------------------------Store the clonal & subclonal CN segments
    simulation$sample$segmented_copynumber <- segmented_copynumber
    return(simulation)
}
