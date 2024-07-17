library(dplyr)

# Function to calculate Subchallenge 1A score
setwd("/Users/khanhngocdinh/Desktop/DREAM-test")

# Read the CSV files
true_cellularity <- read.csv("DREAM_sample_information.csv", stringsAsFactors = FALSE)
predicted_cellularity <- read.csv("Parameters_DREAM_DECODE.csv", stringsAsFactors = FALSE)

# Merge the two data frames based on the Sample column
merged_data <- merge(true_cellularity, predicted_cellularity, by = "Sample")

# Modify the function to calculate Subchallenge 1A score
calculate_subchallenge_1A_score <- function(p, cf1, cf2, cf3, cf4) {
    # Check for NA values first
    if (is.na(p)) {
        return(NA)
    }

    # Determine 'c' based on the availability of cluster frequencies
    c <- if (!is.na(cf4)) {
        cf4
    } else if (!is.na(cf3)) {
        cf3
    } else if (!is.na(cf2)) {
        cf2
    } else {
        cf1
    }

    if (p < 0 || p > 1 || c < 0 || c > 1) {
        stop("p (true cellularity) must be between 0 and 1 and c (predicted cellularity) must be determined and between 0 and 1.")
    }

    a <- c * 2

    score <- 1 - abs(p - a)
    return(score)
}

# Apply the modified function to each row in the merged data frame
scores_vector <- mapply(calculate_subchallenge_1A_score, merged_data$Purity, merged_data$Cluster_frequency_1, merged_data$Cluster_frequency_2, merged_data$Cluster_frequency_3, merged_data$Cluster_frequency_4)

if (length(scores_vector) != nrow(merged_data)) {
    stop("The calculated scores vector length does not match the number of rows in the merged data.")
}

merged_data$Score_1A <- scores_vector

# Select relevant columns to output
output_data <- merged_data[, c("Sample", "Score_1A")]

# Write the output to a new CSV file
write.csv(output_data, "Calculated_Scores_1A.csv", row.names = FALSE)
