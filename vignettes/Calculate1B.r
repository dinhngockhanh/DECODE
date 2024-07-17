library(dplyr)

# Function to calculate Subchallenge 1A score
setwd("/Users/khanhngocdinh/Desktop/DREAM-test")

# Read the CSV files
true_subclonal_lineage <- read.csv("DREAM_sample_information.csv", stringsAsFactors = FALSE)
predicted_subclonal_lineage <- read.csv("Parameters_DREAM_DECODE.csv", stringsAsFactors = FALSE)

# Merge the two data frames based on the Sample column
merged_data <- merge(true_subclonal_lineage, predicted_subclonal_lineage, by = "Sample")

# Define the function to calculate Subchallenge 1B score
calculate_subchallenge_1B_score <- function(L, k) {
    # Check for NA values first
    if (is.na(L) || is.na(k)) {
        return(NA)
    }

    # Ensure L is at least 1
    if (L < 1) {
        stop("L (true number of subclonal lineages) must be at least 1.")
    }

    # Calculate the absolute difference d, with the constraint that d cannot exceed L + 1
    d <- min(abs(k - L), L + 1)

    # Calculate the Subchallenge 1B score
    score <- (L - d + 1) / (L + 1)
    return(score)
}

# Apply the function to each row in the merged data frame
scores_vector <- mapply(calculate_subchallenge_1B_score, merged_data$Subclone_count, merged_data$Cluster_count)

if (length(scores_vector) != nrow(merged_data)) {
    stop("The calculated scores vector length does not match the number of rows in the merged data.")
}

merged_data$Score_1B <- scores_vector

# Select relevant columns to output
output_data <- merged_data[, c("Sample", "Score_1B")]

# Write the output to a new CSV file
write.csv(output_data, "Calculated_Scores_1B.csv", row.names = FALSE)
