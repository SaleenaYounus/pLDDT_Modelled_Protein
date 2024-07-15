###Modelled protein pLDDT calculated based on atom

install.packages("bio3d")
library(bio3d)
library(dplyr)
pdb <- read.pdb("/home/user/Downloads/model.pdb")

# Extract the B-factor column
plddt_scores <- pdb$atom$b

# Print the first few pLDDT scores to check
head(plddt_scores)

# Verify the range of pLDDT scores
range_plddt <- range(plddt_scores, na.rm = TRUE)
cat("Range of pLDDT scores:", range_plddt, "\n")

# Check if scaling is needed (if scores are in a different range, e.g., 0-1 instead of 0-100)
if (max(plddt_scores, na.rm = TRUE) <= 1) {
  cat("Scaling pLDDT scores by 100\n")
  plddt_scores <- plddt_scores * 100
}


# Calculate summary statistics
summary_stats <- data.frame(
  Mean = mean(plddt_scores, na.rm = TRUE),
  Median = median(plddt_scores, na.rm = TRUE),
  Standard_Deviation = sd(plddt_scores, na.rm = TRUE),
  Minimum = min(plddt_scores, na.rm = TRUE),
  Maximum = max(plddt_scores, na.rm = TRUE),
  # Calculate the Average pLDDT
  average_plddt <- mean(plddt_scores, na.rm = TRUE)

)
# Print the summary statistics
print(summary_stats)

# Create a data frame for ggplot2
plddt_data <- data.frame(PLDDT = plddt_scores)

# Plot histogram
ggplot(plddt_data, aes(x = PLDDT)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black") +
  labs(title = "Distribution of pLDDT Scores",
       x = "pLDDT Score",
       y = "Frequency") +
  theme_minimal()

# Plot density plot
ggplot(plddt_data, aes(x = PLDDT)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Density Plot of pLDDT Scores",
       x = "pLDDT Score",
       y = "Density") +
  theme_minimal()

# Define a threshold for low confidence
low_confidence_threshold <- 50

# Find residues with pLDDT scores below the threshold
low_confidence_residues <- which(plddt_scores < low_confidence_threshold)

# Print the indices and scores of low confidence residues
low_confidence_info <- data.frame(
  Residue_Index = low_confidence_residues,
  pLDDT_Score = plddt_scores[low_confidence_residues]
)

print(low_confidence_info)
