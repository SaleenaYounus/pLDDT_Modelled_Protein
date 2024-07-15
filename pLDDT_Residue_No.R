####Modelled protein pLDDT calculated based on residue number
library(bio3d)
library(dplyr)
# Read the PDB file
pdb <- read.pdb("/home/user/Downloads/model.pdb")
# Extract pLDDT scores (B-factor column)
plddt_scores <- pdb$atom$b

# Check if scaling is needed (if scores are in the range of 0-1 instead of 0-100)
if (max(plddt_scores, na.rm = TRUE) <= 1) {
  cat("Scaling pLDDT scores by 100\n")
  plddt_scores <- plddt_scores * 100
}

# Extract residue numbers (sequence positions)
residue_numbers <- pdb$atom$resno

# Aggregate scores by residue number
residue_plddt <- tapply(plddt_scores, residue_numbers, mean, na.rm = TRUE)

# Calculate summary statistics
summary_stats <- data.frame(
  Mean = mean(residue_plddt),
  Median = median(residue_plddt),
  Standard_Deviation = sd(residue_plddt),
  Minimum = min(residue_plddt),
  Maximum = max(residue_plddt),
  Average_pLDDT = mean(residue_plddt)
)

# Print summary statistics
print(summary_stats)


# Create a plot of residue-based pLDDT scores
plot(residue_plddt, type = "l", col = "blue",
     main = "Residue-based pLDDT Scores",
     xlab = "Residue Number", ylab = "Average pLDDT Score")
# Add mean pLDDT line
abline(h = mean(residue_plddt, na.rm = TRUE), col = "red", lty = 2)




# Calculate the plot dimensions
plot_dim <- par("usr")

# Print the plot dimensions
print(plot_dim)

# Add legend on the right-middle of the plot
legend_x <- plot_dim[2] + 0.05 * diff(plot_dim[1:2])  # Adjust the x-coordinate as needed
legend_y <- mean(plot_dim[3:4])  # Y-coordinate for the middle of the plot
legend("right", legend = c("pLDDT Scores", "Mean pLDDT"),
       col = c("blue", "red"), lty = c(1, 2),
       bg = "white", box.lwd = 0.1, x.intersp = 0.5, y.intersp = 1,
       xjust = 1, yjust = 0.5, inset = c(0, 0))



#####or use this for plotting the legend

# Add legend at the top
legend("topright", legend = c("pLDDT Scores", "Mean pLDDT"),
       col = c("blue", "red"), lty = c(1, 2),
       bg = "white", box.lwd = 0.5)




