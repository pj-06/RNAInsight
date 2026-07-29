############################################################
# RNAInsight
# Step 3 : Exploratory Data Analysis (EDA)
############################################################

library(ggplot2)
library(dplyr)
library(corrplot)
setwd("C:/Users/Admin/Desktop/PDS/RNAInsight")
############################################################
# Load Dataset
############################################################

data <- read.csv(
  "C:/Users/Admin/Desktop/PDS/RNAInsight/data/processed/cleaned_dataset.csv",
  stringsAsFactors = FALSE
)

# Convert Target Family to factor
data$Target_Family <- as.factor(data$Target_Family)

############################################################
# Create Output Folders
############################################################

dir.create("results", showWarnings = FALSE)
dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

############################################################
# Dataset Overview
############################################################

cat("Rows :", nrow(data), "\n")
cat("Columns :", ncol(data), "\n\n")

summary(data)

############################################################
# Target Family Distribution
############################################################

p1 <- ggplot(data, aes(x = Target_Family)) +
  geom_bar(fill = "steelblue") +
  labs(
    title = "RNA Target Family Distribution",
    x = "Target Family",
    y = "Number of Molecules"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)

ggsave(
  "results/plots/target_family_distribution.png",
  plot = p1,
  width = 8,
  height = 5,
  dpi = 300
)

############################################################
# Molecular Weight Distribution
############################################################

p2 <- ggplot(data, aes(MolecularWeight)) +
  geom_histogram(bins = 30, fill = "steelblue") +
  labs(title = "Distribution of Molecular Weight") +
  theme_bw()

print(p2)

ggsave(
  "results/plots/molecular_weight_distribution.png",
  plot = p2,
  width = 8,
  height = 5,
  dpi = 300
)

############################################################
# LogP Distribution
############################################################

p3 <- ggplot(data, aes(LogP)) +
  geom_histogram(bins = 30, fill = "steelblue") +
  labs(title = "Distribution of LogP") +
  theme_bw()

print(p3)

ggsave(
  "results/plots/logp_distribution.png",
  plot = p3,
  width = 8,
  height = 5,
  dpi = 300
)

############################################################
# TPSA Distribution
############################################################

p4 <- ggplot(data, aes(TPSA)) +
  geom_histogram(bins = 30, fill = "steelblue") +
  labs(title = "Distribution of TPSA") +
  theme_bw()

print(p4)

ggsave(
  "results/plots/tpsa_distribution.png",
  plot = p4,
  width = 8,
  height = 5,
  dpi = 300
)

############################################################
# Molecular Weight by RNA Family
############################################################

p5 <- ggplot(data,
             aes(Target_Family,
                 MolecularWeight)) +
  geom_boxplot(fill = "lightblue") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p5)

ggsave(
  "results/plots/mw_by_family.png",
  plot = p5,
  width = 9,
  height = 6,
  dpi = 300
)

############################################################
# LogP by RNA Family
############################################################

p6 <- ggplot(data,
             aes(Target_Family,
                 LogP)) +
  geom_boxplot(fill = "lightgreen") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p6)

ggsave(
  "results/plots/logp_by_family.png",
  plot = p6,
  width = 9,
  height = 6,
  dpi = 300
)

############################################################
# TPSA by RNA Family
############################################################

p7 <- ggplot(data,
             aes(Target_Family,
                 TPSA)) +
  geom_boxplot(fill = "orange") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p7)

ggsave(
  "results/plots/tpsa_by_family.png",
  plot = p7,
  width = 9,
  height = 6,
  dpi = 300
)

############################################################
# Correlation Matrix
############################################################

descriptor_data <- data %>%
  select(
    MolecularWeight,
    ExactMolWt,
    LogP,
    MolMR,
    TPSA,
    HBA,
    HBD,
    RotatableBonds,
    RingCount,
    NumAromaticRings,
    NumAliphaticRings,
    HeavyAtomCount,
    NumHeteroAtoms,
    FractionCSP3
  )

corr_matrix <- cor(descriptor_data)

png(
  "results/plots/correlation_heatmap.png",
  width = 1200,
  height = 1000,
  res = 150
)

corrplot(
  corr_matrix,
  method = "color",
  type = "upper"
)

dev.off()

# Display in RStudio
corrplot(
  corr_matrix,
  method = "color",
  type = "upper"
)

############################################################
# Pair Plot
############################################################

png(
  "results/plots/pair_plot.png",
  width = 1400,
  height = 1200,
  res = 150
)

pairs(descriptor_data[,1:6])

dev.off()

# Display in RStudio
pairs(descriptor_data[,1:6])

############################################################
# Summary Statistics by Target Family
############################################################

summary_stats <- data %>%
  group_by(Target_Family) %>%
  summarise(
    Mean_MW = mean(MolecularWeight),
    Mean_LogP = mean(LogP),
    Mean_TPSA = mean(TPSA),
    Mean_HBA = mean(HBA),
    .groups = "drop"
  )

print(summary_stats)

write.csv(
  summary_stats,
  "results/tables/summary_statistics.csv",
  row.names = FALSE
)

############################################################
# Save Final Dataset
############################################################

write.csv(
  data,
  "data/processed/cleaned_dataset.csv",
  row.names = FALSE
)

cat("\n========================================\n")
cat("EDA Completed Successfully!\n")
cat("Plots saved in : results/plots/\n")
cat("Tables saved in : results/tables/\n")
cat("========================================\n")