############################################################
# RNAInsight
# Step 10 : Formal Statistical Analysis & Chemical Space Visualization
############################################################
#
# Purpose:
#   Fill the "Statistical Analysis" rubric item, which is distinct from
#   both EDA (03_eda.R, descriptive only) and Model Development
#   (06-09, predictive only). Here we formally test whether each
#   molecular descriptor differs significantly across the 11 RNA
#   target families (inferential statistics), and visualize the
#   overall chemical space with PCA -- independent of any ML model.
#
# Methods:
#   - Kruskal-Wallis test per descriptor (non-parametric alternative to
#     one-way ANOVA; appropriate here because class sizes are small/
#     imbalanced and descriptor distributions are not reliably normal).
#   - Epsilon-squared effect size for each test (how much of the
#     variance in a descriptor is explained by target family).
#   - Pairwise Wilcoxon post-hoc tests (BH-corrected) for the most
#     significant descriptors, to see *which* families actually differ.
#   - PCA on the scaled 15-descriptor space, colored by target family,
#     to visualize how separable the classes are geometrically.
#
# Input  : data/processed/model_ready.csv
#          results/tables/target_family_label_map.csv
# Output : results/tables/statistical_tests_kruskal_wallis.csv
#          results/tables/pca_variance_explained.csv
#          results/plots/statistical/*.png
############################################################

library(ggplot2)
library(dplyr)

dir.create("results/plots/statistical", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

############################################################
# Load Data
############################################################

data <- read.csv("data/processed/model_ready.csv", stringsAsFactors = FALSE)
label_map <- read.csv("results/tables/target_family_label_map.csv", stringsAsFactors = FALSE)

# Restore readable family names for plots/tables (data currently uses
# make.names()-safe labels, e.g. "G.Quadruplex" instead of "G-Quadruplex")
name_lookup <- setNames(label_map$Original, label_map$Safe)
data$Target_Family_Label <- name_lookup[data$Target_Family]
data$Target_Family_Label <- factor(data$Target_Family_Label, levels = label_map$Original)

feature_cols <- setdiff(names(data), c("Target_Family", "Target_Family_Label"))

cat("========================================\n")
cat("RNAInsight - Statistical Analysis\n")
cat("========================================\n\n")
cat("Descriptors tested:", length(feature_cols), "\n")
cat("Families compared :", nlevels(data$Target_Family_Label), "\n\n")

############################################################
# Kruskal-Wallis Test + Epsilon-Squared Effect Size per Descriptor
############################################################
# Epsilon-squared = (H - k + 1) / (n - k), a common effect-size
# measure for the Kruskal-Wallis statistic (0 = no effect, 1 = max).

n_total <- nrow(data)
n_groups <- nlevels(data$Target_Family_Label)

kw_results <- data.frame(
  Descriptor = character(), H_statistic = numeric(),
  df = integer(), p_value = numeric(), epsilon_squared = numeric(),
  stringsAsFactors = FALSE
)

for (feat in feature_cols) {
  test <- kruskal.test(data[[feat]] ~ data$Target_Family_Label)
  eps_sq <- (test$statistic - n_groups + 1) / (n_total - n_groups)
  kw_results <- rbind(kw_results, data.frame(
    Descriptor = feat,
    H_statistic = round(unname(test$statistic), 3),
    df = unname(test$parameter),
    p_value = test$p.value,
    epsilon_squared = round(unname(eps_sq), 4)
  ))
}

# Benjamini-Hochberg correction across all 15 descriptor tests
kw_results$p_adjusted <- round(p.adjust(kw_results$p_value, method = "BH"), 6)
kw_results$p_value <- round(kw_results$p_value, 6)
kw_results$Significant <- ifelse(kw_results$p_adjusted < 0.05, "Yes", "No")
kw_results <- kw_results[order(-kw_results$epsilon_squared), ]

cat("Kruskal-Wallis results (sorted by effect size):\n")
print(kw_results, row.names = FALSE)

write.csv(kw_results, "results/tables/statistical_tests_kruskal_wallis.csv", row.names = FALSE)
cat("\nSaved: results/tables/statistical_tests_kruskal_wallis.csv\n")

############################################################
# Effect Size Ranking Plot
############################################################

p_effect <- ggplot(kw_results, aes(x = reorder(Descriptor, epsilon_squared), y = epsilon_squared, fill = Significant)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Yes" = "#1f77b4", "No" = "#cccccc")) +
  labs(
    title = "Effect Size of Each Descriptor on RNA Target Family",
    subtitle = "Kruskal-Wallis epsilon-squared (BH-adjusted significance at p < 0.05)",
    x = "Descriptor", y = "Epsilon-squared (effect size)", fill = "Significant"
  ) +
  theme_bw()

ggsave("results/plots/statistical/effect_size_ranking.png", p_effect, width = 8, height = 6, dpi = 300)
cat("Saved: results/plots/statistical/effect_size_ranking.png\n")

############################################################
# Annotated Boxplots for Top Significant Descriptors
############################################################

top_features <- head(kw_results$Descriptor[kw_results$Significant == "Yes"], 4)

for (feat in top_features) {
  p_value <- kw_results$p_adjusted[kw_results$Descriptor == feat]
  p_label <- if (p_value < 0.001) "p < 0.001" else paste0("p = ", round(p_value, 4))

  p_box <- ggplot(data, aes(x = Target_Family_Label, y = .data[[feat]], fill = Target_Family_Label)) +
    geom_boxplot(outlier.alpha = 0.4, show.legend = FALSE) +
    labs(
      title = paste0(feat, " by RNA Target Family"),
      subtitle = paste("Kruskal-Wallis (BH-adjusted):", p_label),
      x = "RNA Target Family", y = feat
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  safe_name <- gsub("[^A-Za-z0-9]", "_", feat)
  out_path <- paste0("results/plots/statistical/boxplot_", safe_name, ".png")
  ggsave(out_path, p_box, width = 8, height = 6, dpi = 300)
  cat("Saved:", out_path, "\n")
}

############################################################
# Pairwise Post-Hoc Tests (Wilcoxon, BH-corrected) for the
# Single Strongest Descriptor
############################################################

strongest_feature <- kw_results$Descriptor[1]
cat("\nPost-hoc pairwise Wilcoxon tests for strongest descriptor:", strongest_feature, "\n")

posthoc <- pairwise.wilcox.test(
  data[[strongest_feature]], data$Target_Family_Label,
  p.adjust.method = "BH", exact = FALSE
)
print(posthoc)

posthoc_df <- as.data.frame(as.table(posthoc$p.value))
names(posthoc_df) <- c("Family_A", "Family_B", "p_adjusted")
posthoc_df <- posthoc_df[!is.na(posthoc_df$p_adjusted), ]
posthoc_df$p_adjusted <- round(posthoc_df$p_adjusted, 6)
posthoc_df$Significant <- ifelse(posthoc_df$p_adjusted < 0.05, "Yes", "No")

write.csv(
  posthoc_df,
  paste0("results/tables/posthoc_", gsub("[^A-Za-z0-9]", "_", strongest_feature), ".csv"),
  row.names = FALSE
)
cat("Saved: results/tables/posthoc_", gsub("[^A-Za-z0-9]", "_", strongest_feature), ".csv\n", sep = "")

############################################################
# PCA: Whole-Dataset Chemical Space, Colored by Target Family
############################################################
# Independent of any ML model -- this shows geometrically how
# separable (or overlapping) the 11 target families are using
# only the 15 descriptors, which is the honest visual explanation
# for why classification accuracy tops out where it does.

pca_data <- data[, feature_cols]
pca_fit <- prcomp(pca_data, center = TRUE, scale. = TRUE)

var_explained <- round(100 * (pca_fit$sdev^2) / sum(pca_fit$sdev^2), 2)

pca_var_table <- data.frame(
  PC = paste0("PC", seq_along(var_explained)),
  Variance_Explained_Pct = var_explained,
  Cumulative_Pct = round(cumsum(var_explained), 2)
)
write.csv(pca_var_table, "results/tables/pca_variance_explained.csv", row.names = FALSE)
cat("\nSaved: results/tables/pca_variance_explained.csv\n")
cat("PC1 + PC2 explain", round(var_explained[1] + var_explained[2], 1), "% of total variance\n")

pca_scores <- as.data.frame(pca_fit$x[, 1:2])
pca_scores$Target_Family_Label <- data$Target_Family_Label

p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Target_Family_Label)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(
    title = "PCA of Molecular Descriptor Space",
    subtitle = paste0("PC1 (", var_explained[1], "%) vs PC2 (", var_explained[2], "%) -- colored by true RNA target family"),
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)"),
    color = "RNA Target Family"
  ) +
  theme_bw()

ggsave("results/plots/statistical/pca_scatter.png", p_pca, width = 9, height = 6, dpi = 300)
cat("Saved: results/plots/statistical/pca_scatter.png\n")

# Scree plot
scree_df <- pca_var_table[1:min(10, nrow(pca_var_table)), ]
p_scree <- ggplot(scree_df, aes(x = reorder(PC, -Variance_Explained_Pct), y = Variance_Explained_Pct)) +
  geom_col(fill = "#1f77b4") +
  geom_line(aes(group = 1), color = "firebrick") +
  geom_point(color = "firebrick") +
  labs(title = "PCA Scree Plot", x = "Principal Component", y = "% Variance Explained") +
  theme_bw()

ggsave("results/plots/statistical/pca_scree.png", p_scree, width = 7, height = 5, dpi = 300)
cat("Saved: results/plots/statistical/pca_scree.png\n")

cat("\n========================================\n")
cat("Statistical analysis complete.\n")
cat("Significant descriptors (BH-adjusted p < 0.05):", sum(kw_results$Significant == "Yes"), "of", nrow(kw_results), "\n")
cat("========================================\n")
