############################################################
# RNAInsight
# Step 9 : Comparative Visualizations (DA2)
############################################################
#
# Purpose:
#   Produce the visual evidence for model comparison:
#     - Confusion matrix heatmaps (top 3 models)
#     - One-vs-rest multiclass ROC curves (best model)
#     - One-vs-rest Precision-Recall curves (best model)
#     - Feature importance bar charts (tree/ensemble models)
#     - CV Accuracy distribution across all 15 models
#
# Input  : results/models/tuned_models.rds, test_data.rds,
#          resamples.rds
#          results/tables/model_comparison.csv
# Output : results/plots/da2/*.png
############################################################

required_packages <- c("ggplot2", "pROC", "PRROC", "caret")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing missing package:", pkg, "\n")
    install.packages(pkg, dependencies = TRUE)
  }
}

library(ggplot2)
library(caret)
library(dplyr)

dir.create("results/plots/da2", recursive = TRUE, showWarnings = FALSE)

############################################################
# Load Everything
############################################################

final_models <- readRDS("results/models/tuned_models.rds")
test_data    <- readRDS("results/models/test_data.rds")
resamps      <- readRDS("results/models/resamples.rds")
leaderboard  <- read.csv("results/tables/model_comparison.csv", stringsAsFactors = FALSE)
label_map    <- read.csv("results/tables/target_family_label_map.csv", stringsAsFactors = FALSE)

# Class labels were sanitized (make.names) for caret compatibility in
# 04_feature_engineering.R. Use the original, readable names on axes
# and legends here; all computation still uses the safe labels.
pretty_label <- function(x) {
  idx <- match(x, label_map$Safe)
  ifelse(is.na(idx), x, label_map$Original[idx])
}

cat("========================================\n")
cat("RNAInsight - Comparative Visualizations\n")
cat("========================================\n\n")

top3 <- head(leaderboard$Model, 3)
best_model_name <- leaderboard$Model[1]
cat("Top 3 models :", paste(top3, collapse = ", "), "\n")
cat("Best model   :", best_model_name, "\n\n")

############################################################
# 1. Confusion Matrix Heatmaps (top 3 models)
############################################################

for (model_name in top3) {

  fit <- final_models[[model_name]]
  pred <- predict(fit, newdata = test_data)
  cm <- confusionMatrix(pred, test_data$Target_Family)
  cm_df <- as.data.frame(cm$table)
  cm_df$Reference <- pretty_label(as.character(cm_df$Reference))
  cm_df$Prediction <- pretty_label(as.character(cm_df$Prediction))

  p <- ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Freq), size = 3) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(
      title = paste("Confusion Matrix -", model_name),
      x = "Actual Target Family", y = "Predicted Target Family"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  safe_name <- gsub("[^A-Za-z0-9]+", "_", model_name)
  out_path <- paste0("results/plots/da2/confusion_matrix_", safe_name, ".png")

  ggsave(out_path, plot = p, width = 8, height = 7, dpi = 300)
  cat("Saved:", out_path, "\n")
}

############################################################
# 2. One-vs-Rest ROC Curves (best model)
############################################################

best_fit <- final_models[[best_model_name]]
prob_matrix <- predict(best_fit, newdata = test_data, type = "prob")
truth <- test_data$Target_Family
classes <- levels(truth)

roc_path <- "results/plots/da2/roc_curves_best_model.png"
png(roc_path, width = 1400, height = 1200, res = 150)

plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "False Positive Rate", ylab = "True Positive Rate",
     main = paste("One-vs-Rest ROC Curves -", best_model_name))
abline(a = 0, b = 1, lty = 2, col = "gray60")

colors <- rainbow(length(classes))
auc_values <- c()

for (i in seq_along(classes)) {
  cls <- classes[i]
  binary_truth <- as.numeric(truth == cls)

  if (length(unique(binary_truth)) < 2 || !(cls %in% colnames(prob_matrix))) next

  roc_obj <- pROC::roc(binary_truth, prob_matrix[[cls]], quiet = TRUE)
  lines(1 - roc_obj$specificities, roc_obj$sensitivities, col = colors[i], lwd = 2)
  auc_values[cls] <- as.numeric(pROC::auc(roc_obj))
}

legend("bottomright",
       legend = paste0(pretty_label(names(auc_values)), " (AUC=", round(auc_values, 3), ")"),
       col = colors[seq_along(auc_values)], lwd = 2, cex = 0.65, bty = "n")

dev.off()
cat("Saved:", roc_path, "\n")

write.csv(
  data.frame(Class = pretty_label(names(auc_values)), AUC = round(auc_values, 4)),
  "results/tables/roc_auc_by_class.csv",
  row.names = FALSE
)

############################################################
# 3. One-vs-Rest Precision-Recall Curves (best model)
############################################################

if (requireNamespace("PRROC", quietly = TRUE)) {

  pr_path <- "results/plots/da2/precision_recall_curves_best_model.png"
  png(pr_path, width = 1400, height = 1200, res = 150)

  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Recall", ylab = "Precision",
       main = paste("One-vs-Rest Precision-Recall Curves -", best_model_name))

  for (i in seq_along(classes)) {
    cls <- classes[i]
    binary_truth <- as.numeric(truth == cls)

    if (length(unique(binary_truth)) < 2 || !(cls %in% colnames(prob_matrix))) next

    pr_obj <- PRROC::pr.curve(
      scores.class0 = prob_matrix[[cls]][binary_truth == 1],
      scores.class1 = prob_matrix[[cls]][binary_truth == 0],
      curve = TRUE
    )
    lines(pr_obj$curve[, 1], pr_obj$curve[, 2], col = colors[i], lwd = 2)
  }

  legend("bottomleft", legend = pretty_label(classes), col = colors, lwd = 2, cex = 0.65, bty = "n")
  dev.off()
  cat("Saved:", pr_path, "\n")
} else {
  cat("PRROC not available - skipping precision-recall curves.\n")
}

############################################################
# 4. Feature Importance (tree / ensemble models)
############################################################

importance_candidates <- intersect(
  c("Random Forest", "XGBoost", "Gradient Boosting Machine", "Decision Tree (C5.0)", "AdaBoost (M1)"),
  names(final_models)
)

for (model_name in importance_candidates) {

  fit <- final_models[[model_name]]
  imp <- tryCatch(varImp(fit), error = function(e) NULL)
  if (is.null(imp)) next

  imp_df <- imp$importance
  imp_df$Feature <- rownames(imp_df)

  # Some caret importance objects are per-class; collapse to a single
  # "Overall" score by averaging across columns if needed.
  score_cols <- setdiff(names(imp_df), "Feature")
  imp_df$Overall <- rowMeans(imp_df[, score_cols, drop = FALSE], na.rm = TRUE)

  imp_df <- imp_df %>% arrange(desc(Overall)) %>% head(15)

  p <- ggplot(imp_df, aes(x = reorder(Feature, Overall), y = Overall)) +
    geom_col(fill = "darkorange") +
    coord_flip() +
    labs(
      title = paste("Feature Importance -", model_name),
      x = "Feature", y = "Importance"
    ) +
    theme_bw()

  safe_name <- gsub("[^A-Za-z0-9]+", "_", model_name)
  out_path <- paste0("results/plots/da2/feature_importance_", safe_name, ".png")

  ggsave(out_path, plot = p, width = 8, height = 6, dpi = 300)
  cat("Saved:", out_path, "\n")
}

############################################################
# 5. CV Accuracy Distribution Across All Models
############################################################

bwplot_path <- "results/plots/da2/cv_accuracy_comparison.png"
png(bwplot_path, width = 1600, height = 1000, res = 150)
print(bwplot(resamps, metric = "Accuracy", main = "Cross-Validated Accuracy by Model"))
dev.off()
cat("Saved:", bwplot_path, "\n")

############################################################
# 6. Leaderboard Bar Chart
############################################################

leaderboard$Model <- factor(leaderboard$Model, levels = rev(leaderboard$Model))

p_leader <- ggplot(leaderboard, aes(x = Model, y = Accuracy)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Test-Set Accuracy by Model",
    x = "", y = "Accuracy"
  ) +
  theme_bw()

ggsave("results/plots/da2/leaderboard_accuracy.png", plot = p_leader, width = 8, height = 7, dpi = 300)
cat("Saved: results/plots/da2/leaderboard_accuracy.png\n")

cat("\n========================================\n")
cat("Comparative visualizations complete!\n")
cat("All plots saved in: results/plots/da2/\n")
cat("========================================\n")
