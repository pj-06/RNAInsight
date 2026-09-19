############################################################
# RNAInsight
# Step 8 : Comparative Performance Analysis (DA2)
############################################################
#
# Purpose:
#   Evaluate all 15 trained/tuned models on the held-out test
#   set with a common set of metrics (Accuracy, Kappa, macro
#   Precision/Recall/F1), then run a statistical comparison
#   across their cross-validation resamples to see whether
#   differences in Accuracy are meaningful or within noise.
#
# Input  : results/models/tuned_models.rds
#          results/models/test_data.rds
# Output : results/tables/model_comparison.csv
#          results/tables/resamples_summary.txt
#          results/tables/pairwise_accuracy_pvalues.csv
#          results/models/resamples.rds
############################################################

# ---- Package bootstrap -------------------------------------------------
required_packages <- c("caret", "dplyr")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing missing package:", pkg, "\n")
    install.packages(pkg, dependencies = TRUE)
  }
}

library(caret)
library(dplyr)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

############################################################
# Load Models & Test Data
############################################################

final_models <- readRDS("results/models/tuned_models.rds")
test_data    <- readRDS("results/models/test_data.rds")

cat("========================================\n")
cat("RNAInsight - Comparative Performance Analysis\n")
cat("========================================\n\n")
cat("Models to compare:", length(final_models), "\n")
cat("Test set size    :", nrow(test_data), "\n\n")

############################################################
# 1. Test-Set Metrics per Model
############################################################

leaderboard <- data.frame(
  Model = character(), Accuracy = numeric(), Kappa = numeric(),
  Macro_Precision = numeric(), Macro_Recall = numeric(), Macro_F1 = numeric(),
  stringsAsFactors = FALSE
)

for (model_name in names(final_models)) {

  fit <- final_models[[model_name]]
  pred <- tryCatch(predict(fit, newdata = test_data), error = function(e) NULL)

  if (is.null(pred)) {
    cat("Skipping", model_name, "- prediction failed\n")
    next
  }

  cm <- confusionMatrix(pred, test_data$Target_Family)

  # byClass is a matrix when there are >2 classes (one row per class)
  by_class <- cm$byClass
  macro_precision <- mean(by_class[, "Precision"], na.rm = TRUE)
  macro_recall    <- mean(by_class[, "Recall"], na.rm = TRUE)
  macro_f1        <- mean(by_class[, "F1"], na.rm = TRUE)

  leaderboard <- rbind(leaderboard, data.frame(
    Model = model_name,
    Accuracy = round(unname(cm$overall["Accuracy"]), 4),
    Kappa = round(unname(cm$overall["Kappa"]), 4),
    Macro_Precision = round(macro_precision, 4),
    Macro_Recall = round(macro_recall, 4),
    Macro_F1 = round(macro_f1, 4)
  ))
}

leaderboard <- leaderboard %>% arrange(desc(Accuracy))

cat("Test-set leaderboard:\n")
print(leaderboard)

write.csv(
  leaderboard,
  "results/tables/model_comparison.csv",
  row.names = FALSE
)

############################################################
# 2. Statistical Comparison Across CV Resamples
############################################################
# All models share the exact same 15 CV resamples (see
# 06_model_training.R), so this is a fair paired comparison of
# cross-validated Accuracy/Kappa, not just a single test-set
# snapshot.

resamps <- resamples(final_models)

cat("\nCross-validated performance across resamples:\n")
resamps_summary <- summary(resamps)
print(resamps_summary)

sink("results/tables/resamples_summary.txt")
cat("RNAInsight - Cross-Validated Resamples Summary\n")
cat("================================================\n\n")
print(resamps_summary)
sink()

# Pairwise significance of Accuracy differences (paired t-tests).
# Wrapped defensively: this is a "nice to have" statistical add-on,
# and the exact internal structure of summary(diff.resamples) has
# varied slightly across caret versions -- if it doesn't match here,
# the leaderboard/resamples outputs above are unaffected.
pvalue_table <- tryCatch({
  diffs <- diff(resamps, metric = "Accuracy")
  diffs_summary <- summary(diffs)
  tbl <- diffs_summary$table[["Accuracy"]]
  if (is.null(tbl)) tbl <- diffs_summary$table[[1]]

  pt <- as.data.frame(as.table(tbl))
  names(pt) <- c("Model_A", "Model_B", "p_value")
  pt[!is.na(pt$p_value), ]
}, error = function(e) {
  cat("Pairwise significance table skipped (", conditionMessage(e), ")\n")
  data.frame(Model_A = character(), Model_B = character(), p_value = numeric())
})

write.csv(
  pvalue_table,
  "results/tables/pairwise_accuracy_pvalues.csv",
  row.names = FALSE
)

saveRDS(resamps, "results/models/resamples.rds")

cat("\n========================================\n")
cat("Comparative analysis complete.\n")
cat("Best model by test Accuracy :", leaderboard$Model[1],
    "(", leaderboard$Accuracy[1], ")\n")
cat("Saved: results/tables/model_comparison.csv\n")
cat("Saved: results/tables/resamples_summary.txt\n")
cat("Saved: results/tables/pairwise_accuracy_pvalues.csv\n")
cat("========================================\n")
