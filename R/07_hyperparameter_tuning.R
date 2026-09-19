############################################################
# RNAInsight
# Step 7 : Hyperparameter Tuning & Model Optimization (DA2)
############################################################
#
# Purpose:
#   06_model_training.R trained 15 algorithms with a modest
#   automatic tuneLength. Here we take the top 3 models by
#   baseline cross-validated Accuracy and re-tune them with a
#   wider, method-specific hyperparameter grid, then report the
#   before/after improvement.
#
# Input  : results/models/trained_models.rds
#          results/models/train_data.rds, test_data.rds
# Output : results/models/tuned_models.rds
#          results/tables/hyperparameter_tuning_report.csv
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
# Load Baseline Models
############################################################

trained_models <- readRDS("results/models/trained_models.rds")
train_data     <- readRDS("results/models/train_data.rds")
test_data      <- readRDS("results/models/test_data.rds")

cat("========================================\n")
cat("RNAInsight - Hyperparameter Tuning\n")
cat("========================================\n\n")

baseline_acc <- sapply(trained_models, function(fit) max(fit$results$Accuracy))
baseline_acc <- sort(baseline_acc, decreasing = TRUE)

cat("Baseline CV Accuracy (from 06_model_training.R):\n")
print(round(baseline_acc, 4))

top_models <- names(baseline_acc)[1:min(3, length(baseline_acc))]
cat("\nSelected for deeper tuning:", paste(top_models, collapse = ", "), "\n\n")

############################################################
# Method-specific "wide" tuning grids
############################################################

n_features <- ncol(train_data) - 1

get_wide_grid <- function(method) {
  switch(method,
    "rf" = expand.grid(mtry = unique(pmax(1, round(seq(1, n_features, length.out = 6))))),
    "xgbTree" = expand.grid(
      nrounds = c(50, 100, 150),
      max_depth = c(2, 4, 6),
      eta = c(0.05, 0.1, 0.3),
      gamma = 0,
      colsample_bytree = 0.8,
      min_child_weight = 1,
      subsample = 0.8
    ),
    "gbm" = expand.grid(
      n.trees = c(50, 100, 150),
      interaction.depth = c(1, 3, 5),
      shrinkage = c(0.01, 0.1),
      n.minobsinnode = 10
    ),
    "svmRadial" = expand.grid(
      sigma = c(0.01, 0.05, 0.1, 0.2),
      C = c(0.25, 0.5, 1, 2, 4, 8)
    ),
    "svmLinear" = expand.grid(C = c(0.01, 0.05, 0.25, 0.5, 1, 2, 4, 8)),
    "knn" = expand.grid(k = seq(3, 25, by = 2)),
    "glmnet" = expand.grid(
      alpha = c(0, 0.25, 0.5, 0.75, 1),
      lambda = 10 ^ seq(-4, -1, length.out = 6)
    ),
    "nnet" = expand.grid(
      size = c(3, 5, 7, 9),
      decay = c(0, 0.01, 0.1, 0.5)
    ),
    "C5.0" = expand.grid(
      trials = c(1, 10, 20, 30),
      model = "tree",
      winnow = FALSE
    ),
    "AdaBoost.M1" = expand.grid(
      mfinal = c(50, 100, 150),
      maxdepth = c(1, 3, 5),
      coeflearn = "Breiman"
    ),
    "rpart" = expand.grid(cp = seq(0.001, 0.05, length.out = 10)),
    "multinom" = expand.grid(decay = c(0, 0.001, 0.01, 0.1, 0.5)),
    "naive_bayes" = expand.grid(laplace = c(0, 0.5, 1), usekernel = c(TRUE, FALSE), adjust = 1),
    NULL # lda / treebag have no tunable hyperparameters in caret
  )
}

############################################################
# Shared CV control (SAME fold indices as 06_model_training.R,
# so tuned-vs-baseline and model-vs-model comparisons are fair)
############################################################

cv_folds <- readRDS("results/models/cv_folds.rds")

ctrl <- trainControl(
  method = "cv",
  index = cv_folds,
  classProbs = TRUE,
  savePredictions = "final",
  verboseIter = FALSE
)

############################################################
# Re-tune Top Models
############################################################

tuned_models <- list()
tuning_report <- data.frame(
  Model = character(), Method = character(),
  Baseline_Accuracy = numeric(), Tuned_Accuracy = numeric(),
  Improvement = numeric(), Best_Params = character(),
  stringsAsFactors = FALSE
)

for (model_name in top_models) {

  method <- trained_models[[model_name]]$method
  grid <- get_wide_grid(method)

  cat("----------------------------------------\n")
  cat("Tuning:", model_name, "(", method, ")\n")

  if (is.null(grid)) {
    cat("  No tunable hyperparameters for this method - keeping baseline fit.\n")
    tuned_models[[model_name]] <- trained_models[[model_name]]
    tuning_report <- rbind(tuning_report, data.frame(
      Model = model_name, Method = method,
      Baseline_Accuracy = round(baseline_acc[[model_name]], 4),
      Tuned_Accuracy = round(baseline_acc[[model_name]], 4),
      Improvement = 0,
      Best_Params = "n/a (no hyperparameters)"
    ))
    next
  }

  # Random Forest needs importance = TRUE at fit time for varImp() to
  # work later in 09_visualizations.R.
  extra_args <- if (method == "rf") list(importance = TRUE) else list()

  fit <- tryCatch({
    args <- c(
      list(
        Target_Family ~ .,
        data = train_data,
        method = method,
        trControl = ctrl,
        preProcess = c("center", "scale"),
        tuneGrid = grid,
        metric = "Accuracy"
      ),
      extra_args
    )
    do.call(caret::train, args)
  }, error = function(e) {
    cat("  FAILED:", conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(fit)) {
    tuned_acc <- max(fit$results$Accuracy)
    tuned_models[[model_name]] <- fit

    cat("  Baseline Accuracy:", round(baseline_acc[[model_name]], 4), "\n")
    cat("  Tuned Accuracy   :", round(tuned_acc, 4), "\n")
    cat("  Best params      :\n")
    print(fit$bestTune)

    tuning_report <- rbind(tuning_report, data.frame(
      Model = model_name, Method = method,
      Baseline_Accuracy = round(baseline_acc[[model_name]], 4),
      Tuned_Accuracy = round(tuned_acc, 4),
      Improvement = round(tuned_acc - baseline_acc[[model_name]], 4),
      Best_Params = paste(names(fit$bestTune), unlist(fit$bestTune), sep = "=", collapse = ", ")
    ))
  }
}

############################################################
# Merge tuned models back into the full model set
############################################################

final_models <- trained_models
for (model_name in names(tuned_models)) {
  final_models[[model_name]] <- tuned_models[[model_name]]
}

saveRDS(final_models, "results/models/tuned_models.rds")

write.csv(
  tuning_report,
  "results/tables/hyperparameter_tuning_report.csv",
  row.names = FALSE
)

cat("\n========================================\n")
cat("Hyperparameter tuning complete.\n")
print(tuning_report)
cat("\nSaved: results/models/tuned_models.rds\n")
cat("Saved: results/tables/hyperparameter_tuning_report.csv\n")
cat("========================================\n")
