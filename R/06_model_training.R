############################################################
# RNAInsight
# Step 6 : Model Development - 15 ML/DL Algorithms (DA2)
############################################################
#
# Purpose:
#   Train 15 classification algorithms (covering linear,
#   probabilistic, instance-based, tree, ensemble, boosting,
#   margin-based and neural-network families) on the same
#   stratified train/test split and the same 5x3 repeated-CV
#   scheme, so results are directly comparable in
#   08_comparative_analysis.R.
#
# Input  : data/processed/model_ready.csv (pulled from the DB
#          by 05_database_connectivity.R)
# Output : results/models/trained_models.rds
#          results/models/test_predictions.rds
#          results/tables/model_training_log.csv
############################################################

# ---- Package bootstrap -------------------------------------------------
required_packages <- c(
  "caret", "dplyr",
  "nnet", "MASS", "naivebayes", "kernlab", "class",
  "rpart", "C50", "ipred", "e1071", "plyr",
  "randomForest", "gbm", "xgboost", "adabag", "glmnet"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing missing package:", pkg, "\n")
    install.packages(pkg, dependencies = TRUE)
  }
}

library(caret)
library(dplyr)

############################################################
# Create Output Folders
############################################################

dir.create("results/models", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

############################################################
# Load Dataset (retrieved from the database)
############################################################

data <- read.csv(
  "data/processed/model_ready.csv",
  stringsAsFactors = FALSE
)

data$Target_Family <- as.factor(data$Target_Family)

cat("========================================\n")
cat("RNAInsight - Model Development (15 algorithms)\n")
cat("========================================\n\n")
cat("Rows :", nrow(data), " | Features :", ncol(data) - 1, "\n")
cat("Classes:\n")
print(table(data$Target_Family))

############################################################
# Stratified Train / Test Split (80 / 20)
############################################################

set.seed(42)
train_idx <- createDataPartition(data$Target_Family, p = 0.8, list = FALSE)

train_data <- data[train_idx, ]
test_data  <- data[-train_idx, ]

cat("\nTrain rows:", nrow(train_data), " | Test rows:", nrow(test_data), "\n\n")

############################################################
# Shared Cross-Validation Control
############################################################
# Fixed fold indices (5-fold x 3 repeats = 15 resamples) are
# generated once and reused for every model below. This is what
# makes the resamples()-based statistical comparison in
# 08_comparative_analysis.R valid -- every model is evaluated on
# the exact same train/validation splits, not just the same
# scheme with different random folds.

set.seed(42)
cv_folds <- createMultiFolds(train_data$Target_Family, k = 5, times = 3)

ctrl <- trainControl(
  method = "cv",
  index = cv_folds,
  classProbs = TRUE,
  savePredictions = "final",
  verboseIter = FALSE
)

saveRDS(cv_folds, "results/models/cv_folds.rds")

############################################################
# 15 Model Specifications
############################################################
# tuneLength controls how many hyperparameter combinations
# caret tries automatically for each algorithm during this
# baseline pass. 07_hyperparameter_tuning.R goes deeper on the
# top performers.

# XGBoost gets an explicit, small grid instead of caret's automatic
# tuneLength grid (which can silently balloon into dozens of
# combinations for xgbTree) and is forced to single-threaded mode --
# the xgboost R package has a known hang/deadlock issue with its
# internal OpenMP multi-threading inside RStudio on Windows.
xgb_grid <- expand.grid(
  nrounds = c(50, 100),
  max_depth = c(3, 6),
  eta = c(0.1, 0.3),
  gamma = 0,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.8
)

model_specs <- list(
  list(name = "Multinomial Logistic Regression", method = "multinom",     tuneLength = 3, extra = list(trace = FALSE)),
  list(name = "Linear Discriminant Analysis",     method = "lda",         tuneLength = 1, extra = list()),
  list(name = "Naive Bayes",                      method = "naive_bayes", tuneLength = 3, extra = list()),
  list(name = "K-Nearest Neighbors",               method = "knn",        tuneLength = 5, extra = list()),
  list(name = "Decision Tree (CART)",              method = "rpart",      tuneLength = 5, extra = list()),
  list(name = "Decision Tree (C5.0)",               method = "C5.0",      tuneLength = 3, extra = list()),
  list(name = "Bagged CART",                        method = "treebag",   tuneLength = 1, extra = list()),
  list(name = "Random Forest",                      method = "rf",        tuneLength = 3, extra = list(importance = TRUE)),
  list(name = "Gradient Boosting Machine",          method = "gbm",       tuneLength = 3, extra = list(verbose = FALSE)),
  list(name = "XGBoost",                            method = "xgbTree",   tuneLength = NULL, grid = xgb_grid, extra = list(verbose = 0, nthread = 1)),
  #list(name = "AdaBoost (M1)",                      method = "AdaBoost.M1", tuneLength = 3, extra = list()),
  list(name = "SVM (Linear Kernel)",                method = "svmLinear", tuneLength = 3, extra = list()),
  list(name = "SVM (Radial Kernel)",                method = "svmRadial", tuneLength = 5, extra = list()),
  list(name = "Regularized Multinomial (glmnet)",   method = "glmnet",    tuneLength = 5, extra = list()),
  list(name = "Neural Network (nnet)",              method = "nnet",      tuneLength = 5, extra = list(trace = FALSE))
)

############################################################
# Train All Models
############################################################

trained_models <- list()
training_log <- data.frame(
  Model = character(), Method = character(),
  Status = character(), Seconds = numeric(),
  Message = character(), stringsAsFactors = FALSE
)

for (spec in model_specs) {

  cat("----------------------------------------\n")
  cat("Training:", spec$name, "(", spec$method, ")\n")

  start_time <- Sys.time()

  fit <- tryCatch({
    base_args <- list(
      Target_Family ~ .,
      data = train_data,
      method = spec$method,
      trControl = ctrl,
      preProcess = c("center", "scale"),
      metric = "Accuracy"
    )

    # Use an explicit tuneGrid when the spec provides one (e.g. XGBoost),
    # otherwise fall back to caret's automatic tuneLength search.
    if (!is.null(spec$grid)) {
      base_args$tuneGrid <- spec$grid
    } else {
      base_args$tuneLength <- spec$tuneLength
    }

    args <- c(base_args, spec$extra)
    do.call(caret::train, args)
  }, error = function(e) {
    cat("  FAILED:", conditionMessage(e), "\n")
    return(NULL)
  })

  elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)

  if (!is.null(fit)) {
    trained_models[[spec$name]] <- fit
    cat("  OK (", elapsed, "sec ) - best Accuracy:",
        round(max(fit$results$Accuracy), 4), "\n")
    training_log <- rbind(training_log, data.frame(
      Model = spec$name, Method = spec$method,
      Status = "OK", Seconds = elapsed, Message = ""
    ))
  } else {
    training_log <- rbind(training_log, data.frame(
      Model = spec$name, Method = spec$method,
      Status = "FAILED", Seconds = elapsed,
      Message = "See console output above"
    ))
  }
}

cat("\n========================================\n")
cat("Trained", length(trained_models), "of", length(model_specs), "models successfully.\n")
cat("========================================\n\n")

############################################################
# Test-Set Predictions (for all successful models)
############################################################

test_predictions <- list(
  truth = test_data$Target_Family
)

for (model_name in names(trained_models)) {
  fit <- trained_models[[model_name]]

  pred_class <- tryCatch(predict(fit, newdata = test_data), error = function(e) NULL)
  pred_prob  <- tryCatch(predict(fit, newdata = test_data, type = "prob"), error = function(e) NULL)

  test_predictions[[model_name]] <- list(
    class = pred_class,
    prob  = pred_prob
  )
}

############################################################
# Save Everything
############################################################

saveRDS(trained_models, "results/models/trained_models.rds")
saveRDS(test_predictions, "results/models/test_predictions.rds")
saveRDS(test_data, "results/models/test_data.rds")
saveRDS(train_data, "results/models/train_data.rds")

write.csv(
  training_log,
  "results/tables/model_training_log.csv",
  row.names = FALSE
)

cat("Saved: results/models/trained_models.rds\n")
cat("Saved: results/models/test_predictions.rds\n")
cat("Saved: results/tables/model_training_log.csv\n")
cat("\nModel development stage complete.\n")
