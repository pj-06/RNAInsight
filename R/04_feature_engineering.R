############################################################
# RNAInsight
# Step 4 : Feature Engineering & Feature Selection (DA2)
############################################################
#
# Purpose:
#   1. Derive new chemically-meaningful descriptors from the
#      existing raw descriptors (feature engineering).
#   2. Remove redundant / low-value features (feature selection)
#      before the data is loaded into the database and modelled.
#
# Input  : data/processed/cleaned_dataset.csv
# Output : data/processed/model_dataset.csv   (final modelling table)
#          results/tables/feature_selection_report.csv
############################################################

# ---- Package bootstrap -------------------------------------------------
required_packages <- c("dplyr", "caret", "randomForest")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing missing package:", pkg, "\n")
    install.packages(pkg, dependencies = TRUE)
  }
}

library(dplyr)
library(caret)

############################################################
# Load Dataset
############################################################

data <- read.csv(
  "data/processed/cleaned_dataset.csv",
  stringsAsFactors = FALSE
)

data$Target_Family <- as.factor(data$Target_Family)

cat("========================================\n")
cat("RNAInsight - Feature Engineering & Selection\n")
cat("========================================\n\n")

cat("Rows :", nrow(data), "\n")
cat("Columns :", ncol(data), "\n\n")

############################################################
# Create Output Folders
############################################################

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

############################################################
# 0. Sanitize Class Labels
############################################################
# Several caret methods (gbm, multinom, xgbTree, ...) call
# make.names() internally on the class labels when classProbs = TRUE.
# Labels like "G-Quadruplex" or "HIV TAR" (hyphens/spaces) can then
# cause mismatched column names between predictions and truth. We
# sanitize once here, up front, and every downstream script (DB,
# training, tuning, comparison, visualization) inherits the safe
# labels automatically because they all read from files this script
# writes. The mapping back to the original, readable names is saved
# for use in the final report.

original_levels <- levels(data$Target_Family)
safe_levels <- make.names(original_levels)
levels(data$Target_Family) <- safe_levels

label_map <- data.frame(Original = original_levels, Safe = safe_levels)
write.csv(label_map, "results/tables/target_family_label_map.csv", row.names = FALSE)

cat("Class labels sanitized for R/caret compatibility (see\n")
cat("results/tables/target_family_label_map.csv for the original names):\n")
print(label_map)
cat("\n")

############################################################
# 1. Feature Engineering
############################################################
# New descriptors derived from the base RDKit-style features.
# These follow standard cheminformatics conventions and are
# computed directly from existing columns (no external lookup,
# no leakage of the target label).

data <- data %>%
  mutate(
    # Total hydrogen bonding capacity
    HBondTotal = HBA + HBD,

    # Lipinski "Rule of Five" violation count (drug-likeness)
    LipinskiViolations =
      as.integer(MolecularWeight > 500) +
      as.integer(LogP > 5) +
      as.integer(HBA > 10) +
      as.integer(HBD > 5),

    # Fraction of rings that are aromatic (0 when no rings at all)
    AromaticRingFraction = ifelse(RingCount == 0, 0, NumAromaticRings / RingCount),

    # Molecular flexibility: rotatable bonds per heavy atom
    FlexibilityIndex = RotatableBonds / HeavyAtomCount,

    # Heteroatom density
    HeteroatomFraction = NumHeteroAtoms / HeavyAtomCount,

    # Polarity normalised by size
    TPSA_per_MW = TPSA / MolecularWeight
  )

cat("Engineered 6 new features:\n")
cat(" - HBondTotal, LipinskiViolations, AromaticRingFraction,\n")
cat("   FlexibilityIndex, HeteroatomFraction, TPSA_per_MW\n\n")

############################################################
# 2. Feature Selection
############################################################

# Candidate numeric predictors (raw + engineered), target excluded
candidate_features <- c(
  "MolecularWeight", "ExactMolWt", "LogP", "MolMR", "TPSA",
  "HBA", "HBD", "RotatableBonds", "RingCount",
  "NumAromaticRings", "NumAliphaticRings", "HeavyAtomCount",
  "NumHeteroAtoms", "FractionCSP3",
  "HBondTotal", "LipinskiViolations", "AromaticRingFraction",
  "FlexibilityIndex", "HeteroatomFraction", "TPSA_per_MW"
)

feature_matrix <- data[, candidate_features]

# --- 2a. Near-zero-variance filter ---
nzv <- nearZeroVar(feature_matrix, saveMetrics = TRUE)
nzv_drop <- rownames(nzv)[nzv$nzv]

cat("Near-zero-variance features flagged for removal:\n")
if (length(nzv_drop) == 0) {
  cat(" (none)\n")
} else {
  print(nzv_drop)
}

feature_matrix <- feature_matrix[, !(names(feature_matrix) %in% nzv_drop)]

# --- 2b. Collinearity filter (|r| > 0.90) ---
corr_matrix <- cor(feature_matrix, use = "complete.obs")
high_corr <- findCorrelation(corr_matrix, cutoff = 0.90, names = TRUE)

cat("\nHighly-correlated features flagged for removal (|r| > 0.90):\n")
if (length(high_corr) == 0) {
  cat(" (none)\n")
} else {
  print(high_corr)
}

selected_features <- setdiff(names(feature_matrix), high_corr)

cat("\nFinal selected feature set (", length(selected_features), " features):\n", sep = "")
print(selected_features)

# --- 2c. Quick Random-Forest importance ranking (for reporting only) ---
set.seed(42)
rf_rank <- caret::train(
  x = data[, selected_features],
  y = data$Target_Family,
  method = "rf",
  ntree = 200,
  importance = TRUE,
  trControl = trainControl(method = "none"),
  tuneGrid = data.frame(mtry = floor(sqrt(length(selected_features))))
)

importance_df <- varImp(rf_rank)$importance
importance_df$Feature <- rownames(importance_df)

# varImp() can return either a single "Overall" column or one column
# per class, depending on the importance type computed -- handle both.
if (!"Overall" %in% names(importance_df)) {
  score_cols <- setdiff(names(importance_df), "Feature")
  importance_df$Overall <- rowMeans(importance_df[, score_cols, drop = FALSE], na.rm = TRUE)
}

importance_df <- importance_df %>% arrange(desc(Overall))

cat("\nFeature importance ranking (mean decrease, RF filter model):\n")
print(importance_df)

############################################################
# 3. Save Outputs
############################################################

model_dataset <- data[, c("Name", "SMILES", "Target_Family", selected_features)]

write.csv(
  model_dataset,
  "data/processed/model_dataset.csv",
  row.names = FALSE
)

selection_report <- data.frame(
  Feature = candidate_features,
  Dropped_NZV = candidate_features %in% nzv_drop,
  Dropped_Collinear = candidate_features %in% high_corr,
  Kept = candidate_features %in% selected_features
)

write.csv(
  selection_report,
  "results/tables/feature_selection_report.csv",
  row.names = FALSE
)

write.csv(
  importance_df,
  "results/tables/feature_importance_filter.csv",
  row.names = FALSE
)

cat("\n========================================\n")
cat("Feature engineering & selection complete!\n")
cat("Modelling dataset  : data/processed/model_dataset.csv\n")
cat("Selection report   : results/tables/feature_selection_report.csv\n")
cat("Importance ranking : results/tables/feature_importance_filter.csv\n")
cat("========================================\n")
