# DA2 Progress Report - RNAInsight

**Course project:** PDS Project Execution Plan (Dr. Vijayalakshmi A, VIT Chennai)
**Assignment:** Digital Assignment 2 - Model Development, Database Connectivity & Comparative Analysis
**Submission date:** 17 September 2026
**Reported completion:** 75%

## 1. Recap: What DA1 Covered

- Python phase (`scripts/01`-`05`): parsed the raw R-BIND Excel file, standardized target
  labels, generated 14 molecular descriptors per molecule/target pair, cleaned and validated
  the data down to a final set of 750 records across 11 RNA target families.
- R phase (`R/03_eda.R`): exploratory data analysis - distributions, boxplots by family,
  correlation heatmap, pair plot, summary statistics table.

## 2. DA2 Scope and How It Maps to the Rubric

| Rubric item | Marks | Script(s) |
|---|---|---|
| Feature engineering and feature selection | 1 | `R/04_feature_engineering.R` |
| Database connectivity and data retrieval | 2 | `R/05_database_connectivity.R` |
| Implementation of 10-15 ML/DL algorithms | 3 | `R/06_model_training.R` |
| Hyperparameter tuning and model optimization | 1 | `R/07_hyperparameter_tuning.R` |
| Comparative performance analysis with evaluation metrics | 1 | `R/08_comparative_analysis.R` |
| Comparative visualizations | 1 | `R/09_visualizations.R` |
| Progress demonstration and documentation | 1 | This document + updated README |

## 3. Feature Engineering & Selection

Six descriptors were derived from the existing 14 raw features, following standard
cheminformatics conventions (no external data, no target leakage):

- `HBondTotal` = HBA + HBD
- `LipinskiViolations` - count of Rule-of-Five violations (drug-likeness)
- `AromaticRingFraction` = aromatic rings / total rings
- `FlexibilityIndex` = rotatable bonds / heavy atom count
- `HeteroatomFraction` = heteroatom count / heavy atom count
- `TPSA_per_MW` - polarity normalized by molecular size

Selection then removes:
- Near-zero-variance features (`caret::nearZeroVar`)
- Highly collinear features at `|r| > 0.90` (`caret::findCorrelation`) - notably
  `ExactMolWt` is redundant with `MolecularWeight` (r ≈ 0.9999), so one is dropped.

The final feature set and a Random-Forest-based importance ranking are written to
`results/tables/feature_selection_report.csv` and `feature_importance_filter.csv`.

## 4. Database Connectivity

A normalized SQLite schema (`data/db/rnainsight.sqlite`) replaces the flat CSV as the data
source for modelling:

- `molecules(mol_id, Name, SMILES, Target_Family)`
- `descriptors(mol_id, <selected descriptor columns>)`

`R/05_database_connectivity.R` connects via `DBI`/`RSQLite`, writes both tables, runs a
filtered `SELECT ... JOIN`, an aggregated `GROUP BY` query (mean descriptors per RNA family,
computed inside the database), and pulls the full joined table back into R as
`data/processed/model_ready.csv` - this is the file every downstream script consumes, so the
ML pipeline is genuinely fed from the database. The same `DBI` code works unchanged against
MySQL/PostgreSQL by swapping the driver package.

## 5. Model Development - 15 Algorithms

All models are trained with `caret::train()` on an identical stratified 80/20 split and an
identical fixed 5-fold x 3-repeat cross-validation scheme (same fold indices for every model,
saved to `results/models/cv_folds.rds`), so results are directly comparable:

1. Multinomial Logistic Regression (`multinom`)
2. Linear Discriminant Analysis (`lda`)
3. Naive Bayes (`naive_bayes`)
4. K-Nearest Neighbors (`knn`)
5. Decision Tree - CART (`rpart`)
6. Decision Tree - C5.0 (`C5.0`)
7. Bagged CART (`treebag`)
8. Random Forest (`rf`)
9. Gradient Boosting Machine (`gbm`)
10. XGBoost (`xgbTree`)
11. AdaBoost.M1 (`AdaBoost.M1`)
12. SVM - Linear Kernel (`svmLinear`)
13. SVM - Radial Kernel (`svmRadial`)
14. Regularized Multinomial / Elastic Net (`glmnet`)
15. Neural Network - single hidden layer (`nnet`)

Each is wrapped in `tryCatch` so a missing package or a single failing model does not stop
the run; `results/tables/model_training_log.csv` records status and training time per model.

## 6. Hyperparameter Tuning

The top 3 models by baseline CV Accuracy are re-tuned with a wider, method-specific grid
(e.g. `mtry` for Random Forest, `nrounds`/`max_depth`/`eta` for XGBoost, `sigma`/`C` for SVM
radial). `results/tables/hyperparameter_tuning_report.csv` reports baseline vs. tuned Accuracy
and the improvement for each.

## 7. Comparative Performance Analysis

`results/tables/model_comparison.csv` ranks all 15 models by test-set Accuracy, Kappa, and
macro-averaged Precision/Recall/F1 (macro-averaging matters here because the classes are
imbalanced, from 288 "Other" down to 8 "Viral RNA"). `caret::resamples()` + paired t-tests
(`diff()`) on the shared CV folds test whether Accuracy differences between models are
statistically meaningful (`results/tables/pairwise_accuracy_pvalues.csv`).

## 8. Comparative Visualizations

Saved under `results/plots/da2/`:
- Confusion matrix heatmaps for the top 3 models
- One-vs-rest ROC curves with per-class AUC (best model)
- One-vs-rest Precision-Recall curves (best model)
- Feature importance bar charts for the tree/ensemble models
- CV Accuracy distribution across all 15 models (box-and-whisker)
- Test-set Accuracy leaderboard bar chart

## 9. Known Limitations

- Class imbalance is severe (8 to 288 records per class); the smallest classes (`Viral RNA`,
  `Ribozyme`, `miRNA`) have only 1-3 samples in the 20% test split, so their individual
  Precision/Recall/AUC estimates are noisy. This is disclosed rather than hidden, and macro
  averages plus the confusion matrices make the effect visible.
- "Other" is a catch-all label (no specific assay target reported) rather than a true RNA
  family, which likely caps achievable accuracy on that class.

## 10. What's Left for the Remaining 25%

- Final model selection and justification (from the tuned leaderboard).
- Consolidating results into the final written report / slide deck.
- Optional: a small Shiny dashboard for interactive exploration (mentioned as optional in the
  course's mandatory repo structure).
- Final live demonstration.
