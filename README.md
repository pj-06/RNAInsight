RNAInsight: Machine Learning-Based RNA Target Classification using Molecular Descriptors

The goal is:

Can we predict which RNA target family a small molecule is likely to bind to using only its molecular descriptors?

This is a classification problem.

1. Raw R-BIND Excel Dataset
2. Python Phase (Data Preparation)
3. Cleaned Dataset (750 molecules)
4. R Phase (EDA + Statistics + ML + Dashboard)

## Dataset

- Source: R-BIND (RNA-targeted BIoactive ligaNd Database) supplementary Excel file.
- After cleaning/deduplication: 750 molecule-target records, 11 RNA target family classes
  (Other, Repeat RNA, G-Quadruplex, HIV TAR, mRNA, Riboswitch, HCV IRES, tRNA, miRNA, Ribozyme, Viral RNA).
- Features: 14 molecular descriptors (MolecularWeight, LogP, TPSA, HBA/HBD, ring counts, etc.).

## Project Structure

```
DS_project/
├── scripts/            Python data-preparation phase (01-05)
├── R/                  R analysis phase
│   ├── 03_eda.R                    DA1: exploratory data analysis
│   ├── 04_feature_engineering.R    DA2: feature engineering + selection
│   ├── 05_database_connectivity.R  DA2: SQLite schema, load, retrieval
│   ├── 06_model_training.R         DA2: 15 ML/DL algorithms (baseline)
│   ├── 07_hyperparameter_tuning.R  DA2: tuning of top-3 models
│   ├── 08_comparative_analysis.R   DA2: metrics + statistical comparison
│   └── 09_visualizations.R         DA2: ROC / PR / confusion / importance plots
├── data/
│   ├── raw/             Original R-BIND Excel file
│   ├── processed/       Cleaned + engineered + DB-retrieved CSVs
│   └── db/              SQLite database (rnainsight.sqlite)
├── results/
│   ├── plots/            DA1 EDA plots
│   ├── plots/da2/        DA2 comparative visualizations
│   ├── tables/           Summary stats, feature selection, model comparison
│   └── models/           Serialized trained/tuned model objects (.rds)
└── docs/
    └── DA2_progress_report.md
```

## How to Run (DA2)

Open `DS_project.Rproj` in RStudio, then run in order from the project root:

```r
source("R/04_feature_engineering.R")
source("R/05_database_connectivity.R")
source("R/06_model_training.R")        # trains 15 algorithms (can take 15-30 min)
source("R/07_hyperparameter_tuning.R")
source("R/08_comparative_analysis.R")
source("R/09_visualizations.R")
```

Each script installs any R packages it needs on first run (`caret`, `DBI`, `RSQLite`,
`pROC`, `PRROC`, `xgboost`, `randomForest`, `gbm`, `C50`, `adabag`, `glmnet`, `nnet`, etc.).

See `docs/DA2_progress_report.md` for methodology, results and current completion status.
