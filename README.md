# RNAInsight

**Predicting the RNA target family of small molecules from their molecular descriptors.**

RNAInsight takes RNA-binding small molecules (ligands) from the R-BIND database, describes each one with
14 RDKit molecular descriptors, and trains classifiers that predict which of 11 RNA target families the
ligand binds. An R Shiny dashboard presents the results and lets you explore every ligand in 3D/2D,
showing how each descriptor is calculated atom by atom.

> Course project, Programming for Data Science, VIT Chennai.

---

## Highlights

- **Data:** 750 molecule-target records (152 unique ligands) across 11 RNA target families, from the R-BIND supplementary dataset.
- **Features:** 14 RDKit descriptors, plus 6 engineered in R. 15 features kept after near-zero-variance and correlation filtering.
- **Storage:** a normalised SQLite database (`data/db/rnainsight.sqlite`) feeds the modelling step.
- **Models:** C5.0, Random Forest, Gradient Boosting, k-Nearest Neighbours and Multinomial Logistic Regression, all trained with `caret` on the same stratified 80/20 split and the same 5-fold x 3 cross-validation folds.
- **Dashboard:** an interactive R Shiny app with results charts and a Ligand Structure Explorer.

## Results

Accuracy on the 145 held-out test records, with the mean accuracy over the 15 cross-validation folds:

| Model | Test accuracy | CV accuracy |
|---|---|---|
| Decision Tree (C5.0) | 0.628 | 0.545 |
| Gradient Boosting Machine | 0.614 | 0.561 |
| Random Forest | 0.600 | 0.577 |
| k-Nearest Neighbours | 0.566 | 0.538 |
| Multinomial Logistic Regression | 0.510 | 0.489 |

Most errors come from the catch-all *Other* class. The same ligand often appears under several RNA
families, so identical descriptor rows can carry different labels.

## Pipeline

```
R-BIND Excel ─► Python + RDKit ─► Cleaning ─► Feature engineering ─► SQLite ─► caret models ─► Shiny dashboard
 data/raw        scripts/01-03     scripts/04-05   R/03-04            R/05      R/06-09        dashboard/
```

| Step | Script | Output |
|---|---|---|
| Dataset audit | `scripts/01_dataset_audit.py` | `results/tables/all_targets.csv` |
| Target standardisation | `scripts/02_target_standardization.py` | `data/processed/target_standardized.csv` |
| Descriptor generation (RDKit) | `scripts/03_descriptor_generation.py` | `data/processed/descriptor_dataset.csv` |
| Cleaning | `scripts/04_data_cleaning.py` | `data/processed/cleaned_dataset.csv` |
| Validation | `scripts/05_dataset_validation.py` | `results/reports/validation_report.txt` |
| Structures for the dashboard | `scripts/06_structure_export.py` | `dashboard/data/molecules.json` |
| Exploratory data analysis | `R/03_eda.R` | `results/figures/eda/`, `results/tables/summary_statistics.csv` |
| Feature engineering and selection | `R/04_feature_engineering.R` | `data/processed/model_dataset.csv` |
| Database connectivity (SQLite) | `R/05_database_connectivity.R` | `data/db/rnainsight.sqlite`, `data/processed/model_ready.csv` |
| Model training | `R/06_model_training.R` | `results/models/`, `results/tables/model_training_log.csv` |
| Hyperparameter tuning | `R/07_hyperparameter_tuning.R` | `results/models/tuned_models.rds` |
| Comparative analysis | `R/08_comparative_analysis.R` | `results/tables/model_comparison.csv` |
| Visualisations | `R/09_visualizations.R` | `results/figures/models/` |

## Repository structure

```
RNAInsight/
├── README.md
├── RNAInsight.Rproj          RStudio project (open this)
├── requirements.txt          Python dependencies
├── scripts/                  Python: data preparation and descriptor generation (01-06)
├── R/                        R: EDA, feature engineering, database, modelling (03-09)
├── dashboard/                R Shiny dashboard (see dashboard/README.md)
│   ├── app.R
│   ├── R/                    dashboard modules
│   └── data/molecules.json   ligand structures + per-atom descriptor breakdowns
├── data/
│   ├── raw/                  original R-BIND Excel file
│   ├── processed/            cleaned, engineered and database-retrieved datasets
│   └── db/                   SQLite database
├── results/
│   ├── figures/              EDA and model figures
│   ├── models/               trained caret models, train/test split, CV folds (.rds)
│   ├── tables/               summary statistics, feature selection, model results
│   └── reports/              dataset validation report
└── docs/
    ├── DA2_progress_report.md
    └── presentations/        review slide decks
```

## Getting started

### Requirements

- **R** 4.2 or later (RStudio recommended)
- **Python** 3.10 or later, needed only to re-run the data-preparation scripts

### 1. Run the dashboard

The repository already contains every processed dataset and trained model, so the dashboard runs straight
away. Open `RNAInsight.Rproj` in RStudio and run:

```r
shiny::runApp("dashboard")
```

Or from a terminal in the project folder:

```bash
Rscript -e "shiny::runApp('dashboard', launch.browser = TRUE)"
```

Missing R packages (`shiny`, `bslib`, `plotly`, `ggplot2`, `DT`, `jsonlite`, `r3dmol`, `caret`,
`randomForest`, `gbm`, `C50`, `nnet`) are installed automatically on the first run.

### 2. Reproduce the pipeline (optional)

Python data preparation:

```bash
python -m venv .venv
.venv\Scripts\activate          # Windows  (macOS/Linux: source .venv/bin/activate)
pip install -r requirements.txt

python scripts/01_dataset_audit.py
python scripts/02_target_standardization.py
python scripts/03_descriptor_generation.py
python scripts/04_data_cleaning.py
python scripts/05_dataset_validation.py
python scripts/06_structure_export.py
```

R analysis and modelling, run from the project root in this order:

```r
source("R/03_eda.R")
source("R/04_feature_engineering.R")
source("R/05_database_connectivity.R")
source("R/06_model_training.R")
source("R/07_hyperparameter_tuning.R")
source("R/08_comparative_analysis.R")
source("R/09_visualizations.R")
```

Each R script installs the packages it needs on first run.

## Dashboard

| Section | Contents |
|---|---|
| Overview | Key figures, the pipeline, records vs unique ligands, and ligands shared between families |
| Data Exploration | Descriptor distributions and an interactive ligand scatter plot |
| Feature Engineering | Feature selection (correlation filter) and model-based feature importance |
| Model Comparison | Cross-validated vs test accuracy, and the leaderboard |
| Class Diagnostics | Confusion matrix of each model |
| Ligand Structure Explorer | 3D (r3dmol) and 2D (ggplot2) structure of every ligand; click a descriptor to see which atoms, bonds or rings it counts and how its value is calculated |

The dashboard is written entirely in R. See [`dashboard/README.md`](dashboard/README.md) for details.

## Data source

R-BIND (RNA-targeted BIoactive ligaNd Database), supplementary dataset
`NIHMS1824842-supplement-RBIND_Info_A.xlsx`. The 14 descriptors are computed with RDKit from each
ligand's SMILES; the dataset gives each ligand's RNA target family, not the RNA structure itself.
