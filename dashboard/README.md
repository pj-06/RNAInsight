# RNAInsight Dashboard

An R Shiny dashboard with two sections:

1. **Results & Visualizations**: the pipeline's results as interactive charts (ggplot2, plotly, DT), in five tabs. Each chart's dropdowns sit in its own card header:
   - **Overview**: KPIs, the pipeline, records vs unique ligands, and ligands shared between RNA families.
   - **Data Exploration**: descriptor distribution and a ligand scatter (click a point to open it in the explorer). One family can be highlighted against the rest.
   - **Feature Engineering**: why 15 of 20 features were kept, model-based feature importance, and the six engineered features.
   - **Model Comparison**: cross-validated vs test accuracy for every model, and the leaderboard.
   - **Class Diagnostics**: the confusion matrix of any model (counts or row %).
2. **Ligand Structure Explorer**: the 3D structure (r3dmol) and 2D structure (ggplot2) of every ligand. Pick any of the 14 RDKit descriptors to light up the atoms, bonds or rings it counts and get a step-by-step calculation. Every breakdown adds up to the value stored in `data/processed/cleaned_dataset.csv`. The explorer also shows where the ligand sits in the dataset, the RNA targets it was tested on, and the most similar ligands.

The dashboard is written entirely in R: Shiny, bslib, ggplot2, plotly, DT and r3dmol. Styles are kept in `R/08_styles.R`.

> The 14 descriptors are computed from the **ligand's** SMILES (`scripts/03_descriptor_generation.py`).
> The dataset has no RNA structures, only the RNA target family/assay label, so the explorer
> visualises the small molecule that binds the RNA.

## Run it

From the project root (`DS_project/`):

```r
shiny::runApp("dashboard")
```

or from Command Prompt:

```cmd
cd C:\Users\priya\Documents\DS_project
Rscript -e "shiny::runApp('dashboard', launch.browser = TRUE)"
```

Missing R packages are installed on the first run. The first start computes all metrics from `results/`
(a few seconds) and caches them in `dashboard/cache/`; the cache is rebuilt automatically whenever a model
or data file changes.

## Files

| File | Purpose |
|---|---|
| `app.R` | Entry point: page layout, navigation, shared ligand selection |
| `R/01_theme.R` | Colour tokens, `theme_rnainsight()` for ggplot2, plotly styling |
| `R/02_data.R` | Loads pipeline outputs (models from `results/models`), computes accuracies and confusion matrices, caching |
| `R/03_plots.R` | Every chart of the Results section (pure functions, usable from the console) |
| `R/04_mod_results.R` | Results & Visualizations module (UI + server) |
| `R/05_mod_explorer.R` | Ligand Structure Explorer module (UI + server) |
| `R/07_structure.R` | Descriptor highlights, step-by-step calculations, 2D drawing (ggplot2) and 3D viewer (r3dmol) |
| `R/08_styles.R` | Page styles (dark theme, cards, explorer layout) |
| `data/molecules.json` | Ligand structures + per-atom descriptor breakdowns, produced by `scripts/06_structure_export.py` |

To rebuild any chart outside the app (e.g. for the report):

```r
setwd("dashboard")
for (f in list.files("R", full.names = TRUE)) source(f)
d <- load_dashboard_data("..", ".")
plot_confusion(d, "Decision Tree (C5.0)", normalise = TRUE)
```
