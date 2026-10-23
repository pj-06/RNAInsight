############################################################
# RNAInsight - Interactive Dashboard (Review 3)
############################################################
#
# Two sections:
#   1. Results & Visualizations  - EDA, feature engineering,
#      model comparison (accuracy) and confusion matrices, built
#      from the outputs of scripts/01-06 and R/03-06.
#   2. Ligand Structure Explorer - interactive 3D (r3dmol) and
#      2D (ggplot2) structure of every ligand, with each of the
#      14 RDKit descriptors explained atom by atom.
#
# The dashboard is written entirely in R (Shiny, bslib, ggplot2,
# plotly, DT, r3dmol); styles live in R/08_styles.R.
#
# Run from the project root (DS_project/):
#   shiny::runApp("dashboard")
# or open dashboard/app.R in RStudio and click "Run App".
#
# Ligand structures are read from data/molecules.json, a data file
# produced by the project's data-preparation step
# (scripts/06_structure_export.py); the dashboard only reads it.
############################################################

# ---- Package bootstrap (same convention as R/04-09) -------------------
required_packages <- c(
  "shiny", "bslib", "plotly", "ggplot2", "DT", "jsonlite", "htmltools", "r3dmol",
  "caret", "randomForest", "gbm", "C50", "nnet"
)
# Rscript has no CRAN mirror set by default, so set one when missing.
if (is.null(getOption("repos")) || getOption("repos")[["CRAN"]] %in% c("@CRAN@", "")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
  install.packages(missing_packages, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(plotly)
  library(ggplot2)
})

# ---- Locations ----------------------------------------------------------
# shiny::runApp() sets the working directory to this folder, so the
# project root (with data/ and results/) is its parent.
APP_DIR <- normalizePath(getwd(), winslash = "/")
PROJECT_ROOT <- normalizePath(file.path(APP_DIR, ".."), winslash = "/")

for (f in sort(list.files(file.path(APP_DIR, "R"), pattern = "\\.R$", full.names = TRUE))) {
  source(f, local = TRUE, encoding = "UTF-8")
}

# ---- Data (computed once, cached in dashboard/cache) --------------------
DASH <- load_dashboard_data(PROJECT_ROOT, APP_DIR)

DEFAULT_LIGAND <- local({
  preferred <- DASH$ligands$key[DASH$ligands$Name == "Ribocil-C"]
  if (length(preferred)) preferred[1] else DASH$ligands$key[which.max(DASH$ligands$N_Records)]
})

# ---- UI -----------------------------------------------------------------
ui <- bslib::page_navbar(
  id = "main_nav",
  title = htmltools::div(
    class = "ri-brand",
    htmltools::span(class = "ri-brand-mark", "RNA"),
    htmltools::span(class = "ri-brand-name", "Insight")
  ),
  window_title = "RNAInsight Dashboard",
  theme = ri_bs_theme(),
  fillable = FALSE,
  header = htmltools::tags$head(ri_style_tag()),
  bslib::nav_panel(
    "Results & Visualizations", value = "results", icon = NULL,
    htmltools::div(class = "ri-page", results_ui("results", DASH))
  ),
  bslib::nav_panel(
    "Ligand Structure Explorer", value = "explorer",
    htmltools::div(class = "ri-page ri-page-explorer", explorer_ui("explorer", DASH, DEFAULT_LIGAND))
  ),
  bslib::nav_spacer(),
  bslib::nav_item(htmltools::tags$a(
    class = "ri-nav-link", href = "https://github.com/pj-06/RNAInsight", target = "_blank", "GitHub ↗"
  ))
)

# ---- Server -------------------------------------------------------------
server <- function(input, output, session) {
  selected_ligand <- shiny::reactiveVal(DEFAULT_LIGAND)

  open_ligand <- function(key) {
    if (!is.na(key) && key %in% DASH$ligands$key) {
      selected_ligand(key)
      bslib::nav_select("main_nav", "explorer")
    }
  }

  results_server("results", DASH, open_ligand)
  explorer_server("explorer", DASH, selected_ligand)
}

shiny::shinyApp(ui, server)
