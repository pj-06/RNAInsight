############################################################
# RNAInsight Dashboard
# 02_data.R : load project outputs and pre-compute metrics
############################################################
#
# Everything the dashboard shows is read from files the
# existing pipeline already produced (scripts/ and R/ steps
# 03-09) - nothing is re-trained here:
#
#   data/processed/cleaned_dataset.csv     750 records, 14 descriptors
#   data/processed/model_ready.csv         modelling table (from SQLite)
#   results/tables/*.csv                   selection / importance tables
#   results/models/*.rds                   trained models + test split
#   dashboard/data/molecules.json          ligand structures + per-atom
#                                          descriptor breakdowns
#                                          (scripts/06_structure_export.py)
#
# Metrics are computed once and cached in dashboard/cache/ so the
# app starts quickly; the cache is rebuilt automatically when any
# input file is newer than it.
############################################################

#' The 14 RDKit descriptors generated in 03_descriptor_generation.py.
RAW_DESCRIPTORS <- c(
  "MolecularWeight", "ExactMolWt", "LogP", "MolMR", "TPSA",
  "HBA", "HBD", "RotatableBonds", "RingCount", "NumAromaticRings",
  "NumAliphaticRings", "HeavyAtomCount", "NumHeteroAtoms", "FractionCSP3"
)

#' The 6 features engineered in R/04_feature_engineering.R.
ENGINEERED_FEATURES <- c(
  "HBondTotal", "LipinskiViolations", "AromaticRingFraction",
  "FlexibilityIndex", "HeteroatomFraction", "TPSA_per_MW"
)

#' Human-readable descriptor labels (with units where they exist).
FEATURE_LABELS <- c(
  MolecularWeight      = "Molecular weight (g/mol)",
  ExactMolWt           = "Exact (monoisotopic) mass (Da)",
  LogP                 = "LogP (Crippen)",
  MolMR                = "Molar refractivity",
  TPSA                 = "TPSA (Å²)",
  HBA                  = "H-bond acceptors",
  HBD                  = "H-bond donors",
  RotatableBonds       = "Rotatable bonds",
  RingCount            = "Rings",
  NumAromaticRings     = "Aromatic rings",
  NumAliphaticRings    = "Aliphatic rings",
  HeavyAtomCount       = "Heavy atoms",
  NumHeteroAtoms       = "Heteroatoms",
  FractionCSP3         = "Fraction sp³ carbons",
  HBondTotal           = "H-bond total (HBA + HBD)",
  LipinskiViolations   = "Lipinski violations",
  AromaticRingFraction = "Aromatic ring fraction",
  FlexibilityIndex     = "Flexibility index",
  HeteroatomFraction   = "Heteroatom fraction",
  TPSA_per_MW          = "TPSA per MW"
)

#' Return a readable label for a feature name.
#'
#' @param x Character vector of feature names.
#' @return Character vector of labels (falls back to the name).
feature_label <- function(x) {
  out <- unname(FEATURE_LABELS[x])
  ifelse(is.na(out), x, out)
}

#' Add the 6 engineered features exactly as 04_feature_engineering.R does.
#'
#' @param data Data frame with the 14 raw descriptors.
#' @return The data frame with the engineered columns appended.
add_engineered_features <- function(data) {
  data$HBondTotal <- data$HBA + data$HBD
  data$LipinskiViolations <-
    as.integer(data$MolecularWeight > 500) +
    as.integer(data$LogP > 5) +
    as.integer(data$HBA > 10) +
    as.integer(data$HBD > 5)
  data$AromaticRingFraction <- ifelse(data$RingCount == 0, 0, data$NumAromaticRings / data$RingCount)
  data$FlexibilityIndex <- data$RotatableBonds / data$HeavyAtomCount
  data$HeteroatomFraction <- data$NumHeteroAtoms / data$HeavyAtomCount
  data$TPSA_per_MW <- data$TPSA / data$MolecularWeight
  data
}

# ----------------------------------------------------------
# Classification metrics
# ----------------------------------------------------------

#' Evaluate one fitted caret model on the held-out test set.
#'
#' @param fit A caret train object.
#' @param test_data Test data frame (Target_Family + predictors).
#' @param pretty Function mapping safe class labels to readable names.
#' @return list(overall = data.frame(Accuracy), confusion = data.frame)
#'         or NULL when the model cannot predict (missing package).
evaluate_model <- function(fit, test_data, pretty) {
  truth <- test_data$Target_Family
  pred <- tryCatch(predict(fit, newdata = test_data), error = function(e) NULL)
  if (is.null(pred)) {
    return(NULL)
  }
  pred <- factor(pred, levels = levels(truth))

  cm <- caret::confusionMatrix(pred, truth)
  tab <- as.data.frame(cm$table, stringsAsFactors = FALSE)
  names(tab) <- c("Predicted", "Actual", "Count")
  tab$Predicted <- pretty(tab$Predicted)
  tab$Actual <- pretty(tab$Actual)

  list(overall = data.frame(Accuracy = unname(cm$overall["Accuracy"])), confusion = tab)
}

#' Evaluate every model on the held-out test set.
#'
#' @param models Named list of caret train objects.
#' @param test_data Test data frame.
#' @param pretty Label prettifier.
#' @return Named list of per-model evaluations, plus their CV resamples.
evaluate_experiment <- function(models, test_data, pretty) {
  models <- models[intersect(MODEL_ORDER, names(models))]
  lapply(stats::setNames(names(models), names(models)), function(name) {
    fit <- models[[name]]
    eval <- evaluate_model(fit, test_data, pretty)
    if (is.null(eval)) {
      return(NULL)
    }
    eval$resample <- fit$resample
    eval
  })
}

#' Model-based variable importance for the tree / ensemble models.
#'
#' @param models Named list of caret train objects.
#' @param pretty Label prettifier (for per-class columns).
#' @return Long data frame: Model, Feature, Class, Importance.
compute_importance <- function(models, pretty) {
  keep <- intersect(c("Random Forest", "Gradient Boosting Machine", "Decision Tree (C5.0)"), names(models))
  # caret's varImp() for gbm calls gbm::relative.influence() unqualified,
  # so the package has to be attached, not just installed.
  if ("Gradient Boosting Machine" %in% keep && requireNamespace("gbm", quietly = TRUE)) {
    suppressPackageStartupMessages(library(gbm))
  }
  out <- lapply(keep, function(name) {
    imp <- tryCatch(caret::varImp(models[[name]])$importance, error = function(e) NULL)
    if (is.null(imp)) {
      return(NULL)
    }
    imp$Feature <- rownames(imp)
    score_cols <- setdiff(names(imp), "Feature")
    overall <- rowMeans(imp[, score_cols, drop = FALSE], na.rm = TRUE)
    data.frame(Model = name, Feature = imp$Feature, Importance = overall)
  })
  do.call(rbind, out)
}

# ----------------------------------------------------------
# Ligand-level tables (152 unique structures)
# ----------------------------------------------------------

#' Load the per-ligand structure records exported by Python.
#'
#' @param path Path to molecules.json.
#' @return list(json = parsed list, table = one row per ligand)
load_ligands <- function(path) {
  mols <- jsonlite::fromJSON(path, simplifyVector = FALSE)

  table <- do.call(rbind, lapply(mols, function(m) {
    families <- vapply(m$records, `[[`, character(1), "family")
    rows <- vapply(m$records, `[[`, numeric(1), "row")
    specific <- setdiff(families, "Other")
    primary <- if (length(specific)) names(sort(table(specific), decreasing = TRUE))[1] else "Other"
    desc <- as.data.frame(lapply(m$descriptors, as.numeric))
    data.frame(
      key = m$key,
      Name = m$name,
      SMILES = m$smiles,
      Formula = m$formula,
      Primary_Family = primary,
      Families = paste(sort(unique(families)), collapse = ", "),
      N_Families = length(unique(families)),
      N_Records = length(rows),
      First_Row = min(rows),
      Rows = I(list(rows)),
      desc,
      stringsAsFactors = FALSE
    )
  }))
  table <- add_engineered_features(table)

  list(json = mols, table = table)
}

#' How many ligands are shared between each pair of RNA families.
#'
#' @param cleaned Cleaned 750-record data with Target_Family and SMILES.
#' @return Long data frame: Family_A, Family_B, Shared.
family_overlap <- function(cleaned) {
  fams <- sort(unique(cleaned$Target_Family))
  sets <- lapply(stats::setNames(fams, fams), function(f) unique(cleaned$SMILES[cleaned$Target_Family == f]))
  grid <- expand.grid(Family_A = fams, Family_B = fams, stringsAsFactors = FALSE)
  grid$Shared <- mapply(function(a, b) length(intersect(sets[[a]], sets[[b]])), grid$Family_A, grid$Family_B)
  grid
}

# ----------------------------------------------------------
# Main loader (with cache)
# ----------------------------------------------------------

#' Folder (relative to the project root) holding the models shown in the
#' dashboard, as saved by R/06_model_training.R.
MODEL_DIR <- "results/models"

#' Files the dashboard depends on (relative to the project root).
dashboard_inputs <- function(root, app_dir) {
  c(
    file.path(root, "data/processed/cleaned_dataset.csv"),
    file.path(root, "data/processed/model_ready.csv"),
    file.path(root, "results/tables/target_family_label_map.csv"),
    file.path(root, "results/tables/feature_selection_report.csv"),
    file.path(root, "results/tables/feature_importance_filter.csv"),
    file.path(root, MODEL_DIR, "trained_models.rds"),
    file.path(root, MODEL_DIR, "test_data.rds"),
    file.path(app_dir, "data/molecules.json")
  )
}

#' Build every table and metric the dashboard displays.
#'
#' @param root Project root folder (DS_project).
#' @param app_dir Dashboard folder.
#' @return A named list consumed by the UI modules.
build_dashboard_data <- function(root, app_dir) {
  inputs <- dashboard_inputs(root, app_dir)
  missing <- inputs[!file.exists(inputs)]
  if (length(missing)) {
    stop(
      "The dashboard needs these pipeline outputs, which were not found:\n  ",
      paste(missing, collapse = "\n  "),
      "\nRun scripts/06_structure_export.py and R/04-06 first (see dashboard/README.md).",
      call. = FALSE
    )
  }

  label_map <- utils::read.csv(inputs[3], stringsAsFactors = FALSE)
  pretty <- function(x) {
    x <- as.character(x)
    idx <- match(x, label_map$Safe)
    ifelse(is.na(idx), x, label_map$Original[idx])
  }

  cleaned <- utils::read.csv(inputs[1], stringsAsFactors = FALSE)
  cleaned <- add_engineered_features(cleaned)
  cleaned$Assay <- trimws(cleaned$Assay.Target.s.)
  family_counts <- sort(table(cleaned$Target_Family), decreasing = TRUE)
  cleaned$Target_Family <- factor(cleaned$Target_Family, levels = names(family_counts))


  models <- readRDS(inputs[6])
  test_data <- readRDS(inputs[7])

  # KNN breaks voting ties at random; a fixed seed keeps the metrics
  # reproducible from one dashboard build to the next.
  set.seed(42)
  evaluation <- evaluate_experiment(models, test_data, pretty)

  leaderboard <- do.call(rbind, lapply(names(evaluation), function(m) {
    ev <- evaluation[[m]]
    if (is.null(ev)) return(NULL)
    cbind(Model = m, ev$overall, CV_Accuracy = mean(ev$resample$Accuracy),
          CV_Accuracy_SD = stats::sd(ev$resample$Accuracy))
  }))

  ligands <- load_ligands(inputs[8])
  molecules <- lapply(ligands$json, prepare_molecule)
  names(molecules) <- vapply(molecules, function(m) as.character(m$key), character(1))
  test_rows <- as.integer(rownames(test_data))
  ligands$table$In_Test <- vapply(ligands$table$Rows, function(r) any(r %in% test_rows), logical(1))
  ligands$table$Test_Records <- vapply(ligands$table$Rows, function(r) sum(r %in% test_rows), numeric(1))

  list(
    built_at = Sys.time(),
    families = names(family_counts),
    family_counts = family_counts,
    cleaned = cleaned,
    label_map = label_map,
    selection = utils::read.csv(inputs[4], stringsAsFactors = FALSE),
    evaluation = evaluation,
    leaderboard = leaderboard,
    importance = compute_importance(models, pretty),
    ligands = ligands$table,
    molecules = molecules,
    overlap = family_overlap(cleaned),
    test_size = nrow(test_data),
    train_size = nrow(cleaned) - nrow(test_data)
  )
}

#' Load the dashboard data, rebuilding the cache only when needed.
#'
#' @param root Project root folder.
#' @param app_dir Dashboard folder.
#' @return The list produced by build_dashboard_data().
load_dashboard_data <- function(root, app_dir) {
  cache_dir <- file.path(app_dir, "cache")
  cache_file <- file.path(cache_dir, "dashboard_data.rds")
  inputs <- dashboard_inputs(root, app_dir)
  inputs <- inputs[file.exists(inputs)]

  cache_ok <- file.exists(cache_file) &&
    length(inputs) > 0 &&
    all(file.mtime(inputs) <= file.mtime(cache_file)) &&
    file.mtime(cache_file) >= file.mtime(file.path(app_dir, "R", "02_data.R"))

  if (cache_ok) {
    cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    if (!is.null(cached)) {
      return(cached)
    }
  }

  message("RNAInsight dashboard: computing metrics from pipeline outputs (first run only)...")
  data <- build_dashboard_data(root, app_dir)
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  try(saveRDS(data, cache_file), silent = TRUE)
  data
}
