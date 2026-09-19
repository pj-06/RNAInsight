############################################################
# RNAInsight
# Step 5 : Database Connectivity & Data Retrieval (DA2)
############################################################
#
# Purpose:
#   Persist the engineered/selected dataset into a relational
#   database (SQLite) using a small normalised schema, then
#   demonstrate SQL-based retrieval (filter, join, aggregate)
#   and read the result back into R as the modelling dataframe.
#
#   SQLite is used because it is serverless and file-based, so
#   the project stays fully reproducible for anyone who clones
#   the repo (no DB server to install/configure). The same DBI
#   code would work unchanged against MySQL/PostgreSQL by
#   swapping the driver (RMySQL / RPostgres).
#
# Input  : data/processed/model_dataset.csv
# Output : data/db/rnainsight.sqlite
#          results/tables/db_query_demo.csv
#          data/processed/model_ready.csv (data pulled back from the DB)
############################################################

# ---- Package bootstrap -------------------------------------------------
required_packages <- c("DBI", "RSQLite", "dplyr")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing missing package:", pkg, "\n")
    install.packages(pkg, dependencies = TRUE)
  }
}

library(DBI)
library(RSQLite)
library(dplyr)

############################################################
# Create Output Folders
############################################################

dir.create("data/db", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

############################################################
# Load Dataset
############################################################

data <- read.csv(
  "data/processed/model_dataset.csv",
  stringsAsFactors = FALSE
)

data$mol_id <- seq_len(nrow(data))

descriptor_cols <- setdiff(names(data), c("mol_id", "Name", "SMILES", "Target_Family"))

molecules_tbl <- data %>%
  select(mol_id, Name, SMILES, Target_Family)

descriptors_tbl <- data %>%
  select(mol_id, all_of(descriptor_cols))

cat("========================================\n")
cat("RNAInsight - Database Connectivity\n")
cat("========================================\n\n")
cat("molecules table   :", nrow(molecules_tbl), "rows\n")
cat("descriptors table :", nrow(descriptors_tbl), "rows,", length(descriptor_cols), "descriptor columns\n\n")

############################################################
# 1. Connect & Write Schema
############################################################

db_path <- "data/db/rnainsight.sqlite"

# Best-effort cleanup only -- dbWriteTable(overwrite = TRUE) below already
# guarantees fresh tables, so a locked file here (e.g. Windows still
# holding a handle from a previous run) is not fatal.
if (file.exists(db_path)) {
  try(suppressWarnings(file.remove(db_path)), silent = TRUE)
}

con <- dbConnect(RSQLite::SQLite(), db_path)

dbWriteTable(con, "molecules", molecules_tbl, overwrite = TRUE)
dbWriteTable(con, "descriptors", descriptors_tbl, overwrite = TRUE)

# Helpful index for the join used below
dbExecute(con, "CREATE INDEX idx_descriptors_mol_id ON descriptors(mol_id);")

cat("Connected to SQLite database:", db_path, "\n")
cat("Tables created: molecules, descriptors\n\n")

############################################################
# 2. Retrieval Demo (filter + join + aggregate)
############################################################
# The exact set of surviving descriptor columns depends on what
# 04_feature_engineering.R's collinearity filter dropped (e.g.
# MolecularWeight is often removed as redundant with HeavyAtomCount).
# Rather than hardcode column names that might not exist, pick a
# small display set from whatever actually survived.

display_priority <- c("LogP", "TPSA", "MolMR", "RotatableBonds", "HBA", "HBD")
display_cols <- intersect(display_priority, descriptor_cols)
if (length(display_cols) < 2) {
  display_cols <- head(descriptor_cols, 3)
}
display_cols <- head(display_cols, 3)
order_col <- display_cols[1]

cat("Using descriptor columns for the demo queries:", paste(display_cols, collapse = ", "), "\n\n")

# 2a. Simple filtered retrieval
select_cols_sql <- paste(paste0("d.", display_cols), collapse = ", ")

riboswitch_query <- dbGetQuery(con, sprintf("
  SELECT m.Name, m.Target_Family, %s
  FROM molecules m
  JOIN descriptors d ON m.mol_id = d.mol_id
  WHERE m.Target_Family = 'Riboswitch'
  ORDER BY d.%s DESC
  LIMIT 10;
", select_cols_sql, order_col))

cat("Top 10 Riboswitch-binding molecules (by", order_col, "):\n")
print(riboswitch_query)

# 2b. Aggregate retrieval: mean descriptors per target family, via SQL
agg_cols_sql <- paste(
  sprintf("ROUND(AVG(d.%s), 2) AS avg_%s", display_cols, tolower(display_cols)),
  collapse = ",\n    "
)

family_summary_query <- dbGetQuery(con, sprintf("
  SELECT
    m.Target_Family,
    COUNT(*) AS n_molecules,
    %s
  FROM molecules m
  JOIN descriptors d ON m.mol_id = d.mol_id
  GROUP BY m.Target_Family
  ORDER BY n_molecules DESC;
", agg_cols_sql))

cat("\nPer-family summary computed inside the database:\n")
print(family_summary_query)

write.csv(
  family_summary_query,
  "results/tables/db_query_demo.csv",
  row.names = FALSE
)

############################################################
# 3. Full Retrieval for Modelling
############################################################
# Pull the full joined table back into R -- this is the
# dataframe that 06_model_training.R consumes, proving the
# ML pipeline is fed via the database rather than the raw CSV.

model_ready <- dbGetQuery(con, "
  SELECT m.Target_Family, d.*
  FROM molecules m
  JOIN descriptors d ON m.mol_id = d.mol_id;
")

model_ready$mol_id <- NULL

write.csv(
  model_ready,
  "data/processed/model_ready.csv",
  row.names = FALSE
)

cat("\nFull modelling table retrieved from DB:", nrow(model_ready), "rows,", ncol(model_ready), "columns\n")
cat("Saved to: data/processed/model_ready.csv\n")

############################################################
# 4. Disconnect
############################################################

dbDisconnect(con)

cat("\n========================================\n")
cat("Database connectivity demo complete!\n")
cat("Database file : ", db_path, "\n")
cat("========================================\n")
