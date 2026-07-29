import pandas as pd
from pathlib import Path

# ---------------------------------
# Paths
# ---------------------------------

BASE_DIR = Path(__file__).resolve().parent.parent

input_file = BASE_DIR / "data" / "processed" / "cleaned_dataset.csv"
results_file = BASE_DIR / "results" / "validation_report.txt"

# ---------------------------------
# Load Dataset
# ---------------------------------

df = pd.read_csv(input_file)

print("=" * 50)
print("DATASET VALIDATION REPORT")
print("=" * 50)

print(f"\nRows: {df.shape[0]}")
print(f"Columns: {df.shape[1]}")

# ---------------------------------
# Missing Values
# ---------------------------------

print("\nMissing Values:")
print(df.isnull().sum())

# ---------------------------------
# Duplicate Rows
# ---------------------------------

duplicates = df.duplicated().sum()
print(f"\nDuplicate Rows: {duplicates}")

# ---------------------------------
# Target Family Distribution
# ---------------------------------

print("\nTarget Family Distribution:")
print(df["Target_Family"].value_counts())

# ---------------------------------
# Summary Statistics
# ---------------------------------

print("\nSummary Statistics:")
print(df.describe())

# ---------------------------------
# Save Report
# ---------------------------------

with open(results_file, "w") as f:

    f.write("DATASET VALIDATION REPORT\n")
    f.write("=" * 50 + "\n\n")

    f.write(f"Rows: {df.shape[0]}\n")
    f.write(f"Columns: {df.shape[1]}\n\n")

    f.write("Missing Values\n")
    f.write(df.isnull().sum().to_string())

    f.write("\n\nDuplicate Rows\n")
    f.write(str(duplicates))

    f.write("\n\nTarget Family Distribution\n")
    f.write(df["Target_Family"].value_counts().to_string())

    f.write("\n\nSummary Statistics\n")
    f.write(df.describe().to_string())

print("\nValidation report saved successfully!")
print(results_file)