import pandas as pd
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
file = BASE_DIR / "data" / "raw" / "NIHMS1824842-supplement-RBIND_Info_A.xlsx"

df = pd.read_excel(file, header=1)

# Remove extra spaces from column names
df.columns = df.columns.str.strip()

# Remove empty columns
df = df.dropna(axis=1, how="all")
print("=" * 60)
print("DATASET AUDIT")
print("=" * 60)

print(f"\nDataset Shape : {df.shape}")

print(f"\nUnique Ligands : {df['Name'].nunique()}")

print(f"Unique SMILES : {df['SMILES'].nunique()}")

print(f"Unique Assay Targets : {df['Assay Target(s)'].nunique()}")

print("\nTop 20 Assay Targets\n")
print(df["Assay Target(s)"].value_counts().head(20))

print("\nMissing Values\n")
print(df.isnull().sum())

# Show all unique assay targets
targets = sorted(df["Assay Target(s)"].dropna().unique())

print("\nAll Unique Assay Targets:\n")

for target in targets:
    print(target)

pd.DataFrame(targets, columns=["Target"]).to_csv(
    "results/all_targets.csv",
    index=False
)