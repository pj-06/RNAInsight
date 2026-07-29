import pandas as pd
from pathlib import Path

# Paths
BASE_DIR = Path(__file__).resolve().parent.parent

input_file = BASE_DIR / "data" / "processed" / "descriptor_dataset.csv"
output_file = BASE_DIR / "data" / "processed" / "cleaned_dataset.csv"

# Load dataset
df = pd.read_csv(input_file)

print("Original rows:", len(df))

# Count invalid SMILES
invalid = df["MolecularWeight"].isna().sum()
print("Invalid SMILES:", invalid)

# Remove rows with missing descriptor values
df = df[df["MolecularWeight"].notna()]
before = len(df)

df = df.drop_duplicates()

after = len(df)

print(f"Duplicate rows removed: {before - after}")

print("Remaining rows:", len(df))

# Save cleaned dataset
df.to_csv(output_file, index=False)

print("\nCleaned dataset saved successfully!")
print(f"Saved to: {output_file}")