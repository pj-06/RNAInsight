import pandas as pd
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent

file = BASE_DIR / "data" / "raw" / "NIHMS1824842-supplement-RBIND_Info_A.xlsx"

df = pd.read_excel(file, header=1)

df.columns = df.columns.str.strip()
df = df.dropna(axis=1, how="all")

df = df[[
    "#",
    "Name",
    "SMILES",
    "Assay Target(s)"
]]

# Create new column
df["Target_Family"] = "Other"

# ------------------------
# Existing mapping rules
# ------------------------

df.loc[
    df["Assay Target(s)"].str.contains("HIV", case=False, na=False),
    "Target_Family"
] = "HIV TAR"

df.loc[
    df["Assay Target(s)"].str.contains("HCV", case=False, na=False),
    "Target_Family"
] = "HCV IRES"

df.loc[
    df["Assay Target(s)"].str.contains("Riboswitch", case=False, na=False),
    "Target_Family"
] = "Riboswitch"

df.loc[
    df["Assay Target(s)"].str.contains("miRNA|miR", case=False, na=False),
    "Target_Family"
] = "miRNA"

df.loc[
    df["Assay Target(s)"].str.contains("G-quadruplex|G4", case=False, na=False),
    "Target_Family"
] = "G-Quadruplex"

# ------------------------
# ADD THE NEW RULES HERE
# ------------------------

# Repeat Expansion RNAs
df.loc[
    df["Assay Target(s)"].str.contains(
        r"CUG|CGG|CCUG|CAG", case=False, na=False
    ),
    "Target_Family"
] = "Repeat RNA"

# mRNA
df.loc[
    df["Assay Target(s)"].str.contains(
        r"mRNA|pre-mRNA|SMN2|MAPT|ADAM10|RpoH|Bcl",
        case=False,
        na=False
    ),
    "Target_Family"
] = "mRNA"

# Viral RNA
df.loc[
    df["Assay Target(s)"].str.contains(
        r"Influenza|SARS",
        case=False,
        na=False
    ),
    "Target_Family"
] = "Viral RNA"

# Ribozymes
df.loc[
    df["Assay Target(s)"].str.contains(
        r"Intron|ribozyme",
        case=False,
        na=False
    ),
    "Target_Family"
] = "Ribozyme"

# tRNA
df.loc[
    df["Assay Target(s)"].str.contains(
        r"tRNA",
        case=False,
        na=False
    ),
    "Target_Family"
] = "tRNA"

print(df["Target_Family"].value_counts())

output = BASE_DIR / "data" / "processed" / "target_standardized.csv"

df.to_csv(output, index=False)

print("Target standardized dataset saved!")
