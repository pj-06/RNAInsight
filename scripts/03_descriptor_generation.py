import pandas as pd
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors

# --------------------------------------------------
# Paths
# --------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent.parent

input_file = BASE_DIR / "data" / "processed" / "target_standardized.csv"
output_file = BASE_DIR / "data" / "processed" / "descriptor_dataset.csv"

# --------------------------------------------------
# Load Dataset
# --------------------------------------------------

df = pd.read_csv(input_file)

print("Dataset Loaded!")
print("Rows:", len(df))

# --------------------------------------------------
# Function to calculate molecular descriptors
# --------------------------------------------------

def calculate_descriptors(smiles):

    # Handle missing SMILES
    if pd.isna(smiles) or smiles == "":
        return pd.Series({
            "MolecularWeight": None,
            "ExactMolWt": None,
            "LogP": None,
            "MolMR": None,
            "TPSA": None,
            "HBA": None,
            "HBD": None,
            "RotatableBonds": None,
            "RingCount": None,
            "NumAromaticRings": None,
            "NumAliphaticRings": None,
            "HeavyAtomCount": None,
            "NumHeteroAtoms": None,
            "FractionCSP3": None
        })

    # Convert SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)

    # Handle invalid SMILES
    if mol is None:
        return pd.Series({
            "MolecularWeight": None,
            "ExactMolWt": None,
            "LogP": None,
            "MolMR": None,
            "TPSA": None,
            "HBA": None,
            "HBD": None,
            "RotatableBonds": None,
            "RingCount": None,
            "NumAromaticRings": None,
            "NumAliphaticRings": None,
            "HeavyAtomCount": None,
            "NumHeteroAtoms": None,
            "FractionCSP3": None
        })

    # Calculate descriptors
    return pd.Series({

        "MolecularWeight": round(Descriptors.MolWt(mol), 3),

        "ExactMolWt": round(Descriptors.ExactMolWt(mol), 3),

        "LogP": round(Descriptors.MolLogP(mol), 3),

        "MolMR": round(Descriptors.MolMR(mol), 3),

        "TPSA": round(Descriptors.TPSA(mol), 3),

        "HBA": Descriptors.NumHAcceptors(mol),

        "HBD": Descriptors.NumHDonors(mol),

        "RotatableBonds": Descriptors.NumRotatableBonds(mol),

        "RingCount": Descriptors.RingCount(mol),

        "NumAromaticRings": Descriptors.NumAromaticRings(mol),

        "NumAliphaticRings": Descriptors.NumAliphaticRings(mol),

        "HeavyAtomCount": Descriptors.HeavyAtomCount(mol),

        "NumHeteroAtoms": Descriptors.NumHeteroatoms(mol),

        "FractionCSP3": round(Descriptors.FractionCSP3(mol), 3)

    })


# --------------------------------------------------
# Generate Descriptors
# --------------------------------------------------

descriptor_df = df["SMILES"].apply(calculate_descriptors)

# Merge descriptors with original dataset
final_df = pd.concat([df, descriptor_df], axis=1)

# --------------------------------------------------
# Display Summary
# --------------------------------------------------

invalid_smiles = final_df["MolecularWeight"].isna().sum()

print("\nInvalid SMILES:", invalid_smiles)

print("\nTotal Columns:", len(final_df.columns))

print("\nPreview:")
print(final_df.head())

# --------------------------------------------------
# Save Dataset
# --------------------------------------------------

final_df.to_csv(output_file, index=False)

print("\nDescriptor dataset created successfully!")
print(f"Saved to: {output_file}")