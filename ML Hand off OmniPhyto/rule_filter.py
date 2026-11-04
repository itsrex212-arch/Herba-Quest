"""
Rule-based drug-likeness feature extraction.
Input: Compounds CSV file
Output: rules_features.csv (saved in script directory)
"""

import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors, QED
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

# Define input and output file paths
HERE = Path(__file__).parent
SAMPLE_FILE = "neem_leave compound.csv"
in_path = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / SAMPLE_FILE

if not in_path.exists():
    sys.exit(f"STOP: cannot find {in_path}\n"
             f"Put the CSV in {HERE} or pass its path when running.")

# Load dataset and validate columns and uniqueness
df = pd.read_csv(in_path)

required = ["compound_id", "coconut_id", "name", "source_plant",
            "smiles_canonical", "inchikey", "mol_weight"]

if "name" not in df.columns and "name_x" in df.columns:
    df = df.rename(columns={"name_x": "name"})

missing = [c for c in required if c not in df.columns]
if missing:
    sys.exit(f"STOP: missing columns {missing}.")

df = df[required].copy()

if df["compound_id"].duplicated().any():
    sys.exit("STOP: compound_id is not unique.")

print(f"Loaded {len(df)} compounds from {in_path.name}")

# Initialize PAINS filter catalog
params = FilterCatalogParams()
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
pains_catalog = FilterCatalog(params)

# Compute molecular descriptors and rule features per molecule
rows = []
for _, r in df.iterrows():
    mol = Chem.MolFromSmiles(str(r["smiles_canonical"]))
    if mol is None:
        rows.append({"compound_id": r["compound_id"], "parse_ok": False})
        continue

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = rdMolDescriptors.CalcNumLipinskiHBD(mol)
    hba = rdMolDescriptors.CalcNumLipinskiHBA(mol)
    hbd_rdkit = Lipinski.NumHDonors(mol)
    hba_rdkit = Lipinski.NumHAcceptors(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)

    lipinski_violations = int(mw > 500) + int(logp > 5) + int(hbd > 5) + int(hba > 10)
    veber_pass = (rotb <= 10) and (tpsa <= 140)

    hit = pains_catalog.GetFirstMatch(mol)
    pains_flag = hit is not None
    pains_name = hit.GetDescription() if hit is not None else ""

    rows.append({
        "compound_id": r["compound_id"],
        "parse_ok": True,
        "mw_calc": round(mw, 2),
        "logp": round(logp, 2),
        "hbd": hbd,
        "hba": hba,
        "hbd_rdkit": hbd_rdkit,
        "hba_rdkit": hba_rdkit,
        "tpsa": round(tpsa, 1),
        "rot_bonds": rotb,
        "lipinski_violations": lipinski_violations,
        "veber_pass": veber_pass,
        "qed": round(QED.qed(mol), 3),
        "pains_flag": pains_flag,
        "pains_name": pains_name,
    })

feat = pd.DataFrame(rows)

# Merge feature data with original dataframe and export to CSV
out = df.merge(feat, on="compound_id", how="left")
out_path = HERE / "rules_features.csv"
out.to_csv(out_path, index=False)

# Print execution summary statistics
ok = out[out["parse_ok"] == True]
print("\n=== Sanity report ===")
print(f"Could not parse SMILES : {(out['parse_ok'] != True).sum()}")
print("Lipinski violations    :")
print(ok["lipinski_violations"].value_counts().sort_index().to_string())
print(f"Lipinski HBA/HBD differ from RDKit default for {(ok['hba'] != ok['hba_rdkit']).sum() + (ok['hbd'] != ok['hbd_rdkit']).sum()} counts")
print(f"Veber pass             : {ok['veber_pass'].sum()} / {len(ok)}")
print(f"PAINS flagged          : {ok['pains_flag'].sum()} / {len(ok)}")
print(f"QED  min/median/max    : {ok['qed'].min()} / {ok['qed'].median()} / {ok['qed'].max()}")

diff = (ok["mw_calc"] - pd.to_numeric(ok["mol_weight"], errors="coerce")).abs()
print(f"MW mismatch > 1 Da     : {(diff > 1).sum()} compounds")
print(f"\nSaved {out_path}")