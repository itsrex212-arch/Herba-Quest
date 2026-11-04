"""
Score compounds using extracted features and trained models to produce ADMET evaluation files.
"""

from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")
HERE = Path(__file__).parent

W_QED = 0.4
W_HERG = 0.3
W_HEPATO = 0.3
PAINS_PENALTY = 0.10

MAX_LIPINSKI = 1
MAX_HERG_RISK = 0.70
MAX_HEPATO_RISK = 0.70
REQUIRE_VEBER = False

MAX_PASS = 300
MIN_PASS_TARGET = 100

UPPERCASE_HEADERS_IN_HANDOFF = False

OUT_COLS = ["compound_id", "qed", "lipinski_violations", "herg_risk",
            "hepatotox_risk", "admet_score", "pass_flag"]
MODEL_FILES = {"herg_risk": HERE / "models" / "herg_model.joblib",
               "hepatotox_risk": HERE / "models" / "hepatotox_model.joblib"}


def as_bool(series):
    """Convert CSV boolean text representations to actual boolean values."""
    return series.astype(str).str.strip().str.lower().eq("true")


def fingerprints(smiles, bits, radius):
    """Generate Morgan fingerprint arrays for a list of SMILES strings."""
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=bits)
    rows = []
    for s in smiles:
        mol = Chem.MolFromSmiles(str(s))
        arr = np.zeros((bits,), dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(gen.GetFingerprint(mol), arr)
        rows.append(arr)
    return np.array(rows)


rules_path = HERE / "rules_features.csv"
if not rules_path.exists():
    raise SystemExit("STOP: rules_features.csv not found.")
df = pd.read_csv(rules_path)
n_in = len(df)

ok = as_bool(df["parse_ok"])
pains = as_bool(df["pains_flag"])
veber = as_bool(df["veber_pass"])
print(f"Loaded {n_in} compounds; {int(ok.sum())} readable, {int((~ok).sum())} unreadable")

for col, path in MODEL_FILES.items():
    if not path.exists():
        raise SystemExit(f"STOP: {path.name} not found.")
    bundle = joblib.load(path)

    df[col] = np.nan
    X = fingerprints(df.loc[ok, "smiles_canonical"], bundle["fp_bits"], bundle["radius"])
    df.loc[ok, col] = bundle["model"].predict_proba(X)[:, 1]
    print(f"{col}: trained on {bundle.get('dataset', 'unknown dataset')}, model = {bundle['model_name']}, "
          f"scaffold-CV AUC = {bundle['scaffold_cv_auc']:.3f}")

df["admet_score"] = (W_QED * df["qed"]
                     + W_HERG * (1 - df["herg_risk"])
                     + W_HEPATO * (1 - df["hepatotox_risk"])
                     - PAINS_PENALTY * pains.astype(float)).clip(0, 1)

fails = pd.DataFrame({
    "unreadable_smiles": ~ok,
    f"lipinski_over_{MAX_LIPINSKI}": ok & (df["lipinski_violations"] > MAX_LIPINSKI),
    "herg_risk_high": ok & (df["herg_risk"] >= MAX_HERG_RISK),
    "hepatotox_risk_high": ok & (df["hepatotox_risk"] >= MAX_HEPATO_RISK),
    "veber_fail": ok & ~veber if REQUIRE_VEBER else pd.Series(False, index=df.index),
})
passes_rules = ~fails.any(axis=1)
fail_reason = fails.apply(lambda r: ";".join(r.index[r.values]), axis=1)

pass_flag = passes_rules.copy()
if MAX_PASS is not None and passes_rules.sum() > MAX_PASS:
    keep_idx = df[passes_rules].nlargest(MAX_PASS, "admet_score").index
    capped = passes_rules & ~df.index.isin(keep_idx)
    pass_flag = passes_rules & ~capped
    fail_reason.loc[capped] = "ranked_below_cap"
    print(f"\n{int(passes_rules.sum())} passed the rules; kept the top {MAX_PASS} by admet_score.")

df["pass_flag"] = pass_flag
df["fail_reason"] = fail_reason

out = df[OUT_COLS].copy()
out["lipinski_violations"] = out["lipinski_violations"].astype("Int64")
for c in ["qed", "herg_risk", "hepatotox_risk", "admet_score"]:
    out[c] = out[c].round(3)
upper = out.rename(columns=str.upper)
upper.to_csv(HERE / "compounds_scored_UPPERCASE.csv", index=False)
(upper if UPPERCASE_HEADERS_IN_HANDOFF else out).to_csv(HERE / "compounds_scored.csv", index=False)

detail_cols = (["compound_id", "name", "source_plant", "mw_calc", "logp", "hbd", "hba", "tpsa",
                "rot_bonds", "veber_pass", "pains_flag", "pains_name"] + OUT_COLS[1:] + ["fail_reason"])
df[[c for c in detail_cols if c in df.columns]].round(3).to_csv(
    HERE / "compounds_scored_detail.csv", index=False)

assert list(out.columns) == OUT_COLS, "Column names/order do not match"
assert len(out) == n_in, "Row count changed"
assert out["compound_id"].is_unique, "compound_id must be unique"
assert out["pass_flag"].dtype == bool, "pass_flag must be boolean"

n_pass = int(out["pass_flag"].sum())
print("\n=== Summary ===")
print(f"Compounds in            : {n_in}")
print(f"pass_flag = True        : {n_pass}")
print("Why compounds failed    :")
reasons = fail_reason[~pass_flag].str.split(";").explode().value_counts()
print(reasons.to_string() if len(reasons) else "  (none)")
print(f"herg_risk median        : {df['herg_risk'].median():.3f}")
print(f"hepatotox_risk median   : {df['hepatotox_risk'].median():.3f}")
print(f"admet_score min/med/max : {df['admet_score'].min():.3f} / "
      f"{df['admet_score'].median():.3f} / {df['admet_score'].max():.3f}")
if n_pass < MIN_PASS_TARGET:
    print(f"\nNOTE: only {n_pass} passed.")
print("\nSaved compounds_scored.csv (handoff), compounds_scored_UPPERCASE.csv and compounds_scored_detail.csv")