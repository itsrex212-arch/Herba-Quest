"""
Train and evaluate hERG and liver-toxicity classification models using TDC datasets.
"""

from pathlib import Path

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedGroupKFold
from tdc.single_pred import Tox

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).parent
DATA_DIR = HERE / "data"
MODEL_DIR = HERE / "models"
DATA_DIR.mkdir(exist_ok=True)
MODEL_DIR.mkdir(exist_ok=True)

SEED = 42
N_FOLDS = 5
FP_BITS = 2048
RADIUS = 2

HERG_DATASET = "herg" if "herg" else "hERG"
TASKS = {"herg": HERG_DATASET, "hepatotox": "DILI"}

fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=RADIUS, fpSize=FP_BITS)


def featurize(smiles_list):
    """Generate Morgan fingerprint matrix, kept indices, and Murcko scaffolds."""
    X, kept, scaffolds = [], [], []
    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            continue
        arr = np.zeros((FP_BITS,), dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(fp_gen.GetFingerprint(mol), arr)
        X.append(arr)
        kept.append(i)
        scaffolds.append(MurckoScaffold.MurckoScaffoldSmiles(mol=mol))
    return np.array(X), kept, scaffolds


def make_model(name):
    """Instantiate a classification model by name."""
    if name == "dummy_baseline":
        return DummyClassifier(strategy="prior")
    if name == "logistic_regression":
        return LogisticRegression(max_iter=2000, class_weight="balanced")
    if name == "random_forest":
        return RandomForestClassifier(n_estimators=300, class_weight="balanced",
                                     random_state=SEED, n_jobs=2)
    raise ValueError(name)


def scaffold_cv(model_name, X, y, groups):
    """Perform 5-fold stratified group cross-validation using scaffolds."""
    folds = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED).split(X, y, groups)
    aucs = []
    oof = np.full(len(y), np.nan)
    for train_idx, test_idx in folds:
        model = make_model(model_name)
        model.fit(X[train_idx], y[train_idx])
        p = model.predict_proba(X[test_idx])[:, 1]
        oof[test_idx] = p
        if len(np.unique(y[test_idx])) == 2:
            aucs.append(roc_auc_score(y[test_idx], p))
    return aucs, oof


results = []
for task, tdc_name in TASKS.items():
    print("\n" + "=" * 60)
    print(f"TASK: {task}   (TDC dataset: {tdc_name})")
    print("=" * 60)

    df = Tox(name=tdc_name, path=str(DATA_DIR)).get_data()
    df.to_csv(DATA_DIR / f"{task}_dataset.csv", index=False)
    X, kept, scaffolds = featurize(df["Drug"].tolist())
    y = df["Y"].astype(int).values[kept]
    groups = pd.factorize(pd.Series(scaffolds))[0]

    assert set(np.unique(y)) <= {0, 1}, "Labels are not 0/1"
    print(f"Molecules used   : {len(y)}  (skipped {len(df) - len(y)} unreadable)")
    print(f"Risky (label 1)  : {y.sum()}  ({100 * y.mean():.1f}%)")

    print(f"\n  SCAFFOLD {N_FOLDS}-fold CV, ROC-AUC (mean +/- std)")
    best = {"name": None, "auc": -1, "oof": None}
    for name in ["dummy_baseline", "logistic_regression", "random_forest"]:
        aucs, oof = scaffold_cv(name, X, y, groups)
        m, s = float(np.mean(aucs)), float(np.std(aucs))
        print(f"    {name:20s} {m:.3f} +/- {s:.3f}   (folds used: {len(aucs)})")
        results.append({"task": task, "dataset": tdc_name, "model": name, "auc_mean": round(m, 3),
                        "auc_std": round(s, 3), "n_molecules": len(y)})
        if name != "dummy_baseline" and m > best["auc"]:
            best = {"name": name, "auc": m, "oof": oof}

    ok = ~np.isnan(best["oof"])
    pooled_auc = roc_auc_score(y[ok], best["oof"][ok])
    fpr, tpr, _ = roc_curve(y[ok], best["oof"][ok])
    plt.figure(figsize=(4.5, 4.5))
    plt.plot(fpr, tpr, label=f"{best['name']} (pooled AUC {pooled_auc:.2f})")
    plt.plot([0, 1], [0, 1], "--", color="grey", label="coin flip")
    plt.xlabel("False positive rate")
    plt.ylabel("True positive rate")
    plt.title(f"{task}: ROC, scaffold cross-validation")
    plt.legend(loc="lower right", fontsize=8)
    plt.tight_layout()
    plt.savefig(HERE / f"roc_{task}.png", dpi=150)
    plt.close()

    final = make_model(best["name"])
    final.fit(X, y)
    joblib.dump({"model": final, "fp_bits": FP_BITS, "radius": RADIUS, "task": task, "dataset": tdc_name,
                 "model_name": best["name"], "scaffold_cv_auc": best["auc"]},
                MODEL_DIR / f"{task}_model.joblib")
    print(f"\n  -> chosen model for {task}: {best['name']}  (saved to models/{task}_model.joblib)")

    dummy = [r for r in results if r["task"] == task and r["model"] == "dummy_baseline"][0]
    assert abs(dummy["auc_mean"] - 0.5) < 0.05, "Dummy baseline is not ~0.5"

pd.DataFrame(results).to_csv(HERE / "cv_results.csv", index=False)
print("\nSaved cv_results.csv, roc_*.png and models/*.joblib")