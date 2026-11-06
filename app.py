"""
OmniPhyto results app  (self-contained: this one file + optional assets/ folder).

Run from the repo root:   streamlit run app.py

HOW THIS FILE IS BUILT (3 parts)
  PART 1  imports, constants, registries, loader, helpers, sections 0 and 10, router  (this file)
  PART 2  sections 1-5  -> paste at the line marked  "# ===== PART 2 GOES HERE ====="
  PART 3  sections 6-9  -> paste at the line marked  "# ===== PART 3 GOES HERE ====="
The router at the bottom already knows all 11 sections, so the app runs after every part.
Sections that are not pasted yet show a small "not added yet" note instead of crashing.

RULES THIS APP FOLLOWS
  * It only DISPLAYS files. It never invents data and never computes a composite score.
  * A lower (more negative) Vina score is better. Scores are predictions, in kcal/mol.
  * Missing data shows "-" or a WAITING card, never 0.
  * A broken file shows a short message and the rest of the app carries on.
"""
from __future__ import annotations

import fnmatch
import html
import io
import json
import os
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

try:  # RDKit is optional: without it we show the SMILES text instead of a drawing.
    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D

    HAVE_RDKIT = True
except Exception:  # pragma: no cover
    HAVE_RDKIT = False

st.set_page_config(page_title="OmniPhyto", page_icon="🍃", layout="wide")

# =============================================================================
# 1. CONSTANTS
# =============================================================================
ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
INCOMING = DATA / "incoming"

SEED = 42                       # fixed random seed (used by any chart jitter / sampling)
NO_VALID_POSE_ABOVE = -2.0      # a Vina score above this is not a real pose (kcal/mol)
NOISE_KCAL = 0.5                # scores closer than this are inside the search noise
EXHAUSTIVENESS_SCREEN = 4       # lowest search effort, from the guide / files
N_STEREO_NOT_DOCKED = 62        # from target_card.json (used only if the card lacks it)

ROLE_COLORS = {
    "Chemist": "#2E7D4F",
    "Bioinformatics": "#1F5FA8",
    "ML Engineer": "#D9822B",
    "Developer": "#7B4FA3",
    "Project Manager": "#5F6B73",
    "Everyone": "#5F6B73",
}
ROLE_DOT = {  # emoji dots for the sidebar radio list (plain text only there)
    "Chemist": "🟢", "Bioinformatics": "🔵", "ML Engineer": "🟠",
    "Developer": "🟣", "Project Manager": "⚪", "Everyone": "⚪",
}

# The four proteins. These are the verified project facts.
PROTEINS = {
    "dpp4":   {"name": "DPP-4",  "pdb": "3D4L", "chains": "A",       "ligand": "605", "role": "disease target"},
    "herg":   {"name": "hERG",   "pdb": "8ZYO", "chains": "A, B, C, D", "ligand": "XB7", "role": "safety protein"},
    "cyp3a4": {"name": "CYP3A4", "pdb": "2V0M", "chains": "A",       "ligand": "KLN", "role": "safety protein"},
    "cyp2d6": {"name": "CYP2D6", "pdb": "3TDA", "chains": "A",       "ligand": "PN0", "role": "safety protein"},
}
ANTITARGET_KEY = {"hERG": "herg", "CYP3A4": "cyp3a4", "CYP2D6": "cyp2d6"}

# Used ONLY for the target card on Home / Section 1 when target_card.json lacks a field.
TARGET_FALLBACK = {
    "disease": "Type 2 diabetes mellitus (T2DM)",
    "target_name": "Human dipeptidyl peptidase IV / CD26 (DPP-4)",
    "pdb_id": "3D4L",
    "resolution": 2.0,
    "method": "X-ray diffraction",
    "chain": "A",
    "ref_code": "605",
    "organism": "Homo sapiens",
}

STORY = (
    "OmniPhyto is a hackathon project from a 5-person team. We take chemicals found in neem leaves "
    "(Azadirachta indica), filter them for drug-likeness and safety risk, and dock them on the computer "
    "to ONE disease protein (DPP-4, for type 2 diabetes) and to three safety proteins (hERG, CYP3A4, CYP2D6). "
    "Then we rank them. Everything here is a computer prediction that needs lab confirmation."
)

# =============================================================================
# 2. REGISTRIES  (STEPS = who does what; FILES = which file feeds which section)
# =============================================================================
STEPS = [
    {"id": 1, "name": "Target & library",        "owner": "Chemist",        "day": "1-2", "needs": ["target_card", "master"]},
    {"id": 2, "name": "Structures & validation", "owner": "Bioinformatics", "day": "1-2", "needs": ["receptors", "redocking"]},
    {"id": 3, "name": "ADMET filter",            "owner": "ML Engineer",    "day": "3",   "needs": ["ml_detail"]},
    {"id": 4, "name": "Docking: DPP-4",          "owner": "Bioinformatics", "day": "3-4", "needs": ["docking_target", "poses_target"]},
    {"id": 5, "name": "Safety panel",            "owner": "Bioinformatics", "day": "5",   "needs": ["docking_antitargets"]},
    {"id": 6, "name": "Final ranking",           "owner": "ML Engineer",    "day": "5",   "needs": ["ranking"]},
    {"id": 7, "name": "Simulation (MD)",         "owner": "Bioinformatics", "day": "5-6", "needs": ["md_summary"]},
]
STEP_BY_ID = {s["id"]: s for s in STEPS}

SCORE_ALIASES = ["vina_affinity_kcal_mol", "affinity_kcal_mol", "docking_score", "vina_score", "score_kcal"]
ID_ALIASES = ["compound", "compoundid", "compound id", "cid"]

# kind: csv | json | folder.  patterns = (folder relative to repo, file-name pattern); newest file wins.
# feeds = the section numbers/names this file fills.  canonical = name used when saving an upload.
FILES: dict[str, dict] = {
    "target_card": dict(
        label="target_card.json", kind="json", owner="Chemist", day="1-2",
        patterns=[("data/incoming", "target_card*.json")], canonical="target_card.json",
        dest="data/incoming", feeds=["0 Home", "1 Target & library"]),
    "master": dict(
        label="neem_compounds_with_DOIs.csv (master table)", kind="csv", owner="Chemist", day="1-2",
        patterns=[("data/incoming", "neem_compounds_with_DOIs*.csv")], canonical="neem_compounds_with_DOIs.csv",
        dest="data/incoming", required=["compound_id", "smiles_canonical"], unique=["compound_id"],
        aliases={"compound_id": ID_ALIASES, "smiles_canonical": ["canonical_smiles", "smiles"]},
        feeds=["1 Target & library", "8 Compound profile"]),
    "ml_detail": dict(
        label="*scored_detail*.csv (ML table)", kind="csv", owner="ML Engineer", day="3",
        patterns=[("data/incoming", "*scored_detail*.csv")], canonical="compounds_scored_detail.csv",
        dest="data/incoming", required=["compound_id", "pass_flag"], unique=["compound_id"],
        numeric=["mw_calc", "logp", "hbd", "hba", "tpsa", "rot_bonds", "qed", "lipinski_violations", "admet_score"],
        aliases={"compound_id": ID_ALIASES},
        feeds=["0 Home", "3 ADMET filter", "4 Docking: DPP-4", "8 Compound profile"]),
    "ml_compact": dict(
        label="*scored.csv (ML compact table, fallback)", kind="csv", owner="ML Engineer", day="3",
        patterns=[("data/incoming", "*scored.csv")], canonical="compounds_scored.csv",
        dest="data/incoming", required=["compound_id", "pass_flag"], unique=["compound_id"],
        numeric=["qed", "lipinski_violations", "admet_score"], aliases={"compound_id": ID_ALIASES},
        feeds=["3 ADMET filter (only if the detail table is missing)"]),
    "ligand_prep_log": dict(
        label="ligand_prep_log.csv", kind="csv", owner="Developer", day="3",
        patterns=[("data/incoming", "ligand_prep_log*.csv")], canonical="ligand_prep_log.csv",
        dest="data/incoming", required=[], feeds=["10 Data inbox"]),
    "docking_target": dict(
        label="docking_target.csv", kind="csv", owner="Bioinformatics", day="3-4",
        patterns=[("data/incoming", "docking_target*.csv")], canonical="docking_target.csv",
        dest="data/incoming", required=["compound_id", "vina_score_kcal"], unique=["compound_id"],
        numeric=["vina_score_kcal", "heavy_atoms", "ligand_efficiency"],
        aliases={"compound_id": ID_ALIASES, "vina_score_kcal": SCORE_ALIASES, "pose_file": ["pose", "pose_path"]},
        feeds=["0 Home", "4 Docking: DPP-4", "8 Compound profile"]),
    "docking_antitargets": dict(
        label="docking_antitargets.csv", kind="csv", owner="Bioinformatics", day="5",
        patterns=[("data/incoming", "docking_antitargets*.csv")], canonical="docking_antitargets.csv",
        dest="data/incoming", required=["compound_id", "antitarget", "vina_score_kcal"],
        unique=["compound_id", "antitarget"],
        numeric=["vina_score_kcal", "target_score_kcal", "delta_vs_target_kcal"],
        aliases={"compound_id": ID_ALIASES, "vina_score_kcal": SCORE_ALIASES,
                 "antitarget": ["protein", "receptor", "target_protein"]},
        feeds=["5 Safety panel", "6 Final ranking", "8 Compound profile"]),
    "summary_report": dict(
        label="summary_report.csv (top-10 report)", kind="csv", owner="Bioinformatics", day="5",
        patterns=[("data/incoming", "summary_report*.csv")], canonical="summary_report.csv",
        dest="data/incoming", required=["compound_id", "vina_score_kcal"], unique=["compound_id"],
        numeric=["vina_score_kcal", "heavy_atoms", "ligand_efficiency", "CYP2D6_score_kcal", "CYP3A4_score_kcal",
                 "hERG_score_kcal", "worst_delta_vs_target_kcal", "n_hbonds", "n_hydrophobic",
                 "reference_score_kcal", "exhaustiveness"],
        aliases={"compound_id": ID_ALIASES, "vina_score_kcal": SCORE_ALIASES},
        feeds=["0 Home", "2 Structures & validation", "4 Docking: DPP-4", "5 Safety panel", "8 Compound profile"]),
    "ranking": dict(
        label="ranking.csv (final ranking)", kind="csv", owner="ML Engineer", day="5",
        patterns=[("data/incoming", "ranking*.csv")], canonical="ranking.csv",
        dest="data/incoming", required=["compound_id", "composite_score"], unique=["compound_id"],
        numeric=["rank", "efficacy_norm", "selectivity_norm", "admet_score", "composite_score"],
        aliases={"compound_id": ID_ALIASES, "rank": ["final_rank", "rank_final"],
                 "composite_score": ["composite", "final_score"], "efficacy_norm": ["efficacy"],
                 "selectivity_norm": ["selectivity"], "admet_score": ["admet"]},
        feeds=["6 Final ranking", "8 Compound profile", "9 Evidence & methods"]),
    "redocking": dict(
        label="redocking_results.csv / redock_summary.csv", kind="csv", owner="Bioinformatics", day="1-2",
        patterns=[("data/incoming", "redocking_results*.csv"), ("data/incoming", "redock_summary*.csv"),
                  ("data/validation", "redock_summary*.csv")],
        canonical="redocking_results.csv", dest="data/incoming",
        required=["receptor", "verdict"], unique=["receptor"],
        numeric=["top1_rmsd", "best_rmsd", "ref_score"],
        aliases={"receptor": ["protein"], "verdict": ["status"], "top1_rmsd": ["rmsd_top_pose"],
                 "best_rmsd": ["rmsd_best_of_9"], "ref_score": ["reference_score_kcal"]},
        feeds=["0 Home", "2 Structures & validation", "5 Safety panel"]),
    "redock_rmsd": dict(
        label="redock_rmsd*.csv", kind="csv", owner="Bioinformatics", day="1-2",
        patterns=[("data/validation", "redock_rmsd*.csv")], canonical="redock_rmsd.csv",
        dest="data/validation", required=[], feeds=["2 Structures & validation"]),
    "md_summary": dict(
        label="md_report/md_summary.csv", kind="csv", owner="Bioinformatics", day="5-6",
        patterns=[("data/incoming/md_report", "md_summary*.csv")], canonical="md_summary.csv",
        dest="data/incoming/md_report", required=["compound_id", "label"],
        numeric=["length_ns", "ligand_rmsd_mean_A", "ligand_rmsd_sd_A", "protein_rmsd_mean_A", "hbond_occupancy_pct"],
        aliases={"compound_id": ID_ALIASES}, feeds=["7 Simulation (MD)"]),
    # ---- folders of structure files ----
    "receptors": dict(
        label="data/receptors/*.pdb (raw RCSB structures)", kind="folder", owner="Bioinformatics", day="1-2",
        patterns=[("data/receptors", "*.pdb")], expected=4, dest="data/receptors",
        feeds=["2 Structures & validation", "4 Docking: DPP-4", "8 Compound profile"]),
    "receptors_prepared": dict(
        label="data/prepared/receptors/* (prepared receptors, *_vina.conf)", kind="folder", owner="Bioinformatics", day="1-2",
        patterns=[("data/prepared/receptors", "*")], expected=None, dest="data/prepared/receptors",
        feeds=["2 Structures & validation"]),
    "ligands_prepared": dict(
        label="data/prepared/ligands/*.pdbqt", kind="folder", owner="Developer", day="3",
        patterns=[("data/prepared/ligands", "*.pdbqt")], expected=None, dest="data/prepared/ligands",
        feeds=["10 Data inbox"]),
    "poses_target": dict(
        label="results/target/poses/*.pdbqt (DPP-4 poses)", kind="folder", owner="Bioinformatics", day="3-4",
        patterns=[("data/incoming/results/target/poses", "*.pdbqt")], expected=238,
        dest="data/incoming/results/target/poses", feeds=["4 Docking: DPP-4", "8 Compound profile"]),
    "poses_herg": dict(
        label="results/herg/poses/*.pdbqt", kind="folder", owner="Bioinformatics", day="5",
        patterns=[("data/incoming/results/herg/poses", "*.pdbqt")], expected=10,
        dest="data/incoming/results/herg/poses", feeds=["5 Safety panel"]),
    "poses_cyp3a4": dict(
        label="results/cyp3a4/poses/*.pdbqt", kind="folder", owner="Bioinformatics", day="5",
        patterns=[("data/incoming/results/cyp3a4/poses", "*.pdbqt")], expected=10,
        dest="data/incoming/results/cyp3a4/poses", feeds=["5 Safety panel"]),
    "poses_cyp2d6": dict(
        label="results/cyp2d6/poses/*.pdbqt", kind="folder", owner="Bioinformatics", day="5",
        patterns=[("data/incoming/results/cyp2d6/poses", "*.pdbqt")], expected=10,
        dest="data/incoming/results/cyp2d6/poses", feeds=["5 Safety panel"]),
    "overlays": dict(
        label="validation/*_overlay_complex.pdb", kind="folder", owner="Bioinformatics", day="1-2",
        patterns=[("data/incoming/validation", "*_overlay_complex.pdb")], expected=4,
        dest="data/incoming/validation", feeds=["2 Structures & validation"]),
    "plip": dict(
        label="data/plip/OPxxxx_complex.pdb", kind="folder", owner="Bioinformatics", day="5",
        patterns=[("data/plip", "OP*_complex.pdb")], expected=10, dest="data/plip",
        feeds=["4 Docking: DPP-4", "8 Compound profile"]),
    "validation_files": dict(
        label="data/validation/* (redocking details)", kind="folder", owner="Bioinformatics", day="1-2",
        patterns=[("data/validation", "*")], expected=None, dest="data/validation",
        feeds=["2 Structures & validation"]),
    "md_plots": dict(
        label="md_report/*.png (MD plots)", kind="folder", owner="Bioinformatics", day="5-6",
        patterns=[("data/incoming/md_report", "*.png")], expected=None, dest="data/incoming/md_report",
        feeds=["7 Simulation (MD)"]),
}
# Files we know about but that no section needs (so they do not appear as "not used").
EXTRA_RECOGNISED = [("data/incoming/md_report", "README.md")]


# =============================================================================
# 3. SMALL HELPERS (text, numbers, paths)
# =============================================================================
def esc(x) -> str:
    """HTML-escape any text that came from a file."""
    return html.escape("" if x is None else str(x), quote=True)


def is_missing(x) -> bool:
    try:
        return x is None or (isinstance(x, float) and np.isnan(x)) or pd.isna(x)
    except (TypeError, ValueError):
        return False


def fmt(x, nd: int = 2, unit: str = "") -> str:
    """Number -> text. Missing is '-', never 0."""
    if is_missing(x):
        return "-"
    try:
        s = f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return esc(x)
    return f"{s} {unit}".strip()


def is_valid_pose(score) -> bool:
    """A score above NO_VALID_POSE_ABOVE (e.g. +4.678) is not a real pose."""
    return (not is_missing(score)) and float(score) <= NO_VALID_POSE_ABOVE


def safe_path(rel, base: Path = INCOMING):
    """Resolve rel inside base only. Returns a Path to an existing file, or None."""
    try:
        p = (base / str(rel).replace("\\", "/")).resolve()
        p.relative_to(base.resolve())
    except Exception:
        return None
    return p if p.is_file() else None


def read_text_safe(rel_from_root: str, max_mb: float = 8.0):
    """Read a repo file (path given relative to the repo root). None if missing / unsafe / too big."""
    p = safe_path(rel_from_root, ROOT)
    if p is None or p.stat().st_size > max_mb * 1_000_000:
        return None
    try:
        return p.read_text(errors="replace")
    except OSError:
        return None


def split_list(s, sep=";") -> list[str]:
    if is_missing(s):
        return []
    return [t.strip() for t in str(s).split(sep) if t.strip()]


def doi_list(s) -> list[str]:
    """Pipe-separated DOIs -> unique list (case-insensitive). The first spelling wins."""
    seen, out = set(), []
    for d in split_list(s, "|"):
        k = d.lower()
        if k not in seen:
            seen.add(k)
            out.append(d)
    return out


def doi_links_html(s) -> str:
    from urllib.parse import quote
    links = [f'<a href="https://doi.org/{quote(d, safe="/:()-._;")}" target="_blank" rel="noopener">{esc(d)}</a>'
             for d in doi_list(s)]
    return " · ".join(links) if links else "-"


def split_pdbqt_models(text: str):
    """Multi-model PDBQT -> list of (model_text, vina_score or None). MODEL 1 is the best pose."""
    models, cur, score = [], [], None
    for line in text.splitlines():
        if line.startswith("MODEL"):
            cur, score = [], None
        elif line.startswith("ENDMDL"):
            models.append(("\n".join(cur) + "\n", score))
            cur = []
        else:
            if line.startswith("REMARK VINA RESULT") and score is None:
                try:
                    score = float(line.split(":")[1].split()[0])
                except (IndexError, ValueError):
                    pass
            cur.append(line)
    if cur and not models:  # single-model file without MODEL/ENDMDL
        models.append(("\n".join(cur) + "\n", None))
    return models


# =============================================================================
# 4. THE LOADER  (find files by pattern, check them, cache by modification time)
# =============================================================================
def find_files(patterns) -> list[Path]:
    """Case-insensitive pattern search. Newest file first (newest wins)."""
    found: dict[str, Path] = {}
    for folder, pat in patterns:
        base = ROOT / folder
        if not base.is_dir():
            continue
        try:
            with os.scandir(base) as it:
                for entry in it:
                    if entry.is_file() and not entry.name.startswith(".") \
                            and fnmatch.fnmatch(entry.name.lower(), pat.lower()):
                        found[str(Path(entry.path))] = Path(entry.path)
        except OSError:
            continue
    return sorted(found.values(), key=lambda p: (p.stat().st_mtime, p.name), reverse=True)


def rel(p: Path) -> str:
    try:
        return str(p.relative_to(ROOT)).replace("\\", "/")
    except ValueError:
        return str(p)


def canon_columns(df: pd.DataFrame, aliases: dict | None, required: list[str], numeric: list[str]) -> pd.DataFrame:
    """Rename columns to the contract names: exact (any case), then aliases."""
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    wanted = list(dict.fromkeys(list(required) + list(numeric) + list((aliases or {}).keys())))
    lower = {c.lower(): c for c in df.columns}
    ren = {}
    for canon in wanted:
        if canon in df.columns:
            continue
        if canon.lower() in lower:
            ren[lower[canon.lower()]] = canon
            continue
        for al in (aliases or {}).get(canon, []):
            if al.lower() in lower and lower[al.lower()] not in ren:
                ren[lower[al.lower()]] = canon
                break
    return df.rename(columns=ren)


def _read_csv_any(src) -> pd.DataFrame:
    def attempt(**kw):
        if hasattr(src, "seek"):
            src.seek(0)
        return pd.read_csv(src, **kw)
    try:
        return attempt()
    except UnicodeDecodeError:
        return attempt(encoding="latin-1")
    except pd.errors.ParserError:
        return attempt(sep=None, engine="python")


def table_from(key: str, src, label: str):
    """Read + check one table against its contract.
    Returns (df or None, problems list, state: ok | warning | invalid)."""
    spec = FILES[key]
    problems: list[str] = []
    try:
        df = _read_csv_any(src)
    except Exception as exc:  # unreadable file must never crash the app
        return None, [f"{label}: could not be read ({type(exc).__name__}: {str(exc)[:80]})"], "invalid"
    required, numeric = spec.get("required", []), spec.get("numeric", [])
    df = canon_columns(df, spec.get("aliases"), required, numeric)
    if df.empty:
        return None, [f"{label}: the file has no rows"], "invalid"
    missing = [c for c in required if c not in df.columns]
    if missing:
        return None, [f"{label}: missing column(s) {', '.join(missing)}"], "invalid"
    if "compound_id" in df.columns:
        df["compound_id"] = df["compound_id"].astype("string").str.strip().str.upper()
        bad = df["compound_id"].isna() | (df["compound_id"] == "")
        if bad.any():
            problems.append(f"{int(bad.sum())} row(s) without compound_id were ignored")
            df = df[~bad]
    for c in numeric:
        if c in df.columns:
            before = df[c].notna().sum()
            df[c] = pd.to_numeric(df[c], errors="coerce")
            lost = int(before - df[c].notna().sum())
            if lost:
                problems.append(f"{lost} non-numeric value(s) in column {c}")
    if "pass_flag" in df.columns:
        df["pass_flag"] = df["pass_flag"].astype("string").str.strip().str.upper().map(
            {"TRUE": True, "FALSE": False, "1": True, "0": False, "YES": True, "NO": False}).astype("boolean")
        if df["pass_flag"].isna().any():
            problems.append(f"{int(df['pass_flag'].isna().sum())} row(s) with unreadable pass_flag")
    uniq = [c for c in spec.get("unique", []) if c in df.columns]
    if uniq and df.duplicated(uniq).any():
        problems.append(f"{int(df.duplicated(uniq).sum())} duplicate row(s) for {'+'.join(uniq)}; first one kept")
        df = df.drop_duplicates(uniq)
    df = _post_process(key, df, problems)
    problems = [f"{label}: {p}" for p in problems]
    return df.reset_index(drop=True), problems, ("warning" if problems else "ok")


def _post_process(key: str, df: pd.DataFrame, problems: list[str]) -> pd.DataFrame:
    """Small, safe tidy-ups. No science happens here."""
    if key == "docking_target":
        df = df.sort_values("vina_score_kcal", ascending=True, kind="stable").reset_index(drop=True)
        df["rank"] = np.arange(1, len(df) + 1)                       # shown by the app (1 = most negative)
        df["valid_pose"] = df["vina_score_kcal"].map(is_valid_pose)
        if df["vina_score_kcal"].isna().any():
            problems.append(f"{int(df['vina_score_kcal'].isna().sum())} compound(s) without a score")
    elif key == "docking_antitargets":
        def norm(a):
            s = re.sub(r"[^a-z0-9]", "", str(a).lower())
            return {"herg": "hERG", "cyp3a4": "CYP3A4", "cyp2d6": "CYP2D6"}.get(s, str(a).strip())
        df["antitarget"] = df["antitarget"].map(norm)
        df["valid_pose"] = df["vina_score_kcal"].map(is_valid_pose)
    elif key == "redocking":
        df["receptor"] = df["receptor"].astype(str).map(lambda s: re.sub(r"[^a-z0-9]", "", s.lower()))
        df["verdict"] = df["verdict"].astype(str).str.strip().str.upper()
    elif key == "ranking":
        if "rank" not in df.columns:
            df = df.sort_values("composite_score", ascending=False, kind="stable").reset_index(drop=True)
            df["rank"] = np.arange(1, len(df) + 1)
            problems.append("no rank column; rank taken from the order of composite_score (shown by the app)")
        else:
            df = df.sort_values("rank", kind="stable").reset_index(drop=True)
    elif key == "summary_report" and "rank" in df.columns:
        df = df.sort_values("rank", kind="stable").reset_index(drop=True)
    elif key == "ml_detail":
        for c in ("name", "source_plant", "fail_reason"):
            if c in df.columns:
                df[c] = df[c].astype("string").str.strip().replace("", pd.NA)
    return df


def _signature():
    """A fingerprint of every file in data/ (path, modified time, size). Changes -> cache refreshes."""
    items = []
    if DATA.is_dir():
        for p in DATA.rglob("*"):
            try:
                if p.is_file():
                    s = p.stat()
                    items.append((rel(p), s.st_mtime_ns, s.st_size))
            except OSError:
                continue
    return tuple(sorted(items))


def _load_all_uncached():
    data: dict = {}
    status: dict = {}
    recognised: set[str] = set()
    for key, spec in FILES.items():
        entry = {"key": key, "label": spec["label"], "feeds": spec["feeds"], "kind": spec["kind"],
                 "owner": spec.get("owner", ""), "day": spec.get("day", ""), "state": "missing",
                 "file": "", "rows": None, "expected": spec.get("expected"), "problems": []}
        try:
            paths = find_files(spec["patterns"])
        except Exception as exc:
            paths = []
            entry["problems"].append(f"search failed: {exc}")
        recognised.update(rel(p) for p in paths)
        if paths:
            newest = paths[0]
            entry["file"] = newest.name if spec["kind"] != "folder" else ""
            if spec["kind"] == "folder":
                data[key] = [rel(p) for p in paths]
                entry["rows"] = len(paths)
            elif spec["kind"] == "json":
                try:
                    obj = json.loads(newest.read_text(encoding="utf-8", errors="replace"))
                    if isinstance(obj, dict) and obj:
                        data[key], entry["state"] = obj, "ok"
                    else:
                        entry["state"], entry["problems"] = "invalid", [f"{newest.name}: not a JSON object"]
                except Exception as exc:
                    entry["state"], entry["problems"] = "invalid", [f"{newest.name}: bad JSON ({str(exc)[:60]})"]
            else:
                try:
                    df, probs, state = table_from(key, str(newest), newest.name)
                except Exception as exc:  # last safety net
                    df, probs, state = None, [f"{newest.name}: unexpected problem ({str(exc)[:80]})"], "invalid"
                entry["state"], entry["problems"] = state, probs
                if df is not None:
                    data[key] = df
                    entry["rows"] = len(df)
            if len(paths) > 1 and spec["kind"] != "folder":
                entry["problems"].append(f"{len(paths)} files match; the newest one ({newest.name}) is used")
        status[key] = entry

    # folders: ok when the expected number of files is there
    dt = data.get("docking_target")
    at = data.get("docking_antitargets")
    for key, spec in FILES.items():
        if spec["kind"] != "folder":
            continue
        n, entry = entry_n(status[key]), status[key]
        expected = spec.get("expected")
        if key == "poses_target" and dt is not None:
            expected = len(dt)
        if key.startswith("poses_") and key != "poses_target" and at is not None:
            prot = {v: k for k, v in ANTITARGET_KEY.items()}.get(key.split("_", 1)[1])
            if prot:
                expected = int((at["antitarget"] == prot).sum()) or expected
        entry["expected"] = expected
        if n == 0:
            entry["state"] = "missing"
        elif expected and n < expected:
            entry["state"] = "warning"
            entry["problems"].append(f"{n} of {expected} files found")
        else:
            entry["state"] = "ok"

    # names: ML table first, then the master table's own name column
    names: dict[str, str] = {}
    for src_key, col in (("ml_detail", "name"), ("ml_compact", "name"), ("master", "name_x"), ("master", "name")):
        df = data.get(src_key)
        if df is not None and col in df.columns:
            for cid, nm in zip(df["compound_id"], df[col]):
                if not is_missing(nm) and str(nm).strip() and cid not in names:
                    names[cid] = str(nm).strip()
    data["names"] = names

    # files in data/ that no section uses
    for folder, pat in EXTRA_RECOGNISED:
        recognised.update(rel(p) for p in find_files([(folder, pat)]))
    unrec = []
    if DATA.is_dir():
        for p in DATA.rglob("*"):
            if p.is_file() and not any(part.startswith(".") for part in p.relative_to(DATA).parts):
                if rel(p) not in recognised:
                    unrec.append(rel(p))
    data["unrecognised"] = sorted(unrec)
    return data, status


def entry_n(entry) -> int:
    return int(entry["rows"] or 0)


@st.cache_data(show_spinner=False)
def _load_all_cached(sig):
    return _load_all_uncached()


def load_all():
    """Find and read every file ONCE per change. Returns (data, status)."""
    try:
        return _load_all_cached(_signature())
    except Exception as exc:  # never crash the app
        st.error(f"Could not read the data folder: {exc}")
        return {"names": {}, "unrecognised": []}, {k: {"key": k, "label": v["label"], "feeds": v["feeds"],
                "kind": v["kind"], "owner": v.get("owner", ""), "day": v.get("day", ""), "state": "missing",
                "file": "", "rows": None, "expected": None, "problems": []} for k, v in FILES.items()}


def get_data():
    return load_all()


def usable(status, key) -> bool:
    return status.get(key, {}).get("state") in ("ok", "warning")


def step_state(status, step) -> str:
    """READY if all files ok; PARTIAL if some usable; WAITING if none."""
    states = [status[k]["state"] for k in step["needs"]]
    if all(s == "ok" for s in states):
        return "READY"
    if any(s in ("ok", "warning") for s in states):
        return "PARTIAL"
    return "WAITING"


def short(s, n: int = 42) -> str:
    s = str(s)
    return s if len(s) <= n else s[: n - 1] + "…"


def display_name(data, cid) -> str:
    nm = data.get("names", {}).get(cid)
    return nm if nm else f"{cid} (unnamed)"


def get_ml(data):
    """ML table: the detail file, else the compact one. Returns (df or None, which)."""
    if data.get("ml_detail") is not None:
        return data["ml_detail"], "detail"
    if data.get("ml_compact") is not None:
        return data["ml_compact"], "compact"
    return None, None


def get_redock(data) -> dict:
    """Redocking facts per protein key (dpp4/herg/cyp3a4/cyp2d6).
    Main source = the redocking file; DPP-4 can fall back to the target card / summary report."""
    out: dict[str, dict] = {}
    rd = data.get("redocking")
    if rd is not None:
        for _, r in rd.iterrows():
            out[r["receptor"]] = {"verdict": r.get("verdict"), "top1_rmsd": r.get("top1_rmsd"),
                                  "best_rmsd": r.get("best_rmsd"), "ref_score": r.get("ref_score"),
                                  "ligand": r.get("ligand"), "source": "redocking file"}
    if "dpp4" not in out:
        card = (data.get("target_card") or {}).get("docking", {}).get("reference_ligand_redocking", {})
        sr = data.get("summary_report")
        verdict = card.get("reported_verdict")
        ref = card.get("reference_docking_score_kcal")
        if sr is not None and "redock_verdict" in sr.columns and sr["redock_verdict"].notna().any():
            verdict = verdict or sr["redock_verdict"].dropna().iloc[0]
        if sr is not None and "reference_score_kcal" in sr.columns and sr["reference_score_kcal"].notna().any():
            ref = ref if ref is not None else float(sr["reference_score_kcal"].dropna().iloc[0])
        if verdict or ref is not None:
            out["dpp4"] = {"verdict": str(verdict).upper() if verdict else None, "top1_rmsd": np.nan,
                           "best_rmsd": np.nan, "ref_score": ref, "ligand": "605",
                           "source": "target card / summary report"}
    return out


def card_value(data, path: str, fallback=None):
    """Read 'target.name' style paths from target_card.json; fall back to the verified facts."""
    cur = data.get("target_card") or {}
    for part in path.split("."):
        if isinstance(cur, dict) and part in cur and cur[part] not in (None, ""):
            cur = cur[part]
        else:
            return fallback
    return cur


# =============================================================================
# 5. LOOK AND FEEL: theme, chips, headers, waiting card
# =============================================================================
def inject_css(dark: bool):
    bg, fg, card, line, muted = ("#0e1117", "#e8eaed", "#171b22", "#2b323d", "#9aa3b0") if dark else \
                                ("#ffffff", "#1d2430", "#f5f7fa", "#d9dee5", "#5b6572")
    st.markdown(f"""<style>
:root {{ --op-bg:{bg}; --op-fg:{fg}; --op-card:{card}; --op-line:{line}; --op-muted:{muted}; }}
.stApp, [data-testid="stAppViewContainer"], [data-testid="stHeader"] {{ background:var(--op-bg); color:var(--op-fg); }}
[data-testid="stSidebar"] {{ background:var(--op-card); }}
.stApp, .stApp p, .stApp li, .stApp label, .stApp h1, .stApp h2, .stApp h3, .stApp h4 {{ color:var(--op-fg); }}
[data-testid="stCaptionContainer"], .stApp small {{ font-size:13px !important; color:var(--op-muted) !important; }}
.op-chip {{ display:inline-block; padding:2px 9px; border-radius:999px; font-size:13px; font-weight:600;
  line-height:1.5; margin:1px 4px 1px 0; white-space:nowrap; border:1.5px solid transparent; }}
.op-green {{ background:#1f8a4c; color:#fff; }} .op-amber {{ background:#e0a100; color:#1d2430; }}
.op-red {{ background:#c0392b; color:#fff; }} .op-grey {{ background:#8a939e; color:#fff; }}
.op-blue {{ background:#2f6fd0; color:#fff; }} .op-outline {{ background:transparent; color:var(--op-fg); border-color:var(--op-muted); }}
.op-card {{ background:var(--op-card); border:1px solid var(--op-line); border-radius:10px; padding:12px 14px; height:100%; }}
.op-metric-l {{ font-size:13px; color:var(--op-muted); }} .op-metric-v {{ font-size:24px; font-weight:700; line-height:1.25; }}
.op-metric-s {{ font-size:13px; color:var(--op-muted); }}
.op-wait {{ background:var(--op-card); border:1px solid var(--op-line); border-left:6px solid #2f6fd0; border-radius:8px;
  padding:14px 16px; margin:8px 0; font-size:15px; }}
.op-hdr {{ display:flex; align-items:center; gap:14px; margin:4px 0 2px 0; }}
.op-step {{ min-width:48px; height:48px; border-radius:50%; color:#fff; font-weight:700; display:flex;
  align-items:center; justify-content:center; font-size:20px; }}
.op-hdr h2 {{ margin:0; padding:0; font-size:26px; }} .op-owner {{ font-size:14px; font-weight:600; }}
.op-blurb {{ font-size:16px; margin:2px 0 10px 0; }}
.op-foot {{ border-top:1px solid var(--op-line); margin-top:22px; padding-top:6px; font-size:13px; color:var(--op-muted); }}
.op-kv td {{ padding:3px 12px 3px 0; vertical-align:top; font-size:14px; }} .op-kv td:first-child {{ color:var(--op-muted); white-space:nowrap; }}
a {{ color:#4c8dff; }}
</style>""", unsafe_allow_html=True)


CHIP_KIND = {
    "PASS": "green", "READY": "green", "OK": "green", "FOUND": "green",
    "PARTIAL": "amber", "FLAGGED": "amber", "LOW CONFIDENCE": "amber", "PROBLEMS": "amber",
    "FAIL": "red", "INVALID": "red",
    "NO VALID POSE": "grey", "NOT TESTED": "outline", "NOT FLAGGED": "outline",
    "WAITING": "blue", "PENDING": "blue", "MISSING": "blue",
}


def chip(text, kind: str | None = None, color: str | None = None) -> str:
    """One status chip as HTML. Always has TEXT, never colour alone."""
    t = str(text)
    k = kind or CHIP_KIND.get(t.upper(), "outline")
    style = f' style="background:{color};color:#fff"' if color else ""
    return f'<span class="op-chip op-{k}"{style}>{esc(t)}</span>'


def verdict_chip(v) -> str:
    return chip(v if not is_missing(v) and str(v).strip() else "NOT TESTED")


def show_html(s: str):
    st.markdown(s, unsafe_allow_html=True)


def go(section_name: str):
    """Button callback: jump to a section."""
    st.session_state["section"] = section_name


def open_profile(cid: str, section_name: str = "8. Compound profile"):
    """Button callback: open the Compound profile of one compound."""
    st.session_state["profile_cid"] = cid
    st.session_state["section"] = section_name


def section_header(sec: dict, status):
    color = ROLE_COLORS.get(sec["owner"], "#5F6B73")
    state = section_state(sec, status)
    badge = f'{sec["step"]}' if sec.get("step") else ("★" if sec["id"] == 0 else str(sec["id"]))
    step_txt = f"Step {sec['step']} of the relay · " if sec.get("step") else ""
    show_html(
        f'<div class="op-hdr"><div class="op-step" style="background:{color}">{esc(badge)}</div>'
        f'<div style="flex:1"><h2>{esc(sec["name"])}</h2>'
        f'<div class="op-owner" style="color:{color}">● {esc(step_txt)}{esc(sec["owner"])}</div></div>'
        f'<div>{chip(state)}</div></div><p class="op-blurb">{esc(sec["blurb"])}</p>')


def section_footer(text: str):
    show_html(f'<div class="op-foot">{esc(text)}</div>')


def waiting_card(step: dict, status: dict, only: list[str] | None = None):
    """The standard WAITING card. Lists the missing / invalid files of the step."""
    keys = only or step["needs"]
    lines = []
    for k in keys:
        s = status[k]
        if s["state"] in ("missing", "invalid"):
            lines.append(esc(FILES[k]["label"]))
    what = ", ".join(lines) if lines else "files for this step"
    show_html(f'<div class="op-wait"><b>Waiting for {what}</b> &nbsp;|&nbsp; from: <b>{esc(step["owner"])}</b> '
              f'&nbsp;|&nbsp; expected: <b>Day {esc(step["day"])}</b> &nbsp;|&nbsp; '
              f'Drop it in <code>data/incoming/</code> and refresh.</div>')
    for k in keys:
        if status[k]["state"] == "invalid":
            for p in status[k]["problems"]:
                st.warning(p)


def problem_notes(status: dict, keys: list[str], limit: int = 3):
    """Short, specific messages for files that loaded with problems."""
    msgs = [p for k in keys if status[k]["state"] == "warning" for p in status[k]["problems"]]
    for m in msgs[:limit]:
        st.warning(m)
    if len(msgs) > limit:
        st.caption(f"+ {len(msgs) - limit} more note(s) in 10. Data inbox.")


def metric_card(label: str, value: str | None, sub: str = ""):
    v = f'<div class="op-metric-v">{esc(value)}</div>' if value else f'<div class="op-metric-v">{chip("WAITING")}</div>'
    show_html(f'<div class="op-card"><div class="op-metric-l">{esc(label)}</div>{v}<div class="op-metric-s">{esc(sub)}</div></div>')


# ---- 2D structure and 3D viewer -------------------------------------------------
@st.cache_data(show_spinner=False)
def mol_svg(smiles: str, w: int = 260, h: int = 200):
    """2D drawing as SVG text, or None if RDKit is not installed / the SMILES cannot be read."""
    if not HAVE_RDKIT or is_missing(smiles):
        return None
    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            return None
        d = rdMolDraw2D.MolDraw2DSVG(w, h)
        d.drawOptions().clearBackground = True
        d.DrawMolecule(mol)
        d.FinishDrawing()
        return d.GetDrawingText()
    except Exception:
        return None


def show_structure(smiles, w: int = 260, h: int = 200):
    svg = mol_svg(smiles, w, h)
    if svg:
        st.image(svg)
    elif not is_missing(smiles):
        st.caption("2D drawing needs RDKit. SMILES text:")
        st.code(str(smiles), language=None)
    else:
        st.caption("No structure available.")


_VIEWER = """<div style="position:relative;width:100%;height:__H__px;font-family:sans-serif">
<div id="v" style="width:100%;height:100%;position:relative"></div>
<div id="leg" style="position:absolute;left:8px;bottom:8px;font:13px sans-serif;padding:4px 8px;border-radius:6px;background:__LBG__;color:__LFG__;z-index:5"></div>
<div id="msg" style="position:absolute;left:12px;top:12px;font:14px sans-serif;color:__LFG__;z-index:5"></div></div>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<script>
const P = __PAYLOAD__;
(function(){
  const msg = document.getElementById('msg');
  if (typeof $3Dmol === 'undefined') { msg.textContent = '3D viewer could not load (it needs internet access to 3dmol.org).'; return; }
  try {
    const v = $3Dmol.createViewer(document.getElementById('v'), {backgroundColor: P.bg});
    const models = [];
    P.layers.forEach(function(L){
      const m = v.addModel(L.data, L.fmt); models.push(m);
      L.styles.forEach(function(S){ m.setStyle(S.sel, S.style); });
    });
    if (P.around) {
      const sel = {within: {distance: P.around.distance, sel: {model: models[P.around.layer]}}, byres: true, chain: 'A', hetflag: false};
      v.addStyle(sel, {stick: {radius: 0.12, colorscheme: 'whiteCarbon'}});
      v.addResLabels(sel, {font: 'sans-serif', fontSize: 11, fontColor: P.labelfg, showBackground: true, backgroundColor: P.labelbg, backgroundOpacity: 0.7});
    }
    v.zoomTo({model: models[P.zoom_layer]});
    v.render();
    document.getElementById('leg').innerHTML = P.legend.map(function(x){return '<span style="color:'+x[1]+'">&#9632;</span> '+x[0];}).join(' &nbsp; ');
  } catch (e) { msg.textContent = '3D viewer problem: ' + e; }
})();
</script>"""


def viewer(layers: list[dict], legend: list[tuple[str, str]], dark: bool, height: int = 520,
           around: dict | None = None, zoom_layer: int = 0):
    """Embed 3Dmol.js. layers = [{"data": text, "fmt": "pdb"|"pdbqt", "styles": [{"sel": {...}, "style": {...}}]}].
    around = {"layer": index, "distance": 4} labels the protein residues near that layer."""
    payload = {"layers": layers, "legend": legend, "around": around, "zoom_layer": zoom_layer,
               "bg": "#0e1117" if dark else "#ffffff",
               "labelfg": "#ffffff" if dark else "#111111", "labelbg": "#2b323d" if dark else "#e9edf2"}
    page = (_VIEWER.replace("__H__", str(height)).replace("__PAYLOAD__", json.dumps(payload).replace("</", "<\\/"))
            .replace("__LBG__", "#171b22" if dark else "#f1f3f6").replace("__LFG__", "#e8eaed" if dark else "#1d2430"))
    components.html(page, height=height + 6)


PROTEIN_STYLE = {"sel": {"chain": "A", "hetflag": False}, "style": {"cartoon": {"color": "#9aa0a6"}}}


def ligand_style(color: str, sel: dict | None = None) -> dict:
    return {"sel": sel or {}, "style": {"stick": {"radius": 0.22, "color": color}}}


def is_dark() -> bool:
    return bool(st.session_state.get("dark", False))


# =============================================================================
# 6. SECTION REGISTRY  (order = sidebar order)
# =============================================================================
SECTIONS = [
    dict(id=0, name="Home", owner="Everyone", step=None, func="section_home",
         blurb="The project in one page: what we did, what we found, and what is still waiting."),
    dict(id=1, name="Target & library", owner="Chemist", step=1, func="section_target",
         blurb="Picks the disease protein and the neem compounds we will test."),
    dict(id=2, name="Structures & validation", owner="Bioinformatics", step=2, func="section_structures",
         blurb="Prepares the proteins and checks that our docking can put a known drug back where the crystal had it."),
    dict(id=3, name="ADMET filter", owner="ML Engineer", step=3, func="section_admet",
         blurb="Filters the compounds for drug-likeness and predicted risk before any docking."),
    dict(id=4, name="Docking: DPP-4", owner="Bioinformatics", step=4, func="section_docking",
         blurb="Docks every compound that passed the filter into DPP-4 and scores how well it fits."),
    dict(id=5, name="Safety panel", owner="Bioinformatics", step=5, func="section_safety",
         blurb="Docks the best compounds into three safety proteins and compares each with its own known drug."),
    dict(id=6, name="Final ranking", owner="ML Engineer", step=6, func="section_ranking",
         blurb="Combines fit, safety docking and drug-likeness into one ranked list."),
    dict(id=7, name="Simulation (MD)", owner="Bioinformatics", step=7, func="section_md",
         blurb="Runs a computer simulation of the top compound to see if its pose stays stable over time."),
    dict(id=8, name="Compound profile", owner="Everyone", step=None, func="section_profile",
         blurb="Everything we know about one compound, with its 3D pose in the protein."),
    dict(id=9, name="Evidence & methods", owner="Project Manager", step=None, func="section_evidence",
         blurb="How the project meets the judging criteria, its limits, its methods and a glossary."),
    dict(id=10, name="Data inbox", owner="Developer", step=None, func="section_inbox",
         blurb="Shows which files the app found, which section each one feeds, and lets you add new files."),
]
SECTION_NAMES = [f"{s['id']}. {s['name']}" if s["id"] else "0. Home" for s in SECTIONS]
SECTION_BY_NAME = dict(zip(SECTION_NAMES, SECTIONS))
NAME_BY_STEP = {s["step"]: n for n, s in SECTION_BY_NAME.items() if s["step"]}


def section_state(sec: dict, status: dict) -> str:
    if sec["step"]:
        return step_state(status, STEP_BY_ID[sec["step"]])
    if sec["id"] == 0:
        states = [step_state(status, s) for s in STEPS]
        return "READY" if all(x == "READY" for x in states) else ("WAITING" if all(x == "WAITING" for x in states) else "PARTIAL")
    if sec["id"] == 8:
        return {"ok": "READY", "warning": "PARTIAL"}.get(status["docking_target"]["state"], "WAITING")
    return "READY"


# =============================================================================
# 7. ENGINE PLACEHOLDERS  (later versions plug real engines in here)
# =============================================================================
def docking_engine_run(*args, **kwargs):
    """PLACEHOLDER for a real docking engine (AutoDock Vina). NOT CONNECTED.
    When connected it must WRITE the same files the sections already read, so no UI change is needed:
      data/incoming/docking_target.csv        compound_id, vina_score_kcal, pose_file, heavy_atoms, ligand_efficiency
      data/incoming/docking_antitargets.csv   compound_id, antitarget, vina_score_kcal, pose_file, target_score_kcal,
                                              delta_vs_target_kcal, stronger_than_target
      data/incoming/results/<protein>/poses/<compound_id>_out.pdbqt     (multi-model, MODEL 1 = best pose)
    """
    raise NotImplementedError("Docking engine is not connected in this version.")


def md_engine_run(*args, **kwargs):
    """PLACEHOLDER for a GROMACS run. NOT CONNECTED.
    When connected it must WRITE:
      data/incoming/md_report/md_summary.csv   compound_id, label (candidate|control), length_ns, ligand_rmsd_mean_A,
                                               ligand_rmsd_sd_A, protein_rmsd_mean_A, hbond_occupancy_pct, notes
      data/incoming/md_report/<compound_id>_rmsd.png, <compound_id>_rmsf.png, <compound_id>_hbonds.png
    """
    raise NotImplementedError("GROMACS engine is not connected in this version.")


# =============================================================================
# 8. SECTION 0: HOME
# =============================================================================
def section_home():
    data, status = get_data()
    sec = SECTION_BY_NAME["0. Home"]
    section_header(sec, status)
    st.write(STORY)

    # (b) target card -------------------------------------------------------
    st.subheader("The target")
    rows = [
        ("Disease", card_value(data, "disease", TARGET_FALLBACK["disease"])),
        ("Target protein", card_value(data, "target.name", TARGET_FALLBACK["target_name"])),
        ("PDB structure", f'{card_value(data, "target.pdb_id", TARGET_FALLBACK["pdb_id"])} · '
                          f'{card_value(data, "target.resolution_angstrom", TARGET_FALLBACK["resolution"])} Å · '
                          f'{card_value(data, "target.method", TARGET_FALLBACK["method"])} · '
                          f'chain {card_value(data, "target.chain_used", TARGET_FALLBACK["chain"])}'),
        ("Known drug in the crystal", f'{card_value(data, "target.reference_ligand_code", TARGET_FALLBACK["ref_code"])}'
                                      f' ({card_value(data, "target.reference_ligand_name", "name not in target card")})'),
        ("Plant", f'{card_value(data, "plant_source.scientific_name", "Azadirachta indica")} '
                  f'({card_value(data, "plant_source.common_name", "neem")})'),
        ("Why DPP-4?", card_value(data, "biological_rationale", "-")),
    ]
    show_html('<table class="op-kv">' + "".join(f"<tr><td>{esc(k)}</td><td>{esc(v)}</td></tr>" for k, v in rows) + "</table>")
    if data.get("target_card") is None:
        st.caption("target_card.json was not found, so this card shows the verified project facts.")

    # (c) four metric cards --------------------------------------------------
    st.subheader("Key numbers")
    ml, _ = get_ml(data)
    dt = data.get("docking_target")
    rd = get_redock(data)
    n_total = len(ml) if ml is not None else None
    n_pass = int(ml["pass_flag"].fillna(False).sum()) if ml is not None else None
    n_dock = len(dt) if dt is not None else None
    best = dt.iloc[0] if dt is not None else None
    ref = rd.get("dpp4", {}).get("ref_score")
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        metric_card("Compounds that passed the ML filter", f"{n_pass} of {n_total}" if ml is not None else None,
                    "drug-likeness and ADMET risk filter")
    with c2:
        if n_dock is None:
            metric_card("Compounds docked into DPP-4", None)
        else:
            extra = f"{n_dock} of {n_pass} passed" if n_pass else ""
            metric_card("Compounds docked into DPP-4", str(n_dock), extra)
    with c3:
        if best is None:
            metric_card("Best DPP-4 score vs reference drug", None)
        else:
            refs = fmt(ref, 3) if not is_missing(ref) else "-"
            metric_card("Best DPP-4 score vs reference drug", f"{fmt(best['vina_score_kcal'], 3)} vs {refs}",
                        f"kcal/mol, predicted · best = {short(display_name(data, best['compound_id']))} ({best['compound_id']}); "
                        "reference = known drug 605 redocked")
    with c4:
        v = rd.get("dpp4", {}).get("verdict")
        r = rd.get("dpp4", {}).get("top1_rmsd")
        if not v:
            metric_card("Redocking check (DPP-4)", None)
        else:
            metric_card("Redocking check (DPP-4)", str(v), f"top-pose RMSD {fmt(r, 2, 'Å')}" if not is_missing(r) else "known drug docked back")

    # (d) relay strip ----------------------------------------------------------
    st.subheader("The relay: who hands what to whom")
    cols = st.columns(len(STEPS))
    for col, step in zip(cols, STEPS):
        color = ROLE_COLORS[step["owner"]]
        state = step_state(status, step)
        with col:
            show_html(f'<div class="op-card" style="border-top:5px solid {color}">'
                      f'<div style="font-size:13px;color:{color};font-weight:700">Step {step["id"]} · {esc(step["owner"])}</div>'
                      f'<div style="font-size:15px;font-weight:600;margin:3px 0 6px 0">{esc(step["name"])}</div>{chip(state)}</div>')
            st.button("Open", key=f"relay_{step['id']}", on_click=go, args=(NAME_BY_STEP[step["id"]],))
    st.caption("The Software Developer connects every step into this app. The Project Manager runs the schedule and the final deck.")

    # (e) honest banner ---------------------------------------------------------
    st.subheader("Read this first")
    msgs = ["These are **computer predictions, not lab measurements**. Nothing here shows that a compound works or is safe in a body."]
    if best is not None and not is_missing(ref):
        if float(best["vina_score_kcal"]) > float(ref):
            msgs.append(f"**No neem compound beats the known drug** in predicted DPP-4 fit: the best neem score is "
                        f"{fmt(best['vina_score_kcal'], 3)} kcal/mol and the known drug scored {fmt(ref, 3)} kcal/mol.")
        else:
            msgs.append("At least one neem compound scores as well as or better than the known drug in this prediction. "
                        "Treat that as a prompt to check the pose, not as proof.")
    h = rd.get("herg", {}).get("verdict")
    if h and str(h).upper() != "PASS":
        msgs.append(f"The hERG redocking check was **{esc(h)}**, so every hERG score is **low confidence**.")
    if dt is not None:
        msgs.append(f"The DPP-4 screen used the lowest search effort (exhaustiveness {EXHAUSTIVENESS_SCREEN}), so scores within about "
                    f"{NOISE_KCAL} kcal/mol of each other are effectively tied.")
    st.warning("\n\n".join(f"- {m}" for m in msgs))
    section_footer("OmniPhyto · EthnoHack 2026 · Natural Products for Drug Discovery · all values are predictions.")


# =============================================================================
# 9. SECTION 10: DATA INBOX
# =============================================================================
STATE_TEXT = {"ok": "✅ found", "warning": "⚠️ found, with problems", "invalid": "❌ invalid", "missing": "⏳ missing"}


def pattern_specificity(pat: str) -> int:
    return len(re.sub(r"[*?\[\]]", "", pat))


def identify_upload(name: str, raw: bytes):
    """Decide what an uploaded file is: by file-name pattern first, then by its columns.
    Returns dict(keys=[candidates], how=str)."""
    lname = name.lower()
    is_table = lname.endswith((".csv", ".json"))
    hits = []
    for key, spec in FILES.items():
        # tables only match table entries; structure files only match folder entries
        if (spec["kind"] == "folder") == is_table:
            continue
        for _, pat in spec["patterns"]:
            if pat.strip("*") == "":      # catch-all pattern such as "*": too vague to decide anything
                continue
            if fnmatch.fnmatch(lname, pat.lower()):
                hits.append((pattern_specificity(pat), key))
                break
    if hits:
        top = max(h[0] for h in hits)
        return {"keys": [k for s, k in hits if s == top], "how": "its file name"}
    # by content
    if lname.endswith(".json"):
        try:
            obj = json.loads(raw.decode("utf-8", errors="replace"))
            if isinstance(obj, dict) and ("target" in obj or "project" in obj):
                return {"keys": ["target_card"], "how": "its contents"}
        except Exception:
            pass
        return {"keys": [], "how": ""}
    if lname.endswith(".csv"):
        best, best_keys = -1, []
        for key, spec in FILES.items():
            if spec["kind"] != "csv" or not spec.get("required"):
                continue
            try:
                df = _read_csv_any(io.BytesIO(raw))
            except Exception:
                return {"keys": [], "how": ""}
            df = canon_columns(df, spec.get("aliases"), spec["required"], spec.get("numeric", []))
            if all(c in df.columns for c in spec["required"]):
                score = len(spec["required"]) + 0.01 * sum(c in df.columns for c in spec.get("numeric", []))
                if score > best + 1e-9:
                    best, best_keys = score, [key]
                elif abs(score - best) < 1e-9:
                    best_keys.append(key)
        return {"keys": best_keys, "how": "its columns"}
    return {"keys": [], "how": ""}


def check_upload(key: str, name: str, raw: bytes):
    """Check the upload against the contract. Returns (rows or None, problems, usable?)."""
    spec = FILES[key]
    if spec["kind"] == "csv":
        df, probs, state = table_from(key, io.BytesIO(raw), name)
        return (len(df) if df is not None else None), probs, state != "invalid"
    if spec["kind"] == "json":
        try:
            obj = json.loads(raw.decode("utf-8", errors="replace"))
            ok = isinstance(obj, dict) and bool(obj)
            return None, ([] if ok else [f"{name}: not a JSON object"]), ok
        except Exception as exc:
            return None, [f"{name}: bad JSON ({str(exc)[:60]})"], False
    text = raw.decode("utf-8", errors="replace")
    n_atoms = sum(1 for ln in text.splitlines() if ln.startswith(("ATOM", "HETATM")))
    probs = [] if n_atoms else [f"{name}: no ATOM/HETATM lines found"]
    return None, probs, n_atoms > 0


def clean_filename(name: str) -> str | None:
    n = Path(name).name
    return n if re.fullmatch(r"[A-Za-z0-9._ \-()]+", n) and not n.startswith(".") else None


def section_inbox():
    data, status = get_data()
    sec = SECTION_BY_NAME["10. Data inbox"]
    section_header(sec, status)

    # (a) status table ---------------------------------------------------------
    st.subheader("What the app found")
    rows = []
    for key, s in status.items():
        if s["kind"] == "folder":
            exp = f" of {s['expected']}" if s.get("expected") else ""
            count = f"{s['rows'] or 0}{exp} files"
            found = "-"
        else:
            count = f"{s['rows']} rows" if s["rows"] is not None else ("JSON card" if s["state"] == "ok" and s["kind"] == "json" else "-")
            found = s["file"] or "-"
        rows.append({"File needed": FILES[key]["label"], "State": STATE_TEXT.get(s["state"], s["state"]),
                     "Found as": found, "Size": count, "Problems": "; ".join(s["problems"]) or "none",
                     "Feeds section": ", ".join(s["feeds"]), "From": s["owner"], "Expected day": s["day"]})
    st.dataframe(pd.DataFrame(rows), hide_index=True)

    unrec = data.get("unrecognised", [])
    with st.expander(f"Files in data/ that no section uses ({len(unrec)})"):
        if unrec:
            for u in unrec[:200]:
                st.text(f"{u}   -> not used by any section")
            if len(unrec) > 200:
                st.caption(f"... and {len(unrec) - 200} more")
        else:
            st.write("None. Every file in data/ feeds a section.")

    # (b) drop a file -----------------------------------------------------------
    st.subheader("Drop a file")
    st.caption("CSV, JSON, PDB or PDBQT. The app works out what it is and which section it will fill.")
    ups = st.file_uploader("Choose one or more files", accept_multiple_files=True,
                           type=["csv", "json", "pdb", "pdbqt"], key="inbox_upload")
    for i, up in enumerate(ups or []):
        raw = up.getvalue()
        ident = identify_upload(up.name, raw)
        keys = ident["keys"]
        with st.container(border=True):
            st.markdown(f"**{esc(up.name)}**")
            if not keys:
                st.warning("The app cannot tell what this file is (name and columns do not match any expected file). It was not saved.")
                continue
            key = keys[0]
            if len(keys) > 1:
                key = st.selectbox("Several places fit. Which one is it?", keys, key=f"pick_{i}",
                                   format_func=lambda k: FILES[k]["label"])
            spec = FILES[key]
            n_rows, probs, ok = check_upload(key, up.name, raw)
            feeds = " and ".join(spec["feeds"])
            size_txt = f" It has {n_rows} rows." if n_rows is not None else ""
            probs_txt = "; ".join(probs) if probs else "none"
            st.markdown(f"This looks like **{esc(key)}** (matched by {ident['how']}). It feeds **Section {esc(feeds)}**."
                        f"{size_txt} Problems: {esc(probs_txt)}.")
            if not ok:
                st.error("This file cannot be used as it is, so saving is switched off.")
                continue
            target_name = clean_filename(up.name)
            if spec["kind"] in ("csv", "json") and not any(
                    fnmatch.fnmatch(up.name.lower(), p.lower()) for _, p in spec["patterns"]):
                target_name = spec["canonical"]
                st.caption(f"The name does not match the pattern, so it will be saved as {target_name}.")
            if target_name is None:
                st.error("This file name has unusual characters. Please rename it and try again.")
                continue
            dest = f'{spec["dest"]}/{target_name}'
            if st.button(f"Save into {spec['dest']}/", key=f"save_{i}_{up.name}"):
                try:
                    out = ROOT / dest
                    out.parent.mkdir(parents=True, exist_ok=True)
                    out.write_bytes(raw)
                    st.cache_data.clear()
                    st.success(f"Saved {dest}. Refresh the page (or switch section) and it will fill in.")
                except OSError as exc:
                    st.error(f"Could not save: {exc}")
    st.info("On a hosted website, a saved file may disappear when the app restarts. To keep it permanently, "
            "commit it to GitHub in `data/incoming/`.")

    # (c) engines ---------------------------------------------------------------
    with st.expander("Advanced: engines (later versions)"):
        c1, c2 = st.columns(2)
        with c1:
            st.button("Docking engine: not connected", disabled=True, key="eng_dock")
            st.caption("Will write docking_target.csv, docking_antitargets.csv and the pose files (see docking_engine_run in app.py).")
        with c2:
            st.button("GROMACS engine: not connected", disabled=True, key="eng_md")
            st.caption("Will write md_summary.csv and the MD plots (see md_engine_run in app.py).")
    section_footer("Data inbox · refresh the page after adding files; the app reads them again automatically.")


# ===== PART 2 GOES HERE =====  (sections 1-5: section_target, section_structures, section_admet, section_docking, section_safety)

# ===== PART 3 GOES HERE =====  (sections 6-9: section_ranking, section_md, section_profile, section_evidence)


# =============================================================================
# 10. SIDEBAR + ROUTER + MAIN   (keep this block LAST in the file)
# =============================================================================
def not_added_yet(sec: dict, data_status):
    _, status = data_status
    section_header(sec, status)
    st.info(f"The code for '{sec['name']}' has not been pasted into app.py yet. "
            f"It arrives in {'Part 2' if sec['id'] <= 5 else 'Part 3'}.")


def sidebar(status: dict):
    sb = st.sidebar
    sb.markdown("## 🍃 OmniPhyto")
    sb.caption("DPP-4 | PDB 3D4L | type 2 diabetes")
    sb.toggle("Dark theme", key="dark")
    labels = {}
    for name, sec in SECTION_BY_NAME.items():
        labels[name] = f"{ROLE_DOT[sec['owner']]} {name} · {section_state(sec, status)}"
    sb.radio("Sections", SECTION_NAMES, key="section", format_func=lambda n: labels[n])
    sb.markdown(
        "<div style='font-size:13px'>" + " ".join(
            f"<span style='color:{ROLE_COLORS[r]}'>●</span> {r}" for r in
            ["Chemist", "Bioinformatics", "ML Engineer", "Developer", "Project Manager"]) + "</div>",
        unsafe_allow_html=True)


def main():
    st.session_state.setdefault("section", "0. Home")
    st.session_state.setdefault("dark", False)
    data_status = get_data()
    inject_css(is_dark())
    sidebar(data_status[1])
    name = st.session_state["section"]
    sec = SECTION_BY_NAME.get(name, SECTIONS[0])
    fn = globals().get(sec["func"])
    if fn is None:
        not_added_yet(sec, data_status)
        return
    try:
        fn()
    except Exception as exc:  # a broken section must never take the whole app down
        st.error(f"This section hit a problem and was skipped: {type(exc).__name__}: {exc}")
        st.caption("Other sections still work. Check 10. Data inbox for file problems.")


main()
