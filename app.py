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
import altair as alt
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
EXTRA_RECOGNISED = [("data/incoming/md_report", "README.md"), ("data/incoming/md_report", "*_hbonds_detail.csv")]


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
.op-tbl {{ border-collapse:collapse; width:100%; }}
.op-tbl th, .op-tbl td {{ border-bottom:1px solid var(--op-line); padding:6px 8px; font-size:14px; text-align:left; vertical-align:top; }}
.op-tbl th {{ color:var(--op-muted); font-weight:600; }}
</style>""", unsafe_allow_html=True)


CHIP_KIND = {
    "PASS": "green", "READY": "green", "OK": "green", "FOUND": "green",
    "PARTIAL": "amber", "FLAGGED": "amber", "LOW CONFIDENCE": "amber", "PROBLEMS": "amber",
    "FAIL": "red", "INVALID": "red",
    "NO VALID POSE": "grey", "NOT TESTED": "outline", "NOT FLAGGED": "outline",
    "WAITING": "blue", "PENDING": "blue", "MISSING": "blue", "BENCHMARK PENDING": "blue",
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
      const sel = {within: {distance: P.around.distance, sel: Object.assign({model: models[P.around.layer]}, P.around.sel || {})}, byres: true, chain: 'A', hetflag: false};
      v.addStyle(sel, {stick: {radius: 0.12, colorscheme: 'whiteCarbon'}});
      v.addResLabels(sel, {font: 'sans-serif', fontSize: 11, fontColor: P.labelfg, showBackground: true, backgroundColor: P.labelbg, backgroundOpacity: 0.7});
    }
    v.zoomTo(P.zoom_sel || {model: models[P.zoom_layer]});
    v.render();
    document.getElementById('leg').innerHTML = P.legend.map(function(x){return '<span style="color:'+x[1]+'">&#9632;</span> '+x[0];}).join(' &nbsp; ');
  } catch (e) { msg.textContent = '3D viewer problem: ' + e; }
})();
</script>"""


def viewer(layers: list[dict], legend: list[tuple[str, str]], dark: bool, height: int = 520,
           around: dict | None = None, zoom_layer: int = 0, zoom_sel: dict | None = None):
    """Embed 3Dmol.js. layers = [{"data": text, "fmt": "pdb"|"pdbqt", "styles": [{"sel": {...}, "style": {...}}]}].
    around = {"layer": index, "distance": 4, "sel": {...}} labels the protein residues near that ligand;
    zoom_sel (optional) is the atom selection to zoom to."""
    payload = {"layers": layers, "legend": legend, "around": around, "zoom_layer": zoom_layer, "zoom_sel": zoom_sel,
               "bg": "#0e1117" if dark else "#ffffff",
               "labelfg": "#ffffff" if dark else "#111111", "labelbg": "#2b323d" if dark else "#e9edf2"}
    page = (_VIEWER.replace("__H__", str(height)).replace("__PAYLOAD__", json.dumps(payload).replace("</", "<\\/"))
            .replace("__LBG__", "#171b22" if dark else "#f1f3f6").replace("__LFG__", "#e8eaed" if dark else "#1d2430"))
    if hasattr(st, "iframe"):          # newer Streamlit: components.html is being retired
        st.iframe(page, height=height + 6)
    else:
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
# =============================================================================
# PART 2 - SHARED HELPERS FOR SECTIONS 1-9
# =============================================================================
GREEN_POSE = "#2EA043"     # compound pose in the 3D viewer
ORANGE_XTL = "#E8872B"     # known (crystal) drug in the 3D viewer
SAFETY_PROTEINS = ["herg", "cyp3a4", "cyp2d6"]
PROT_LABEL = {k: v["name"] for k, v in PROTEINS.items()}


def style_chart(ch):
    """Make an Altair chart readable in both light and dark theme (text >= 13 px)."""
    fg = "#e8eaed" if is_dark() else "#1d2430"
    return (ch.properties(background="transparent")
            .configure_axis(labelColor=fg, titleColor=fg, gridColor="#8883", labelFontSize=13, titleFontSize=13)
            .configure_legend(labelColor=fg, titleColor=fg, labelFontSize=13, titleFontSize=13)
            .configure_title(color=fg, fontSize=15).configure_view(strokeWidth=0))


def compounds_frame(data):
    """One row per DPP-4-docked compound, with names and ML columns added. None if no docking file."""
    dt = data.get("docking_target")
    if dt is None:
        return None
    df = dt.copy()
    df["name"] = df["compound_id"].map(lambda c: display_name(data, c))
    ml, _ = get_ml(data)
    if ml is not None:
        cols = [c for c in ["source_plant", "qed", "lipinski_violations", "herg_risk", "hepatotox_risk",
                            "admet_score", "pass_flag", "fail_reason"] if c in ml.columns and c not in df.columns]
        df = df.merge(ml[["compound_id"] + cols], on="compound_id", how="left")
    return df


def safety_matrix(data) -> dict:
    """{compound_id: {'herg': score, 'cyp3a4': score, 'cyp2d6': score}}. Antitarget file first,
    summary_report columns as a fallback. Missing = not in the dict (shown as 'not tested')."""
    out: dict[str, dict] = {}
    at = data.get("docking_antitargets")
    if at is not None:
        for _, r in at.iterrows():
            k = ANTITARGET_KEY.get(r["antitarget"])
            if k:
                out.setdefault(r["compound_id"], {})[k] = r["vina_score_kcal"]
    sr = data.get("summary_report")
    if sr is not None:
        for _, r in sr.iterrows():
            for k, col in (("herg", "hERG_score_kcal"), ("cyp3a4", "CYP3A4_score_kcal"), ("cyp2d6", "CYP2D6_score_kcal")):
                if col in sr.columns and not is_missing(r[col]) and k not in out.get(r["compound_id"], {}):
                    out.setdefault(r["compound_id"], {})[k] = r[col]
    return out


def safety_status(score, ref) -> str:
    """Flagged = valid pose AND as good as or better (more negative) than the protein's OWN reference drug."""
    if is_missing(score):
        return "NOT TESTED"
    if not is_valid_pose(score):
        return "NO VALID POSE"
    if is_missing(ref):
        return "BENCHMARK PENDING"
    return "FLAGGED" if float(score) <= float(ref) else "NOT FLAGGED"


def ref_score(redock: dict, prot: str):
    v = redock.get(prot, {}).get("ref_score")
    return None if is_missing(v) else float(v)


def parse_conf(text: str) -> dict:
    out = {}
    for line in text.splitlines():
        m = re.match(r"\s*([A-Za-z_]+)\s*=\s*(-?[\d.]+)", line)
        if m:
            out[m.group(1).lower()] = float(m.group(2))
    return out


def conf_for(data, prot: str):
    """Box centre / size from the prepared receptor's *_vina.conf, or None."""
    pdb = PROTEINS[prot]["pdb"].lower()
    for p in data.get("receptors_prepared", []):
        n = Path(p).name.lower()
        if n.endswith("_vina.conf") and (prot in n or pdb in n):
            txt = read_text_safe(p, 1)
            if txt:
                c = parse_conf(txt)
                if all(k in c for k in ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z")):
                    return c
    return None


def file_for(data, key: str, *needles: str):
    """First repo-relative path in data[key] whose file name contains any needle (case-insensitive)."""
    for p in data.get(key, []):
        n = Path(p).name.lower()
        if any(x.lower() in n for x in needles):
            return p
    return None


def pdb_residue_lines(text: str, resn: str) -> str:
    """Keep only the HETATM lines of one residue name (used to pull the crystal drug 605 'XTL')."""
    keep = [ln for ln in text.splitlines() if ln.startswith(("HETATM", "ATOM")) and ln[17:20].strip() == resn]
    return "\n".join(keep) + "\nEND\n" if keep else ""


def score_background(df: pd.DataFrame, col: str):
    """Pandas Styler with a simple colour scale on the score column (best = green)."""
    s = df[col]
    valid = s[s.map(is_valid_pose)]
    lo, hi = (valid.min(), valid.max()) if len(valid) else (0, 1)

    def colour(v):
        if not is_valid_pose(v) or hi == lo:
            return ""
        t = (float(v) - lo) / (hi - lo)            # 0 = best, 1 = worst
        r, g, b = (int(a + (c - a) * t) for a, c in zip((143, 209, 158), (243, 244, 246)))
        return f"background-color: rgb({r},{g},{b}); color: #1d2430"
    return colour


def profile_button(cid_options: list[str], key: str, label: str = "Open compound profile"):
    """A small 'pick a compound -> open its profile' control."""
    if not cid_options:
        return
    c1, c2 = st.columns([3, 1])
    with c1:
        pick = st.selectbox("Compound", cid_options, key=f"{key}_pick", label_visibility="collapsed")
    with c2:
        st.button(label, key=f"{key}_btn", on_click=open_profile, args=(pick,))


def row_selection(event) -> list[int]:
    try:
        return list(event.selection.rows)
    except Exception:
        return []


# =============================================================================
# SECTION 1: TARGET & LIBRARY
# =============================================================================
def card_rows(data) -> list[tuple[str, str]]:
    g = lambda p, fb="-": card_value(data, p, fb)  # noqa: E731
    rows = [
        ("Disease", g("disease", TARGET_FALLBACK["disease"])),
        ("Target protein", g("target.name", TARGET_FALLBACK["target_name"])),
        ("Organism · class", f'{g("target.organism", TARGET_FALLBACK["organism"])} · {g("target.target_class")}'),
        ("PDB structure", f'{g("target.pdb_id", TARGET_FALLBACK["pdb_id"])} · {g("target.resolution_angstrom", TARGET_FALLBACK["resolution"])} Å · '
                          f'{g("target.method", TARGET_FALLBACK["method"])} · chain {g("target.chain_used", TARGET_FALLBACK["chain"])}'),
        ("Known drug in the crystal", f'{g("target.reference_ligand_code", TARGET_FALLBACK["ref_code"])} · {g("target.reference_ligand_name")}'),
        ("Reference paper", g("target.reference")),
        ("Why this target?", g("biological_rationale")),
        ("What we try to find", g("drug_discovery_rationale")),
        ("Plant", f'{g("plant_source.scientific_name", "Azadirachta indica")} ({g("plant_source.common_name", "neem")}; '
                  f'local name {g("plant_source.local_name")})'),
        ("Compound source", g("library.source")),
        ("Event · track", f'{g("hackathon")} · {g("track")}'),
    ]
    return rows


def section_target():
    data, status = get_data()
    sec = SECTION_BY_NAME["1. Target & library"]
    step = STEP_BY_ID[1]
    section_header(sec, status)
    master, card = data.get("master"), data.get("target_card")
    if master is None and card is None:
        waiting_card(step, status)
        section_footer("Step 1 hands the target card and the compound library to everyone else.")
        return
    problem_notes(status, ["target_card", "master"])
    if master is None:
        waiting_card(step, status, ["master"])

    # --- target card in full
    st.subheader("Target card")
    show_html('<table class="op-kv">' + "".join(f"<tr><td>{esc(k)}</td><td>{esc(v)}</td></tr>" for k, v in card_rows(data)) + "</table>")
    if card is None:
        waiting_card(step, status, ["target_card"])
        st.caption("Showing the verified project facts until target_card.json arrives.")
    panel = (card or {}).get("safety_panel")
    if isinstance(panel, list) and panel:
        st.markdown("**Safety proteins (from the target card)**")
        st.dataframe(pd.DataFrame(panel).rename(columns={"key": "Protein", "pdb_id": "PDB ID", "chains": "Chains",
                     "reference_ligand_code": "Known drug", "note": "Note"}), hide_index=True)

    # --- library counts
    st.subheader("Library")
    ml, _ = get_ml(data)
    lib = (card or {}).get("library", {}) if isinstance((card or {}).get("library"), dict) else {}
    total = len(master) if master is not None else lib.get("initial_compound_records")
    n_pass = int(ml["pass_flag"].fillna(False).sum()) if ml is not None else lib.get("passed_admet_filter")
    c1, c2, c3 = st.columns(3)
    with c1:
        metric_card("Neem compound records", str(total) if total else None, f'source: {esc(lib.get("source", "COCONUT"))}')
    with c2:
        metric_card("Passed the ML filter", str(n_pass) if n_pass is not None else None, "see 3. ADMET filter")
    with c3:
        n_doi = int(master["dois"].map(lambda s: len(doi_list(s)) > 0).sum()) if master is not None and "dois" in master.columns else None
        metric_card("Compounds with at least one DOI", f"{n_doi} of {total}" if n_doi is not None else None, "literature link in the master table")
    if master is None:
        section_footer("Step 1 · Chemist")
        return
    st.caption("Every compound here has a neem source in the master table, but botanical provenance and exact identities "
               "still need confirming against the primary literature.")

    # --- class charts
    st.subheader("What kinds of molecules?")
    cc1, cc2 = st.columns(2)
    for col, label, box in (("chemical_class", "Chemical class (top 10)", cc1), ("np_classifier_class", "NP-Classifier class (top 10)", cc2)):
        with box:
            if col not in master.columns:
                st.caption(f"Column {col} is not in the master table.")
                continue
            vc = master[col].dropna().astype(str).value_counts().head(10).rename_axis("class").reset_index(name="compounds")
            ch = alt.Chart(vc).mark_bar(color=ROLE_COLORS["Chemist"]).encode(
                x=alt.X("compounds:Q", title="compounds"), y=alt.Y("class:N", sort="-x", title=None),
                tooltip=["class", "compounds"]).properties(title=label, height=26 * max(len(vc), 3) + 30)
            st.altair_chart(style_chart(ch))
            st.caption(f"{int(master[col].isna().sum())} compounds have no {col.replace('_', ' ')} in the table.")

    # --- compound table
    st.subheader("Compound table")
    df = master.copy()
    df["name"] = df["compound_id"].map(lambda c: display_name(data, c))
    df["n_dois"] = df["dois"].map(lambda s: len(doi_list(s))) if "dois" in df.columns else 0
    if ml is not None:
        df = df.merge(ml[["compound_id", "pass_flag"]], on="compound_id", how="left")
    q = st.text_input("Search by name, ID or class", key="lib_search")
    prov_opts = sorted(df["provenance_status"].dropna().unique()) if "provenance_status" in df.columns else []
    prov = st.multiselect("Provenance status (from the master table)", prov_opts, default=prov_opts, key="lib_prov") if prov_opts else []
    view = df
    if q:
        blob = view[["compound_id", "name"] + [c for c in ("chemical_class", "np_classifier_class") if c in view.columns]].astype(str).agg(" ".join, axis=1)
        view = view[blob.str.contains(re.escape(q), case=False, na=False)]
    if prov:
        view = view[view["provenance_status"].isin(prov)]
    cols = [c for c in ["compound_id", "name", "molecular_formula", "chemical_class", "np_classifier_class",
                        "provenance_status", "n_dois", "pass_flag"] if c in view.columns]
    view = view[cols].reset_index(drop=True)
    cfg = {"compound_id": st.column_config.TextColumn("Compound ID", pinned=True),
           "name": "Name", "molecular_formula": "Formula", "chemical_class": "Chemical class",
           "np_classifier_class": "NP-Classifier class", "provenance_status": "Provenance status",
           "n_dois": st.column_config.NumberColumn("DOIs", help="Number of different DOIs (case duplicates removed)"),
           "pass_flag": st.column_config.CheckboxColumn("Passed ML filter")}
    st.caption(f"{len(view)} of {len(df)} compounds. Click a row to see its structure and references.")
    ev = st.dataframe(view, hide_index=True, column_config=cfg, on_select="rerun", selection_mode="single-row", key="lib_tbl")
    sel = row_selection(ev)
    if sel and sel[0] < len(view):
        cid = view.iloc[sel[0]]["compound_id"]
        r = master[master["compound_id"] == cid].iloc[0]
        with st.container(border=True):
            a, b = st.columns([1, 2])
            with a:
                show_structure(r.get("smiles_canonical"))
            with b:
                st.markdown(f"**{esc(display_name(data, cid))}** ({esc(cid)})")
                show_html("<b>References (DOI):</b> " + doi_links_html(r.get("dois")))
                st.caption(f'Provenance status: {r.get("provenance_status", "-")}')
                st.button("Open compound profile", key="lib_open", on_click=open_profile, args=(cid,))
    section_footer("Step 1 · Chemist · the library is the list of compounds every later step works on.")


# =============================================================================
# SECTION 2: STRUCTURES & VALIDATION
# =============================================================================
def section_structures():
    data, status = get_data()
    sec = SECTION_BY_NAME["2. Structures & validation"]
    step = STEP_BY_ID[2]
    section_header(sec, status)
    rd = get_redock(data)
    if not rd and not usable(status, "receptors") and not usable(status, "overlays"):
        waiting_card(step, status)
        section_footer("Step 2 · Bioinformatics")
        return
    problem_notes(status, ["redocking", "receptors", "overlays"])
    st.info("**Redocking** = remove the known drug, dock it back, and check it lands where the crystal had it. "
            "Pass = the top pose is under 2 Å from the crystal position.")

    rows = []
    for k, p in PROTEINS.items():
        r = rd.get(k, {})
        conf = conf_for(data, k)
        centre = f'{conf["center_x"]:.1f}, {conf["center_y"]:.1f}, {conf["center_z"]:.1f}' if conf else "-"
        size = f'{conf["size_x"]:.0f} × {conf["size_y"]:.0f} × {conf["size_z"]:.0f}' if conf else "-"
        verdict = r.get("verdict")
        ref_txt = fmt(r.get("ref_score"), 3)
        extra = chip("LOW CONFIDENCE") if verdict and str(verdict).upper() == "PARTIAL" else ""
        rows.append(f"<tr><td><b>{esc(p['name'])}</b><br><span style='font-size:13px'>{esc(p['role'])}</span></td>"
                    f"<td>{esc(p['pdb'])}</td><td>{esc(p['chains'])}</td><td>{esc(r.get('ligand') if not is_missing(r.get('ligand')) else p['ligand'])}</td>"
                    f"<td>{esc(centre)}</td><td>{esc(size)}</td><td>{fmt(r.get('top1_rmsd'), 2)}</td><td>{fmt(r.get('best_rmsd'), 2)}</td>"
                    f"<td>{esc(ref_txt)}</td><td>{verdict_chip(verdict)}{extra}</td></tr>")
    head = ["Protein", "PDB ID", "Chains", "Crystal ligand", "Box centre (Å)", "Box size (Å)", "Top-pose RMSD (Å)",
            "Best RMSD (Å)", "Reference score (kcal/mol)", "Result"]
    show_html('<table class="op-tbl"><tr>' + "".join(f"<th>{h}</th>" for h in head) + "</tr>" + "".join(rows) + "</table>")
    st.caption("Box = the search box drawn around the pocket (from the crystal drug). Reference score = the known drug redocked, in kcal/mol, predicted.")
    if str(rd.get("herg", {}).get("verdict", "")).upper() == "PARTIAL":
        h = rd["herg"]
        st.warning(f"hERG redocking was PARTIAL: the top pose is {fmt(h.get('top1_rmsd'), 2)} Å from the crystal position, but the best of the "
                   f"poses is {fmt(h.get('best_rmsd'), 2)} Å. The right pose was found but not ranked first, so **all hERG scores are low confidence**.")
    if status["receptors"]["state"] != "ok":
        st.caption(f"Raw receptor files found: {status['receptors']['rows'] or 0} of 4 (data/receptors/).")

    rm = data.get("redock_rmsd")
    if rm is not None:
        with st.expander("Detailed redocking RMSD file (data/validation)"):
            st.dataframe(rm, hide_index=True)

    st.subheader("Crystal drug vs redocked pose")
    keys = list(PROTEINS)
    pick = st.selectbox("Protein", keys, format_func=lambda k: f"{PROT_LABEL[k]} ({PROTEINS[k]['pdb']})", key="struct_pick")
    path = file_for(data, "overlays", pick)
    if not path:
        show_html(f'<div class="op-wait">Overlay file for <b>{esc(PROT_LABEL[pick])}</b> not delivered '
                  f'(expected <code>{esc(pick)}_..._overlay_complex.pdb</code> in data/incoming/validation/).</div>')
    else:
        txt = read_text_safe(path, 12)
        if not txt:
            st.warning(f"{Path(path).name}: could not be read or is too large.")
        else:
            viewer([{"data": txt, "fmt": "pdb", "styles": [PROTEIN_STYLE, ligand_style(ORANGE_XTL, {"resn": "XTL"}),
                                                           ligand_style(GREEN_POSE, {"resn": "TOP"})]}],
                   legend=[("Protein", "#9aa0a6"), ("Crystal drug", ORANGE_XTL), ("Redocked pose", GREEN_POSE)],
                   dark=is_dark(), zoom_sel={"resn": ["XTL", "TOP"]})
            st.caption("Orange = known drug as found in the crystal. Green = our redocked top pose. The closer they overlap, the better the docking setup.")
    section_footer("Step 2 · Bioinformatics · a docking setup is only trusted after it can put a known drug back.")


# =============================================================================
# SECTION 3: ADMET FILTER
# =============================================================================
ADMET_HELP = {
    "qed": "QED: a 0-1 score of how drug-like a molecule looks (size, greasiness, etc.). Higher is more drug-like.",
    "lipinski_violations": "Lipinski's rule of five: quick size/greasiness rules for drug-likeness. This counts how many rules a molecule breaks.",
    "herg_risk": "Model output (0-1) for predicted risk of blocking the hERG heart channel. Higher = more concern. A prediction, not a test.",
    "hepatotox_risk": "Model output (0-1) for predicted liver-toxicity concern. Higher = more concern. A prediction, not a test.",
    "admet_score": "Combined drug-likeness and predicted ADMET score from the ML Engineer (0-1).",
    "pass_flag": "True if the compound passed the filter and went on to docking.",
}


def section_admet():
    data, status = get_data()
    sec = SECTION_BY_NAME["3. ADMET filter"]
    step = STEP_BY_ID[3]
    section_header(sec, status)
    ml, which = get_ml(data)
    if ml is None:
        waiting_card(step, status)
        section_footer("Step 3 · ML Engineer")
        return
    problem_notes(status, ["ml_detail"] if which == "detail" else ["ml_compact"])
    if which == "compact":
        st.caption("Only the compact ML table is available, so some columns are missing.")
    n, passed = len(ml), int(ml["pass_flag"].fillna(False).sum())
    c1, c2, c3 = st.columns(3)
    with c1:
        metric_card("Library", str(n), "compounds scored")
    with c2:
        metric_card("Passed the filter", f"{passed} of {n}", "go on to docking")
    with c3:
        metric_card("Did not pass", str(n - passed), "see the reasons below")

    if "fail_reason" in ml.columns:
        failed = ml[ml["pass_flag"] == False]  # noqa: E712
        reasons = failed["fail_reason"].dropna().astype(str).str.split(";").explode().str.strip()
        rc = reasons[reasons != ""].value_counts().rename_axis("reason").reset_index(name="compounds")
        if len(rc):
            st.subheader("Why compounds did not pass")
            ch = alt.Chart(rc).mark_bar(color=ROLE_COLORS["ML Engineer"]).encode(
                x=alt.X("compounds:Q"), y=alt.Y("reason:N", sort="-x", title=None), tooltip=["reason", "compounds"]).properties(height=30 * len(rc) + 30)
            st.altair_chart(style_chart(ch))
            multi = int(failed["fail_reason"].astype(str).str.contains(";").sum())
            st.caption(f"{multi} compound(s) had several reasons at once; each reason is counted, so the bars add up to more than {n - passed}.")

    st.subheader("Browse the table")
    f1, f2, f3, f4 = st.columns(4)
    with f1:
        plants = sorted(ml["source_plant"].dropna().unique()) if "source_plant" in ml.columns else []
        plant = st.multiselect("Plant", plants, default=plants, key="adm_plant") if plants else plants
    with f2:
        outcome = st.radio("Outcome", ["All", "Passed", "Did not pass"], horizontal=True, key="adm_out")
    with f3:
        hmax = st.slider("hERG risk up to", 0.0, 1.0, 1.0, 0.01, key="adm_herg") if "herg_risk" in ml.columns else 1.0
    with f4:
        lmax = st.slider("Liver risk up to", 0.0, 1.0, 1.0, 0.01, key="adm_liv") if "hepatotox_risk" in ml.columns else 1.0
    v = ml.copy()
    if plants and plant:
        v = v[v["source_plant"].isin(plant)]
    if outcome == "Passed":
        v = v[v["pass_flag"] == True]  # noqa: E712
    elif outcome == "Did not pass":
        v = v[v["pass_flag"] == False]  # noqa: E712
    if "herg_risk" in v.columns:
        v = v[v["herg_risk"].fillna(0) <= hmax]
    if "hepatotox_risk" in v.columns:
        v = v[v["hepatotox_risk"].fillna(0) <= lmax]
    v = v.copy()
    v.insert(1, "name_shown", v["compound_id"].map(lambda c: display_name(data, c)))
    cols = [c for c in ["compound_id", "name_shown", "qed", "lipinski_violations", "herg_risk", "hepatotox_risk",
                        "admet_score", "pass_flag", "fail_reason"] if c in v.columns]
    show = v[cols].reset_index(drop=True)
    st.caption(f"{len(show)} of {n} compounds shown.")
    cfg = {"compound_id": st.column_config.TextColumn("Compound ID", pinned=True), "name_shown": "Name",
           "qed": st.column_config.NumberColumn("QED", format="%.2f", help=ADMET_HELP["qed"]),
           "lipinski_violations": st.column_config.NumberColumn("Lipinski violations", format="%d", help=ADMET_HELP["lipinski_violations"]),
           "herg_risk": st.column_config.NumberColumn("hERG risk (0-1)", format="%.2f", help=ADMET_HELP["herg_risk"]),
           "hepatotox_risk": st.column_config.NumberColumn("Liver risk (0-1)", format="%.2f", help=ADMET_HELP["hepatotox_risk"]),
           "admet_score": st.column_config.NumberColumn("ADMET score (0-1)", format="%.2f", help=ADMET_HELP["admet_score"]),
           "pass_flag": st.column_config.CheckboxColumn("Passed", help=ADMET_HELP["pass_flag"]), "fail_reason": "Fail reason"}
    st.dataframe(show, hide_index=True, column_config=cfg)
    st.download_button("Download the filtered table (CSV)", show.to_csv(index=False).encode("utf-8"),
                       file_name="admet_filtered.csv", mime="text/csv", key="adm_dl")
    with st.expander("What do these words mean?"):
        for k, t in ADMET_HELP.items():
            st.markdown(f"- {t}")
    section_footer("Step 3 · ML Engineer · risk numbers are model predictions, not lab results.")


# =============================================================================
# SECTION 4: DOCKING: DPP-4
# =============================================================================
def section_docking():
    data, status = get_data()
    sec = SECTION_BY_NAME["4. Docking: DPP-4"]
    step = STEP_BY_ID[4]
    section_header(sec, status)
    df = compounds_frame(data)
    if df is None:
        waiting_card(step, status)
        section_footer("Step 4 · Bioinformatics")
        return
    problem_notes(status, ["docking_target", "poses_target"])
    ml, _ = get_ml(data)
    n_pass = int(ml["pass_flag"].fillna(False).sum()) if ml is not None else None
    n = len(df)
    note = card_value(data, "library.docking_note", None)
    if n_pass:
        extra = f" ({note})" if note and n_pass - n == N_STEREO_NOT_DOCKED else (f" ({n_pass - n} not docked)" if n_pass != n else "")
        st.markdown(f"**{n} of {n_pass} docked**{extra}. Coverage is shown by the app.")
    else:
        st.markdown(f"**{n} compounds docked.**")
    rd = get_redock(data)
    ref = ref_score(rd, "dpp4")
    sr = data.get("summary_report")
    exh = int(sr["exhaustiveness"].dropna().iloc[0]) if sr is not None and "exhaustiveness" in sr.columns and sr["exhaustiveness"].notna().any() else EXHAUSTIVENESS_SCREEN
    top10 = df.head(10)["vina_score_kcal"]
    spread = float(top10.max() - top10.min()) if len(top10) == 10 else None
    msg = f"The screen used a low search effort (exhaustiveness {exh}), so scores within about {NOISE_KCAL} kcal/mol are inside the noise."
    if spread is not None and spread <= NOISE_KCAL:
        msg += f" The top 10 span only {spread:.2f} kcal/mol, so they are effectively tied."
    st.info(msg)
    if ref is not None and float(df["vina_score_kcal"].min()) > ref:
        st.caption(f"No compound scores better than the known drug 605 redocked ({ref:.3f} kcal/mol).")

    st.subheader("How the scores are spread")
    base = alt.Chart(df).mark_bar(color=ROLE_COLORS["Bioinformatics"]).encode(
        x=alt.X("vina_score_kcal:Q", bin=alt.Bin(step=0.25), title="Vina score (kcal/mol, predicted). More negative = better fit"),
        y=alt.Y("count():Q", title="compounds"))
    layers = base
    if ref is not None:
        rule = alt.Chart(pd.DataFrame({"x": [ref], "label": [f"Known drug 605: {ref:.3f}"]})).mark_rule(color="#c0392b", strokeWidth=3).encode(x="x:Q")
        txt = alt.Chart(pd.DataFrame({"x": [ref], "label": [f"Known drug 605: {ref:.3f}"]})).mark_text(
            align="right", dx=-6, dy=-6, baseline="top", fontSize=13, color="#c0392b").encode(x="x:Q", y=alt.value(4), text="label:N")
        layers = base + rule + txt
    else:
        st.caption("Reference drug score is pending, so no reference line is drawn.")
    st.altair_chart(style_chart(layers.properties(height=260)))

    st.subheader("All docked compounds")
    q = st.text_input("Search by name or ID", key="dock_search")
    view = df
    if q:
        blob = view[["compound_id", "name"]].astype(str).agg(" ".join, axis=1)
        view = view[blob.str.contains(re.escape(q), case=False, na=False)]
    cols = [c for c in ["rank", "compound_id", "name", "vina_score_kcal", "ligand_efficiency", "heavy_atoms", "qed",
                        "admet_score", "herg_risk", "hepatotox_risk", "pass_flag"] if c in view.columns]
    show = view[cols].reset_index(drop=True)
    fmt_map = {c: "{:.2f}" for c in ("vina_score_kcal", "ligand_efficiency", "qed", "admet_score", "herg_risk", "hepatotox_risk") if c in show.columns}
    styler = show.style.format(fmt_map, na_rep="-").map(score_background(show, "vina_score_kcal"), subset=["vina_score_kcal"])
    cfg = {"rank": st.column_config.NumberColumn("Rank", help="Shown by the app: 1 = most negative score"),
           "compound_id": st.column_config.TextColumn("Compound ID", pinned=True), "name": "Name",
           "vina_score_kcal": st.column_config.NumberColumn("Vina score (kcal/mol)", help="Predicted docking score. More negative = better predicted fit.", pinned=True),
           "ligand_efficiency": st.column_config.NumberColumn("Ligand efficiency", help="Score per heavy atom, shown as a positive number. Higher = more efficient."),
           "heavy_atoms": st.column_config.NumberColumn("Heavy atoms", help="Atoms other than hydrogen"),
           "qed": st.column_config.NumberColumn("QED", help=ADMET_HELP["qed"]),
           "admet_score": st.column_config.NumberColumn("ADMET score (0-1)", help=ADMET_HELP["admet_score"]),
           "herg_risk": st.column_config.NumberColumn("hERG risk (0-1)", help=ADMET_HELP["herg_risk"]),
           "hepatotox_risk": st.column_config.NumberColumn("Liver risk (0-1)", help=ADMET_HELP["hepatotox_risk"]),
           "pass_flag": st.column_config.CheckboxColumn("Passed ML filter")}
    st.caption(f"{len(show)} of {n} compounds. Click a column header to sort; click a row to open its profile.")
    ev = st.dataframe(styler, hide_index=True, column_config=cfg, on_select="rerun", selection_mode="single-row", key="dock_tbl")
    sel = row_selection(ev)
    if sel and sel[0] < len(show):
        cid = show.iloc[sel[0]]["compound_id"]
        st.button(f"Open profile of {cid}", key="dock_open", on_click=open_profile, args=(cid,))
    st.download_button("Download the filtered table (CSV)", show.to_csv(index=False).encode("utf-8"),
                       file_name="docking_dpp4_filtered.csv", mime="text/csv", key="dock_dl")
    section_footer("Step 4 · Bioinformatics · Vina scores are predictions, in kcal/mol; lower is better.")


# =============================================================================
# SECTION 5: SAFETY PANEL
# =============================================================================
def section_safety():
    data, status = get_data()
    sec = SECTION_BY_NAME["5. Safety panel"]
    step = STEP_BY_ID[5]
    section_header(sec, status)
    mat = safety_matrix(data)
    if not mat:
        waiting_card(step, status)
        section_footer("Step 5 · Bioinformatics")
        return
    problem_notes(status, ["docking_antitargets"])
    if status["docking_antitargets"]["state"] in ("missing", "invalid"):
        st.warning("docking_antitargets.csv is not usable, so only the safety scores listed in summary_report.csv (top 10) are shown.")
    rd = get_redock(data)
    dt = data.get("docking_target")
    dpp = dict(zip(dt["compound_id"], dt["vina_score_kcal"])) if dt is not None else {}
    sr = data.get("summary_report")
    if sr is not None:
        for cid, sc in zip(sr["compound_id"], sr["vina_score_kcal"]):
            dpp.setdefault(cid, sc)
    order = sorted(mat, key=lambda c: (dpp.get(c, 0) if not is_missing(dpp.get(c)) else 99, c))

    # redocking chips per protein
    cols = st.columns(3)
    for col, k in zip(cols, SAFETY_PROTEINS):
        r = rd.get(k, {})
        ref = ref_score(rd, k)
        v = r.get("verdict")
        chips = verdict_chip(v)
        if v and str(v).upper() != "PASS":
            chips += chip("LOW CONFIDENCE")
        if ref is None:
            chips += chip("BENCHMARK PENDING")
        with col:
            show_html(f'<div class="op-card"><b>{esc(PROT_LABEL[k])}</b> <span style="font-size:13px">(PDB {PROTEINS[k]["pdb"]})</span><br>'
                      f'Redocking: {chips}<br><span style="font-size:13px">Known drug {esc(PROTEINS[k]["ligand"])} redocked: '
                      f'{esc(fmt(ref, 3, "kcal/mol"))}</span></div>')
    st.info("**Flagged** = docks as strongly as the known drug in this protein. It is a prompt to check, not a verdict.")
    st.warning("Scores are not comparable between proteins (big greasy pockets give strong scores to almost anything), "
               "so each protein is compared only with its own known drug.")

    # heatmap data
    records = []
    for cid in order:
        label = f"{cid} · {short(display_name(data, cid), 22)}"
        records.append({"compound": label, "protein": "DPP-4", "score": dpp.get(cid), "state": "target", "delta": np.nan})
        for k in SAFETY_PROTEINS:
            sc = mat[cid].get(k)
            ref = ref_score(rd, k)
            stt = safety_status(sc, ref)
            state = {"NOT TESTED": "missing", "NO VALID POSE": "invalid"}.get(stt, "valid")
            records.append({"compound": label, "protein": PROT_LABEL[k], "score": sc, "state": state,
                            "delta": (float(sc) - ref) if state == "valid" and ref is not None else np.nan,
                            "flag": stt == "FLAGGED"})
    refrow = "Known drug (redocked)"
    for k in ["dpp4"] + SAFETY_PROTEINS:
        records.append({"compound": refrow, "protein": PROT_LABEL[k], "score": ref_score(rd, k), "state": "reference", "delta": np.nan})
    hm = pd.DataFrame(records)
    hm["flag"] = hm.get("flag", False)
    hm["flag"] = hm["flag"].fillna(False).astype(bool)

    def label(r):
        if r["state"] == "invalid":
            return "no valid pose"
        if is_missing(r["score"]):
            return "-"
        return f"{r['score']:.2f}" + (" ⚑" if r["flag"] else "")
    hm["label"] = hm.apply(label, axis=1)
    y_order = order_labels = list(dict.fromkeys(hm["compound"]))
    x_order = ["DPP-4"] + [PROT_LABEL[k] for k in SAFETY_PROTEINS]
    enc = dict(x=alt.X("protein:N", sort=x_order, title=None, axis=alt.Axis(orient="top", labelAngle=0)),
               y=alt.Y("compound:N", sort=y_order, title=None))
    base = alt.Chart(hm)
    fill_valid = base.transform_filter(alt.datum.state == "valid").mark_rect(stroke="#8a939e", strokeWidth=0.5).encode(
        **enc, color=alt.Color("delta:Q", scale=alt.Scale(domain=[-3, 0, 3], range=["#d9822b", "#f5e6c8", "#e9edf2"], clamp=True),
                               legend=alt.Legend(title="score minus own known drug (kcal/mol)")),
        tooltip=["compound", "protein", alt.Tooltip("score:Q", format=".3f"), alt.Tooltip("delta:Q", format=".2f")])
    fill_grey = base.transform_filter(alt.datum.state == "invalid").mark_rect(color="#8a939e").encode(**enc, tooltip=["compound", "protein"])
    fill_target = base.transform_filter(alt.datum.state == "target").mark_rect(color="#dce6f2").encode(**enc)
    fill_ref = base.transform_filter(alt.datum.state == "reference").mark_rect(color="#2f6fd0", opacity=0.35).encode(**enc)
    flag_box = base.transform_filter(alt.datum.flag == True).mark_rect(fill=None, stroke="#c0392b", strokeWidth=3).encode(**enc)  # noqa: E712
    text_dark = base.transform_filter("datum.state != 'invalid'").mark_text(fontSize=13, color="#1d2430").encode(**enc, text="label:N")
    text_light = base.transform_filter(alt.datum.state == "invalid").mark_text(fontSize=13, color="#ffffff").encode(**enc, text="label:N")
    chart = (fill_target + fill_ref + fill_valid + fill_grey + flag_box + text_dark + text_light).properties(
        width=140 * len(x_order), height=30 * len(y_order), title="Vina score (kcal/mol, predicted) per protein")
    st.altair_chart(style_chart(chart))
    st.caption("Amber = this protein's score is as strong as or stronger than its own known drug. ⚑ + red outline = flagged. "
               "Grey = no valid pose (score above -2.0 kcal/mol), left out of the colour scale. Blue row = the known drug redocked.")

    # flag table
    st.subheader("Flags per protein")
    rows = []
    for k in SAFETY_PROTEINS:
        ref = ref_score(rd, k)
        for cid in order:
            sc = mat[cid].get(k)
            rows.append({"Protein": PROT_LABEL[k], "Compound ID": cid, "Name": display_name(data, cid),
                         "Vina score (kcal/mol)": None if is_missing(sc) else float(sc),
                         "Known drug (kcal/mol)": ref,
                         "Result": safety_status(sc, ref).capitalize() if safety_status(sc, ref) != "FLAGGED" else "Flagged",
                         "Raw difference vs DPP-4 (not comparable across proteins)":
                             None if is_missing(sc) or is_missing(dpp.get(cid)) or not is_valid_pose(sc) else round(float(sc) - float(dpp[cid]), 2)})
    ft = pd.DataFrame(rows)
    for k in SAFETY_PROTEINS:
        sub = ft[ft["Protein"] == PROT_LABEL[k]]
        valid = sub[sub["Result"].isin(["Flagged", "Not flagged"])]
        flagged = sub[sub["Result"] == "Flagged"]["Compound ID"].tolist()
        closest = valid.sort_values("Vina score (kcal/mol)").head(1)
        low = " · low confidence (redocking " + str(rd.get(k, {}).get("verdict", "?")).lower() + ")" if str(rd.get(k, {}).get("verdict", "PASS")).upper() != "PASS" else ""
        cl = f"; closest: {closest.iloc[0]['Compound ID']} at {closest.iloc[0]['Vina score (kcal/mol)']:.3f} kcal/mol" if len(closest) else ""
        st.markdown(f"**{PROT_LABEL[k]}**: {len(flagged)} flagged" + (f" ({', '.join(flagged)})" if flagged else "") +
                    f" · {len(valid)} of {len(sub)} with a valid pose{cl}{low}")
    pick = st.multiselect("Show proteins", [PROT_LABEL[k] for k in SAFETY_PROTEINS], default=[PROT_LABEL[k] for k in SAFETY_PROTEINS], key="saf_prot")
    st.dataframe(ft[ft["Protein"].isin(pick)], hide_index=True,
                 column_config={"Vina score (kcal/mol)": st.column_config.NumberColumn(format="%.2f"),
                                "Known drug (kcal/mol)": st.column_config.NumberColumn(format="%.2f")})
    profile_button(order, "saf")
    section_footer("Step 5 · Bioinformatics · flags are prompts to check, not verdicts; nothing here says a compound is safe or harmful.")


# ===== PART 3 GOES HERE =====  (sections 6-9: section_ranking, section_md, section_profile, section_evidence)
# =============================================================================
# SECTION 6: FINAL RANKING
# =============================================================================
COMPOSITE_FORMULA = "Composite = 0.5 × efficacy_norm + 0.3 × selectivity_norm + 0.2 × ADMET score"
SCORE_NAME = "Selectivity-adjusted efficacy score (in silico proxy)"


def section_ranking():
    data, status = get_data()
    sec = SECTION_BY_NAME["6. Final ranking"]
    step = STEP_BY_ID[6]
    section_header(sec, status)
    rk = data.get("ranking")
    if rk is None:
        show_html('<div class="op-wait"><b>Pending ML ranking</b> &nbsp;|&nbsp; from: <b>ML Engineer</b> &nbsp;|&nbsp; expected: <b>Day 5</b> '
                  '&nbsp;|&nbsp; Drop it in <code>data/incoming/</code> and refresh.</div>')
        for p in status["ranking"]["problems"]:
            st.warning(p)
        section_footer("Step 6 · ML Engineer · this app never computes a composite score itself.")
        return
    problem_notes(status, ["ranking"])
    st.markdown(f"**{SCORE_NAME}**")
    st.code(COMPOSITE_FORMULA, language=None)
    st.caption("Shown exactly as the ML Engineer delivered it. It is a proxy built from predictions, not a therapeutic index.")

    top = rk.iloc[0]
    cid = top["compound_id"]
    dt = data.get("docking_target")
    dpp = dict(zip(dt["compound_id"], dt["vina_score_kcal"])) if dt is not None else {}
    dpp_rank = dict(zip(dt["compound_id"], dt["rank"])) if dt is not None else {}
    mat = safety_matrix(data)
    rd = get_redock(data)

    msg = (f"**#1 in the final ranking: {display_name(data, cid)} ({cid})**, score {float(top['composite_score']):.3f}.")
    if dt is not None:
        best = dt.iloc[0]["compound_id"]
        if best != cid:
            br = rk[rk["compound_id"] == best]
            where = f"#{int(br.iloc[0]['rank'])} in the final ranking ({float(br.iloc[0]['composite_score']):.3f})" if len(br) else "not in the final ranking"
            msg += (f" Note: the best DPP-4 docking score belongs to **{best}** ({float(dt.iloc[0]['vina_score_kcal']):.3f} kcal/mol), "
                    f"which is {where}. The ranking also weighs the safety-panel scores and drug-likeness.")
    st.success(msg)

    # enriched table
    t = rk.copy()
    t.insert(2, "name", t["compound_id"].map(lambda c: display_name(data, c)))
    t.insert(3, "dpp4_score", t["compound_id"].map(lambda c: dpp.get(c)))
    notes = []
    for c in t["compound_id"]:
        bad = [PROT_LABEL[k] for k in SAFETY_PROTEINS if k in mat.get(c, {}) and safety_status(mat[c][k], ref_score(rd, k)) == "NO VALID POSE"]
        notes.append("no valid pose: " + ", ".join(bad) if bad else "")
    t["safety_note"] = notes
    cfg = {"rank": st.column_config.NumberColumn("Rank", pinned=True),
           "compound_id": st.column_config.TextColumn("Compound ID", pinned=True), "name": "Name",
           "dpp4_score": st.column_config.NumberColumn("DPP-4 score (kcal/mol)", format="%.2f", help="Predicted docking score; more negative = better fit."),
           "efficacy_norm": st.column_config.NumberColumn("Efficacy (0-1)", format="%.2f", help="Normalised DPP-4 docking performance, as delivered."),
           "selectivity_norm": st.column_config.NumberColumn("Selectivity (0-1)", format="%.2f", help="Normalised computational selectivity assessment, as delivered."),
           "admet_score": st.column_config.NumberColumn("ADMET score (0-1)", format="%.2f"),
           "composite_score": st.column_config.NumberColumn("Composite score (0-1)", format="%.3f", help=COMPOSITE_FORMULA),
           "safety_note": "Safety note"}
    if not any(notes):
        t = t.drop(columns=["safety_note"])
    ev = st.dataframe(t, hide_index=True, column_config=cfg, on_select="rerun", selection_mode="single-row", key="rank_tbl")
    sel = row_selection(ev)
    if sel and sel[0] < len(t):
        c2 = t.iloc[sel[0]]["compound_id"]
        st.button(f"Open profile of {c2}", key="rank_open", on_click=open_profile, args=(c2,))
    if any(notes):
        show_html(chip("NO VALID POSE") + " One or more safety scores for this compound were not real poses (score above -2.0 kcal/mol). "
                  "Treat its safety-panel numbers as missing.")

    st.subheader("Composite score by compound")
    ch = alt.Chart(t.assign(label=t["compound_id"] + " · " + t["name"].map(lambda s: short(s, 24)))).mark_bar(color=ROLE_COLORS["ML Engineer"]).encode(
        x=alt.X("composite_score:Q", title="composite score (0-1, higher is ranked higher)", scale=alt.Scale(domain=[0, 1])),
        y=alt.Y("label:N", sort=alt.EncodingSortField("rank", order="ascending"), title=None),
        tooltip=["compound_id", alt.Tooltip("composite_score:Q", format=".3f"), "rank"]).properties(height=30 * len(t) + 30)
    st.altair_chart(style_chart(ch))

    sel_note = card_value(data, "composite_scoring.definitions_checked_against_files.selectivity_norm", None)
    with st.expander("How was selectivity built? (read before quoting this ranking)"):
        st.markdown("The selectivity part is built from raw differences between a safety-protein score and the DPP-4 score. "
                    "Vina scores are **not comparable between different proteins**, so treat this part with care and prefer the "
                    "per-protein comparison with each known drug in 5. Safety panel.")
        if sel_note:
            st.caption(f"Note from the target card: {sel_note}")
    section_footer("Step 6 · ML Engineer · a ranking of predictions, to choose what to test in the lab first.")


# =============================================================================
# SECTION 7: SIMULATION (MD)
# =============================================================================
def section_md():
    data, status = get_data()
    sec = SECTION_BY_NAME["7. Simulation (MD)"]
    step = STEP_BY_ID[7]
    section_header(sec, status)
    md = data.get("md_summary")
    md_dir = ROOT / "data" / "incoming" / "md_report"
    if md is None:
        waiting_card(step, status)
        st.caption("Expected columns: compound_id, label (candidate|control), length_ns, ligand_rmsd_mean_A, ligand_rmsd_sd_A, "
                   "protein_rmsd_mean_A, hbond_occupancy_pct, notes. Plots: <compound_id>_rmsd.png, _rmsf.png, _hbonds.png. "
                   "Optional: <compound_id>_hbonds_detail.csv.")
    else:
        problem_notes(status, ["md_summary"])
        st.markdown("Molecular dynamics results from supplied files; no simulation is run here. "
                    "Candidate and control (the known crystal drug) side by side. Lower RMSD = steadier in the pocket.")
        cols = [c for c in ["compound_id", "label", "length_ns", "ligand_rmsd_mean_A", "ligand_rmsd_sd_A",
                            "protein_rmsd_mean_A", "hbond_occupancy_pct", "notes"] if c in md.columns]
        cfg = {"compound_id": "Compound ID", "label": "Role", "length_ns": st.column_config.NumberColumn("Length (ns)", format="%.1f"),
               "ligand_rmsd_mean_A": st.column_config.NumberColumn("Ligand RMSD mean (Å)", format="%.2f"),
               "ligand_rmsd_sd_A": st.column_config.NumberColumn("Ligand RMSD sd (Å)", format="%.2f"),
               "protein_rmsd_mean_A": st.column_config.NumberColumn("Protein RMSD mean (Å)", format="%.2f"),
               "hbond_occupancy_pct": st.column_config.NumberColumn("H-bond occupancy (%)", format="%.1f"), "notes": "Notes"}
        st.dataframe(md[cols], hide_index=True, column_config=cfg)
        st.caption("Ligand RMSD is measured against the first production frame after fitting on the protein backbone. "
                   "A flat, low curve (like the crystal-drug control) means the pose held; a rising curve means it moved. "
                   "H-bond occupancy = share of frames with at least one ligand-protein hydrogen bond.")

        ids = list(dict.fromkeys(md["compound_id"].astype(str)))
        control = [str(c) for c, lab in zip(md["compound_id"], md["label"]) if str(lab).strip().lower() == "control"]
        default = (control[:1] + [i for i in ids if i not in control])[:3]
        picked = st.multiselect("Compounds to compare side by side", ids, default=default, key="md_compare")
        plots = [("rmsd", "RMSD over time"), ("hbonds", "Ligand-protein H-bonds"), ("rmsf", "Protein flexibility (C-alpha RMSF)")]
        if picked:
            for col, cid in zip(st.columns(len(picked)), picked):
                row = md[md["compound_id"].astype(str) == cid].iloc[0]
                with col:
                    st.markdown(f"**{cid}** · {row['label']}")
                    for kind, title in plots:
                        p = file_for(data, "md_plots", f"{cid}_{kind}.png")
                        path = safe_path(Path(p).name, ROOT / Path(p).parent) if p else None
                        if path:
                            st.image(str(path), caption=title)
                        else:
                            st.caption(f"{title}: plot not delivered")
                    detail = safe_path(f"{cid}_hbonds_detail.csv", md_dir)
                    if detail is None and md_dir.is_dir():  # compound ids are upper-cased on load; match file names ignoring case
                        want = f"{cid}_hbonds_detail.csv".lower()
                        detail = next((q for q in md_dir.glob("*_hbonds_detail.csv") if q.name.lower() == want), None)
                    if detail is not None:
                        with st.expander("H-bond partners"):
                            try:
                                st.dataframe(pd.read_csv(detail), hide_index=True)
                            except (OSError, ValueError, pd.errors.ParserError) as exc:
                                st.caption(f"Could not read {detail.name}: {exc}")
        st.caption("A simulation of one compound is a plausibility check on its pose, not proof of binding. "
                   "Report previews do not validate simulation workflow or interpret missing values as results.")
    st.button("Run simulation (engine not connected)", disabled=True, key="md_run")
    section_footer("Step 7 · Bioinformatics · the simulation engine connects in a later version (see md_engine_run in app.py).")


# =============================================================================
# SECTION 8: COMPOUND PROFILE
# =============================================================================
def _chips(items, color):
    return "".join(chip(i, color=color) for i in items) if items else '<span style="font-size:13px">none listed</span>'


def _pose_layers(data, cid: str, pose_n: int, show_xtl: bool):
    """Build the 3D viewer layers. Returns (layers, legend, around, zoom_sel, message)."""
    legend = [("Protein", "#9aa0a6"), ("This compound", GREEN_POSE)]
    layers, around, zoom_sel, msg = [], None, None, ""
    plip = file_for(data, "plip", f"{cid}_")
    pose_rel = None
    dt = data.get("docking_target")
    if dt is not None and "pose_file" in dt.columns:
        pf = dt.loc[dt["compound_id"] == cid, "pose_file"]
        if len(pf) and not is_missing(pf.iloc[0]):
            p = safe_path(pf.iloc[0])
            pose_rel = rel(p) if p else None
    if not pose_rel:
        pose_rel = file_for(data, "poses_target", f"{cid}_")
    if pose_n == 1 and plip and (txt := read_text_safe(plip, 12)):
        layers.append({"data": txt, "fmt": "pdb", "styles": [PROTEIN_STYLE, ligand_style(GREEN_POSE, {"resn": "LIG"})]})
        around, zoom_sel = {"layer": 0, "distance": 4, "sel": {"resn": "LIG"}}, {"resn": "LIG"}
    else:
        rec = file_for(data, "receptors", "3d4l")
        rtxt = read_text_safe(rec, 12) if rec else None
        ptxt = read_text_safe(pose_rel, 12) if pose_rel else None
        if not pose_rel or not ptxt:
            return [], legend, None, None, "pose file not delivered"
        models = split_pdbqt_models(ptxt)
        if pose_n > len(models):
            return [], legend, None, None, f"this pose file has only {len(models)} model(s)"
        if not rtxt:
            return [], legend, None, None, "protein structure (data/receptors/3D4L.pdb) not delivered"
        layers.append({"data": rtxt, "fmt": "pdb", "styles": [PROTEIN_STYLE]})
        layers.append({"data": models[pose_n - 1][0], "fmt": "pdbqt", "styles": [ligand_style(GREEN_POSE)]})
        around = {"layer": 1, "distance": 4, "sel": {}}
    if show_xtl:
        ov = file_for(data, "overlays", "dpp4")
        otxt = read_text_safe(ov, 12) if ov else None
        xtl = pdb_residue_lines(otxt, "XTL") if otxt else ""
        if xtl:
            layers.append({"data": xtl, "fmt": "pdb", "styles": [ligand_style(ORANGE_XTL)]})
            legend.append(("Known drug 605 (crystal)", ORANGE_XTL))
        else:
            msg = "known drug 605 file not available (needs the DPP-4 overlay file)"
    return layers, legend, around, zoom_sel, msg


def section_profile():
    data, status = get_data()
    sec = SECTION_BY_NAME["8. Compound profile"]
    step = STEP_BY_ID[4]
    section_header(sec, status)
    dt = data.get("docking_target")
    master, ml = data.get("master"), get_ml(data)[0]
    pool = list(dt["compound_id"]) if dt is not None else (list(master["compound_id"]) if master is not None else (list(ml["compound_id"]) if ml is not None else []))
    if not pool:
        waiting_card(step, status, ["docking_target"])
        section_footer("Profile · everyone")
        return
    cur = st.session_state.get("profile_cid")
    if cur and cur not in pool:
        pool = pool + [cur]
    if cur not in pool:
        st.session_state["profile_cid"] = pool[0]            # default = rank 1
    cid = st.selectbox("Compound", pool, key="profile_cid",
                       format_func=lambda c: f"{c} · {short(display_name(data, c), 40)}")
    problem_notes(status, ["docking_target", "summary_report"])

    mrow = master[master["compound_id"] == cid].iloc[0] if master is not None and (master["compound_id"] == cid).any() else None
    drow = dt[dt["compound_id"] == cid].iloc[0] if dt is not None and (dt["compound_id"] == cid).any() else None
    mlrow = ml[ml["compound_id"] == cid].iloc[0] if ml is not None and (ml["compound_id"] == cid).any() else None
    srow = data["summary_report"][data["summary_report"]["compound_id"] == cid] if data.get("summary_report") is not None else None
    srow = srow.iloc[0] if srow is not None and len(srow) else None

    a, b = st.columns([1, 2])
    with a:
        show_structure(mrow.get("smiles_canonical") if mrow is not None else None)
    with b:
        st.markdown(f"### {esc(display_name(data, cid))}")
        plant = mlrow.get("source_plant") if mlrow is not None else None
        kv = [("Compound ID", cid), ("Plant", plant if not is_missing(plant) else "-"),
              ("Formula", mrow.get("molecular_formula") if mrow is not None else "-"),
              ("Chemical class", mrow.get("chemical_class") if mrow is not None else "-"),
              ("NP-Classifier class", mrow.get("np_classifier_class") if mrow is not None else "-"),
              ("Provenance status", mrow.get("provenance_status") if mrow is not None else "-")]
        show_html('<table class="op-kv">' + "".join(f"<tr><td>{esc(k)}</td><td>{esc('-' if is_missing(v) else v)}</td></tr>" for k, v in kv) + "</table>")
        show_html("<b>References (DOI):</b> " + doi_links_html(mrow.get("dois") if mrow is not None else None))

    # docking + ranking + ML
    st.subheader("Scores")
    c1, c2, c3 = st.columns(3)
    rd = get_redock(data)
    with c1:
        if drow is None:
            metric_card("DPP-4 score", None)
            st.caption("This compound was not docked into DPP-4.")
        else:
            metric_card("DPP-4 score (kcal/mol, predicted)", f"{float(drow['vina_score_kcal']):.3f}",
                        f"rank {int(drow['rank'])} of {len(dt)} (shown by the app) · ligand efficiency {fmt(drow.get('ligand_efficiency'), 3)}")
    with c2:
        rk = data.get("ranking")
        r = rk[rk["compound_id"] == cid] if rk is not None else None
        if r is not None and len(r):
            metric_card(SCORE_NAME, f"{float(r.iloc[0]['composite_score']):.3f}", f"final rank {int(r.iloc[0]['rank'])} of {len(rk)}")
        else:
            metric_card("Final ranking", None, "not in the final ranking" if rk is not None else "")
    with c3:
        if mlrow is None:
            metric_card("ML filter", None)
        else:
            metric_card("ML filter", "Passed" if mlrow.get("pass_flag") is True else "Did not pass",
                        f"QED {fmt(mlrow.get('qed'), 2)} · ADMET {fmt(mlrow.get('admet_score'), 2)} · Lipinski violations {fmt(mlrow.get('lipinski_violations'), 0)}")
    if mlrow is not None:
        st.caption(f"Predicted risk (0-1): hERG {fmt(mlrow.get('herg_risk'), 2)} · liver {fmt(mlrow.get('hepatotox_risk'), 2)}"
                   + (f" · fail reason: {mlrow.get('fail_reason')}" if not is_missing(mlrow.get("fail_reason")) else ""))

    st.subheader("Safety panel")
    mat = safety_matrix(data).get(cid, {})
    cols = st.columns(3)
    for col, k in zip(cols, SAFETY_PROTEINS):
        sc, ref = mat.get(k), ref_score(rd, k)
        stt = safety_status(sc, ref)
        low = chip("LOW CONFIDENCE") if str(rd.get(k, {}).get("verdict", "PASS")).upper() != "PASS" else ""
        with col:
            show_html(f'<div class="op-card"><b>{esc(PROT_LABEL[k])}</b><br><span style="font-size:20px;font-weight:700">'
                      f'{esc(fmt(sc, 3) if not is_missing(sc) and is_valid_pose(sc) else "-")}</span> '
                      f'<span style="font-size:13px">kcal/mol</span><br>{chip(stt)}{low}<br>'
                      f'<span style="font-size:13px">own known drug: {esc(fmt(ref, 3))}</span></div>')
    if not mat:
        st.caption("This compound was not sent to the safety panel (only the top 10 were).")

    st.subheader("Interactions with DPP-4")
    if srow is None:
        st.info("Interaction analysis was only done for the top 10 compounds.")
    else:
        for label, col, color in (("Hydrogen bonds", "hbond_residues", "#1F5FA8"), ("Hydrophobic", "hydrophobic_residues", "#5F6B73"),
                                  ("Salt bridges", "salt_bridge_residues", "#D9822B"), ("Pi interactions", "pi_residues", "#2E7D4F")):
            if col in srow.index:
                items = [i.replace("(A)", "").strip() for i in split_list(srow[col])]
                show_html(f"<b>{label}:</b> " + _chips(items, color))
        st.caption("Residue names are shown as in the report (e.g. SER209 = serine at position 209).")

    # 3D viewer
    st.subheader("3D view")
    pose_options = [1, 2, 3, 4, 5]
    ptxt = None
    pf = file_for(data, "poses_target", f"{cid}_")
    if pf:
        ptxt = read_text_safe(pf, 8)
        if ptxt:
            n = len(split_pdbqt_models(ptxt))
            pose_options = list(range(1, min(max(n, 1), 5) + 1))
    v1, v2 = st.columns([1, 1])
    with v1:
        pose_n = st.radio("Pose (1 = best)", pose_options, horizontal=True, key=f"pose_{cid}")
    with v2:
        show_xtl = st.toggle("Show known drug 605 (orange)", value=True, key="show_xtl")
    layers, legend, around, zoom_sel, msg = _pose_layers(data, cid, pose_n, show_xtl)
    if not layers:
        show_html(f'<div class="op-wait">3D view unavailable: <b>{esc(msg or "pose file not delivered")}</b>.</div>')
    else:
        viewer(layers, legend, dark=is_dark(), around=around, zoom_sel=zoom_sel)
        st.caption("Grey = DPP-4 protein, green = this compound, orange = known drug 605. Labelled residues are within 4 Å of the compound.")
        if msg:
            st.caption(msg)
    if pf and ptxt:
        scores = [s for _, s in split_pdbqt_models(ptxt)]
        if any(s is not None for s in scores):
            st.caption("Vina score per pose (kcal/mol): " + ", ".join(f"{i}: {s:.2f}" for i, s in enumerate(scores[:5], 1) if s is not None))
    section_footer("Profile · everything above is predicted, not measured.")


# =============================================================================
# SECTION 9: EVIDENCE & METHODS
# =============================================================================
GLOSSARY = {
    "Protein": "A tiny machine in the body. Many diseases involve a protein that is working wrongly.",
    "Target": "The protein we want a medicine to block or switch off (here DPP-4).",
    "Anti-target (safety protein)": "A protein we do not want a medicine to hit, because that can cause side effects. Ours are hERG, CYP3A4 and CYP2D6.",
    "hERG": "A channel in heart cells. Blocking it can cause dangerous heart rhythm problems.",
    "CYP450 (CYP3A4, CYP2D6)": "Liver enzymes that break down most medicines. Blocking them can cause liver trouble or harmful interactions with other drugs.",
    "Ligand": "A chemical molecule once it is being fitted into a protein.",
    "PDB ID": "A 4-character code for a 3D protein structure in the free public PDB library (rcsb.org), for example 3D4L.",
    "SMILES": "A way of writing a chemical structure as one line of text. For example, CCO is ethanol.",
    "Docking": "The computer tries to fit a chemical into a protein pocket and scores the fit (AutoDock Vina).",
    "Vina score": "The docking score in kcal/mol. More negative is better (a stronger predicted fit). A prediction, not a lab measurement.",
    "Pose": "The position and angle in which the chemical sits in the pocket.",
    "Redocking": "Take the known drug out of its protein, dock it back in, and check the computer puts it where it was found.",
    "RMSD": "How far apart two positions are, in Å (angstrom, a tiny unit of length). Under 2 Å means close enough.",
    "ADMET": "How the body Absorbs, Distributes, Metabolizes and Excretes a chemical, and how Toxic it is.",
    "QED": "A 0-1 score of how drug-like a molecule looks.",
    "Lipinski's rule": "Quick size and greasiness rules for drug-likeness. We count how many a molecule breaks.",
    "PAINS": "Compounds that look active in many tests for the wrong reasons. We flag them.",
    "Ligand efficiency": "Docking score per heavy atom, shown as a positive number (higher = more efficient), so small compounds are not punished.",
    "MD (molecular dynamics)": "A computer simulation of how atoms move over time, used to see whether the best compound stays in the pocket (GROMACS is the software).",
    "Composite score": "One number that combines efficacy, selectivity and drug-likeness, used to rank compounds. Here: a selectivity-adjusted efficacy score (in silico proxy), not a therapeutic index.",
}


def section_evidence():
    data, status = get_data()
    sec = SECTION_BY_NAME["9. Evidence & methods"]
    section_header(sec, status)
    st.info("**Predictions, not lab measurements.** Docking scores and risk numbers are computer predictions. "
            "They help decide what to test first; they do not show that any compound works or is safe.")

    rd = get_redock(data)
    red_txt = ", ".join(f"{PROT_LABEL[k]} {rd[k]['verdict']}" for k in PROTEINS if k in rd and rd[k].get("verdict")) or "redocking results pending"
    criteria = [
        ("Clearly defined problem", "Natural products are usually ranked by potency alone and safety is looked at late. OmniPhyto adds an early safety check "
                                    "for neem compounds against DPP-4, a type 2 diabetes target.", "1. Target & library"),
        ("Computational / data-driven approach", "Open data (COCONUT, PDB), an ML drug-likeness and ADMET filter, AutoDock Vina docking into DPP-4 "
                                                 "and three safety proteins.", "3. ADMET filter"),
        ("Working prototype", "This app: every step of the pipeline appears as a section and fills in when its file is added.", "10. Data inbox"),
        ("Evidence supporting the approach", f"Redocking of the known drug in each protein: {red_txt}. The limits of that evidence are listed below.", "2. Structures & validation"),
        ("Clear documentation", "Methods, limits and a plain-English glossary on this page, plus the team guide.", "9. Evidence & methods"),
        ("Real-world relevance", "A short list of neem compounds to test next in the lab (DPP-4 enzyme assay, then selectivity and toxicity checks).", "6. Final ranking"),
    ]
    st.subheader("The six judging criteria")
    for i, (crit, shows, target) in enumerate(criteria):
        c1, c2, c3 = st.columns([2, 5, 2])
        c1.markdown(f"**{crit}**")
        c2.write(shows)
        if target != "9. Evidence & methods":
            c3.button(f"See: {target}", key=f"crit_{i}", on_click=go, args=(target,))
        else:
            c3.caption("You are here")

    st.subheader("Limits you should know about")
    n_pass = int(get_ml(data)[0]["pass_flag"].fillna(False).sum()) if get_ml(data)[0] is not None else None
    n_dock = len(data["docking_target"]) if data.get("docking_target") is not None else None
    stereo = f"{n_pass - n_dock} stereoisomer duplicates of docked structures were not docked." if n_pass and n_dock and n_pass > n_dock else "Stereoisomer duplicates of docked structures were not docked."
    limits = [
        "Vina scores are predictions. They do not show real binding or enzyme inhibition.",
        "Scores are not comparable between proteins, so each safety protein is compared only with its own known drug.",
        "hERG docking is a flag, not a verdict, and its redocking was only PARTIAL, so hERG scores are low confidence.",
        f"The screen used the lowest search effort (exhaustiveness {EXHAUSTIVENESS_SCREEN}); scores within about {NOISE_KCAL} kcal/mol are inside the noise.",
        "No neem compound scores better than the known drug 605 in predicted DPP-4 fit." if (data.get("docking_target") is not None and not is_missing(ref_score(rd, "dpp4")) and float(data["docking_target"]["vina_score_kcal"].min()) > ref_score(rd, "dpp4")) else "Compare every neem score with the known drug 605 before drawing conclusions.",
        stereo,
        "ADMET risk numbers are model predictions and need experimental validation.",
        "Botanical provenance and exact compound identities need confirming against the primary literature.",
        "Everything needs laboratory confirmation. No experimental DPP-4 inhibition results are reported.",
    ]
    st.warning("\n".join(f"- {x}" for x in limits))
    nxt = (data.get("target_card") or {}).get("proposed_next_steps")
    if isinstance(nxt, list) and nxt:
        with st.expander("Proposed next steps (from the target card)"):
            for x in nxt:
                st.markdown(f"- {x}")

    st.subheader("Methods in short")
    st.markdown(
        "- **Library:** neem (*Azadirachta indica*) compounds from COCONUT, with literature DOIs.\n"
        "- **Filter:** drug-likeness (RDKit: Lipinski, Veber, QED, PAINS) plus ML risk scores for hERG and liver toxicity.\n"
        "- **Docking:** AutoDock Vina 1.2.x with ligands prepared by Meeko and receptors prepared with PDBFixer; the search box comes from the crystal drug.\n"
        f"- **Search effort:** exhaustiveness {EXHAUSTIVENESS_SCREEN} for the screen; higher for the redocking checks (as recorded in the files).\n"
        "- **Interactions:** PLIP, for the top 10 compounds.\n"
        "- **Safety panel:** hERG (8ZYO), CYP3A4 (2V0M), CYP2D6 (3TDA); each compared with its own redocked known drug.")

    st.subheader("Glossary")
    with st.expander("Plain-English glossary (20 terms)", expanded=False):
        for k, v in GLOSSARY.items():
            st.markdown(f"**{k}** — {v}")

    st.subheader("Evidence pack for judges")
    buf = io.BytesIO()
    parts = []
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for key, name in (("ranking", "ranking.csv"), ("docking_antitargets", "safety_panel.csv")):
            if data.get(key) is not None:
                z.writestr(name, data[key].to_csv(index=False)); parts.append(name)
        if rd:
            rt = pd.DataFrame([{"receptor": k, "ligand": v.get("ligand"), "verdict": v.get("verdict"), "top1_rmsd_A": v.get("top1_rmsd"),
                                "best_rmsd_A": v.get("best_rmsd"), "ref_score_kcal_mol": v.get("ref_score")} for k, v in rd.items()])
            z.writestr("redocking.csv", rt.to_csv(index=False)); parts.append("redocking.csv")
        z.writestr("README.txt", "OmniPhyto evidence pack. All values are computer predictions, not lab measurements.\n"
                                 "Vina scores are in kcal/mol (lower = better predicted fit) and are not comparable between proteins.\n")
    if parts:
        st.download_button("Download evidence pack (ZIP)", buf.getvalue(), file_name="omniphyto_evidence_pack.zip",
                           mime="application/zip", key="evi_dl")
        st.caption("Contains: " + ", ".join(parts))
    else:
        st.button("Download evidence pack (waiting for files)", disabled=True, key="evi_dl_off")
    section_footer("Evidence & methods · Project Manager · calibrated honesty beats inflated claims.")



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


# =============================================================================
# 11. SAGE & CREAM THEME   (ADDED block. It restyles the app and changes no existing code.)
#     It wraps inject_css(), so main() picks it up on its own. Delete this block to go back.
# =============================================================================
SAGE_PALETTES = {
    False: {  # light: cream page, sage sidebar, orange call-to-action
        "BG": "#f7f0e5", "CARD": "#fffaf1", "LINE": "#dfd4c0", "INK": "#2b3a2e", "MUTED": "#5d6b5f",
        "SIDE1": "#a9bba7", "SIDE2": "#94a894", "SIDEINK": "#1f2d22",
        "BLOB1": "rgba(135,153,137,.30)", "BLOB2": "rgba(231,151,107,.16)",
        "ACCENT": "#c2561a", "ACCENTHOVER": "#a94a15", "ACCENTINK": "#ffffff", "LINK": "#2f6b4f",
        "PILL": "rgba(255,255,255,.35)", "PILLON": "#fffaf1", "SAGE": "#59725d", "SHADOW": "rgba(60,80,60,.10)",
    },
    True: {  # dark: deep forest green with the same sage and orange accents
        "BG": "#17211b", "CARD": "#1f2d25", "LINE": "#33473a", "INK": "#f1efe4", "MUTED": "#a9b8ab",
        "SIDE1": "#22332a", "SIDE2": "#1a261f", "SIDEINK": "#eaf0e4",
        "BLOB1": "rgba(111,148,120,.16)", "BLOB2": "rgba(232,121,63,.08)",
        "ACCENT": "#f08a4b", "ACCENTHOVER": "#f39a63", "ACCENTINK": "#1b130c", "LINK": "#9fd3a9",
        "PILL": "rgba(255,255,255,.07)", "PILLON": "#2f4538", "SAGE": "#7fa088", "SHADOW": "rgba(0,0,0,.30)",
    },
}

SAGE_CSS = """<style>
@import url('https://fonts.googleapis.com/css2?family=Nunito:wght@500;700;800;900&display=swap');
.stApp, [data-testid="stAppViewContainer"] {
  background:
    radial-gradient(circle at 96% -4%, __BLOB1__ 0, __BLOB1__ 190px, transparent 191px),
    radial-gradient(ellipse at -6% 106%, __BLOB2__ 0, transparent 46%),
    __BG__ !important; background-attachment: fixed !important; color: __INK__; }
[data-testid="stHeader"] { background: transparent !important; }
.stApp h1, .stApp h2, .stApp h3, .stApp h4, .stApp h5, .stApp p, .stApp li, .stApp label,
.stApp button, .stApp [data-testid="stMarkdownContainer"], .stApp [data-testid="stCaptionContainer"] {
  font-family: 'Nunito', 'Segoe UI', system-ui, sans-serif; }
.stApp h1, .stApp h2, .stApp h3, .stApp h4 { color: __INK__; font-weight: 800; letter-spacing: -0.015em; }
.stApp h1 { font-weight: 900; line-height: 1.12; }
.stApp, .stApp p, .stApp li, .stApp label { color: __INK__; }
[data-testid="stCaptionContainer"], .stApp small { color: __MUTED__ !important; }
.stApp a { color: __LINK__ !important; text-decoration: underline; text-underline-offset: 2px; }

[data-testid="stSidebar"] { background: linear-gradient(180deg, __SIDE1__ 0%, __SIDE2__ 100%) !important;
  border-right: 0; border-radius: 0 30px 30px 0; box-shadow: 8px 0 28px __SHADOW__; }
[data-testid="stSidebar"] h1, [data-testid="stSidebar"] h2, [data-testid="stSidebar"] h3, [data-testid="stSidebar"] p,
[data-testid="stSidebar"] label, [data-testid="stSidebar"] [data-testid="stCaptionContainer"] { color: __SIDEINK__ !important; }
[data-testid="stSidebar"] [data-testid="stRadio"] label { background: __PILL__; border: 1px solid transparent;
  border-radius: 999px; padding: 6px 14px; margin: 3px 0; transition: background .15s; }
[data-testid="stSidebar"] [data-testid="stRadio"] label:hover { background: __PILLON__; }
[data-testid="stSidebar"] [data-testid="stRadio"] label:has(input:checked) { background: __PILLON__;
  border-color: __SAGE__; font-weight: 800; box-shadow: 0 4px 12px __SHADOW__; }

.stApp .stButton > button, .stApp .stDownloadButton > button, .stApp [data-testid^="stBaseButton"] {
  border-radius: 999px; font-weight: 800; padding: .5rem 1.4rem; border: 1.5px solid __SAGE__;
  background: __CARD__; color: __INK__; }
.stApp [data-testid="stBaseButton-primary"], .stApp .stButton > button[kind="primary"] {
  background: __ACCENT__ !important; color: __ACCENTINK__ !important; border-color: __ACCENT__ !important;
  box-shadow: 0 8px 18px rgba(194,86,26,.28); }
.stApp [data-testid="stBaseButton-primary"]:hover { background: __ACCENTHOVER__ !important; }
.stApp [data-testid="stBaseButton-primary"] p, .stApp [data-testid="stBaseButton-primary"] span { color: __ACCENTINK__ !important; }

.op-card, .op-wait { background: __CARD__; border: 1px solid __LINE__; border-radius: 22px;
  box-shadow: 0 10px 26px __SHADOW__; }
.op-wait { border-left: 6px solid __SAGE__; }
.op-metric-v { color: __INK__; font-weight: 800; }
.op-metric-l, .op-metric-s, .op-foot, .op-kv td:first-child, .op-tbl th { color: __MUTED__; }
.op-tbl th, .op-tbl td { border-bottom-color: __LINE__; }
.op-chip { font-weight: 700; }
.op-green { background: #3f7a50; color: #fff; } .op-blue { background: #2f6b7a; color: #fff; }
.op-red { background: #b0392b; color: #fff; } .op-grey { background: #6a756d; color: #fff; }
.op-amber { background: #f2b64a; color: #2b2a1d; } .op-outline { color: __INK__; border-color: __MUTED__; }

[data-testid="stExpander"] { background: __CARD__; border: 1px solid __LINE__ !important; border-radius: 20px; overflow: hidden; }
[data-testid="stAlert"] { border-radius: 18px; }
[data-testid="stDataFrame"] { border: 1px solid __LINE__; border-radius: 18px; overflow: hidden; }
[data-baseweb="input"], [data-baseweb="select"] > div, [data-baseweb="textarea"] { border-radius: 14px !important; }
button[data-baseweb="tab"] { font-weight: 700; }
.block-container { padding-top: 2.2rem; }
</style>"""


def inject_sage_theme(dark: bool) -> None:
    css = SAGE_CSS
    for token, value in SAGE_PALETTES[bool(dark)].items():
        css = css.replace(f"__{token}__", value)
    st.markdown(css, unsafe_allow_html=True)


_inject_css_original = inject_css


def inject_css(dark: bool):  # noqa: F811  (same name on purpose, so main() uses the themed version)
    _inject_css_original(dark)
    inject_sage_theme(dark)


main()
