"""
Generate hERG second opinion reviews and compound rankings from scoring and docking results.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent

SCORED_FILE = "neem leave_scored.csv"
TARGET_FILE = "docking_target.csv"
ANTI_FILE = "docking_antitargets.csv"

NO_POSE_ABOVE = -2.0

TOP_N_HERG = 10
ML_FLAG = 0.50
DOCK_MARGIN = 1.0
ML_HIGH = 0.60
DOCK_HIGH = 2.0

W_EFF, W_SEL, W_ADMET = 0.5, 0.3, 0.2

PREFIX = Path(SCORED_FILE).stem
PREFIX = PREFIX[: -len("_scored")] if PREFIX.endswith("_scored") else PREFIX
RANKING_OUT = f"{PREFIX}_ranking.csv"
EXCLUDED_OUT = f"{PREFIX}_excluded.csv"
HERG_OUT = f"{PREFIX}_herg_second_opinion.csv"


def stop(msg):
    sys.exit("\nSTOP: " + msg)


def load(name, needed):
    p = HERE / name
    if not p.exists():
        stop(f"{name} is not in the folder {HERE}")
    df = pd.read_csv(p, dtype={"compound_id": str}, sep=None, engine="python", encoding="utf-8-sig")
    df.columns = df.columns.str.strip()
    missing = [c for c in needed if c not in df.columns]
    if missing:
        stop(f"{name} has no column {missing}.")
    df["compound_id"] = df["compound_id"].str.strip()
    return df


def check_scores(df, name, key):
    s = pd.to_numeric(df["vina_score_kcal"], errors="coerce")
    if (s > 0).mean() > 0.5:
        stop(f"{name}: most Vina scores are POSITIVE.")
    if df.duplicated(key).any():
        stop(f"{name}: repeated rows found.")
    df["vina_score_kcal"] = s
    bad = s.isna() | (s > NO_POSE_ABOVE)
    n_blank, n_weak = int(s.isna().sum()), int((s > NO_POSE_ABOVE).sum())
    print(f"  {name}: {len(df)} rows; {int((~bad).sum())} valid poses; NO VALID POSE = {int(bad.sum())} "
          f"({n_weak} scored above {NO_POSE_ABOVE}, {n_blank} blank/non-numeric)")
    if n_blank:
        print(f"    blank/non-numeric in: {list(df.loc[s.isna(), 'compound_id'])[:8]}")
    if (s < -14).any():
        print(f"    WARNING: scores below -14.")
    return df, set(df.loc[bad, "compound_id"])


def minmax(s):
    if s.max() == s.min():
        stop("all values are equal, so they cannot be scaled to 0..1.")
    return (s - s.min()) / (s.max() - s.min())


scored = load(SCORED_FILE, ["compound_id", "herg_risk", "admet_score", "pass_flag"])
print("Docking files:")
tgt, tgt_nopose = check_scores(load(TARGET_FILE, ["compound_id", "vina_score_kcal"]), TARGET_FILE, ["compound_id"])
anti = load(ANTI_FILE, ["compound_id", "antitarget", "vina_score_kcal"])
anti["antitarget"] = anti["antitarget"].astype(str).str.strip()
anti, _ = check_scores(anti, ANTI_FILE, ["compound_id", "antitarget"])
anti_all_ids = set(anti["compound_id"])
anti_ok = anti["vina_score_kcal"].notna() & (anti["vina_score_kcal"] <= NO_POSE_ABOVE)
herg_nopose = set(anti.loc[~anti_ok & (anti["antitarget"].str.lower() == "herg"), "compound_id"])
anti = anti[anti_ok].copy()

if scored["compound_id"].duplicated().any():
    stop(f"{SCORED_FILE} has repeated compound_id values.")
scored["pass_flag"] = scored["pass_flag"].astype(str).str.strip().str.lower().eq("true")
unknown = (set(tgt["compound_id"]) | set(anti_all_ids)) - set(scored["compound_id"])
if unknown:
    stop(f"docking files contain compound_ids that are not in {SCORED_FILE}.")
if "MOCK" in " ".join(tgt.get("pose_file", pd.Series(dtype=str)).astype(str).head(50)).upper():
    stop(f"{TARGET_FILE} looks like MOCK data.")
print(f"Loaded: {len(scored)} scored compounds ({int(scored['pass_flag'].sum())} passed), "
      f"{len(tgt)} rows in target file, {anti['compound_id'].nunique()} with valid antitarget pose.")

excluded = [(c, "pass_flag is False") for c in scored.loc[~scored["pass_flag"], "compound_id"]]
docked = scored[scored["pass_flag"]][["compound_id", "herg_risk", "admet_score"]].merge(
    tgt[["compound_id", "vina_score_kcal"]].rename(columns={"vina_score_kcal": "target_kcal"}), on="compound_id", how="left")
docked.loc[docked["compound_id"].isin(tgt_nopose), "target_kcal"] = np.nan
for c in docked.loc[docked["target_kcal"].isna(), "compound_id"]:
    excluded.append((c, f"no valid pose for disease target" if c in tgt_nopose
                        else "no disease-target docking result"))
docked = docked[docked["target_kcal"].notna()].copy()
if len(docked) < 4:
    stop(f"only {len(docked)} passed compound(s) have a valid disease-target pose.")
docked["efficacy_raw"] = -docked["target_kcal"]
docked["efficacy_norm"] = minmax(docked["efficacy_raw"])
eff_lo, eff_hi = docked["efficacy_raw"].min(), docked["efficacy_raw"].max()
print(f"\nDisease-target list: {len(docked)} passed compounds with valid pose.")

idx = anti.groupby("compound_id")["vina_score_kcal"].idxmin()
best = anti.loc[idx, ["compound_id", "antitarget", "vina_score_kcal"]].rename(
    columns={"antitarget": "strongest_antitarget", "vina_score_kcal": "antitarget_kcal"})

herg_rows = anti[anti["antitarget"].str.lower() == "herg"][["compound_id", "vina_score_kcal"]] \
    .rename(columns={"vina_score_kcal": "herg_kcal"})
top = docked.sort_values(["target_kcal", "compound_id"]).head(TOP_N_HERG)[["compound_id", "target_kcal", "herg_risk"]].copy()
top.insert(0, "target_rank", np.arange(1, len(top) + 1))
top = top.merge(herg_rows, on="compound_id", how="left")
gap_raw = top["herg_kcal"] - top["target_kcal"]
top["herg_minus_target_kcal"] = gap_raw.round(2)
has = top["herg_kcal"].notna()
ml_flag = top["herg_risk"] >= ML_FLAG
dock_flag = gap_raw < DOCK_MARGIN
top["ml_says"] = np.where(ml_flag, "flag", "lower risk")
top["docking_says"] = np.where(~has, "no valid result", np.where(dock_flag, "flag", "lower risk"))
none_label = np.where(top["compound_id"].isin(herg_nopose), "no valid hERG pose: ML only", "no hERG docking: ML only")
top["second_opinion"] = np.select([~has, ml_flag & dock_flag, ~ml_flag & ~dock_flag],
                                 [none_label, "agree: flag", "agree: lower risk"], default="disagree: needs human review")
ml_high = top["herg_risk"] >= ML_HIGH
dock_high = gap_raw <= -DOCK_HIGH
top["flag_for_review"] = (ml_flag & dock_flag) | ml_high | dock_high
top["flag_reason"] = np.select(
    [ml_flag & dock_flag, ml_high, dock_high], ["both methods flag it", f"ML risk >= {ML_HIGH}", f"hERG docks stronger than target"], default="")
top["herg_risk"] = top["herg_risk"].round(3)
top[["target_rank", "compound_id", "target_kcal", "herg_kcal", "herg_minus_target_kcal", "herg_risk",
     "ml_says", "docking_says", "second_opinion", "flag_for_review", "flag_reason"]].to_csv(HERE / HERG_OUT, index=False)

print(f"\n=== hERG second opinion: top {len(top)} docking hits ===")
print(top[["target_rank", "compound_id", "target_kcal", "herg_kcal", "herg_minus_target_kcal", "herg_risk",
           "second_opinion", "flag_for_review"]].to_string(index=False))
for label, n in top["second_opinion"].value_counts().items():
    print(f"    {n:3d}  {label}")

df = docked.merge(best, on="compound_id", how="left")
for c in df.loc[df["antitarget_kcal"].isna(), "compound_id"]:
    excluded.append((c, f"no valid antitarget pose" if c in anti_all_ids else "no antitarget docking result"))
df = df[df["antitarget_kcal"].notna()].copy()

assert len(df) + len(excluded) == len(scored), "row count mismatch"
if len(df) < 4:
    stop(f"only {len(df)} compound(s) have everything needed.")
if df["admet_score"].isna().any() or not df["admet_score"].between(0, 1).all():
    stop("admet_score has blanks or values outside 0..1.")

df["selectivity_raw"] = df["antitarget_kcal"] - df["target_kcal"]
df["selectivity_norm"] = minmax(df["selectivity_raw"])
df["composite_score"] = W_EFF * df["efficacy_norm"] + W_SEL * df["selectivity_norm"] + W_ADMET * df["admet_score"]
df["rank"] = df["composite_score"].rank(method="min", ascending=False).astype(int)
df = df.sort_values(["rank", "compound_id"]).reset_index(drop=True)

assert docked.loc[docked["target_kcal"].idxmin(), "efficacy_norm"] == 1, "efficacy sign is backwards"
assert df.loc[df["selectivity_raw"].idxmax(), "selectivity_norm"] == 1, "selectivity sign is backwards"
assert df[["efficacy_norm", "selectivity_norm", "composite_score"]].apply(lambda c: c.between(0, 1).all()).all()

df[["compound_id", "efficacy_norm", "selectivity_norm", "admet_score", "composite_score", "rank"]] \
    .round(6).to_csv(HERE / RANKING_OUT, index=False)
pd.DataFrame(excluded, columns=["compound_id", "reason"]).to_csv(HERE / EXCLUDED_OUT, index=False)

print(f"\n=== Ranking ===\nRanked {len(df)} compounds. Not ranked: {len(excluded)} (see {EXCLUDED_OUT})")
for reason, n in pd.Series([r for _, r in excluded]).value_counts().items():
    print(f"    {n:5d}  {reason}")
print("\nTOP 10")
print(df.head(10)[["rank", "compound_id", "target_kcal", "strongest_antitarget", "antitarget_kcal",
                    "efficacy_norm", "selectivity_norm", "admet_score", "composite_score"]].round(3).to_string(index=False))

slo, shi = df["selectivity_raw"].min(), df["selectivity_raw"].max()
print("\nHAND CHECK (top 3):")
for _, r in df.head(3).iterrows():
    e = (r["efficacy_raw"] - eff_lo) / (eff_hi - eff_lo)
    s = (r["selectivity_raw"] - slo) / (shi - slo)
    c = W_EFF * e + W_SEL * s + W_ADMET * r["admet_score"]
    print(f"  {r['compound_id']}: efficacy ({r['efficacy_raw']:.1f} - {eff_lo:.1f}) / ({eff_hi:.1f} - {eff_lo:.1f}) = {e:.3f} | "
          f"selectivity {r['selectivity_raw']:.1f} -> {s:.3f} | "
          f"composite {W_EFF}*{e:.3f} + {W_SEL}*{s:.3f} + {W_ADMET}*{r['admet_score']:.3f} = {c:.3f} (table {r['composite_score']:.3f})")
    assert abs(c - r["composite_score"]) < 1e-5, "hand calculation does not match table"

X = df[["efficacy_norm", "selectivity_norm", "admet_score"]].to_numpy()
ids = df["compound_id"].to_numpy()
base3 = set(ids[:3])
sets = [(i / 10, j / 10, (10 - i - j) / 10) for i in range(11) for j in range(11 - i)]
sets = [w for w in sets if all(0.1 - 1e-9 <= x <= 0.6 + 1e-9 for x in w)]
count = {c: 0 for c in ids}
same = 0
for w in sets:
    order = np.argsort(-(X @ np.array(w)), kind="stable")
    top3 = set(ids[order[:3]])
    same += top3 == base3
    for c in top3:
        count[c] += 1
print(f"\nWEIGHT CHECK: {len(sets)} reasonable weight sets.")
for c in ids[:3]:
    print(f"  {c} stays in top 3 in {100 * count[c] / len(sets):.0f}% of them")
print(f"  same three compounds are top 3 in {100 * same / len(sets):.0f}% of them")
print(f"\nSaved {HERG_OUT}, {RANKING_OUT} and {EXCLUDED_OUT} in", HERE)