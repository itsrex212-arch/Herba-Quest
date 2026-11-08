# HerbaQuest: neem compounds against DPP-4, screened for safety

**EthnoHack 2026 · Team: HerbaQuest · Challenge: Track 02 – Natural Products for Drug Discovery**

Repository: `makeathon-2026-herbaquest-track-02-natural-products-for-drug-discovery`

HerbaQuest (the OmniPhyto pipeline) screens neem (*Azadirachta indica*) leaf compounds against DPP-4, a validated type 2 diabetes target, and checks every top candidate against three safety proteins before ranking it. A Streamlit web app presents every step, with a built-in guided demo.

> These are computational predictions only. Nothing here is a clinical claim, a measured affinity, or a validated drug lead. Wet-lab confirmation is required.

**Quick links:** [Demo video](https://drive.google.com/file/d/1z15GzfTW1MqKvd9rE-A__TIVJWt1mtJ0/view?usp=sharing) · Slide deck: [[add link]](https://drive.google.com/file/d/18bvgKlF9fZf_KsIwqIzolrPGjZkZv0KY/view?usp=drive_link) · Live app: [[add URL]](https://omniphyto.streamlit.app/)

## 1. Project description

### Chosen challenge
Track 02: Natural Products for Drug Discovery.

### The problem
589 million adults (aged 20-79) live with diabetes and 252 million of them are undiagnosed (IDF Diabetes Atlas, 11th edition, 2025). DPP-4 inhibitors are an established treatment, but a compound that binds DPP-4 is useless if it also blocks hERG (heart rhythm) or CYP3A4 and CYP2D6 (drug metabolism). Many leads fail late for exactly this reason. Our question: can compounds from neem block DPP-4 and pass an early safety check?

### Our solution
1. **Filter.** 437 neem-leaf compounds are scored for drug-likeness (Lipinski, QED, PAINS) and by machine-learning hERG and liver-injury risk models. 300 pass (at most one Lipinski violation, hERG and liver risk both below 0.70).
2. **Dock.** 238 of those are docked on DPP-4 (PDB 3D4L). The top 10 are also docked on hERG (8ZYO), CYP3A4 (2V0M) and CYP2D6 (3TDA). Each protocol is validated by redocking the crystal ligand.
3. **Rank.** Composite = 0.5 efficacy + 0.3 selectivity + 0.2 ADMET. Weight sensitivity was checked over 27 weight sets.
4. **Simulate.** 10 ns molecular dynamics on the two leads (OP0048, OP0145), the crystal ligand 605 (control) and OP0437 (negative control).

### Key results
| Check | Result |
|---|---|
| Redocking RMSD (top pose) | DPP-4 0.33 Å, CYP3A4 0.70 Å, CYP2D6 1.48 Å pass; hERG 3.95 Å (best of 9 poses 0.86 Å), partial |
| Best DPP-4 docking score | OP0048, −9.385 kcal/mol (reference inhibitor 605: −10.077) |
| Combined ranking | OP0145 (0.916), OP0048 (0.859), OP0385 (0.836) |
| MD ligand RMSD, 10 ns | 605: 0.92 Å · OP0048: 3.68 Å · OP0437: 4.39 Å · OP0145: 5.90 Å |
| MD H-bond frames | OP0048 99.8% (ARG358) · OP0145 51.1% · 605 42% · OP0437 0% |

How to read these: no neem compound beats the known inhibitor in predicted DPP-4 fit. In the simulation, none of the neem candidates held its pose as tightly as the crystal ligand. OP0048 stayed closest to its starting pose among the candidates and kept hydrogen bonds almost throughout. The top-ranked OP0145 moved more than the negative control, so its docked pose may not be stable, which tempers its rank.

Honest limits: the safety panel is a flag, not a verdict (18 of 30 safety dockings score stronger than the DPP-4 score, and hERG redocking is only partial); the DPP-4 screen used low search effort (exhaustiveness 4), so scores within about 0.5 kcal/mol are effectively tied; one 10 ns replica per system; chain A only for DPP-4; no binding free energy; the ML models were trained on drug data, not plant compounds.

## 2. Technologies

- **Languages:** Python 3 (notebooks, scripts, app); JavaScript (3Dmol.js, embedded in the app for the 3D viewer).
- **Docking and structure:** AutoDock Vina 1.2.5, PDBFixer, RDKit, Meeko, Open Babel, Dimorphite-DL, PLIP, spyrmsd.
- **Molecular dynamics:** OpenMM 8.6 (Amber14, GAFF 2.11, TIP3P-FB), MDTraj.
- **Machine learning:** scikit-learn (logistic regression for hERG, random forest for liver injury), Morgan fingerprints, Therapeutics Data Commons datasets.
- **App:** Streamlit, pandas, numpy, Altair, RDKit (2D structures), 3Dmol.js (3D viewer with PLIP contacts).
- **Compute:** Google Colab and Kaggle GPU notebooks.

## 3. Repository layout

```
├── app.py                    HerbaQuest dashboard (Streamlit): intro + guided demo, ranking, 3D viewer, safety panel, MD page
├── requirements.txt          app dependencies
├── packages.txt              system packages for Streamlit Community Cloud
├── Codes/
│   ├── Scoring/              ML filter: rule filter, model training, scoring, datasets, trained models
│   ├── Docking&MD/           OmniPhyto_docking.ipynb (Vina docking) and OmniPhyto_MD_Kaggle.ipynb (OpenMM MD)
│   └── Ranking/              ranking_compounds.py and the ranking output
├── ML Hand off OmniPhyto/    ML hand-off package (data, models, ROC plots, scripts, scored tables)
├── data/
│   ├── incoming/             result tables the app reads (docking, safety docking, ranking, redocking, target card),
│   │   ├── results/          docking poses (PDBQT) for DPP-4 and the safety proteins
│   │   ├── validation/       crystal-vs-redocked overlays
│   │   └── md_report/        MD summary table, plots, H-bond tables
│   ├── receptors/            raw PDB structures (3D4L, 8ZYO, 2V0M, 3TDA)
│   ├── prepared/             prepared receptors and ligands (PDBQT), Vina configs
│   ├── plip/                 protein-ligand complexes for the top 10 compounds (PLIP contacts)
│   ├── validation/           redocking inputs, poses and RMSD tables
│   └── MANIFEST.md           what each data folder holds
├── md_report/                copy of the MD plots and H-bond files
├── notebooks/                HerbaQuest.ipynb (Colab notebook, outputs removed)
├── omniphyto/                helper package from the first prototype (not imported by the current app)
├── engines/                  placeholders for connecting a live docking or MD engine later
├── tests/                    smoke tests for the dashboard
└── assets/                   app styling
```

## 4. Execution instructions

### A. Run the app (Windows, VS Code, or any system with Python 3.10+)
1. Clone the repository and open the folder in VS Code, then open a new terminal.
2. Install dependencies: `python -m pip install -r requirements.txt`
3. Start the app: `python -m streamlit run app.py`
4. The app opens on **What we built**. Click **Start the guided tour** (or switch on **Guided demo** in the sidebar) for a step-by-step walk-through with pointers and plain-English interpretation of each result.
5. Sections follow the project's steps: Target & library, Structures & validation, ADMET filter, Docking: DPP-4, Safety panel, Final ranking, Simulation (MD), Compound profile, Evidence & methods, Data inbox.

Notes:
- The app only displays results from the files in `data/`. It does not run docking or MD. Add or update files in `data/incoming/` (or use **Data inbox** to check which section a file feeds) and refresh.
- The first load takes about 10 seconds and the 3D viewers can take about 5 seconds to appear. The sidebar switch **Fast 3D view** sends only the pocket area to the browser to keep it quick; switch it off to see the whole protein.
- 2D structures need RDKit (in `requirements.txt`). Without it the app shows the SMILES text instead.
- The 3D viewer loads 3Dmol.js from the internet, so it needs a connection.

### B. Reproduce docking and ranking (Colab)
1. Open `Codes/Docking&MD/OmniPhyto_docking.ipynb` in Google Colab.
2. Run the cells top to bottom. Key settings: Vina 1.2.5, 22 Å box centred on ligand 605, chain A only, pH 7.4. Exhaustiveness is 4 for the 238-compound screen and 16 for redocking and safety docking.
3. Outputs: redocking results, target docking table, safety docking table. The composite ranking is produced by `Codes/Ranking/ranking_compounds.py`.

### C. Reproduce the ML filter
1. Use `Codes/Scoring/`: `train_models.py` trains the hERG and liver-injury models from the datasets in `Codes/Scoring/data/`, `rule_filter.py` applies the drug-likeness rules, and `score_compounds.py` scores the compound table.
2. Output: the scored tables with `pass_flag` (300 of 437 pass).

### D. Reproduce the molecular dynamics (Kaggle or Colab GPU)
1. Open `Codes/Docking&MD/OmniPhyto_MD_Kaggle.ipynb` and select a GPU runtime.
2. Run the setup and production cells. Settings: OpenMM 8.6, amber14-all, TIP3P-FB, GAFF 2.11 with AM1-BCC charges, 0.15 M NaCl, 300 K, 1 bar, 4 fs steps, 10 ns, first 10% of frames dropped for means.
3. Outputs: `md_summary.csv`, `<id>_rmsd.png`, `<id>_rmsf.png`, `<id>_hbonds.png`, `<id>_hbonds_detail.csv`, saved in `data/incoming/md_report/`.

### E. Tests
```
python -m unittest discover -s tests -v
```
The smoke tests load every section of the dashboard against the data in `data/` and walk through the whole guided demo. They check that nothing crashes. They are not scientific validation.

## 5. API keys and secrets

The pipeline and app need **no API keys**. The only credential ever used is a GitHub personal access token in the Colab notebook that copies results from Google Drive to this repository.

- Create it at GitHub > Settings > Developer settings > Personal access tokens, with the minimum scope needed for this repository only.
- In Colab, open the key icon (Secrets), add a secret named `GITHUB_TOKEN`, and enable notebook access. The notebook reads it with `google.colab.userdata.get("GITHUB_TOKEN")`.
- Never paste a token into a notebook cell, a commit or an issue. If one is exposed, revoke it immediately.
- For local use, put secrets in a `.env` file that is listed in `.gitignore`:

```
.env
*.env
config/secrets.json
```

## 6. Presentation

- **Demo video (max 10 minutes):** https://drive.google.com/file/d/1z15GzfTW1MqKvd9rE-A__TIVJWt1mtJ0/view?usp=sharing
- **Slide deck (max 10 slides):** [add link]
- **Live app:** [add URL]
- **Screenshots:** <!-- add images to docs/ and uncomment: docs/screenshot_ranking.png, docs/screenshot_3d.png, docs/screenshot_md.png -->

## 7. Development practices

- Results were added through pull requests from the `omniphyto-data` branch into `main`.
- Large trajectories and raw simulation files are kept out of the repository (only summary tables and plots are stored); link to external storage instead.
- All team members are added as collaborators with Write or Admin access.
- No credentials are stored in the repository.

## 8. Team

HerbaQuest, EthnoHack 2026. Structural biology and bioinformatics: Goodness. [Add the other members and their roles.]

## 9. Limitations and scope

This is a research prototype. It is not for clinical decisions or patient care and is not validated for hospital use. Pending work: wet-lab enzyme assay and hERG test, repeat MD seeds, binding free-energy calculations, re-docking the top hits at higher exhaustiveness.

## 10. References

Full list with DOIs is on the last slide of the deck. Main sources: IDF Diabetes Atlas 11th ed. (2025); PDB 3D4L, 8ZYO, 2V0M, 3TDA primary citations; Trott and Olson 2010 and Eberhardt et al. 2021 (AutoDock Vina); Adasme et al. 2021 (PLIP); Eastman et al. 2024 (OpenMM 8); Maier et al. 2015 (ff14SB); Wang et al. 2004 (GAFF); McGibbon et al. 2015 (MDTraj).
