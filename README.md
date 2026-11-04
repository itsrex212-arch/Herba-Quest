# OmniPhyto, non-demo project workspace

OmniPhyto accepts user-provided compound CSV data and can prepare ligand PDBQT files when its chemistry dependencies are installed. It does not provide clinical decision support and is not validated for hospital use.

## Start on Windows

1. Open this folder in VS Code.
2. Select Terminal > New Terminal.
3. Install app dependencies: `python -m pip install -r requirements.txt`.
4. Start the app: `python -m streamlit run app.py`.
5. In Upload compounds, select the actual CSV, choose its ID and SMILES headers, review validation errors, and save it.
6. In Settings, enter project/target/protocol information. Leave box fields empty until your team provides verified coordinates. A saved receptor path or coordinates do not establish scientific validity.
7. In Run Pipeline, click Run ligand preparation only when the required tools are available. The action executes implemented ligand-prep steps only.

The app works as an intake/settings/results interface without chemistry dependencies. Real preparation needs RDKit, Meeko and Open Babel installed in the same environment. The UI checks those tools before enabling the actual operation.

## Compound CSV

Required: a unique nonempty identifier column and a nonempty SMILES column. Exact column names are chosen during upload. Upload size is limited to 5 MB; encoding must be UTF-8. CSV validation is not molecular validation. RDKit later checks whether each SMILES can be parsed and embedded.

## Actual workflow boundaries

The implemented ligand prep uses RDKit to parse SMILES, add hydrogens, build a 3D conformer and minimize it with MMFF or UFF, Meeko to write PDBQT, and Open Babel to parse/round-trip validate the resulting file. Per-compound successes and errors are recorded in `data/output/manifest.csv`. Existing output files are not silently replaced.

This app workflow does not perform docking or scoring. Receptor, target-specific box, chosen docking program/version, scoring/reporting rule, and domain review have not been provided. Do not present this project as hospital-ready or use it for patient care. It is a research software prototype pending scientific validation, security/privacy review, reproducible evaluation, and appropriate institutional/regulatory review.

## Tests

From this folder run:

```powershell
python -m unittest discover -s tests -v
python -m compileall -q app.py omniphyto tests
```

Tests use temporary fixtures and mocked external tools. They do not count as scientific validation.
