# OmniPhyto data manifest
_generated 2026-11-07 12:50 by the Drive-to-GitHub notebook_

## Delivered (files per folder)

- `data/incoming`: 30 files
- `data/incoming/results`: 268 files
- `data/plip`: 10 files
- `data/prepared`: 256 files
- `data/receptors`: 4 files
- `data/validation`: 16 files
- `md_report/605_hbonds.png`: 1 files
- `md_report/605_hbonds_detail.csv`: 1 files
- `md_report/605_rmsd.png`: 1 files
- `md_report/605_rmsf.png`: 1 files
- `md_report/OP0048_hbonds.png`: 1 files
- `md_report/OP0048_hbonds_detail.csv`: 1 files
- `md_report/OP0048_rmsd.png`: 1 files
- `md_report/OP0048_rmsf.png`: 1 files
- `md_report/OP0145_hbonds.png`: 1 files
- `md_report/OP0145_hbonds_detail.csv`: 1 files
- `md_report/OP0145_rmsd.png`: 1 files
- `md_report/OP0145_rmsf.png`: 1 files
- `md_report/OP0437_hbonds.png`: 1 files
- `md_report/OP0437_rmsd.png`: 1 files
- `md_report/OP0437_rmsf.png`: 1 files
- `md_report/md_summary.csv`: 1 files

## Expected, not delivered yet

- `data/incoming/ranking.csv`: full ranking of all docked compounds from the ML engineer (only the top 10 so far)
- `data/incoming/md_report/`: MD results (delivered if `omniphyto_run/md_report/` exists on Drive; see the README in that folder)
- poses for compounds docked later: copy the new `results/` folder into `data/incoming/`

## Where things are

- `data/incoming/`: read by the app (tables, poses, target card)
- `data/receptors/`: raw RCSB structures for the 3D viewer
- `data/prepared/`: receptors and ligands as docked (PDBQT, box, config) for reproducibility
- `data/validation/`: redocking inputs, poses and RMSD tables
- `data/plip/`: interaction complexes for the top compounds
- `notebooks/`: the notebooks that produced the results
