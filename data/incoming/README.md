# data/incoming: the drop folder

Copy result files here and refresh the app. File names are matched by pattern; the newest file wins.

| Pattern | What it is |
|---|---|
| docking_target*.csv | DPP-4 docking, all docked compounds |
| docking_antitargets*.csv | safety docking (hERG, CYP3A4, CYP2D6) |
| summary_report*.csv | top-compound report (scores, residues, redocking verdicts) |
| redocking*.csv | redocking RMSD per protein |
| *scored_detail*.csv, *scored.csv | ML engineer tables |
| *compounds_with_DOIs*.csv | master compound table (SMILES, DOIs) |
| ranking*.csv | final ML ranking (arrives later) |
| target_card.json | disease, target, PDB, reference ligand, safety panel |
| results/<target|herg|cyp3a4|cyp2d6>/poses/*.pdbqt | pose files named in the pose_file column |
| md_report/ | GROMACS outputs (later) |

Receptor structures (3D4L.pdb, 8ZYO.pdb, 2V0M.pdb, 3TDA.pdb) go in data/receptors/.
Delivery is by git: the sync notebook commits files here, and the hosted app redeploys with them.
