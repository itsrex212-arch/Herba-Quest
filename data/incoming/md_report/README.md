# data/incoming/md_report: GROMACS results (expected, not delivered yet)

Put here, when the MD run is finished:

- `md_summary.csv` with the columns: compound_id, label (candidate|control), length_ns, ligand_rmsd_mean_A, ligand_rmsd_sd_A, protein_rmsd_mean_A, hbond_occupancy_pct, notes
- plots named `<compound_id>_rmsd.png`, `<compound_id>_rmsf.png`, `<compound_id>_hbonds.png`

`candidate` = the compound; `control` = the crystal drug (ligand 605) simulated the same way. The MD page lights up when these files appear.
