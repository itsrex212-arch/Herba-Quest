"""Dependency-light tests; chemistry/Vina integration requires external tools."""
import csv
import subprocess
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from omniphyto.docking import DockingError, VinaConfig, _read_score, dock_ligand
from omniphyto.ligand_prep import LigandPreparationError, dependency_status, prepare_smiles
from omniphyto.pipeline import _safe_name, run_pipeline
from omniphyto.results import available_pose_files, load_results


class PipelineTests(unittest.TestCase):
    def test_safe_name(self):
        self.assertEqual(_safe_name("sample / compound 1"), "sample_compound_1")
        self.assertEqual(_safe_name("  "), "unnamed")

    def test_missing_dependencies_reported(self):
        with patch("omniphyto.ligand_prep.dependency_status",
                   return_value={"rdkit": False, "meeko": True, "obabel": True}):
            with self.assertRaisesRegex(LigandPreparationError, "rdkit"):
                prepare_smiles("CCO", "unused.pdbqt")

    def make_input(self, root):
        source = root / "input.csv"
        source.write_text("compound_id,smiles\nFAKE_A,CCO\n", encoding="utf-8")
        receptor = root / "receptor.pdbqt"
        receptor.write_text("FAKE receptor placeholder\n", encoding="utf-8")
        def fake_prepare(smiles, path, **kwargs):
            Path(path).write_text("FAKE ligand output, not real PDBQT\n", encoding="utf-8")
            return Path(path)
        return source, receptor, fake_prepare

    def read_row(self, manifest):
        with manifest.open(newline="", encoding="utf-8") as handle:
            return next(csv.DictReader(handle))

    def test_pipeline_without_docking_is_explicit(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            source, _, fake_prepare = self.make_input(root)
            with patch("omniphyto.pipeline.prepare_smiles", side_effect=fake_prepare):
                self.assertEqual(run_pipeline(source, root / "out", id_column="compound_id",
                                              smiles_column="smiles", seed=1, overwrite=False), (1, 1))
            row = self.read_row(root / "out" / "manifest.csv")
            self.assertEqual(row["docking_status"], "not_requested")
            self.assertEqual(row["score"], "")
            self.assertEqual(row["ligand_pdbqt"], "FAKE_A.pdbqt")

    def test_pipeline_docking_success_is_recorded(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            source, receptor, fake_prepare = self.make_input(root)
            cfg = VinaConfig(receptor, (1, 2, 3), (10, 10, 10), "vina")
            with patch("omniphyto.pipeline.prepare_smiles", side_effect=fake_prepare), \
                 patch("omniphyto.pipeline.VinaConfig.validate", return_value="/fake/vina"), \
                 patch("omniphyto.pipeline.dock_ligand", return_value=("-7.2", "mocked executable")):
                run_pipeline(source, root / "out", id_column="compound_id", smiles_column="smiles",
                             seed=1, overwrite=False, vina_config=cfg)
            row = self.read_row(root / "out" / "manifest.csv")
            self.assertEqual(row["docking_status"], "completed")
            self.assertEqual(row["score"], "-7.2")
            self.assertTrue(row["docking_log"].endswith("_vina.log"))

    def test_pipeline_docking_failure_keeps_score_blank(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            source, receptor, fake_prepare = self.make_input(root)
            cfg = VinaConfig(receptor, (1, 2, 3), (10, 10, 10), "vina")
            with patch("omniphyto.pipeline.prepare_smiles", side_effect=fake_prepare), \
                 patch("omniphyto.pipeline.VinaConfig.validate", return_value="/fake/vina"), \
                 patch("omniphyto.pipeline.dock_ligand", return_value=(None, "Vina failed")):
                run_pipeline(source, root / "out", id_column="compound_id", smiles_column="smiles",
                             seed=1, overwrite=False, vina_config=cfg)
            row = self.read_row(root / "out" / "manifest.csv")
            self.assertEqual(row["docking_status"], "failed")
            self.assertEqual(row["score"], "")

    def test_config_rejects_invalid_inputs(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            receptor = root / "receptor.pdbqt"
            receptor.write_text("placeholder\n", encoding="utf-8")
            with self.assertRaisesRegex(DockingError, "greater than zero"):
                VinaConfig(receptor, (0, 0, 0), (10, 0, 10)).validate()
            with self.assertRaisesRegex(DockingError, "does not exist"):
                VinaConfig(root / "missing.pdbqt", (0, 0, 0), (10, 10, 10)).validate()

    def test_score_parser_requires_vina_record(self):
        self.assertEqual(_read_score("REMARK VINA RESULT: -8.125 0.0 0.0\n"), "-8.125")
        self.assertIsNone(_read_score("FAKE score -8.1\n"))

    def test_dock_invocation_safe_and_captures_log(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            receptor, ligand = root / "receptor.pdbqt", root / "ligand.pdbqt"
            receptor.write_text("placeholder\n", encoding="utf-8")
            ligand.write_text("placeholder\n", encoding="utf-8")
            output, log = root / "out.pdbqt", root / "vina.log"
            config = VinaConfig(receptor, (1, 2, 3), (10, 11, 12), "vina")
            def fake_run(command, **kwargs):
                self.assertFalse(kwargs["shell"])
                output.write_text("REMARK VINA RESULT: -6.5 0 0\n", encoding="utf-8")
                return subprocess.CompletedProcess(command, 0, "ok", "")
            with patch("omniphyto.docking.shutil.which", return_value="/fake/vina"), \
                 patch("omniphyto.docking.subprocess.run", side_effect=fake_run):
                score, _ = dock_ligand(ligand, config, output, log)
            self.assertEqual(score, "-6.5")
            self.assertIn("exit_code: 0", log.read_text(encoding="utf-8"))

    def test_vina_nonzero_and_timeout_fail_without_score(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            receptor, ligand = root / "receptor.pdbqt", root / "ligand.pdbqt"
            receptor.write_text("placeholder\n", encoding="utf-8")
            ligand.write_text("placeholder\n", encoding="utf-8")
            config = VinaConfig(receptor, (1, 2, 3), (10, 10, 10), "vina")
            with patch("omniphyto.docking.shutil.which", return_value="/fake/vina"), \
                 patch("omniphyto.docking.subprocess.run",
                       return_value=subprocess.CompletedProcess([], 2, "", "bad input")):
                score, message = dock_ligand(ligand, config, root / "out.pdbqt", root / "nonzero.log")
            self.assertIsNone(score)
            self.assertIn("code 2", message)
            with patch("omniphyto.docking.shutil.which", return_value="/fake/vina"), \
                 patch("omniphyto.docking.subprocess.run", side_effect=subprocess.TimeoutExpired("vina", 1)):
                score, message = dock_ligand(ligand, config, root / "out2.pdbqt", root / "timeout.log")
            self.assertIsNone(score)
            self.assertIn("timed out", message)
            self.assertIn("timeout", (root / "timeout.log").read_text(encoding="utf-8"))

    def test_dependency_status_reports_expected_tools(self):
        self.assertEqual(set(dependency_status()), {"rdkit", "meeko", "obabel"})


class ResultsTests(unittest.TestCase):
    def test_missing_results_are_empty(self):
        with TemporaryDirectory() as temp:
            self.assertEqual(load_results(Path(temp)), ([], None))

    def test_manifest_hides_scores_unless_docking_completed(self):
        with TemporaryDirectory() as temp:
            output = Path(temp)
            (output / "manifest.csv").write_text(
                "compound_id,docking_status,score\nA,not_requested,-9.2\nB,completed,-7.4\n",
                encoding="utf-8")
            rows, source = load_results(output)
            self.assertEqual(source, "manifest.csv")
            self.assertEqual([row["score"] for row in rows], ["", "-7.4"])

    def test_ranking_csv_selected_before_manifest(self):
        with TemporaryDirectory() as temp:
            output = Path(temp)
            (output / "manifest.csv").write_text("compound_id\nmanifest\n", encoding="utf-8")
            (output / "ranking.csv").write_text("compound_id,rank\nA,1\n", encoding="utf-8")
            rows, source = load_results(output)
            self.assertEqual(source, "ranking.csv")
            self.assertEqual(rows[0]["rank"], "1")

    def test_pose_picker_only_lists_docked_or_pose_folder_files(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            output, poses = root / "data" / "output", root / "data" / "poses"
            output.mkdir(parents=True)
            poses.mkdir(parents=True)
            ligand = output / "ligand.pdbqt"
            docked = output / "compound_docked.pdbqt"
            pose = poses / "saved_pose.pdb"
            for path in (ligand, docked, pose):
                path.write_text("fixture\n", encoding="utf-8")
            found = available_pose_files(root)
            self.assertEqual(set(found), {pose.resolve()})
            (output / "manifest.csv").write_text(
                "compound_id,docking_status,docking_output\nA,completed,compound_docked.pdbqt\n",
                encoding="utf-8")
            found = available_pose_files(root)
            self.assertEqual(set(found), {docked.resolve(), pose.resolve()})


if __name__ == "__main__":
    unittest.main()
