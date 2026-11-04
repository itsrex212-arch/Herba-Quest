"""Tests for the no-demo dashboard, upload schema, settings and fail-safe UI."""
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from omniphyto.dashboard import filter_rows, load_dashboard_data, sort_rows
from omniphyto.input_preview import validate_compound_csv
from omniphyto.settings import DEFAULTS, load_settings, save_settings, validate_settings


class DashboardTests(unittest.TestCase):
    def test_no_sample_fallback(self):
        with TemporaryDirectory() as temp:
            rows, source, is_demo = load_dashboard_data(Path(temp))
            self.assertEqual((rows, source, is_demo), ([], "", False))

    def test_manifest_is_loaded_from_real_output(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            output = root / "data" / "output"
            output.mkdir(parents=True)
            (output / "manifest.csv").write_text("compound_id,ligand_status\nX,failed\n", encoding="utf-8")
            self.assertEqual(load_dashboard_data(root), ([{"compound_id": "X", "ligand_status": "failed"}], "manifest.csv", False))

    def test_filter_and_sort(self):
        rows = [{"id": "A", "value": "10", "group": "leaf"}, {"id": "B", "value": "2", "group": "root"}]
        self.assertEqual([row["id"] for row in filter_rows(rows, query="LEAF")], ["A"])
        self.assertEqual([row["id"] for row in filter_rows(rows, field="group", value="root")], ["B"])
        self.assertEqual([row["id"] for row in sort_rows(rows, "value")], ["B", "A"])

    def test_app_loads_without_demo_widgets_or_exceptions(self):
        from streamlit.testing.v1 import AppTest
        app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "app.py"), default_timeout=30).run()
        self.assertEqual(len(app.exception), 0, [str(exc.value) for exc in app.exception])
        self.assertEqual(len(app.number_input), 0)
        rendered_markdown = "\n".join(item.value for item in app.markdown)
        self.assertIn("Welcome to OmniPhyto", rendered_markdown)
        self.assertIn("No compound records or preparation outputs yet", rendered_markdown)
        self.assertIn("Pipeline Progress", rendered_markdown)
        self.assertIn("Molecular View", rendered_markdown)
        self.assertIn("Candidates and Results", rendered_markdown)
        self.assertIn("Project Insights", rendered_markdown)

    def test_settings_and_results_pages_load(self):
        from streamlit.testing.v1 import AppTest
        path = str(Path(__file__).resolve().parents[1] / "app.py")
        for target in ("Settings", "Results", "Ligand prep", "Upload compounds", "Run Pipeline"):
            with self.subTest(page=target):
                app = AppTest.from_file(path, default_timeout=30).run()
                app.radio[0].set_value(target).run()
                self.assertEqual(len(app.exception), 0, [str(exc.value) for exc in app.exception])


class InputSchemaTests(unittest.TestCase):
    def test_valid_upload_schema(self):
        headers, rows = validate_compound_csv(b"compound_id,smiles\nC1,CCO\n")
        self.assertEqual(headers, ["compound_id", "smiles"])
        self.assertEqual(rows[0]["compound_id"], "C1")

    def test_missing_or_blank_required_values_fail(self):
        for body in (b"compound_id,name\nC1,ethanol\n", b"compound_id,smiles\n,CCO\n", b"compound_id,smiles\nC1,\n"):
            with self.subTest(body=body), self.assertRaises(ValueError):
                validate_compound_csv(body)

    def test_duplicate_ids_and_row_shape_fail(self):
        for body in (b"compound_id,smiles\nC1,CCO\nC1,CCC\n", b"compound_id,smiles\nC1,CCO,extra\n"):
            with self.subTest(body=body), self.assertRaises(ValueError):
                validate_compound_csv(body)


class SettingsTests(unittest.TestCase):
    def test_blank_configuration_has_no_scientific_defaults(self):
        self.assertEqual(DEFAULTS["target_name"], "")
        self.assertEqual(validate_settings(DEFAULTS), ["Project name is required", "Target name is required", "Protocol is required"])

    def test_box_requires_six_finite_values_and_positive_sizes(self):
        config = {**DEFAULTS, "project_name": "Study", "target_name": "Target", "protocol": "Protocol", "box_center_x": "1"}
        self.assertTrue(any("all six" in error for error in validate_settings(config)))
        config = {**config, "box_center_x": "1", "box_center_y": "2", "box_center_z": "3", "box_size_x": "10", "box_size_y": "0", "box_size_z": "10"}
        self.assertTrue(any("greater than zero" in error for error in validate_settings(config)))

    def test_settings_round_trip(self):
        with TemporaryDirectory() as temp:
            path = Path(temp) / "settings.json"
            config = {**DEFAULTS, "project_name": "Study", "target_name": "Target", "protocol": "team protocol"}
            save_settings(path, config)
            self.assertEqual(load_settings(path), config)

    def test_docking_guard_requires_actual_receptor_and_box(self):
        config = {**DEFAULTS, "project_name": "Study", "target_name": "Target", "protocol": "Protocol"}
        errors = validate_settings(config, require_receptor=True)
        self.assertTrue(any("Receptor PDBQT" in error for error in errors))
        self.assertTrue(any("box" in error.lower() for error in errors))


if __name__ == "__main__":
    unittest.main()
