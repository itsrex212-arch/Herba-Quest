"""Streamlit integration checks for the non-demo workflow."""
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from streamlit.testing.v1 import AppTest


class DashboardInterfaceTests(unittest.TestCase):
    def test_pages_load_without_exceptions_and_have_no_demo_widgets(self):
        app_path = Path(__file__).resolve().parents[1] / "app.py"
        for page_name in ("Overview", "Upload compounds", "Run Pipeline", "Ligand prep", "Results", "Settings"):
            with self.subTest(page=page_name), TemporaryDirectory() as temp:
                at = AppTest.from_file(str(app_path)).run()
                self.assertFalse(at.exception)
                at.radio[0].set_value(page_name).run()
                self.assertFalse(at.exception)
                visible = " ".join(x.value for x in at.title) + " " + " ".join(x.value for x in at.info)
                self.assertNotIn("DEMO_", visible)

    def test_overview_demo_toggle_shows_tagged_editable_demo_and_defaults_off(self):
        app_path = Path(__file__).resolve().parents[1] / "app.py"
        at = AppTest.from_file(str(app_path)).run()
        self.assertFalse(at.exception)
        self.assertTrue(any("DEMO dashboard" in x.label for x in at.sidebar.toggle))
        rendered = " ".join(x.value for x in at.markdown)
        self.assertNotIn("DEMO DATA ONLY", rendered)
        at.sidebar.toggle[0].set_value(True).run()
        self.assertFalse(at.exception)
        rendered = " ".join(x.value for x in at.markdown)
        self.assertIn("DEMO DATA ONLY", rendered)
        self.assertIn("DEMO PLACEHOLDER", rendered)
        self.assertIn("Top-ranked candidates", rendered)
        self.assertTrue(any("Edit demo dashboard model" in x.label for x in at.expander))
        self.assertTrue(any("DEMO ONLY" in x.value for x in at.markdown))

    def test_demo_mode_does_not_replace_real_zero_input_metric(self):
        app_path = Path(__file__).resolve().parents[1] / "app.py"
        at = AppTest.from_file(str(app_path)).run()
        self.assertFalse(at.exception)
        self.assertTrue(any("COMPOUNDS LOADED" in x.value and "<div class=\"metric-value\">0</div>" in x.value for x in at.markdown))
        at.sidebar.toggle[0].set_value(True).run()
        self.assertFalse(at.exception)
        at.sidebar.toggle[0].set_value(False).run()
        self.assertFalse(at.exception)
        rendered = " ".join(x.value for x in at.markdown)
        self.assertNotIn("DEMO DATA ONLY", rendered)
        self.assertIn("COMPOUNDS LOADED", rendered)

    def test_results_page_has_honest_empty_state_without_outputs(self):
        app_path = Path(__file__).resolve().parents[1] / "app.py"
        at = AppTest.from_file(str(app_path)).run()
        at.radio[0].set_value("Results").run()
        messages = " ".join(x.value for x in at.info)
        self.assertIn("No validated ligand-preparation outputs", messages)


if __name__ == "__main__":
    unittest.main()
