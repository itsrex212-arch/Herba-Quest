"""Smoke test for the HerbaQuest dashboard (app.py).

Loads every section headless against the real files in data/ and walks through the whole guided demo.
It checks that nothing crashes and that no section shows an error. It is NOT a scientific validation.
Run:  python -m unittest discover -s tests -v
"""
import unittest
from pathlib import Path

from streamlit.testing.v1 import AppTest

APP = str(Path(__file__).resolve().parents[1] / "app.py")
SECTIONS = ["▶ What we built", "0. Home", "1. Target & library", "2. Structures & validation", "3. ADMET filter",
            "4. Docking: DPP-4", "5. Safety panel", "6. Final ranking", "7. Simulation (MD)", "8. Compound profile",
            "9. Evidence & methods", "10. Data inbox"]


class DashboardSmokeTests(unittest.TestCase):
    def test_every_section_loads_without_errors(self):
        at = AppTest.from_file(APP, default_timeout=180)
        at.run()
        self.assertFalse(at.exception, "landing page crashed")
        for name in SECTIONS:
            with self.subTest(section=name):
                at.session_state["section"] = name
                at.run()
                self.assertFalse(at.exception, f"{name} raised an exception")
                self.assertFalse(at.error, f"{name} shows an error box")

    def test_guided_demo_walks_through_every_stop(self):
        at = AppTest.from_file(APP, default_timeout=180)
        at.run()
        [b for b in at.button if b.key == "intro_start"][0].click()
        at.run()
        stops = 1
        while True:
            self.assertFalse(at.exception, f"demo stop {stops} crashed")
            nxt = [b for b in at.button if b.key == "demo_next"]
            if not nxt:
                break
            nxt[0].click()
            at.run()
            stops += 1
        self.assertGreaterEqual(stops, 8)
        text = "\n".join(m.value for m in at.markdown)
        self.assertIn("op-demo-t", text, "demo panel was not shown")


if __name__ == "__main__":
    unittest.main()
