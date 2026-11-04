"""Tests for safe, non-persisting CSV previews and no bundled fictional compounds."""
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from omniphyto.input_preview import preview_csv


class CsvPreviewTests(unittest.TestCase):
    def test_preview_parses_headers_without_writing(self):
        content = b"sample_key,canonical_smiles\nA,CCO\n"
        with TemporaryDirectory() as temp:
            before = set(Path(temp).iterdir())
            headers, rows = preview_csv(content)
            self.assertEqual(headers, ["sample_key", "canonical_smiles"])
            self.assertEqual(rows, [{"sample_key": "A", "canonical_smiles": "CCO"}])
            self.assertEqual(set(Path(temp).iterdir()), before)

    def test_preview_rejects_empty_oversize_and_invalid_utf8(self):
        with self.assertRaisesRegex(ValueError, "empty"):
            preview_csv(b"")
        with self.assertRaisesRegex(ValueError, "exceeds"):
            preview_csv(b"a\n1\n", max_bytes=2)
        with self.assertRaisesRegex(ValueError, "UTF-8"):
            preview_csv(b"\xff")

    def test_preview_rejects_bad_headers_and_extra_fields(self):
        for content in (b"\nvalue\n", b"id,id\nA,A\n", b"id,smiles\nA,CCO,extra\n"):
            with self.subTest(content=content), self.assertRaises(ValueError):
                preview_csv(content)

    def test_preview_caps_rows(self):
        content = ("id,smiles\n" + "".join(f"{i},C\n" for i in range(30))).encode()
        headers, rows = preview_csv(content)
        self.assertEqual(headers, ["id", "smiles"])
        self.assertEqual(len(rows), 20)

    def test_sample_folder_contains_no_compound_fixtures(self):
        project_root = Path(__file__).resolve().parents[1]
        sample_dir = project_root / "data" / "sample"
        self.assertEqual(list(sample_dir.glob("*.csv")), [])
        self.assertTrue((sample_dir / "README.md").is_file())


if __name__ == "__main__":
    unittest.main()
