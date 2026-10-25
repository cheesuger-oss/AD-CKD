import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).with_name("phase4_freeze_immune_programs.py")


class ImmuneProgramFreezeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location("phase4_immune", SCRIPT)
        cls.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.module)

    def test_six_go_programs_freeze_with_expected_counts(self):
        root = Path(__file__).parents[2]
        db = root / "project/Rlib-4.5/org.Hs.eg.db/extdata/org.Hs.eg.sqlite"
        with tempfile.TemporaryDirectory() as tmp:
            rows, summary = self.module.freeze_programs(db, Path(tmp) / "programs.tsv", Path(tmp) / "audit.md")
        self.assertEqual(len(summary), 6)
        self.assertEqual({item["go_id"]: item["gene_n"] for item in summary}, {
            "GO:0042116": 115,
            "GO:0042119": 48,
            "GO:0042110": 607,
            "GO:0042113": 291,
            "GO:0030101": 107,
            "GO:0036336": 33,
        })
        self.assertEqual(len(rows), len({(r["program_id"], r["human_gene"]) for r in rows}))
        self.assertTrue(all(r["interpretation"] == "activation_program_not_cell_fraction" for r in rows))


if __name__ == "__main__":
    unittest.main()
