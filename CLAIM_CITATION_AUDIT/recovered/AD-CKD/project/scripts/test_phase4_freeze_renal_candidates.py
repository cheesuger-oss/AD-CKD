import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).with_name("phase4_freeze_renal_candidates.py")


class Phase4RenalCandidateFreezeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location("phase4_freeze", SCRIPT)
        cls.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.module)

    def test_real_inputs_produce_renal_only_unique_candidates(self):
        root = Path(__file__).parents[2]
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "candidates.tsv"
            audit = Path(tmp) / "audit.md"
            rows = self.module.freeze_candidates(
                root / "project/results/geo/GSE146638_DEG.csv",
                root / "project/results/geo/ortholog_rat_human_1to1.csv",
                out,
                audit,
            )
            self.assertEqual(len(rows), 1285)
            self.assertEqual(len({r["human_gene"] for r in rows}), 1285)
            self.assertTrue(all(r["candidate_origin"] == "renal_only" for r in rows))
            self.assertTrue(all(r["aad_label_used"] == "false" for r in rows))
            self.assertTrue(all(float(r["renal_adj_p"]) < 0.05 for r in rows))
            self.assertTrue(all(abs(float(r["renal_log_fc"])) > 0.5 for r in rows))

    def test_conflicting_duplicate_human_gene_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            deg = tmp / "deg.csv"
            mapping = tmp / "map.csv"
            with deg.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["gene_symbol", "logFC", "adj.P.Val"])
                writer.writeheader()
                writer.writerows([
                    {"gene_symbol": "RatA", "logFC": "1", "adj.P.Val": "0.01"},
                    {"gene_symbol": "RatB", "logFC": "-1", "adj.P.Val": "0.01"},
                ])
            with mapping.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["source_symbol", "human_symbol", "orthology_type", "orthology_confidence"])
                writer.writeheader()
                writer.writerows([
                    {"source_symbol": "RatA", "human_symbol": "HUMAN1", "orthology_type": "ortholog_one2one", "orthology_confidence": "1"},
                    {"source_symbol": "RatB", "human_symbol": "HUMAN1", "orthology_type": "ortholog_one2one", "orthology_confidence": "1"},
                ])
            with self.assertRaises(ValueError):
                self.module.freeze_candidates(deg, mapping, tmp / "out.tsv", tmp / "audit.md")


if __name__ == "__main__":
    unittest.main()
