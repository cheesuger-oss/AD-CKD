import csv
import math
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "geo" / "human_validation"
EXPECTED_GO = {
    "GO:0030198", "GO:0001503", "GO:0006939",
    "GO:0048659", "GO:0048660", "GO:0048661",
}


class HumanExpressionValidationOutputs(unittest.TestCase):
    def test_frozen_family_is_complete_for_all_datasets(self):
        with (OUT / "human_validation_six_pathway_GSEA.csv").open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 18)
        for dataset in ("GSE153434", "GSE38750", "GSE38751"):
            observed = {r["GO_ID"] for r in rows if r["dataset"] == dataset}
            self.assertEqual(observed, EXPECTED_GO)

    def test_probabilities_and_effects_are_finite(self):
        with (OUT / "human_validation_six_pathway_GSEA.csv").open() as handle:
            rows = list(csv.DictReader(handle))
        for row in rows:
            for key in ("ES", "NES", "P_value", "BH_q"):
                self.assertTrue(math.isfinite(float(row[key])))
            self.assertGreaterEqual(float(row["P_value"]), 0)
            self.assertLessEqual(float(row["P_value"]), 1)
            self.assertGreaterEqual(float(row["BH_q"]), float(row["P_value"]))
            self.assertLessEqual(float(row["BH_q"]), 1)

    def test_sample_counts_and_required_figures(self):
        expected = {"GSE153434": (10, 10), "GSE38750": (19, 15), "GSE38751": (5, 5)}
        for accession, (controls, cases) in expected.items():
            with (OUT / f"{accession}_sample_metadata.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(sum(r["group"] == "Control" for r in rows), controls)
            self.assertEqual(sum(r["group"] == "Case" for r in rows), cases)
            for suffix in ("PCA.png", "sample_correlation.png", "volcano.png", "top40_heatmap.png", "six_pathway_GSEA.png"):
                path = OUT / "figures" / f"{accession}_{suffix}"
                self.assertTrue(path.exists() and path.stat().st_size > 1000, path)


if __name__ == "__main__":
    unittest.main()
