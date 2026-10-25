import csv
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[2]
OUT = ROOT / "project/results/phase4/analysis"


class Phase4AnalysisOutputTests(unittest.TestCase):
    def test_immune_family_is_complete_for_primary_datasets(self):
        with (OUT / "immune_program_contrasts.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        by_dataset = {}
        for row in rows:
            by_dataset.setdefault(row["dataset"], []).append(row)
            self.assertGreaterEqual(float(row["observed_gene_n"]), 10)
            self.assertGreaterEqual(float(row["q_value"]), 0)
            self.assertLessEqual(float(row["q_value"]), 1)
        self.assertEqual(len(by_dataset["GSE146638"]), 6)
        self.assertEqual(len(by_dataset["GSE153434"]), 6)
        self.assertEqual(len(by_dataset["GSE38750"]), 6)

    def test_primary_oof_predictions_cover_all_samples_once(self):
        with (OUT / "primary_oof_predictions.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 20)
        self.assertEqual(len({row["sample"] for row in rows}), 20)
        self.assertEqual({row["group"] for row in rows}, {"Control", "Case"})
        self.assertTrue(all(0 <= float(row["probability"]) <= 1 for row in rows))

    def test_model_summary_and_external_cohort_accounting(self):
        summary = json.loads((OUT / "model_summary.json").read_text())
        self.assertEqual(summary["permutations"], 1000)
        self.assertEqual(summary["stability_resamples"], 500)
        self.assertEqual(summary["oof_prediction_n"], 20)
        with (OUT / "external_validation_results.tsv").open() as handle:
            external = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([row["dataset"] for row in external], ["GSE52093", "GSE98770", "GSE147026", "GSE190635"])
        self.assertEqual(next(row for row in external if row["dataset"] == "GSE147026")["status"], "NOT_READY")

    def test_required_figures_exist(self):
        required = [
            "immune_program_heatmap.png",
            "immune_effect_forest.png",
            "primary_oof_roc.png",
            "primary_oof_pr.png",
            "primary_calibration.png",
            "feature_stability.png",
            "locked_signature_coefficients.png",
            "external_validation_forest.png",
        ]
        for name in required:
            path = OUT / "figures" / name
            self.assertTrue(path.exists(), name)
            self.assertGreater(path.stat().st_size, 10_000, name)


if __name__ == "__main__":
    unittest.main()
