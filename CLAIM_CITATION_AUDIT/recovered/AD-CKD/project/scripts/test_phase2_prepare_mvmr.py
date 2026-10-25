import tempfile
import unittest
from pathlib import Path

from phase2_prepare_mvmr import build_joint_clump_input, export_clumped_targets, map_sbp_to_reference


class MapSbpToReferenceTests(unittest.TestCase):
    """Catches coordinate-only rsID mapping and loss of strand complements."""

    def test_requires_allele_compatible_reference_variant(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bim = root / "eur.bim"
            sbp = root / "sbp.tsv"
            out = root / "mapped.tsv"
            bim.write_text("1\trs_good\t0\t100\tA\tG\n1\trs_wrong\t0\t200\tA\tC\n")
            sbp.write_text(
                "chr\tpos\tref\talt\taf_EUR\tbeta_EUR\tse_EUR\tneglog10_pval_EUR\n"
                "1\t100\tA\tG\t0.2\t0.1\t0.01\t9\n"
                "1\t200\tA\tG\t0.3\t0.2\t0.02\t10\n"
            )

            summary = map_sbp_to_reference(sbp, bim, out)

            self.assertEqual(summary, {"input_rows": 2, "mapped_rows": 1, "ambiguous_rows": 0})
            self.assertEqual(out.read_text().splitlines()[1].split("\t")[-1], "rs_good")

    def test_accepts_unique_strand_complement(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bim = root / "eur.bim"
            sbp = root / "sbp.tsv"
            out = root / "mapped.tsv"
            bim.write_text("2\trs_comp\t0\t300\tT\tC\n")
            sbp.write_text(
                "chr\tpos\tref\talt\taf_EUR\tbeta_EUR\tse_EUR\tneglog10_pval_EUR\n"
                "2\t300\tA\tG\t0.4\t-0.2\t0.03\t12\n"
            )

            summary = map_sbp_to_reference(sbp, bim, out)

            self.assertEqual(summary["mapped_rows"], 1)
            self.assertEqual(out.read_text().splitlines()[1].split("\t")[-1], "rs_comp")


class BuildJointClumpInputTests(unittest.TestCase):
    """Catches exposure-specific clumping instead of clumping the instrument union."""

    def test_uses_best_p_value_once_per_rsid(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            egfr = root / "egfr.tsv"
            sbp = root / "sbp.tsv"
            out = root / "joint.tsv"
            egfr.write_text(
                "chr\tpos\tP.value.GC\tRSID\n"
                "1\t100\t1e-9\trs_shared\n"
                "1\t200\t2e-8\trs_egfr\n"
            )
            sbp.write_text(
                "chr\tpos\tneglog10_pval_EUR\tRSID\n"
                "1\t100\t12\trs_shared\n"
                "2\t300\t10\trs_sbp\n"
            )

            summary = build_joint_clump_input(egfr, sbp, out)

            self.assertEqual(summary, {"egfr_candidates": 2, "sbp_candidates": 2, "union_variants": 3})
            rows = {line.split("\t")[0]: line.split("\t") for line in out.read_text().splitlines()[1:]}
            self.assertEqual(set(rows), {"rs_shared", "rs_egfr", "rs_sbp"})
            self.assertEqual(float(rows["rs_shared"][1]), 1e-12)

    def test_exports_clumped_ids_and_coordinates(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            clumped = root / "joint.clumped"
            joint = root / "joint.tsv"
            ids = root / "ids.txt"
            positions = root / "positions.tsv"
            clumped.write_text(" CHR F SNP BP P TOTAL\n 1 1 rs2 200 1e-10 3\n 2 1 rs3 300 2e-9 2\n")
            joint.write_text("SNP\tP\tCHR\tBP\nrs2\t1e-10\t1\t200\nrs3\t2e-9\t2\t300\n")

            count = export_clumped_targets(clumped, joint, ids, positions)

            self.assertEqual(count, 2)
            self.assertEqual(ids.read_text(), "rs2\nrs3\n")
            self.assertEqual(positions.read_text(), "chr\tpos\n1\t200\n2\t300\n")


if __name__ == "__main__":
    unittest.main()
