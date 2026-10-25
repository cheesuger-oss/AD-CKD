import gzip
import tempfile
import unittest
from pathlib import Path

from phase2_extract_gwas_rows import extract_rows, extract_rows_by_positions


class ExtractRowsTests(unittest.TestCase):
    """Catches matching the wrong SNP column or losing the source header."""

    def test_extracts_requested_rsids_and_preserves_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.tsv.gz"
            wanted = root / "wanted.txt"
            output = root / "subset.tsv"
            with gzip.open(source, "wt") as handle:
                handle.write("chr\tpos\trsids\tbeta\n")
                handle.write("1\t10\trs1\t0.1\n")
                handle.write("1\t20\trs2\t0.2\n")
                handle.write("1\t30\trs3,rs33\t0.3\n")
            wanted.write_text("rs2\nrs33\n")

            summary = extract_rows(source, wanted, output, "rsids")

            self.assertEqual(summary["rows_scanned"], 3)
            self.assertEqual(summary["rows_matched"], 2)
            self.assertEqual(
                output.read_text().splitlines(),
                ["chr\tpos\trsids\tbeta", "1\t20\trs2\t0.2", "1\t30\trs3,rs33\t0.3"],
            )

    def test_rejects_missing_snp_column(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.tsv.gz"
            wanted = root / "wanted.txt"
            with gzip.open(source, "wt") as handle:
                handle.write("chr\tpos\tbeta\n1\t10\t0.1\n")
            wanted.write_text("rs1\n")

            with self.assertRaisesRegex(ValueError, "SNP column"):
                extract_rows(source, wanted, root / "out.tsv", "rsids")

    def test_extracts_requested_chromosome_positions(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.tsv.gz"
            wanted = root / "positions.tsv"
            output = root / "subset.tsv"
            with gzip.open(source, "wt") as handle:
                handle.write("chr\tpos\tref\talt\tbeta\n")
                handle.write("1\t10\tA\tG\t0.1\n")
                handle.write("2\t20\tC\tT\t0.2\n")
            wanted.write_text("chr\tpos\n2\t20\n")

            summary = extract_rows_by_positions(source, wanted, output, "chr", "pos")

            self.assertEqual(summary["rows_matched"], 1)
            self.assertEqual(output.read_text().splitlines()[-1], "2\t20\tC\tT\t0.2")


if __name__ == "__main__":
    unittest.main()
