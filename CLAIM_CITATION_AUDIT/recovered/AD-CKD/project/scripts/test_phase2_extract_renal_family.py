#!/usr/bin/env python3

import csv
import gzip
import math
import tempfile
import unittest
from pathlib import Path

from phase2_extract_renal_family import extract


class RenalFamilyExtractionTest(unittest.TestCase):
    def test_whitespace_or_is_converted_to_log_odds(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "x.gz"
            output = root / "out.tsv"
            with gzip.open(source, "wt") as handle:
                handle.write("Chr Pos RSID A1 A2 Freq OR SE P N\n")
                handle.write("1 10 rs1 A G 0.4 2 0.1 1e-9 1000\n")
                handle.write("1 11 rs2 C T 0.2 1.1 0.2 1e-3 1000\n")
            row = {
                "exposure_id": "x", "trait": "X", "scale": "odds_ratio", "status": "VERIFIED",
                "source_file": str(source), "delimiter": "whitespace", "chr_col": "Chr", "pos_col": "Pos",
                "snp_col": "RSID", "ea_col": "A1", "oa_col": "A2", "eaf_col": "Freq",
                "effect_col": "OR", "se_col": "SE", "p_col": "P", "n_col": "N",
            }
            summary = extract(row, output)
            self.assertEqual(summary["unique_significant_rsids"], 1)
            with output.open(newline="") as handle:
                result = next(csv.DictReader(handle, delimiter="\t"))
            self.assertAlmostEqual(float(result["beta"]), math.log(2))
            self.assertEqual(result["SNP"], "rs1")


if __name__ == "__main__":
    unittest.main()
