import unittest

from extract_coloc_regions import WindowIndex, harmonize, parse_values


class RegionalExtractionTests(unittest.TestCase):
    def test_parse_values_rejects_collapsed_missing_whitespace_column(self):
        header = "Chr Pos_b37 RSID Allele1 Allele2 Freq1 Effect StdErr P-value n_total_sum"
        malformed = "1 2556125  t c 0.3260 -0.0035140 0.0021112 0.09603 547349"
        self.assertEqual(len(header.split()), 10)
        self.assertIsNone(parse_values(malformed, None, 10))

    def test_window_index_returns_nonoverlapping_window(self):
        idx = WindowIndex([("1", 100, 200, "w1"), ("1", 300, 400, "w2")])
        self.assertEqual(idx.find("1", 150), "w1")
        self.assertEqual(idx.find("1", 400), "w2")
        self.assertIsNone(idx.find("1", 250))

    def test_harmonize_flips_swapped_aad_effect(self):
        renal = {"rsid": "rs1", "a1": "A", "a2": "G", "beta": "0.2", "se": "0.1", "p": "0.05", "freq": "0.3", "n": "1000"}
        aad = {"a1": "G", "a2": "A", "beta": "0.4", "se": "0.2", "p": "0.04", "freq": "0.7"}
        out = harmonize(renal, aad)
        self.assertAlmostEqual(out["aad_beta"], -0.4)
        self.assertAlmostEqual(out["aad_freq"], 0.3)

    def test_harmonize_rejects_incompatible_pair(self):
        renal = {"rsid": "rs1", "a1": "A", "a2": "G", "beta": "0.2", "se": "0.1", "p": "0.05", "freq": "0.3", "n": "1000"}
        aad = {"a1": "A", "a2": "C", "beta": "0.4", "se": "0.2", "p": "0.04", "freq": "0.3"}
        self.assertIsNone(harmonize(renal, aad))


if __name__ == "__main__":
    unittest.main()
