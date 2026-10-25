import unittest

from prepare_coloc_loci import normalize_rsid, parse_p, merge_windows


class ColocLocusTests(unittest.TestCase):
    def test_normalize_rsid_rejects_multi_and_non_rs(self):
        self.assertEqual(normalize_rsid("rs123"), "rs123")
        self.assertIsNone(normalize_rsid("rs1,rs2"))
        self.assertIsNone(normalize_rsid("1:100:A:G"))

    def test_parse_p_allows_underflow_zero_but_rejects_invalid(self):
        self.assertEqual(parse_p("0"), 0.0)
        self.assertIsNone(parse_p("NA"))
        self.assertIsNone(parse_p("1.2"))

    def test_merge_windows_merges_touching_intervals(self):
        rows = [("1", 100, "a"), ("1", 900, "b"), ("2", 500, "c")]
        merged = merge_windows(rows, flank=500)
        self.assertEqual(merged[0][:3], ("1", 1, 1400))
        self.assertEqual(merged[0][3], ["a", "b"])
        self.assertEqual(merged[1], ("2", 1, 1000, ["c"]))


if __name__ == "__main__":
    unittest.main()
