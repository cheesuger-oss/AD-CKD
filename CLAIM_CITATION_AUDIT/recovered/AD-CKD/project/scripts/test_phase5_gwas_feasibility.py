import io
import unittest

from phase5_gwas_feasibility import audit_stream


class TestGwasFeasibility(unittest.TestCase):
    def test_audits_tab_delimited_finngen_rows(self):
        data = io.StringIO(
            "#chrom\tpos\tref\talt\trsids\tpval\tbeta\tsebeta\taf_alt\n"
            "1\t10\tA\tG\trs1\t1e-9\t0.2\t0.1\t0.4\n"
            "1\t20\tC\tT\trs2\t0.2\t-0.1\t0.2\t0.3\n"
        )
        result = audit_stream(data, delimiter="tab")
        self.assertEqual(result["variant_rows"], 2)
        self.assertEqual(result["gws_rows"], 1)
        self.assertEqual(result["rsid_rows"], 2)
        self.assertEqual(result["min_p"], 1e-9)

    def test_audits_whitespace_and_recognizes_p_value_alias(self):
        data = io.StringIO(
            "Chr Pos_b37 RSID Allele1 Allele2 Effect StdErr P-value\n"
            "1 10 rs1 A G 0.2 0.1 5e-8\n"
            "1 20 . C T -0.1 0.2 NA\n"
        )
        result = audit_stream(data, delimiter="whitespace")
        self.assertEqual(result["variant_rows"], 2)
        self.assertEqual(result["valid_p_rows"], 1)
        self.assertEqual(result["gws_rows"], 0)
        self.assertEqual(result["rsid_rows"], 1)

    def test_missing_p_column_is_explicit(self):
        data = io.StringIO("SNP\tbeta\tse\nrs1\t0.2\t0.1\n")
        result = audit_stream(data, delimiter="tab")
        self.assertIsNone(result["p_column"])
        self.assertEqual(result["valid_p_rows"], 0)

    def test_recognizes_harmonised_rsid_column(self):
        data = io.StringIO(
            "hm_variant_id\thm_rsid\tp_value\n1_10_A_G\trs1\t0.01\n"
        )
        result = audit_stream(data, delimiter="tab")
        self.assertEqual(result["rsid_column"], "hm_rsid")
        self.assertEqual(result["rsid_rows"], 1)


if __name__ == "__main__":
    unittest.main()
