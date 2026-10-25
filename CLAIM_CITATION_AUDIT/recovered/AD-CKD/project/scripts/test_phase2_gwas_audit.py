import gzip
import importlib.util
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("phase2_gwas_audit.py")


def load_module():
    spec = importlib.util.spec_from_file_location("phase2_gwas_audit", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Phase2AuditTests(unittest.TestCase):
    def test_missing_primary_files_block_mr(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = module.audit_inputs(root, root)
        self.assertEqual(result["mr_ready"], False)
        self.assertEqual(
            set(result["missing_required"]),
            {"GCST90103634_eGFR", "FinnGen_R13_I9_AORTDIS"},
        )

    def test_sbp_requires_eur_specific_effect_columns(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sbp = root / "continuous-4080-both_sexes-irnt.tsv.bgz"
            with gzip.open(sbp, "wt") as handle:
                handle.write("chr\tpos\tref\talt\tbeta_EUR\tse_EUR\taf_EUR\n")
                handle.write("1\t100\tA\tG\t0.1\t0.02\t0.4\n")
            result = module.audit_inputs(root, root)
        self.assertEqual(result["inputs"]["UKB_4080_SBP"]["field_status"], "PASS")
        self.assertEqual(result["inputs"]["UKB_4080_SBP"]["effect_allele"], "alt")

    def test_unknown_gwas_stays_excluded(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            unknown = root / "GCST90469223.tsv.gz"
            with gzip.open(unknown, "wt") as handle:
                handle.write(
                    "chromosome\tbase_pair_location\teffect_allele\tother_allele\t"
                    "beta\tstandard_error\tp_value\n"
                )
            result = module.audit_inputs(root, root)
        self.assertEqual(result["inputs"]["GCST90469223"]["eligible"], False)
        self.assertEqual(result["inputs"]["GCST90469223"]["trait"], "UNKNOWN")

    def test_accepts_real_ckdgen_and_finngen_header_spellings(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with gzip.open(root / "metal_eGFR_meta_ea1.TBL.map.annot.gc.gz", "wt") as handle:
                handle.write(
                    "MarkerName\tAllele1\tAllele2\tn\tFreq1\tEffect\tStdErr\tP.value\t"
                    "Direction\tmac\tchr\tpos\tP.value.GC\tStdErr.GC\tRSID\n"
                )
            with gzip.open(root / "finngen_R13_I9_AORTDIS.gz", "wt") as handle:
                handle.write(
                    "#chrom\tpos\tref\talt\trsids\tnearest_genes\tpval\tmlogp\tbeta\t"
                    "sebeta\taf_alt\taf_alt_cases\taf_alt_controls\n"
                )
            result = module.audit_inputs(root, root)
        self.assertEqual(result["inputs"]["GCST90103634_eGFR"]["field_status"], "PASS")
        self.assertEqual(result["inputs"]["FinnGen_R13_I9_AORTDIS"]["field_status"], "PASS")
        self.assertEqual(result["mr_ready"], True)


if __name__ == "__main__":
    unittest.main()
