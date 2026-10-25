import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).with_name("phase4_external_aad_audit.py")


class ExternalAADAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location("phase4_external_audit", SCRIPT)
        cls.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.module)

    def test_locked_local_cohorts_have_expected_identity_and_status(self):
        root = Path(__file__).parents[2]
        raw = root / "project/data/geo/phase4/raw"
        with tempfile.TemporaryDirectory() as tmp:
            rows = self.module.audit_external_cohorts(raw, Path(tmp) / "audit.tsv", Path(tmp) / "audit.md")
        by_id = {row["dataset_id"]: row for row in rows}
        self.assertEqual((by_id["GSE52093"]["case_n"], by_id["GSE52093"]["control_n"]), ("7", "5"))
        self.assertEqual((by_id["GSE98770"]["case_n"], by_id["GSE98770"]["control_n"]), ("6", "5"))
        self.assertEqual((by_id["GSE190635"]["case_n"], by_id["GSE190635"]["control_n"]), ("4", "4"))
        self.assertEqual(by_id["GSE147026"]["expression_rows"], "0")
        self.assertEqual(by_id["GSE147026"]["fitness"], "NOT_READY")
        self.assertEqual(by_id["GSE52093"]["execution_order"], "1")
        self.assertEqual(by_id["GSE98770"]["execution_order"], "2")


if __name__ == "__main__":
    unittest.main()
