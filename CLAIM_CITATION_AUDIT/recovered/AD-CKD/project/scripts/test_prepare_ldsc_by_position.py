import importlib.util
from pathlib import Path
import unittest


SCRIPT = Path(__file__).with_name("prepare_ldsc_by_position.py")
SPEC = importlib.util.spec_from_file_location("prepare_ldsc_by_position", SCRIPT)
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


class PositionMatcherTests(unittest.TestCase):
    def test_match_direct_and_swapped_alleles(self):
        candidates = [("rs1", "A", "G"), ("rs2", "C", "T")]
        self.assertEqual(MOD.choose_match(candidates, "A", "G"), ("rs1", "A", "G"))
        self.assertEqual(MOD.choose_match(candidates, "G", "A"), ("rs1", "G", "A"))

    def test_rejects_incompatible_alleles(self):
        self.assertIsNone(MOD.choose_match([("rs1", "A", "G")], "A", "C"))

    def test_rejects_ambiguous_multiple_matches(self):
        candidates = [("rs1", "A", "G"), ("rs2", "A", "G")]
        self.assertIsNone(MOD.choose_match(candidates, "A", "G"))


if __name__ == "__main__":
    unittest.main()
