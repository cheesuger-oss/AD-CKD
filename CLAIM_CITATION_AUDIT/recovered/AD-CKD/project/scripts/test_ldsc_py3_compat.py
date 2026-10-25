import gzip
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


LDSC_ROOT = Path(__file__).resolve().parents[1] / "tools" / "ldsc"
sys.path.insert(0, str(LDSC_ROOT))
SPEC = importlib.util.spec_from_file_location("vendored_munge_sumstats", LDSC_ROOT / "munge_sumstats.py")
MUNGE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MUNGE)


class TestLdscPython3Compatibility(unittest.TestCase):
    def test_read_header_decodes_gzip_text(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tiny.tsv.gz"
            with gzip.open(path, "wt") as handle:
                handle.write("SNP\tA1\tA2\tP\nrs1\tA\tG\t0.5\n")
            self.assertEqual(MUNGE.read_header(str(path)), ["SNP", "A1", "A2", "P"])


if __name__ == "__main__":
    unittest.main()
