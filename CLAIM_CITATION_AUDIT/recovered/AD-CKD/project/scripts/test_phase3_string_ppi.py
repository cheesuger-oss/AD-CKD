#!/usr/bin/env python3

import tempfile
import unittest
from pathlib import Path

from phase3_string_ppi import parse_tsv, write_tsv


class StringPpiHelpersTest(unittest.TestCase):
    def test_parse_tsv(self):
        fields, rows = parse_tsv(b"queryItem\tstringId\nTP53\t9606.ENSP1\n")
        self.assertEqual(fields, ["queryItem", "stringId"])
        self.assertEqual(rows, [{"queryItem": "TP53", "stringId": "9606.ENSP1"}])

    def test_write_tsv_is_deterministic(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "x.tsv"
            write_tsv(path, ["a", "b"], [{"a": "1", "b": "2"}])
            self.assertEqual(path.read_bytes(), b"a\tb\n1\t2\n")


if __name__ == "__main__":
    unittest.main()
