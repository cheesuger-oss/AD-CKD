#!/usr/bin/env python3

import unittest

from phase3_dgidb_query import flatten


class DgidbTest(unittest.TestCase):
    def test_flatten_keeps_zero_interaction_gene(self):
        payload = {"data": {"genes": {"nodes": [{"name": "A", "interactions": []}]}}}
        audit, interactions = flatten(payload, ["A", "B"])
        self.assertEqual([x["interaction_count"] for x in audit], [0, 0])
        self.assertEqual([x["dgidb_mapping_status"] for x in audit], ["mapped", "unmapped"])
        self.assertEqual(interactions, [])


if __name__ == "__main__":
    unittest.main()
