#!/usr/bin/env python3

import unittest

from phase3_network_topology import brandes_betweenness


class TopologyTest(unittest.TestCase):
    def test_three_node_path(self):
        graph = {"A": {"B"}, "B": {"A", "C"}, "C": {"B"}}
        score = brandes_betweenness(graph)
        self.assertAlmostEqual(score["A"], 0.0)
        self.assertAlmostEqual(score["B"], 1.0)
        self.assertAlmostEqual(score["C"], 0.0)

    def test_triangle_has_zero_betweenness(self):
        graph = {"A": {"B", "C"}, "B": {"A", "C"}, "C": {"A", "B"}}
        self.assertEqual(brandes_betweenness(graph), {"A": 0.0, "B": 0.0, "C": 0.0})


if __name__ == "__main__":
    unittest.main()
