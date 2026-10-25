#!/usr/bin/env python3

import unittest

from phase35_virtual_perturbation import components, distances, remove_nodes


class PerturbationTest(unittest.TestCase):
    def setUp(self):
        self.graph = {"A": {"B"}, "B": {"A", "C"}, "C": {"B"}, "D": set()}

    def test_components(self):
        self.assertEqual(sorted(map(len, components(self.graph))), [1, 3])

    def test_distances(self):
        self.assertEqual(distances(self.graph, "A"), {"A": 0, "B": 1, "C": 2})

    def test_remove(self):
        self.assertEqual(remove_nodes(self.graph, {"B"}), {"A": set(), "C": set(), "D": set()})


if __name__ == "__main__":
    unittest.main()
