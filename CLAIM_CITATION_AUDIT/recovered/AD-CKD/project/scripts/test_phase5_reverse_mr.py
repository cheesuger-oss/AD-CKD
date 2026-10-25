import unittest

from phase5_reverse_mr import ivw, bh_adjust


class ReverseMRTests(unittest.TestCase):
    def test_ivw_known_equal_ratios(self):
        b,se,q,qp=ivw([1.0,2.0],[2.0,4.0],[1.0,1.0])
        self.assertAlmostEqual(b,2.0); self.assertAlmostEqual(se,(1/5)**0.5); self.assertAlmostEqual(q,0.0)

    def test_bh_is_monotone_in_rank(self):
        q=bh_adjust([0.01,0.04,0.03])
        self.assertEqual([round(x,3) for x in q],[0.03,0.04,0.04])


if __name__=="__main__": unittest.main()
