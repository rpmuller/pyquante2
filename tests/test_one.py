import unittest
from pyquante2.basis.pgbf import pgbf
from pyquante2.ints.one import S,T,V

class test_one(unittest.TestCase):
    def test_pgbf(self):
        s = pgbf(1.0)
        self.assertAlmostEqual(s(0,0,0),0.7127054703549901)

    def test_overlap(self):
        s = pgbf(1.0)
        self.assertAlmostEqual(S(s,s),1.0)

    def test_kinetic(self):
        s = pgbf(1.0)
        self.assertAlmostEqual(T(s,s),1.5)

    def test_nuclear(self):
        s = pgbf(1.0)
        self.assertAlmostEqual(V(s,s,((0,0,0))),-1.5957691216057328)

    def test_h_atom(self):
        "Single primitive approximation to H atom"
        s = pgbf(0.285)
        self.assertAlmostEqual(T(s,s)+V(s,s,((0,0,0))),-0.424407589178)
