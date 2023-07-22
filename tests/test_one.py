import unittest
from pyquante2.basis.pgbf import pgbf
from pyquante2.ints.one import S,T,V
from pyquante2.cints.one import S as cS
from pyquante2.cints.one import T as cT
from pyquante2.cints.one import V as cV

class test_one(unittest.TestCase):
    def test_pgbf(self):
        s = pgbf(1.0)
        self.assertAlmostEqual(s(0,0,0),0.7127054703549901)

    def test_overlap(self):
        s = pgbf(1.0)
        self.assertAlmostEqual(S(s,s),1.0)
        self.assertAlmostEqual(cS(s,s),1.0)

    def test_kinetic(self):
        s = pgbf(1.0)
        self.assertAlmostEqual(T(s,s),1.5)
        self.assertAlmostEqual(cT(s,s),1.5)

    def test_nuclear(self):
        s = pgbf(1.0)
        self.assertAlmostEqual(V(s,s,((0,0,0))),-1.5957691216057328)
        self.assertAlmostEqual(cV(s,s,((0,0,0))),-1.5957691216057328)

    def test_h_atom(self):
        "Single primitive approximation to H atom"
        s = pgbf(0.285)
        self.assertAlmostEqual(T(s,s)+V(s,s,((0,0,0))),-0.424407589178)
        self.assertAlmostEqual(cT(s,s)+cV(s,s,((0,0,0))),-0.424407589178)
