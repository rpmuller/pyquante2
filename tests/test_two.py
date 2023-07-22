import unittest
from pyquante2.basis.pgbf import pgbf
from pyquante2.ints.two import ERI as coul_huz
from pyquante2.ints.rys import ERI as coul_rys
from pyquante2.ints.hgp import ERI as coul_hgp
from pyquante2.cints.two import ERI as coul_chuz
from pyquante2.cints.hgp import ERI as coul_chgp

s = pgbf(1.0)
px = pgbf(1.0,(0,0,0),(1,0,0))

class test_two(unittest.TestCase):
    def test_huz(self):
        self.assertAlmostEqual(coul_huz(s,s,px,px),0.9403159725793302)
        self.assertEqual(coul_huz(s,s,s,px),0)

    def test_rys(self):
        self.assertAlmostEqual(coul_rys(s,s,px,px),0.9403159725793302)
        self.assertEqual(coul_rys(s,s,s,px),0)

    def test_hgp(self):
        self.assertAlmostEqual(coul_hgp(s,s,px,px),0.9403159725793302)
        self.assertEqual(coul_hgp(s,s,s,px),0)

    def test_chuz(self):
        self.assertAlmostEqual(coul_chuz(s,s,px,px),0.9403159725793302)
        self.assertEqual(coul_chuz(s,s,s,px),0)

    def test_chgp(self):
        self.assertAlmostEqual(coul_chgp(s,s,px,px),0.9403159725793302)
        self.assertEqual(coul_chgp(s,s,s,px),0)

    def test_huz_range(self):
        for r,val in [(0.0, 1.1283791670951362),
                      (1.0, 0.8427007900292186),
                      (2.0, 0.49766113257563993),
                      (3.0, 0.33332596983445223),
                      (4.0, 0.2499999961456855), # coulomb's law hereafter:
                      (5.0, 1/5), 
                      (6.0, 1/6),
                      (7.0, 1/7),
                      (8.0, 1/8),
                      (9.0, 1/9)]:
            s3 = pgbf(1.0, (0.0,0.0,r))
            self.assertAlmostEqual(coul_huz(s,s,s3,s3),val)
        return

    def test_hgp_range(self):
        for r,val in [(0.0, 1.1283791670951362),
                      (1.0, 0.8427007900292186),
                      (2.0, 0.49766113257563993),
                      (3.0, 0.33332596983445223),
                      (4.0, 0.2499999961456855), # coulomb's law hereafter:
                      (5.0, 1/5), 
                      (6.0, 1/6),
                      (7.0, 1/7),
                      (8.0, 1/8),
                      (9.0, 1/9)]:
            s3 = pgbf(1.0, (0.0,0.0,r))
            self.assertAlmostEqual(coul_hgp(s,s,s3,s3),val)
        return

    def test_chuz_range(self):
        for r,val in [(0.0, 1.1283791670951362),
                      (1.0, 0.8427007900292186),
                      (2.0, 0.49766113257563993),
                      (3.0, 0.33332596983445223),
                      (4.0, 0.2499999961456855), # coulomb's law hereafter:
                      (5.0, 1/5), 
                      (6.0, 1/6),
                      (7.0, 1/7),
                      (8.0, 1/8),
                      (9.0, 1/9)]:
            s3 = pgbf(1.0, (0.0,0.0,r))
            self.assertAlmostEqual(coul_chuz(s,s,s3,s3),val)
        return

    def test_chgp_range(self):
        for r,val in [(0.0, 1.1283791670951362),
                      (1.0, 0.8427007900292186),
                      (2.0, 0.49766113257563993),
                      (3.0, 0.33332596983445223),
                      (4.0, 0.2499999961456855), # coulomb's law hereafter:
                      (5.0, 1/5), 
                      (6.0, 1/6),
                      (7.0, 1/7),
                      (8.0, 1/8),
                      (9.0, 1/9)]:
            s3 = pgbf(1.0, (0.0,0.0,r))
            self.assertAlmostEqual(coul_chgp(s,s,s3,s3),val)
        return
