import unittest
from numpy import array
from pyquante2 import pgbf,cgbf
from pyquante2.cints.one import S,T,V
from pyquante2.cints.two import ERI
from pyquante2.cints.two import ERI as ERI_hgp
from pyquante2.cints.hgp import vrr,vrr_recursive
from pyquante2.ints.hgp import vrr as pyvrr
s = pgbf(1)
s2 = cgbf(exps=[1],coefs=[1])
s3 = cgbf((0,0,1),(0,0,0),[1],[1])

class test_cython(unittest.TestCase):
    def test_S(self):
        self.assertAlmostEqual(S(s,s),1)
    def test_T(self):
        self.assertAlmostEqual(T(s,s),1.5)
    def test_V(self):
        self.assertAlmostEqual(V(s,s,(0,0,0)),-1.5957691216)
    def test_ERI(self):
        self.assertAlmostEqual(ERI(s,s,s,s),1.1283791671)
        self.assertAlmostEqual(ERI(s2,s2,s2,s2),1.1283791671)
        self.assertAlmostEqual(ERI(s2,s2,s3,s3),0.84270079)
    def test_ERI_hgp(self):
        self.assertAlmostEqual(ERI_hgp(s,s,s,s),1.1283791671)
        self.assertAlmostEqual(ERI_hgp(s2,s2,s2,s2),1.1283791671)
        self.assertAlmostEqual(ERI_hgp(s2,s2,s3,s3),0.84270079)

    def test_vrr(self):
        zero = array([0,0,0],'d')
        pyval = pyvrr(zero,1.,(0,0,0),1.,
                      zero,1.,1.,
                      zero,1.,(0,0,0),1.,
                      zero,1.,1.,
                      0)
        cval = vrr(0,0,0,1.,0,0,0,1.,
                   0,0,0,1.,1.,
                   0,0,0,1.,0,0,0,1.,
                   0,0,0,1.,1.,
                   0)
        cval2 = vrr_recursive(0,0,0,1.,0,0,0,1.,
                              0,0,0,1.,1.,
                              0,0,0,1.,0,0,0,1.,
                              0,0,0,1.,1.,
                              0)
        self.assertAlmostEqual(cval2,pyval)
        self.assertAlmostEqual(cval,pyval)

if __name__ == '__main__':
    # Manually make the test suite
    testSuite = unittest.TestSuite()
    testSuite.addTests(unittest.makeSuite(test_cython))

    # Run the suite
    unittest.TextTestRunner(verbosity=2).run(testSuite)
    
