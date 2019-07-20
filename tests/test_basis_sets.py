import unittest

import numpy as np
from pyquante2 import basisset,rhf,h2

class test_basissets(unittest.TestCase):
    def test_basissets(self):
        ref = basisset(h2,'sto3g')

        test1 = basisset(h2,'sto-3g')
        test2 = basisset(h2,'sto-3G')
        test3 = basisset(h2,'StO-3g')

        dif = basisset(h2,'6-31g')

        self.assertEqual(str(ref), str(test1))
        self.assertEqual(str(ref), str(test2))
        self.assertEqual(str(ref), str(test3))

        self.assertNotEqual(str(ref), str(dif))#Sanity checks
        self.assertNotEqual(str(test1), str(dif))
        self.assertNotEqual(str(test2), str(dif))
        self.assertNotEqual(str(test3), str(dif))

        return

def runsuite(verbose=True):
    if verbose: verbosity=2
    else: verbosity=1
    suite = unittest.TestLoader().loadTestsFromTestCase(test_basissets)
    unittest.TextTestRunner(verbosity=verbosity).run(suite)
    return

if __name__ == '__main__':
    import sys
    runsuite()
 
