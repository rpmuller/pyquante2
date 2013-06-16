import unittest
import doctest

# This is called "test_all" instead of "test" because I'm hoping that several different test
# modules will be launched from here. Here's an example of how to import other tests:
from test_one import test_one

if __name__ == '__main__':
    # Manually make the test suite
    testSuite = unittest.TestSuite()
    testSuite.addTests(unittest.makeSuite(test_one))

    # Add all doctests
    from pyquante2 import utils
    testSuite.addTest(doctest.DocTestSuite(utils))

    from pyquante2.ints import one,two
    testSuite.addTest(doctest.DocTestSuite(one))
    testSuite.addTest(doctest.DocTestSuite(two))

    from pyquante2.basis import pgbf,cgbf,basisset
    testSuite.addTest(doctest.DocTestSuite(pgbf))
    testSuite.addTest(doctest.DocTestSuite(cgbf))
    testSuite.addTest(doctest.DocTestSuite(basisset))

    from pyquante2.geo import atom, molecule, samples
    testSuite.addTest(doctest.DocTestSuite(atom))
    testSuite.addTest(doctest.DocTestSuite(molecule))
    testSuite.addTest(doctest.DocTestSuite(samples))

    from pyquante2.scf import hamiltonians
    testSuite.addTest(doctest.DocTestSuite(hamiltonians))

    # Run the suite
    unittest.TextTestRunner(verbosity=2).run(testSuite)
