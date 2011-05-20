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
    from pyquante2.ints import one,two
    from pyquante2.basis import pgbf,cgbf
    testSuite.addTest(doctest.DocTestSuite(utils))
    testSuite.addTest(doctest.DocTestSuite(one))
    testSuite.addTest(doctest.DocTestSuite(two))
    testSuite.addTest(doctest.DocTestSuite(pgbf))
    testSuite.addTest(doctest.DocTestSuite(cgbf))

    # Run the suite
    unittest.TextTestRunner(verbosity=2).run(testSuite)
