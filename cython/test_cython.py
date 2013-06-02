import unittest
from pyquante2.cutils import hello

class test_cython(unittest.TestCase):
    def test_hello(self):
        s = hello()
        self.assertEqual(s,"Hello, World!")

if __name__ == '__main__':
    # Manually make the test suite
    testSuite = unittest.TestSuite()
    testSuite.addTests(unittest.makeSuite(test_cython))

    # Run the suite
    unittest.TextTestRunner(verbosity=2).run(testSuite)
    
