# Example of how another test suite can be imported to test_all.py
import unittest

class OtherTest(unittest.TestCase):
    def test_null(self): self.assertTrue(True)
