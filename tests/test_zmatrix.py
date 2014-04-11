import unittest
from pyquante2.geo.zmatrix import parse_zmatrix,z2xyz,cartesians_equal

# Notes:
# The tests for the cartesians_equal is not ideal -- I should just build a PyQuante molecule,
#  and have tests for atoms being AlmostEqual

class test_zmatrix(unittest.TestCase):
    def test_parse(self):
        self.assertEqual(parse_zmatrix("H"),[['H']])
        self.assertEqual(parse_zmatrix("H\nH 1 0.7"),[['H'], ['H', 1, 0.7]])
        self.assertEqual(parse_zmatrix("O\nH 1 1.0\nH 1 1.0 2 90.0"),[['O'],['H',1,1.0],['H',1,1.0,2,90.0]])
        return

    def test_convert(self):
        self.assertTrue(cartesians_equal(z2xyz([['H']]),
                                         [['H',0.0,0.0,0.0]]))
        self.assertTrue(cartesians_equal(z2xyz([['H'], ['H', 1, 0.7]]),
                                         [['H',0.0,0.0,0.0],['H',0.7,0.0,0.0]]))
        self.assertTrue(cartesians_equal(z2xyz([['O'],['H',1,1.0],['H',1,1.0,2,90.0]]),
                                         [['O',0.0,0.0,0.0],['H',1.0,0.0,0.0],['H',0.0,1.0,0.0]]))
        self.assertTrue(cartesians_equal(z2xyz([['H'],['O',1,1],['O',2,1,1,90.],['H',3,1,2,90.,1,180.]]),
                                         [['H',0,0,0],['O',1,0,0],['O',1,1,0],['H',2,1,0]]))
        self.assertTrue(cartesians_equal(z2xyz([['H'],['O',1,1],['O',2,1,1,90.],['H',3,1,2,90.,1,0.]]),
                                         [['H',0,0,0],['O',1,0,0],['O',1,1,0],['H',0,1,0]]))
        self.assertTrue(cartesians_equal(z2xyz([['H'],['O',1,1],['O',2,1,1,90.],['H',3,1,2,90.,1,90.]]),
                                         [['H',0,0,0],['O',1,0,0],['O',1,1,0],['H',1,1,-1]]))
        self.assertTrue(cartesians_equal(z2xyz([['H'],['O',1,1],['O',2,1,1,90.],['H',3,1,2,90.,1,-90.]]),
                                         [['H',0,0,0],['O',1,0,0],['O',1,1,0],['H',1,1,1]]))
        return
        

