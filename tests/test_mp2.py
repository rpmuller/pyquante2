import unittest
from pyquante2 import rhf,basisset,h2,mp2

class test_mp2(unittest.TestCase):
    def test_h2(self):
        bfs = basisset(h2,'6-31g**')
        solver=rhf(h2,bfs)
        solver.converge()
        nvirt = len(bfs)-h2.nocc()
        emp2 = mp2(solver.i2,solver.orbs,solver.orbe,h2.nocc(),len(bfs)-h2.nocc())
        self.assertAlmostEqual(emp2,-0.02632654197486595)
        return

if __name__ == '__main__':
    import sys
    if "-d" in sys.argv:
        debugsuite()
    else:
        runsuite()
