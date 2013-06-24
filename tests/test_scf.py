import unittest
from pyquante2 import molecule,rhf,uhf,h2,h2o,lih,basisset

class test_scf(unittest.TestCase):
    def test_h2(self):
        bfs = basisset(h2,'sto3g')
        solver = rhf(h2,bfs)
        ens = solver.converge()
        self.assertAlmostEqual(solver.energy,-1.117099582955609)

    def test_h2_631(self):
        bfs = basisset(h2,'6-31G**')
        solver = rhf(h2,bfs)
        ens = solver.converge()
        self.assertAlmostEqual(solver.energy,-1.1313335790123258)

    def test_lih(self):
        bfs = basisset(lih,'sto3g')
        solver = rhf(lih,bfs)
        ens = solver.converge()
        self.assertAlmostEqual(solver.energy,-7.8607437)

    def test_lih_averaging(self):
        from pyquante2.scf.iterators import averaging
        bfs = basisset(lih,'sto3g')
        solver = rhf(lih,bfs)
        ens = solver.converge(averaging)
        self.assertAlmostEqual(solver.energy,-7.8607375733271088)

    def test_h4(self):
        h4 = molecule([
            (1,  0.00000000,     0.00000000,     0.36628549),
            (1,  0.00000000,     0.00000000,    -0.36628549),
            (1,  0.00000000,     4.00000000,     0.36628549),
            (1,  0.00000000,     4.00000000,    -0.36628549),
            ],
                      units='Angstrom')
        bfs = basisset(h4,'sto3g')
        solver = rhf(h4,bfs)
        ens = solver.converge()
        self.assertAlmostEqual(solver.energy,-2.234185653441159)
        # This is not quite equal to 2x the h2 energy, but very close

    def test_h2o(self):
        bfs = basisset(h2o,'sto3g')
        solver = rhf(h2o,bfs)
        ens = solver.converge()
        self.assertAlmostEqual(solver.energy,-74.959856073494194,6)

    def test_h2o_averaging(self):
        from pyquante2.scf.iterators import averaging
        bfs = basisset(h2o,'sto3g')
        solver = rhf(h2o,bfs)
        ens = solver.converge(averaging)
        self.assertAlmostEqual(solver.energy,-74.959846854926553,6)

    def test_oh(self):
        from pyquante2 import oh
        bfs = basisset(oh,'sto3g')
        solver = uhf(oh,bfs)
        Es = solver.converge()
        self.assertAlmostEqual(solver.energy,-74.14666861386641,6)

    def test_li(self):
        from pyquante2 import li
        bfs = basisset(li,'sto3g')
        solver = uhf(li,bfs)
        Es = solver.converge()
        self.assertAlmostEqual(solver.energy,-7.2301642412807379)


def runsuite(verbose=True):
    # To use psyco, uncomment this line:
    #import psyco; psyco.full()
    if verbose: verbosity=2
    else: verbosity=1
    # If you want more output, uncomment this line:
    #logging.basicConfig(format="%(message)s",level=logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(test_scf)
    unittest.TextTestRunner(verbosity=verbosity).run(suite)
    # Running without verbosity is equivalent to replacing the above
    # two lines with the following:
    #unittest.main()
    return

def debugsuite():
    import cProfile,pstats
    cProfile.run('runsuite()','prof')
    prof = pstats.Stats('prof')
    prof.strip_dirs().sort_stats('time').print_stats(15)

if __name__ == '__main__':
    import sys
    if "-d" in sys.argv:
        debugsuite()
    else:
        runsuite()
    
    
