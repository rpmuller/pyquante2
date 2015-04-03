import unittest
from pyquante2 import molecule,rhf,uhf,rohf,h2,h2o,lih,basisset

match_digits = 5

class test_scf(unittest.TestCase):
    def test_h2(self):
        bfs = basisset(h2,'sto3g')
        solver = rhf(h2,bfs)
        ens = solver.converge()
        self.assertAlmostEqual(solver.energy,-1.117099582955609,match_digits)

    def test_h2_631(self):
        bfs = basisset(h2,'6-31G**')
        solver = rhf(h2,bfs)
        ens = solver.converge()
        self.assertAlmostEqual(solver.energy,-1.1313335790123258,match_digits)

    def test_lih(self):
        bfs = basisset(lih,'sto3g')
        solver = rhf(lih,bfs)
        ens = solver.converge()
        self.assertAlmostEqual(solver.energy,-7.8607437,match_digits)

    def test_lih_averaging(self):
        from pyquante2.scf.iterators import AveragingIterator
        bfs = basisset(lih,'sto3g')
        solver = rhf(lih,bfs)
        ens = solver.converge(AveragingIterator)
        self.assertAlmostEqual(solver.energy,-7.8607375733271088,match_digits)

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
        self.assertAlmostEqual(solver.energy,-2.234185653441159,match_digits)
        # This is not quite equal to 2x the h2 energy, but very close

    def test_h2o(self):
        bfs = basisset(h2o,'sto3g')
        solver = rhf(h2o,bfs)
        ens = solver.converge()
        self.assertAlmostEqual(solver.energy,-74.959856675848712,match_digits)

    def test_h2o_averaging(self):
        from pyquante2.scf.iterators import AveragingIterator
        bfs = basisset(h2o,'sto3g')
        solver = rhf(h2o,bfs)
        ens = solver.converge(AveragingIterator)
        self.assertAlmostEqual(solver.energy,-74.959847457272502,match_digits)

    def test_oh_uhf(self):
        from pyquante2 import oh
        bfs = basisset(oh,'sto3g')
        solver = uhf(oh,bfs)
        Es = solver.converge()
        self.assertAlmostEqual(solver.energy,-74.3602416207,match_digits)

    def test_li_uhf(self):
        from pyquante2 import li
        bfs = basisset(li,'sto3g')
        solver = uhf(li,bfs)
        Es = solver.converge()
        self.assertAlmostEqual(solver.energy,-7.31552585354,match_digits)

    def test_oh_rohf(self):
        from pyquante2 import oh
        bfs = basisset(oh,'sto3g')
        solver = rohf(oh,bfs)
        Es = solver.converge()
        self.assertAlmostEqual(solver.energy,-74.3591663559,match_digits)

    def test_li_rohf(self):
        from pyquante2 import li
        bfs = basisset(li,'sto3g')
        solver = rohf(li,bfs)
        Es = solver.converge()
        self.assertAlmostEqual(solver.energy,-7.31552591799,match_digits)


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
    
    
