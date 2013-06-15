import unittest
from pyquante2 import molecule,rhf,h2,h2o,lih,basisset
from pyquante2.scf.iterators import simple,averaging

class test_scf(unittest.TestCase):
    def test_h2(self):
        bfs = basisset(h2,'sto3g')
        h2_rhf = rhf(h2,bfs)
        ens = h2_rhf.converge(simple)
        self.assertAlmostEqual(h2_rhf.energy(),-1.117099582955609)

    def test_lih(self):
        bfs = basisset(lih,'sto3g')
        lih_rhf = rhf(lih,bfs)
        #E,c = lih_rhf.iterate()
        ens = lih_rhf.converge(simple)
        self.assertAlmostEqual(lih_rhf.energy(),-7.8607437) # Jaguar energy for 6-31G**

    def test_h4(self):
        h4 = molecule([
            (1,  0.00000000,     0.00000000,     0.36628549),
            (1,  0.00000000,     0.00000000,    -0.36628549),
            (1,  0.00000000,     4.00000000,     0.36628549),
            (1,  0.00000000,     4.00000000,    -0.36628549),
            ],
                      units='Angstrom')
        bfs = basisset(h4,'sto3g')
        h4_rhf = rhf(h4,bfs)
        ens = h4_rhf.converge(simple)
        self.assertAlmostEqual(h4_rhf.energy(),-2.234185653441159)
        # This is not quite equal to 2x the h2 energy, but very close

    # def test_h2o(self):
    #     bfs = basisset(h2o,'sto3g')
    #     h2o_rhf = rhf(h2o,bfs)
    #     ens = h2o_rhf.converge(averaging)
    #     print ens
    #     self.assertAlmostEqual(h2o_rhf.energy(),-1.117099582955609)
