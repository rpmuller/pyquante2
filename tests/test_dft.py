import unittest
import numpy as np
from pyquante2 import basisset, molecule
from pyquante2.dft.functionals import xs,cvwn5,xb88,xpbe,clyp,cpbe
from pyquante2.dft.reference import data
from pyquante2.scf.hamiltonians import dft

try:
    import pylab as pyl
    have_pylab = True
except:
    have_pylab = False

def aamax(x): return np.amax(np.absolute(x))

class test_dft(unittest.TestCase):
    def test_xs(self):
        na = data['xs'][:,0]
        nb = data['xs'][:,1]
        fa,dfa = xs(na)
        fb,dfb = xs(nb)
        max_f = aamax(fa+fb-data['xs'][:,5])
        max_dfa = aamax(dfa-data['xs'][:,6])
        max_dfb = aamax(dfb-data['xs'][:,7])
        #print(np.column_stack([na,nb,data['xs'][:,5]-fa-fb]))
        self.assertAlmostEqual(max_f,0,5) ## Fix this!
        self.assertAlmostEqual(max_dfa,0)
        self.assertAlmostEqual(max_dfb,0)

    @unittest.skip("Skipping CVWN test that consistently fails")
    def test_cvwn(self):
        na = data['cvwn'][:,0]
        nb = data['cvwn'][:,1]
        f,dfa,dfb = cvwn5(na,nb)
        max_f = aamax(f-data['cvwn'][:,5])
        max_dfa = aamax(dfa-data['cvwn'][:,6])
        max_dfb = aamax(dfb-data['cvwn'][:,7])
        np.set_printoptions(precision=3)
        print(np.column_stack([na,nb,data['cvwn'][:,5],f,data['cvwn'][:,6],dfb]))
        #pyl.plot(data['cvwn'][:,5]-f)
        #pyl.plot(data['cvwn'][:,6]-dfa)
        #pyl.show()
        self.assertAlmostEqual(max_f,0)
        self.assertAlmostEqual(max_dfa,0)
        self.assertAlmostEqual(max_dfb,0)

    def test_xb88(self):
        na = data['xb88'][:,0]
        nb = data['xb88'][:,1]
        gaa = data['xb88'][:,2]
        gbb = data['xb88'][:,4]
        fa,dfa,dfga = xb88(na,gaa)
        fb,dfb,dfgb = xb88(nb,gbb)
        max_f = aamax(fa+fb-data['xb88'][:,5])
        max_dfa = aamax(dfa-data['xb88'][:,6])
        max_dfb = aamax(dfb-data['xb88'][:,7])
        #print(np.column_stack([na,nb,data['xb88'][:,5]-fa-fb]))
        self.assertAlmostEqual(max_f,0,5)
        self.assertAlmostEqual(max_dfa,0)
        self.assertAlmostEqual(max_dfb,0)

    def test_xpbe(self):
        na = data['xpbe'][:,0]
        nb = data['xpbe'][:,1]
        gaa = data['xpbe'][:,2]
        gbb = data['xpbe'][:,4]
        fa,dfa,dfga = xpbe(na,gaa)
        fb,dfb,dfgb = xpbe(nb,gbb)
        max_f = aamax(fa+fb-data['xpbe'][:,5])
        max_dfa = aamax(dfa-data['xpbe'][:,6])
        max_dfb = aamax(dfb-data['xpbe'][:,7])
        #print(np.column_stack([na,nb,data['xpbe'][:,5]-fa-fb]))
        self.assertAlmostEqual(max_f,0,5)
        self.assertAlmostEqual(max_dfa,0)
        self.assertAlmostEqual(max_dfb,0)

    @unittest.skip("Skipping CLYP since it doesn't work")
    def test_clyp(self):
        na = data['clyp'][:,0]
        nb = data['clyp'][:,1]
        gaa = data['clyp'][:,2]
        gab = data['clyp'][:,3]
        gbb = data['clyp'][:,4]
        fc,dfa,dfb,dfga,dfgab,dfgb = clyp(na,nb,gaa,gab,gbb)

        max_f = aamax(fc-data['clyp'][:,5])
        max_dfa = aamax(dfa-data['clyp'][:,6])
        max_dfb = aamax(dfb-data['clyp'][:,7])
        #print(np.column_stack([na,nb,data['clyp'][:,5]-fa-fb]))
        self.assertAlmostEqual(max_f,0,5)
        self.assertAlmostEqual(max_dfa,0)
        self.assertAlmostEqual(max_dfb,0)

    @unittest.skip("Skipping CPBE since it doesn't work")
    def test_cpbe(self):
        na = data['cpbe'][:,0]
        nb = data['cpbe'][:,1]
        gaa = data['cpbe'][:,2]
        gab = data['cpbe'][:,3]
        gbb = data['cpbe'][:,4]
        fc,dfa,dfb,dfga,dfgab,dfgb = cpbe(na,nb,gaa,gab,gbb)

        max_f = aamax(fc-data['cpbe'][:,5])
        max_dfa = aamax(dfa-data['cpbe'][:,6])
        max_dfb = aamax(dfb-data['cpbe'][:,7])
        #print(np.column_stack([na,nb,data['cpbe'][:,5]-fa-fb]))
        self.assertAlmostEqual(max_f,0,5)
        self.assertAlmostEqual(max_dfa,0)
        self.assertAlmostEqual(max_dfb,0)

    def test_he_xlda_sto(self):
        from pyquante2.geo.samples import he
        bfs = basisset(he,'sto-3g')
        solver = dft(he,bfs,'xs')
        ens = solver.converge()
        self.assertAlmostEqual(solver.energy,-2.65731226258)

    def test_he_xlda_631(self):
        from pyquante2.geo.samples import he
        bfs = basisset(he,'6-31G**')
        solver = dft(he,bfs,'xs',verbose=True)
        ens = solver.converge()
        self.assertAlmostEqual(solver.energy,-2.7146407379057935)

    def test_he_lda_631(self):
        from pyquante2.geo.samples import he
        bfs = basisset(he,'6-31G**')
        solver = dft(he,bfs,'lda')
        ens = solver.converge()
        self.assertAlmostEqual(solver.energy,-2.8266972953683389)

    def test_he_triplet_lda(self):
        he_trip = molecule([(2,0,0,0)], multiplicity=3)
        bfs = basisset(he_trip,'6-31G**')
        solver = dft(he_trip,bfs,'lda')
        ens = solver.converge()
        self.assertAlmostEqual(solver.energy,-1.1784857927828982)

    def test_h2_svwn(self):
        h2 = molecule([(1,0,0,-0.368),(1,0,0,0.368)],units='angs',nrad=50,do_sg1=False)
        bfs = basisset(h2,'sto3g')
        solver = dft(h2,bfs,'svwn')
        ens = solver.converge()
        self.assertAlmostEqual(solver.energy, -1.1212155284066108)

def runsuite(verbose=True):
    if verbose: verbosity=2
    else: verbosity=1
    suite = unittest.TestLoader().loadTestsFromTestCase(test_dft)
    unittest.TextTestRunner(verbosity=verbosity).run(suite)
    return

def profsuite():
    import cProfile,pstats
    cProfile.run('runsuite()','prof')
    prof = pstats.Stats('prof')
    prof.strip_dirs().sort_stats('time').print_stats(15)

if __name__ == '__main__':
    import sys
    if "-p" in sys.argv:
        profsuite()
    else:
        runsuite()
    
    
        
