"""
This is the simplest form of the SCF procedure. I find the current code with iterators
and hamiltonians can be counterintuitive and hard to play with.
"""
import numpy as np
from pyquante2.geo.samples import h2o,lih
from pyquante2.basis.basisset import basisset
from pyquante2.utils import trace2,geigh,dmat
from pyquante2.ints.integrals import onee_integrals, twoe_integrals
from pyquante2.ints.two import ERI

def scf_simple(geo,basisname='sto3g',maxiter=25,verbose=False):
    bfs = basisset(geo,basisname)
    i1 = onee_integrals(bfs,geo)
    i2 = twoe_integrals(bfs)
    if verbose: print ("S=\n%s" % i1.S)
    h = i1.T + i1.V
    if verbose: print ("h=\n%s" % h)
    if verbose: print ("T=\n%s" % i1.T)
    if verbose: print ("V=\n%s" % i1.V)
    E,U = geigh(h,i1.S)
    if verbose: print ("E=\n%s" % E)
    if verbose: print ("U=\n%s" % U)
    Enuke = geo.nuclear_repulsion()
    nocc = geo.nocc()
    Eold = Energy = 0
    if verbose: print ("2e ints\n%s" % i2)
    for i in xrange(maxiter):
        D = dmat(U,nocc)
        if verbose: print ("D=\n%s" % D)
        Eone = trace2(h,D)
        G = i2.get_2jk(D)
        H = h+G
        Etwo = trace2(H,D)
        E,U = geigh(H,i1.S)
        Energy = Enuke+Eone+Etwo
        print ("HF: %d   %10.4f : %10.4f %10.4f %10.4f" % ((i+1),Energy,Enuke,Eone,Etwo))
        if np.isclose(Energy,Eold):
            break
        Eold = Energy
    else:
        print ("Warning: Maxiter %d hit in scf_simple" % maxiter)
    return Energy,E,U

if __name__ == '__main__':
    scf_simple(lih)
    scf_simple(h2o)
