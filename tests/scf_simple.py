"""
This is the simplest form of the SCF procedure. I find the current code with iterators
and hamiltonians can be counterintuitive and hard to play with.
"""
import numpy as np
from pyquante2.geo.samples import h2o,lih,oh,li
from pyquante2.basis.basisset import basisset
from pyquante2.utils import trace2,geigh,dmat
from pyquante2.ints.integrals import onee_integrals, twoe_integrals
from pyquante2.ints.two import ERI

def scf_simple(geo,basisname='sto3g',maxiter=25,verbose=False):
    if geo.nopen() == 0: return rhf_simple(geo,basisname,maxiter,verbose)
    return uhf_simple(geo,basisname,maxiter,verbose)

def rhf_simple(geo,basisname='sto3g',maxiter=25,verbose=False):
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
        #G = i2.get_2jk(D)
        G = 2*i2.get_j(D)-i2.get_k(D)
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

def uhf_simple(geo,basisname='sto3g',maxiter=25,verbose=False):
    if geo.nopen() == 0: return scf_simple(geo,basisname,maxiter,verbose)
    bfs = basisset(geo,basisname)
    i1 = onee_integrals(bfs,geo)
    i2 = twoe_integrals(bfs)
    h = i1.T + i1.V
    E,U = geigh(h,i1.S)
    Enuke = geo.nuclear_repulsion()
    Eold = Energy = 0
    ca = cb = U
    na,nb = geo.nup(),geo.ndown()

    for i in xrange(maxiter):
        Energy = Enuke
        Da = dmat(ca,na)
        Db = dmat(cb,nb)
        h = i1.T + i1.V
        Energy += trace2(Da+Db,h)/2
        Ja,Ka = i2.get_j(Da),i2.get_k(Da)
        Jb,Kb = i2.get_j(Db),i2.get_k(Db)
        Fa = h + Ja + Jb - Ka
        Fb = h + Ja + Jb - Kb
        orbea,ca = geigh(Fa,i1.S)
        orbeb,cb = geigh(Fb,i1.S)
        Energy += trace2(Fa,Da)/2 + trace2(Fb,Db)/2
        print ("UHF: %d   %10.4f : %10.4f" % ((i+1),Energy,Enuke))
        if np.isclose(Energy,Eold):
            break
        Eold = Energy
    else:
        print ("Warning: Maxiter %d hit in scf_simple" % maxiter)
    return Energy,E,U

def testall():
    print ("LiH energy should be -7.8607")
    scf_simple(lih)
    print ("H2O energy should be -74.9598")
    scf_simple(h2o)
    print ("OH energy should be -73.3602")
    scf_simple(oh)
    print ("Li energy should be -7.3155")
    scf_simple(li)
    

if __name__ == '__main__':
    scf_simple(lih,verbose=True)
