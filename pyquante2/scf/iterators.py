import numpy as np
from pyquante2.utils import dmat

class SCFIterator:
    def __init__(self,H,c=None,tol=1e-5,maxiters=100):
        self.H = H
        self.Eold = 0
        if c is None:
            orbe,self.c = np.linalg.eigh(H.i1.S)
        else:
            self.c = c
        self.maxiters = maxiters
        self.tol = tol
        
        self.converged = False
        self.iterations = 0
        return

    def __iter__(self): return self

    def next(self):
        self.iterations += 1
        if self.iterations > self.maxiters:
            raise StopIteration
        D = dmat(self.c,self.H.geo.nocc())
        self.c = self.H.update(D)
        E = self.H.energy
        if abs(E-self.Eold) < self.tol:
            self.converged = True
            raise StopIteration
        self.Eold = E
        return E

def simple(H,c=None,tol=1e-5,maxiters=100):
    """
    Simple SCF iterator.
    """
    Eold = 0
    if c is None:
        orbe,c = np.linalg.eigh(H.i1.S)
    for i in xrange(maxiters):
        D = dmat(c,H.geo.nocc())
        c = H.update(D)
        E = H.energy
        if abs(E-Eold) < tol:
            break
        Eold = E
        yield E
    return

def usimple(H,c=None,tol=1e-5,maxiters=100):
    Eold = 0
    if c is None:
        orbe,cup = np.linalg.eigh(H.i1.S)
        cdown = cup
    nup,ndown = H.geo.nup(),H.geo.ndown()
    for i in xrange(maxiters):
        Dup = dmat(cup,nup)
        Ddown = dmat(cdown,ndown)
        cup,cdown = H.update(Dup,Ddown)
        E = H.energy
        if abs(E-Eold) < tol:
            break
        Eold = E
        yield E
    return
  
def averaging(H,c=None,fraction=0.5,tol=1e-5,maxiters=100):
    """
    Density matrix averaging iterator.
    """
    Eold = 0
    if c is None:
        orbe,c = np.linalg.eigh(H.i1.S)
    Dold = dmat(c,H.geo.nocc())
    for i in xrange(maxiters):
        D = (1-fraction)*Dold + fraction*dmat(c,H.geo.nocc())
        Dold = D
        c = H.update(D)
        E = H.energy
        if abs(E-Eold) < tol:
            break
        Eold = E
        yield E
    return

if __name__ == '__main__':
    import doctest; doctest.testmod()
