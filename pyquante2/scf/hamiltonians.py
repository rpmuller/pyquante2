from pyquante2.ints.integrals import onee_integrals,twoe_integrals
import numpy as np

class mock:
    """
    A mock Hamiltonian for testing purposes. Implements energy() and update()
    """
    def __init__(self,*args): self.i = 0
    def update(self,*args): self.i += 1
    def energy(self): return 1.0/pow(10,self.i)

class rhf:
    """
    >>> from pyquante2.geo.samples import h
    >>> from pyquante2.basis.basisset import basisset
    >>> from pyquante2.scf.iterators import averaging
    >>> bfs = basisset(h,'sto3g')
    >>> h_rhf = rhf(h,bfs)
    >>> h.converge(averaging)
    """
    def __init__(self,geo,bfs):
        self.geo = geo
        self.bfs = bfs
        self.i1 = onee_integrals(bfs,geo)
        self.i2 = twoe_integrals(bfs)

    def converge(self,iterator): self.scf_energies = list(iterator(self))

    def update(self,c):
        from pyquante2.utils import geigh,dmat,trace2
        S = self.i1.S

        if c is None:
            E,c = np.eigh(S)

        D = dmat(c,geo.nocc)
        H = self.i1.T + self.i1.V
        self.e0 = trace2(H,D)

        J,K = self.i2.jk(D)
        H += 2*J-K
        self.e1 = trace2(H,D)
        E,c = geigh(H,S)
        return c
        
    def energy(self): return self.e0+self.e1

if __name__ == '__main__':
    import doctest; doctest.testmod()
