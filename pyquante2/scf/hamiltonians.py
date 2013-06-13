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
    >>> from pyquante2.geo.samples import h2
    >>> from pyquante2.basis.basisset import basisset
    >>> from pyquante2.scf.iterators import simple,averaging
    >>> bfs = basisset(h2,'sto3g')
    >>> h2_rhf = rhf(h2,bfs)
    >>> ens = h2_rhf.converge(simple)
    >>> round(h2_rhf.energy(),6)
    -1.1171
    """
    def __init__(self,geo,bfs):
        self.geo = geo
        self.bfs = bfs
        self.i1 = onee_integrals(bfs,geo)
        self.i2 = twoe_integrals(bfs)

    def converge(self,iterator,verbose=False):
        return list(iterator(self))

    def update(self,c):
        from pyquante2.utils import geigh,dmat,trace2
        S = self.i1.S

        if c is None:
            E,c = np.linalg.eigh(S)

        self._energy = self.geo.nuclear_repulsion()
        D = dmat(c,self.geo.nocc())
        H = self.i1.T + self.i1.V
        self._energy += trace2(H,D)

        JK = self.i2.jk(D)
        H = H + JK
        self._energy += trace2(H,D)
        E,c = geigh(H,S)
        return c
        
    def energy(self): return self._energy

if __name__ == '__main__':
    import doctest; doctest.testmod()
