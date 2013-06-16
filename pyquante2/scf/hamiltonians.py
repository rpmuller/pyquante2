from pyquante2.ints.integrals import onee_integrals,twoe_integrals
from pyquante2.utils import trace2,geigh
from pyquante2.scf.iterators import simple,usimple
import numpy as np

class rhf:
    """
    >>> from pyquante2.geo.samples import h2
    >>> from pyquante2.basis.basisset import basisset
    >>> from pyquante2.scf.iterators import simple,averaging
    >>> bfs = basisset(h2,'sto3g')
    >>> h2_rhf = rhf(h2,bfs)
    >>> ens = h2_rhf.converge(simple)
    >>> round(h2_rhf.energy,6)
    -1.1171
    """
    def __init__(self,geo,bfs):
        self.geo = geo
        self.bfs = bfs
        self.i1 = onee_integrals(bfs,geo)
        self.i2 = twoe_integrals(bfs)

    def converge(self,iterator=simple,**kwargs):
        return list(iterator(self,**kwargs))

    def update(self,D):
        self.energy = self.geo.nuclear_repulsion()
        H = self.i1.T + self.i1.V
        self.energy += trace2(H,D)

        JK = self.i2.get_2jk(D)
        H = H + JK
        self.energy += trace2(H,D)
        E,c = geigh(H,self.i1.S)
        return c
        
class uhf:
    def __init__(self,geo,bfs):
        self.geo = geo
        self.bfs = bfs
        self.i1 = onee_integrals(bfs,geo)
        self.i2 = twoe_integrals(bfs)

    def converge(self,iterator=usimple,**kwargs):
        return list(iterator(self,**kwargs))

    def update(self,Da,Db):
        self.energy = self.geo.nuclear_repulsion()
        h = self.i1.T + self.i1.V
        self.energy += trace2(Da+Db,h)/2
        Ja,Ka = self.i2.get_j(Da),self.i2.get_k(Da)
        Jb,Kb = self.i2.get_j(Db),self.i2.get_k(Db)
        Fa = h + Ja + Jb - Ka
        Fb = h + Ja + Jb - Kb
        orbea,ca = geigh(Fa,self.i1.S)
        orbeb,cb = geigh(Fb,self.i1.S)
        self.energy += trace2(Fa,Da)/2 + trace2(Fb,Db)/2
        return ca,cb

if __name__ == '__main__':
    import doctest; doctest.testmod()
