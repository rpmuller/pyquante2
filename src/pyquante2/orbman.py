#!/usr/bin/env python
from pyquante2 import molecule
from pyquante2.geo.atom import atom

class orbman(object):
    """\
    >>> from pyquante2 import h2,rhf,basisset
    >>> bfs = basisset(h2,'sto3g')
    >>> solver = rhf(h2,bfs)
    >>> Es = solver.converge()
    >>> orbs = solver.orbs
    >>> o = orbman(orbs,bfs,h2)
    >>> o.norb
    2
    >>> o[0]
       -0.5480  H0 s
       -0.5480  H1 s
    >>> o[1]
       -1.2212  H0 s
        1.2212  H1 s
    """
    def __init__(self,orbs,bfs,geo,cutoff=0.1):
        self.orbs = orbs
        self.tags = bfs.atom_symbols(geo)
        self.cutoff = cutoff
        self.nbf,self.norb = orbs.shape
        return

    def __getitem__(self,i):
        assert 0 <= i < self.norb
        for c,l in zip(self.orbs[:,i],self.tags):
            if abs(c) > self.cutoff:
                print("%10.4f %5s" % (c,l))
        return


    

