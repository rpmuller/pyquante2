"""
>>> from pyquante2.geo.samples import h
>>> from pyquante2.basis.basisset import basisset
>>> bfs = basisset(h,'sto3g')
>>> integrals(bfs)
    array([ 1.76093193])

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 2.0 and later is covered by the GPL
 license. Please see the file LICENSE that is part of this
 distribution. 
"""
try:
    from pyquante2.ctwo import ERI_hgp as ERI
except:
    print "Couldn't find cython int routine"
    from pyquante2.ints.hgp import ERI_hgp as ERI

import numpy as np

class integrals:
    def __init__(self,bfs):
        nbf = self.nbf = len(bfs)
        self.totlen = nbf*(nbf+1)*(nbf*nbf+nbf+2)/8
        self._ints = np.empty(self.totlen,'d')
        for i,j,k,l in iterator(nbf):
            print i,j,k,l,index(i,j,k,l),self.totlen
            self._ints[index(i,j,k,l)] = ERI(bfs[i],bfs[j],bfs[k],bfs[l])
        return
    def __getitem__(self,pos): return self._ints[index(*pos)]
    def __repr__(self): return repr(self._ints)

def iterator(nbf):
    """
    Iterator over n**4 integral indices
    >>> list(iterator(2))
    [(0, 0, 0, 0), (1, 0, 0, 0), (1, 0, 1, 0), (1, 1, 0, 0), (1, 1, 1, 0), (1, 1, 1, 1)]
    """
    for i in xrange(nbf):
        for j in xrange(i+1):
            ij = i*(i+1)/2+j
            for k in xrange(nbf):
                for l in xrange(k+1):
                    kl = k*(k+1)/2+l
                    if ij >= kl:
                        yield i,j,k,l
    return

def index(i,j,k,l):
    """
    Indexing into the integral array
    >>> index(0,0,0,0)
    0
    >>> index(1,0,0,0)
    1
    >>> index(0,1,0,0)
    1
    >>> index(0,0,1,0)
    1
    >>> index(0,0,0,1)
    1
    >>> index(1,0,1,0)
    2
    >>> index(0,1,0,1)
    2
    >>> index(1,1,0,0)
    3
    >>> index(0,0,1,1)
    3
    """
    if i<j: i,j = j,i
    if k<l: k,l = l,k
    ij = i*(i+1)/2+j
    kl = k*(k+1)/2+l
    if ij < kl: ij,kl = kl,ij
    return ij*(ij+1)/2+kl
        

if __name__ == '__main__':
    import doctest; doctest.testmod()
