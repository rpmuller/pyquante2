"""
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

try:
    from pyquante2.cone import S,T,V
except:
    print "Couldn't find cython int routine"
    from pyquante2.ints.one import S,T,V

from pyquante2.utils import pairs
from itertools import product
import numpy as np

# The most time-consuming parts of the code are here. Tread carefully!
class twoe_integrals:
    """
    >>> from pyquante2.geo.samples import h
    >>> from pyquante2.basis.basisset import basisset
    >>> bfs = basisset(h,'sto3g')
    >>> twoe_integrals(bfs)
    array([ 0.77460594])
    """
    def __init__(self,bfs):
        nbf = self.nbf = len(bfs)
        self.totlen = nbf*(nbf+1)*(nbf*nbf+nbf+2)/8
        self._2e_ints = np.empty(self.totlen,'d')
        
        for i,j,k,l in iterator(nbf):
            self._2e_ints[index(i,j,k,l)] = ERI(bfs[i],bfs[j],bfs[k],bfs[l])
        return
    def __getitem__(self,pos): return self._2e_ints[index(*pos)]
    def __repr__(self): return repr(self._2e_ints)

    def fetchjk(self,i,j):
        nbf = self.nbf
        temp = np.empty(nbf**2,'d')
        kl = 0
        for k,l in product(xrange(nbf),repeat=2):
            temp[kl] = 2*self[i,j,k,l]-self[i,k,j,l]
            kl += 1
        return temp

    def jk(self,D):
        nbf = self.nbf
        D1 = np.reshape(D,(nbf*nbf,))
        G = np.empty((nbf,nbf),'d')
        for i,j in pairs(xrange(nbf)):
            temp = self.fetchjk(i,j)
            G[i,j] = G[j,i] = np.dot(D1,temp)
        return G
                    
class onee_integrals:
    """
    >>> from pyquante2.geo.samples import h
    >>> from pyquante2.basis.basisset import basisset
    >>> bfs = basisset(h,'sto3g')
    >>> i1 = onee_integrals(bfs,h)
    >>> i1.S
    array([[ 1.]])
    >>> i1.T
    array([[ 0.76003188]])
    >>> i1.V
    array([[-1.22661373]])
    """
    def __init__(self,bfs,geo):
        nbf = self.nbf = len(bfs)
        self.S = np.empty((nbf,nbf),'d')
        self.T = np.empty((nbf,nbf),'d')
        self.V = np.empty((nbf,nbf),'d')
        for i,j in pairs(xrange(nbf)):
            ibf,jbf = bfs[i],bfs[j]
            self.S[i,j] = self.S[j,i] = S(ibf,jbf)
            self.T[i,j] = self.T[j,i] = T(ibf,jbf)
            self.V[i,j] = self.V[j,i] = sum(V(ibf,jbf,at.r) for at in geo)
        return
        

def iterator(nbf):
    """
    Iterator over n**4 integral indices
    >>> list(iterator(2))
    [(0, 0, 0, 0), (0, 0, 0, 1), (0, 0, 1, 1), (0, 1, 0, 1), (0, 1, 1, 1), (1, 1, 1, 1)]
    """
    for i,j in pairs(xrange(nbf)):
        ij = i*(i+1)/2+j
        for k,l in pairs(xrange(nbf)):
            kl = k*(k+1)/2+l
            if ij <= kl:
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
