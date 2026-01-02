"""
General module for integral generation and access.
"""
import logging

logger = logging.getLogger(__name__)

try:
    from pyquante2.cints.hgp import ERI
except ImportError:
    logger.debug("Couldn't find cython int routine, using pure Python")
    from pyquante2.ints.hgp import ERI

try:
    from pyquante2.cints.one import S,T,V
except ImportError:
    logger.debug("Couldn't find cython one-electron routine, using pure Python")
    from pyquante2.ints.one import S,T,V

from pyquante2.utils import pairs
from itertools import product
import numpy as np

class twoe_integrals(object):
    """
    >>> from pyquante2.geo.samples import h
    >>> from pyquante2.basis.basisset import basisset
    >>> bfs = basisset(h,'sto3g')
    >>> twoe_integrals(bfs)
    array([ 0.77460594])
    """
    def __init__(self,bfs):
        nbf = self.nbf = len(bfs)
        self._2e_ints = np.empty((nbf,nbf,nbf,nbf),'d')
        ints = self._2e_ints
        
        for i,j,k,l in iiterator(nbf):
            ints[i,j,k,l] = ints[j,i,k,l] = ints[i,j,l,k] = ints[j,i,l,k] = \
                            ints[k,l,i,j] = ints[l,k,i,j] = ints[k,l,j,i] = \
                            ints[l,k,j,i] = ERI(bfs[i],bfs[j],bfs[k],bfs[l])
        return
    def __getitem__(self,*args): return self._2e_ints.__getitem__(*args)
    def __repr__(self): return repr(self._2e_ints.ravel())

    def transform(self,c): return np.einsum('aI,bJ,cK,dL,abcd->IJKL',c,c,c,c,self._2e_ints)
    def transform_mp2(self,c,nocc):
        return np.einsum('aI,bJ,cK,dL,abcd->IJKL',c[:,:nocc],c,c[:,:nocc],c,self._2e_ints)


    # This turns out to be slower:
    #def get_j(self,D): return np.einsum('ij,ijkl->kl',D,self._2e_ints)
    def get_j(self,D): return np.einsum('kl,ijkl->ij',D,self._2e_ints)
    def get_k(self,D): return np.einsum('ij,ikjl->kl',D,self._2e_ints)
    def get_2jk(self,D): return 2*self.get_j(D)-self.get_k(D)
                    
class onee_integrals(object):
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
        for i,j in pairs(range(nbf)):
            ibf,jbf = bfs[i],bfs[j]
            self.S[i,j] = self.S[j,i] = S(ibf,jbf)
            self.T[i,j] = self.T[j,i] = T(ibf,jbf)
            self.V[i,j] = self.V[j,i] = sum(at.Z*V(ibf,jbf,at.r) for at in geo)
        return
        

def iiterator(nbf):
    """
    Iterator over n**4 integral indices
    >>> list(iiterator(2))
    [(0, 0, 0, 0), (0, 0, 0, 1), (0, 0, 1, 1), (0, 1, 0, 1), (0, 1, 1, 1), (1, 1, 1, 1)]
    """
    for i,j in pairs(range(nbf)):
        ij = i*(i+1)/2+j
        for k,l in pairs(range(nbf)):
            kl = k*(k+1)/2+l
            if ij <= kl:
                yield i,j,k,l
    return

def iindex(i,j,k,l):
    """
    Indexing into the integral array
    >>> iindex(0,0,0,0)
    0
    >>> iindex(1,0,0,0)
    1
    >>> iindex(0,1,0,0)
    1
    >>> iindex(0,0,1,0)
    1
    >>> iindex(0,0,0,1)
    1
    >>> iindex(1,0,1,0)
    2
    >>> iindex(0,1,0,1)
    2
    >>> iindex(1,1,0,0)
    3
    >>> iindex(0,0,1,1)
    3
    """
    if i<j: i,j = j,i
    if k<l: k,l = l,k
    ij = (i*(i+1))//2+j
    kl = (k*(k+1))//2+l
    if ij < kl: ij,kl = kl,ij
    return (ij*(ij+1))//2+kl
        

if __name__ == '__main__':
    import doctest; doctest.testmod()
