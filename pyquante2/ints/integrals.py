"""
General module for integral generation and access.
"""
try:
    from pyquante2.ctwo import ERI_hgp as ERI
except:
    print("Couldn't find cython int routine")
    from pyquante2.ints.hgp import ERI_hgp as ERI

try:
    from pyquante2.cone import S,T,V
except:
    print("Couldn't find cython int routine")
    from pyquante2.ints.one import S,T,V

from pyquante2.utils import pairs
from itertools import product
import numpy as np

# This is the old part of the code. It has now been replaced with the one
#  below it, which takes 8x as much space, but is significantly faster.
#  It's also more elegant code.
class twoe_integrals_compressed(object):
    """
    >>> from pyquante2.geo.samples import h
    >>> from pyquante2.basis.basisset import basisset
    >>> bfs = basisset(h,'sto3g')
    >>> twoe_integrals_compressed(bfs)
    array([ 0.77460594])
    """
    def __init__(self,bfs):
        nbf = self.nbf = len(bfs)
        self.totlen = nbf*(nbf+1)*(nbf*nbf+nbf+2)/8
        self._2e_ints = np.empty(self.totlen,'d')
        
        for i,j,k,l in iiterator(nbf):
            self._2e_ints[iindex(i,j,k,l)] = ERI(bfs[i],bfs[j],bfs[k],bfs[l])
        return
    def __getitem__(self,pos): return self._2e_ints[iindex(*pos)]
    def __repr__(self): return repr(self._2e_ints)

    def fetch_2jk(self,i,j):
        nbf = self.nbf
        temp = np.empty(nbf**2,'d')
        kl = 0
        for k,l in product(range(nbf),repeat=2):
            temp[kl] = 2*self[i,j,k,l]-self[i,k,j,l]
            kl += 1
        return temp

    def fetch_j(self,i,j):
        nbf = self.nbf
        temp = np.empty(nbf**2,'d')
        kl = 0
        for k,l in product(range(nbf),repeat=2):
            temp[kl] = self[i,j,k,l]
            kl += 1
        return temp

    def fetch_k(self,i,j):
        nbf = self.nbf
        temp = np.empty(nbf**2,'d')
        kl = 0
        for k,l in product(range(nbf),repeat=2):
            temp[kl] = self[i,k,j,l]
            kl += 1
        return temp

    def make_operator(self,D,fetcher):
        nbf = self.nbf
        D1 = np.reshape(D,(nbf*nbf,))
        G = np.empty((nbf,nbf),'d')
        for i,j in pairs(range(nbf)):
            temp = fetcher(i,j) # replace temp with fetcher()
            G[i,j] = G[j,i] = np.dot(D1,temp)
        return G

    def get_2jk(self,D): return self.make_operator(D,self.fetch_2jk)
    def get_j(self,D): return self.make_operator(D,self.fetch_j)
    def get_k(self,D): return self.make_operator(D,self.fetch_k)

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
            self.V[i,j] = self.V[j,i] = sum(at.atno*V(ibf,jbf,at.r) for at in geo)
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
