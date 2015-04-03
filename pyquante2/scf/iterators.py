import numpy as np
from pyquante2.utils import dmat,geigh

class SCFIterator(object):
    def __init__(self,H,c=None,tol=1e-5,maxiters=100):
        self.H = H
        self.Eold = 0
        if c is None:
            orbe,self.c = geigh(H.i1.T+H.i1.V,H.i1.S)
        else:
            self.c = c
        self.maxiters = maxiters
        self.tol = tol
        
        self.converged = False
        self.iterations = 0
        return

    def __iter__(self): return self

    def next(self): return self.__next__()
    def __next__(self):
        self.iterations += 1
        if self.iterations > self.maxiters:
            raise StopIteration
        nc,no = self.H.geo.nclosed(),self.H.geo.nopen()
        if no != 0:
            print ("Warning: nopen != 0 in SCF Iterator: %d" % no )
        D = dmat(self.c,nc,no)
        self.c = self.H.update(D)
        E = self.H.energy
        if abs(E-self.Eold) < self.tol:
            self.converged = True
            raise StopIteration
        self.Eold = E
        return E

class USCFIterator(SCFIterator):
    def __init__(self,H,c=None,tol=1e-5,maxiters=100):
        SCFIterator.__init__(self,H,c,tol,maxiters)
        self.nup,self.ndown = self.H.geo.nup(),self.H.geo.ndown()
        self.cup = self.cdown = self.c
        return

    def __next__(self):
        self.iterations += 1
        if self.iterations > self.maxiters:
            raise StopIteration
        Dup = dmat(self.cup,self.nup)
        Ddown = dmat(self.cdown,self.ndown)
        self.cup,self.cdown = self.H.update(Dup,Ddown)
        E = self.H.energy
        if abs(E-self.Eold) < self.tol:
            self.converged = True
            raise StopIteration
        self.Eold = E
        return E

class ROHFIterator(SCFIterator):
    def __init__(self,H,c=None,tol=1e-5,maxiters=100):
        SCFIterator.__init__(self,H,c,tol,maxiters)
        self.nup,self.ndown = self.H.geo.nup(),self.H.geo.ndown()
        return

    def __next__(self):
        self.iterations += 1
        if self.iterations > self.maxiters:
            raise StopIteration
        Dup = dmat(self.c,self.nup)
        Ddown = dmat(self.c,self.ndown)
        self.c = self.H.update(Dup,Ddown,self.c)
        E = self.H.energy
        if abs(E-self.Eold) < self.tol:
            self.converged = True
            raise StopIteration
        self.Eold = E
        return E

class AveragingIterator(SCFIterator):
    def __init__(self,H,c=None,fraction=0.5,tol=1e-5,maxiters=100):
        SCFIterator.__init__(self,H,c,tol,maxiters)
        self.fraction = fraction
        self.Dold = dmat(self.c,self.H.geo.nocc())
        return

    def __next__(self):
        D = (1-self.fraction)*self.Dold + self.fraction*dmat(self.c,self.H.geo.nocc())
        self.Dold = D
        self.c = self.H.update(D)
        E = self.H.energy
        if abs(E-self.Eold) < self.tol:
            self.converged = True
            raise StopIteration
        self.Eold = E
        return E

if __name__ == '__main__':
    import doctest; doctest.testmod()
