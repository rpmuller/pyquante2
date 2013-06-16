import numpy as np
from pyquante2.utils import dmat

def simple(H,c=None,tol=1e-5,maxiters=100):
    """
    Simple SCF iterator.
    """
    Eold = 0
    orbe,c = np.linalg.eigh(H.i1.S)
    for i in xrange(maxiters):
        D = dmat(c,H.geo.nocc())
        c = H.update(D)
        E = H.energy()
        if abs(E-Eold) < tol:
            break
        Eold = E
        yield E
    return

def averaging(H,c=None,fraction=0.5,tol=1e-5,maxiters=100):
    """
    Simplest possible iterator for SCF.
    >>> list(averaging(mock()))
    [0.1, 0.01, 0.001, 0.0001, 1e-05]
    """
    Eold = 0
    orbe,c = np.linalg.eigh(H.i1.S)
    Dold = None
    for i in xrange(maxiters):
        if Dold is not None:
            D = (1-fraction)*Dold + fraction*dmat(c,H.geo.nocc())
        else:
            D = dmat(c,H.geo.nocc())
        c = H.update(D)
        E = H.energy()
        if abs(E-Eold) < tol:
            break
        Eold = E
        yield E
    return

if __name__ == '__main__':
    import doctest; doctest.testmod()
