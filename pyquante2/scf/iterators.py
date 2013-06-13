from pyquante2.scf.hamiltonians import mock

def simple(H,c=None,tol=1e-5,maxiters=100):
    """
    Simplest possible iterator for SCF.
    >>> list(simple(mock()))
    [0.1, 0.01, 0.001, 0.0001, 1e-05]
    """
    Eold = 0
    for i in xrange(maxiters):
        c = H.update(c)
        E = H.energy()
        if abs(E-Eold) < tol:
            break
        Eold = E
        yield E
    return

def averaging(H,c=None,fraction=0.5,tol=1e-5,maxiters=100):
    """
    Simple orbital averaging.
    >>> list(simple(mock()))
    [0.1, 0.01, 0.001, 0.0001, 1e-05]
    """
    Eold = 0
    for i in xrange(maxiters):
        cnew = H.update(c)
        if c is None:
            c = cnew
        else:
            c = (1-fraction)*c + fraction*cnew
        E = H.energy()
        if abs(E-Eold) < tol:
            break
        Eold = E
        yield E
    return

if __name__ == '__main__':
    import doctest; doctest.testmod()
