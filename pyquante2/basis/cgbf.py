"""\
 cgbf.py Perform basic operations over contracted gaussian basis
  functions. Uses the functions in pgbf.py.

 References:
  OHT = K. O-ohata, H. Taketa, S. Huzinaga. J. Phys. Soc. Jap. 21, 2306 (1966).
  THO = Taketa, Huzinaga, O-ohata, J. Phys. Soc. Jap. 21,2313 (1966).

 This program is part of the PyQuante quantum chemistry program suite

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 See the LICENSE file for licensing information.
"""

class cgbf:
    """
    Class for a contracted Gaussian basis function

    >>> s = cgbf(exps=[1],coefs=[1])
    >>> s(0,0,0)
    0.7127054703549902
    """
    def __init__(self,origin=(0,0,0),powers=(0,0,0),exps=[],coefs=[]):
        assert len(origin)==3
        assert len(powers)==3

        self.origin = origin
        self.powers = powers

        self.pgbfs = []
        self.coefs = []

        for expn,coef in zip(exps,coefs):
            self.add_pgbf(expn,coef)

        if self.pgbfs:
            self.normalize()
        return

    def __getitem__(self,item): return zip(self.coefs,self.pgbfs).__getitem__(item)
    def __call__(self,x,y,z): return sum(c*p(x,y,z) for c,p in self)

    def add_pgbf(self,expn,coef,renormalize=False):
        from pyquante2.basis.pgbf import pgbf

        self.pgbfs.append(pgbf(expn,self.origin,self.powers))
        self.coefs.append(coef)

        if renormalize:
            self.normalize()
        return

    def normalize(self):
        from math import sqrt
        from pyquante2.ints.one import Sc
        Saa = Sc(self,self)
        Saa_sqrt = sqrt(Saa)
        for i in xrange(len(self.coefs)):
            self.coefs[i] /= Saa_sqrt
        return

if __name__ == '__main__':
    import doctest
    doctest.testmod()
