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
    >>> from pyquante2.utils import isnear
    >>> s = cgbf(exps=[1],coefs=[1])
    >>> isnear(s(0,0,0),0.7127054704,1e-9)
    True
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

    def exps(self): return [p.exponent for p in self.pgbfs]

    def __getitem__(self,item): return zip(self.coefs,self.pgbfs).__getitem__(item)
    def __call__(self,x,y,z): return sum(c*p(x,y,z) for c,p in self)
    def as_tuple(self): return (tuple(self.origin),self.powers,zip(self.exps(),self.coefs))
    def __repr__(self): return repr(self.as_tuple())


    def add_pgbf(self,expn,coef,renormalize=False):
        from pyquante2.basis.pgbf import pgbf

        self.pgbfs.append(pgbf(expn,self.origin,self.powers))
        self.coefs.append(coef)

        if renormalize:
            self.normalize()
        return

    def normalize(self):
        from math import sqrt
        Saa = Sc(self,self)
        Saa_sqrt = sqrt(Saa)
        for i in xrange(len(self.coefs)):
            self.coefs[i] /= Saa_sqrt
        return

# Sketch of what could be the contracted version of overlap. Move to one.py when done
def Sc(a,b):
    from pyquante2.ints.one import S,contract
    return contract(S,a,b)
    #return sum(ca*cb*S(pa,pb) for (ca,pa) in a for (cb,pb) in b)

# Alternatively, I could define pgbf.__getitem__  as [(1,pgbf)] and use the above S for all fns.
# Could be appealing to have a single S,T,V, etc.
# Probably dont want to use this trick for ERIs, though.

if __name__ == '__main__':
    import doctest
    doctest.testmod()
