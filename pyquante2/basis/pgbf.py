"""\
  pgbf.py Perform basic operations over primitive
    gaussian basis functions. The equations herein are based upon
    'Gaussian Expansion Methods for Molecular Orbitals.' H. Taketa,
    S. Huzinaga, and K. O-ohata. H. Phys. Soc. Japan, 21, 2313, 1966.
    [THO paper].

  For the purposes of this routine, a gaussian is defined as:

    g(x,y,z) = A*(x^i)*(y^j)*(z^k)*exp{-a*(r-ro)^2}

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 2.0 and later is covered by the GPL
 license. Please see the file LICENSE that is part of this
 distribution. 

"""

import numpy as np
from pyquante2.utils import fact2,norm2

class pgbf:
    """
    Construct a primitive gaussian basis functions.
    >>> s = pgbf(1.0)
    >>> round(s(0,0,0),6)
    0.712705
    >>> px = pgbf(1.0,powers=(1,0,0))
    >>> round(px(0,0,0),6)
    0.0
    """
    contracted = False
    def __init__(self,exponent,origin=(0,0,0),powers=(0,0,0)):
        self.norm = 1
        assert len(origin) == 3
        assert len(powers) == 3
        self.exponent = float(exponent)
        self.origin = np.array(origin,'d')
        self.powers = powers
        self._normalize()

    def __repr__(self): return "pgbf(%f,%s,%s)" % (self.exponent,tuple(self.origin),self.powers)

    def __call__(self,*args):
        "Compute the amplitude of the PGBF at point x,y,z"
        I,J,K = self.powers
        xyz = np.array(args,'d')
        d = xyz-self.origin
        d2 = norm2(d)
        return self.norm*d[0]**I*d[1]**J*d[2]**K*np.exp(-self.exponent*d2)

    def _normalize(self):
        "Normalize basis function. From THO eq. 2.2"
        l,m,n = self.powers
        self.norm = np.sqrt(pow(2,2*(l+m+n)+1.5)*
                            pow(self.exponent,l+m+n+1.5)/
                            fact2(2*l-1)/fact2(2*m-1)/
                            fact2(2*n-1)/pow(np.pi,1.5))
        return

if __name__ == '__main__':
    import doctest
    doctest.testmod()
