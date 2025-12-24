"""\
Basis set constructor

>>> from pyquante2.geo.samples import h
>>> bfs = basisset(h,'sto3g')
>>> bfs
cgbf((0.0, 0.0, 0.0),(0, 0, 0),[3.42525091, 0.62391373, 0.1688554],[0.1543289707029839, 0.5353281424384733, 0.44463454202535485])
>>> bfs.shells
[S shell: [3.42525091, 0.62391373, 0.1688554],[0.15432897, 0.53532814, 0.44463454]]
"""

import numpy as np
from pyquante2 import settings
from pyquante2.basis import data
from pyquante2.basis.cgbf import cgbf
from pyquante2.basis.tools import sym2pow,sym2am,am2pow,am2sym

class shell(object):
    def __init__(self,am,origin,exps,coefs):
        self.am = am
        self.origin = origin
        self.exps = exps
        self.coefs = coefs
        self.bfs = []
        for power in am2pow[am]:
            self.bfs.append(cgbf(origin,power,exps,coefs))
        return

    def __repr__(self):
        return "%s shell: %s,%s" % (am2sym[self.am],self.exps,self.coefs)

class basisset(object):
    def __init__(self,atoms,name='sto3g',**kwargs):
        self.bfs = []
        self.shells = []
        self.name = name.lower()
        basis_data = data.basis[self.name]
        omit_f = kwargs.get('omit_f',settings.omit_f)
        for atom in atoms:
            for sym,prims in basis_data[atom.atno]:
                if omit_f and sym == 'F': continue
                exps = [e for e,c in prims]
                coefs = [c for e,c in prims]
                self.shells.append(shell(sym2am[sym],atom.r,exps,coefs))
                for power in sym2pow[sym]:
                    self.bfs.append(cgbf(atom.r,power,exps,coefs))
        return

    def __getitem__(self,i): return self.bfs.__getitem__(i)
    def __repr__(self): return "\n".join(repr(bf) for bf in self.bfs)
    def __len__(self): return len(self.bfs)
    def mesh(self,points):
        nbf = len(self.bfs)
        ng,sb3 = points.shape
        self.bfmesh = np.empty((ng,nbf),'d')
        for i,bf in enumerate(self.bfs):
            self.bfmesh[:,i] = bf.mesh(points)
        return self.bfmesh

    def atom_symbols(self,geo):
        """\
        >>> from pyquante2.geo.samples import h
        >>> bfs = basisset(h,'sto3g')
        >>> bfs.atom_symbols(h)
        ['H0 s']
        """
        return [bf.atom_symbol(geo) for bf in self.bfs]
            
if __name__ == '__main__':
    import doctest
    doctest.testmod()
