"""\
Basis set constructor

#>>> from pyquante2.geo.samples import h
#>>> print(shell_basisset(h,'sto3g'))
cgbf((0.0, 0.0, 0.0),(0, 0, 0),[3.42525091, 0.62391373, 0.1688554],[0.1543289707029839, 0.5353281424384733, 0.44463454202535485])

"""
import numpy as np
from pyquante2 import settings
from pyquante2.basis import data
from pyquante2.basis.cgbf import cgbf
from pyquante2.basis.tools import sym2pow,sym2am,am2pow

class shell(object):
    def __init__(self,am,origin,exps,coefs):
        self.am = am
        self.origin = origin
        self.exps = exps
        self.coefs = coefs
        self.bfs = []
        for power in am2pow[sym]:
            self.bfs.append(cgbf(atom.r,power,exps,coefs))
        return

class shell_basisset(object):
    def __init__(self,atoms,name='sto3g',**kwargs):
        self.shells = []
        self.name=name
        basis_data = data.basis[name]
        omit_f = kwargs.get('omit_f',settings.omit_f)
        for atom in atoms:
            for sym,prims in basis_data[atom.atno]:
                if omit_f and sym == 'F': continue
                exps = [e for e,c in prims]
                coefs = [c for e,c in prims]
                self.shells.append(shell(sym2am[sym],atom.r,exps,coefs))
        return

    def __len__(self): return sum(len(sh.bfs) for sh in self.shells)
    def bfs(self): return (bf for sh in self.shells for bf in sh.bfs)

