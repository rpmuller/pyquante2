"""
The DFT grids are a little different in pyquante2 from pyquante1. Here we are
only storing the points and the weights, and we will use other data objects to
store, say, the density, the basis functions, or the various gradients at each
point.
"""
import numpy as np
from pyquante2.grid.atomic_grid import atomic_grid

class grid:
    def __init__(self,atoms,**kwargs):
        self.atoms = atoms
        agrids = [atomicgrid(atom,**kwargs) for atom in atoms]
        becke_reweight_atoms(atoms,agrids)
        self.points = np.vstack([agrid.points for agrid in agrids])
        self.ng,sb4 = self.points.shape
        assert sb4==4
        return

    def __len__(self): return self.ng
    def __getitem__(self,*args): self.points.__getitem__(*args)

# These are the functions for the becke projection operator
def fbecke(x,n=3):
    for i in range(n): x = pbecke(x)
    return x
def pbecke(x): return 1.5*x-0.5*pow(x,3)
def sbecke(x,n=3): return 0.5*(1-fbecke(x,n))

def becke_atomic_grid_p(iat,(xp,yp,zp),atoms,**opts):
    do_becke_hetero = opts.get('do_becke_hetero',True)
    nat = len(atoms)
    sprod = 1
    ati = atoms[iat]
    rip2 = dist2(ati.pos(),(xp,yp,zp))
    rip = sqrt(rip2)
    for jat in range(nat):
        if jat == iat: continue
        atj = atoms[jat]
        rjp2 = dist2(atj.pos(),(xp,yp,zp))
        rjp = sqrt(rjp2)
        rij2 = dist2(ati.pos(),atj.pos())
        rij = sqrt(rij2)
        mu = (rip-rjp)/rij
        # Modify mu based on Becke hetero formulas (App A)
        if do_becke_hetero and ati.atno != atj.atno:
            chi = Bragg[ati.atno]/Bragg[atj.atno]
            u = (chi-1.)/(chi+1.)
            a = u/(u*u-1)
            a = min(a,0.5)
            a = max(a,-0.5)
            mu += a*(1-mu*mu)
        sprod *= sbecke(mu)
    return sprod

    
