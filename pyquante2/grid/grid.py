"""
The DFT grids are a little different in pyquante2 from pyquante1. Here we are
only storing the points and the weights, and we will use other data objects to
store, say, the density, the basis functions, or the various gradients at each
point.
"""
import numpy as np

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

class atomic_grid:
    def __init__(self,atom,**kwargs):
        atno,x,y,z,Z = atom.atuple()
        nrad = kwargs.get('nrad',32)
        fineness = kwargs.get('fineness',1)
        
        if kwargs.get('radial','EulerMaclaurin') == 'Legendre':
            grid_params = LegendreGrid(nrad,atno,**kwargs)
        else:
            grid_params = EulerMaclaurinGrid(nrad,atno,**kwargs)
        self.points = []
        for rrad,wrad,nang in grid_params:
            for xang,yang,zang,wang in Lebedev[nang]:
                w = wrad*wang
                self.points.append((rrad*xand+x,rrad*yang+y,rrad*zang+z,w))
        self.points = np.array(self.points,dtype=float)
        return

# The following two routines return [(ri,wi,nangi)] for nrad shells.
# The ri's are properly adjusted to go to the proper distances.
# The wi's are adjusted to only have to be multiplied by wrad from
# the lebedev shell
def EulerMaclaurinGrid(nrad,Z,**opts):
    do_sg1 = opts.get('do_sg1',True)
    nang = opts.get('nang',194)
    radial = EulerMaclaurinRadialGrid(nrad,Z)
    if do_sg1:
        grid = [(r,w,SG1Angs(r,Z)) for r,w in radial]
    else:
        grid = [(r,w,nang) for r,w in radial]
    return grid

def LegendreGrid(nrad,Z,fineness):
    from pyquante2.constants import ang2bohr
    Rmax = 0.5*Bragg[Z]*ang2bohr

    radial = Legendre[nrad]
    grid = []
    for i in range(nrad):
        xrad,wrad = radial[i]
        rrad = BeckeRadMap(xrad,Rmax)
        dr = 2*Rmax/pow(1-xrad,2)
        vol = 4*pi*rrad*rrad*dr
        nangpts = ang_mesh(float(i+1)/nrad,fineness)
        grid.append((rrad,wrad*vol,nangpts))
    return grid
    
def BeckeRadMap(x,Rmax):
    return Rmax*(1.0+x)/(1.0-x)

def ang_mesh(frac,fineness,alevs = None):
    """\
    Determine the number of points in the angular mesh based on
    the fraction of the total radial grid index frac c (0,1).

    You can optionally pass in the number of points for
    the 5 different regions
    """
    if not alevs:
        ang_levels = [
            [ 6, 14, 26, 26, 14], # Coarse
            [ 50, 50,110, 50, 26], # Medium
            [ 50,110,194,110, 50], # Fine
            [194,194,194,194,194]  # ultrafine
            ]
        alevs = ang_levels[fineness]
    nang = alevs[0]
    if frac > 0.4: nang = alevs[1]
    if frac > 0.5: nang = alevs[2]
    if frac > 0.7: nang = alevs[3]
    if frac > 0.8: nang = alevs[4]
    return nang

def EulerMaclaurinRadialGrid(nrad,Z):
    # Radial part of the Gill, Johnson, Pople SG-1 grid
    R = PopleRadii[Z]
    grid = []
    for i in range(1,nrad+1):
        # Changed to include a factor of 4pi
        #w = 2.*pow(R,3)*(nrad+1.)*pow(i,5)*pow(nrad+1-i,-7)
        w = 8.*pi*pow(R,3)*(nrad+1.)*pow(i,5)*pow(nrad+1-i,-7)
        r = R*i*i*pow(nrad+1-i,-2)
        grid.append((r,w))
    return grid

def SG1Angs(r,Z):
    # Gill, Johnson, Pople rules for SG-1 angular densities
    R = PopleRadii[Z]
    if Z in range(1,3): # H-He
        alphas = [0.25,0.5,1.0,4.5]
    elif Z in range(3,11): # Li-Ne
        alphas = [0.1667, 0.500, 0.900, 3.5]
    else: # only fit for Na-Ar
        alphas = [0.1,0.4,0.8,2.5]

    if r < alphas[0]*R: return 6
    elif r < alphas[1]*R: return 38
    elif r < alphas[2]*R: return 86
    elif r < alphas[3]*R: return 194
    return 86

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

    
