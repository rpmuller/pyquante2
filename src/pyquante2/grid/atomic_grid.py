import numpy as np
from pyquante2.grid.lebedev import lebedev

class atomic_grid(object):
    def __init__(self,atom,**kwargs):
        atno,x,y,z = atom.atuple()

        if kwargs.get('radial','EulerMaclaurin') == 'Legendre':
            grid_params = LegendreGrid(atno,**kwargs)
        else:
            grid_params = EulerMaclaurinGrid(atno,**kwargs)
        self.points = []
        for rrad,wrad,nang in grid_params:
            for xang,yang,zang,wang in lebedev[nang]:
                w = wrad*wang
                self.points.append((rrad*xang+x,rrad*yang+y,rrad*zang+z,w))
        self.points = np.array(self.points,dtype=float)
        self.npts = self.points.shape[0]
        return

# The following two routines return [(ri,wi,nangi)] for nrad shells.
# The ri's are properly adjusted to go to the proper distances.
# The wi's are adjusted to only have to be multiplied by wrad from
# the lebedev shell
def EulerMaclaurinGrid(Z,**opts):
    nrad = opts.get('nrad',32)
    do_sg1 = opts.get('do_sg1',True)
    nang = opts.get('nang',194)
    radial = EulerMaclaurinRadialGrid(nrad,Z)
    if do_sg1:
        grid = [(r,w,SG1Angs(r,Z)) for r,w in radial]
    else:
        grid = [(r,w,nang) for r,w in radial]
    return grid

def LegendreGrid(Z,**kwargs):
    from pyquante2.constants import ang2bohr
    from pyquante2.grid.data import Bragg
    from pyquante2.grid.legendre import legendre
    Rmax = 0.5*Bragg[Z]*ang2bohr

    nrad = kwargs.get('nrad',32)
    fineness = kwargs.get('fineness',1)
    radial = legendre[nrad]
    grid = []
    for i in range(nrad):
        xrad,wrad = radial[i]
        rrad = BeckeRadMap(xrad,Rmax)
        dr = 2*Rmax/pow(1-xrad,2)
        vol = 4*np.pi*rrad*rrad*dr
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
    from pyquante2.grid.data import PopleRadii
    # Radial part of the Gill, Johnson, Pople SG-1 grid
    R = PopleRadii[Z]
    grid = []
    for i in range(1,nrad+1):
        # Changed to include a factor of 4pi
        #w = 2.*pow(R,3)*(nrad+1.)*pow(i,5)*pow(nrad+1-i,-7)
        w = 8.*np.pi*pow(R,3)*(nrad+1.)*pow(i,5)*pow(nrad+1-i,-7)
        r = R*i*i*pow(nrad+1-i,-2)
        grid.append((r,w))
    return grid

def SG1Angs(r,Z):
    from pyquante2.grid.data import PopleRadii
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

if __name__ == '__main__':
    import pylab
    nrad = 32
    for atno in [8,1]:
        #print('lgrid:  ')
        #print('elgrid: ',[w for r,w,n in EulerMaclaurinGrid(nrad,atno)])
        pylab.semilogy([w for r,w,n in LegendreGrid(nrad,atno)],label='L%d'%atno)
        pylab.semilogy([w for r,w,n in EulerMaclaurinGrid(nrad,atno)],label='EL%d'%atno)
    pylab.legend(loc='lower right')
    pylab.title("Radial weights for DFT Grids")
    pylab.show()
    
