import numpy as np

from pyquante2 import basisset
from pyquante2.grid.grid import grid
from pyquante2.utils import dmat
from pyquante2.geo.samples import he
from pyquante2.graphics.maya import view_dft_density
#from mayavi import mlab
#from scipy.interpolate import griddata

bfs = basisset(he)
he_grid = grid(he)
he_grid.setbfamps(bfs)
nbf = len(bfs)
C = np.eye(nbf)
D = dmat(C,he.nocc())
view_dft_density(he_grid,D,(-1,1,-1,1,-1,1))
