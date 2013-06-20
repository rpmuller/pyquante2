from pyquante2 import basisset,rhf,h2
from pyquante2.graphics.vtk import vtk_orbital

bfs = basisset(h2,'sto3g')
solver = rhf(h2,bfs)
ens = solver.converge()

# Note: these orbitals are not coming out symmetric. Why not??
print solver
print solver.orbs
vtk_orbital(h2,solver.orbs,bfs)

