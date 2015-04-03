from pyquante2 import molecule,rhf,uhf,h2,lih,basisset
from pyquante2.graphics.maya import view_mol, view_orb, view_bonds

bfs = basisset(h2)
solver = rhf(h2,bfs)
ens = solver.converge()
orbs = solver.orbs

view_mol(h2,doshow=False)
view_bonds(h2,doshow=False)
view_orb(h2,orbs[:,1],bfs,planes=[('x',False)])
