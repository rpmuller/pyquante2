import numpy as np
from pyquante2 import basisset,rhf,h2
from pyquante2.graphics.vtkplot import vtk_orbs
from pyquante2.graphics.lineplot import test_plot_orbs,test_plot_bfs,lineplot_orbs
from pyquante2.graphics.contourplot import test_contour

def test_lineplot_orbs(): return test_plot_orbs()
def test_lineplot_bfs(): return test_plot_bfs()
def contour_orb(): return test_contour(True)
    
def plot_h2_lineplot():
    bfs = basisset(h2,'6-31G')
    solver = rhf(h2,bfs)
    solver.converge()
    points = [(0,0,z) for z in np.linspace(-5,5)]
    lineplot_orbs(points,solver.orbs[:,:2],bfs,True)
    return

def plot_h2_vtk():
    bfs = basisset(h2,'sto3g')
    solver = rhf(h2,bfs)
    ens = solver.converge()

    # Note: these orbitals are not coming out symmetric. Why not??
    print(solver)
    print(solver.orbs)
    vtk_orbs(h2,solver.orbs,bfs,npts=10)

def plot_orbs():
    bfs = basisset(h2,'sto3g')
    orbs = np.array([[1.0,1.0],
                     [1.0,-1.0]],'d')
    
    vtk_orbital(h2,orbs,bfs)
    return

if __name__ == '__main__':
    #plot_h2_lineplot()
    plot_h2_vtk()

