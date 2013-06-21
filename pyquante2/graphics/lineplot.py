import numpy as np
import pylab as pl

def lineplot_orb(points,orbs,bfs,doshow=False):
    zvals = [z for (x,y,z) in points]
    for orb in orbs:
        pl.plot(zvals,orb_on_points(points,orb,bfs))
    #pl.savefig('pyq_orb.png')
    if doshow:
        pl.show()
    return

def bf_on_points(points,bf): return np.array([bf(*point) for point in points],'d')
def orb_on_points(points,orb,bfs):
    f = np.zeros(len(points),'d')
    for c in orb:
        for bf in bfs:
            f = f + c*bf_on_points(points,bf)
    return f

def lineplot_bfs(points,bfs,doshow=False):
    zvals = [z for (x,y,z) in points]
    for bf in bfs:
        pl.plot(zvals,bf_on_points(points,bf))
    if doshow:
        pl.show()
    return
        

def test_plot_bfs():
    from pyquante2 import basisset,h2
    bfs = basisset(h2,'sto3g')
    zvals = np.linspace(-5,5)
    points = [(0,0,z) for z in zvals]
    lineplot_bfs(points,bfs,True)
    return

def test_plot_orbs():
    from pyquante2 import basisset,h2
    bfs = basisset(h2,'sto3g')
    orbs = np.array([[1.0,0.0],
                     [0.0,-1.0]],'d')
    
    zvals = np.linspace(-5,5)
    points = [(0,0,z) for z in zvals]
    lineplot_orb(points,orbs,bfs,True)
    return

if __name__ == '__main__':
    test_plot_orbs()

        
