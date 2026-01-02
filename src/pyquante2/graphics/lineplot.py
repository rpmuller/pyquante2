import numpy as np
try:
    import pylab as pl
except ImportError:
    pass

def lineplot_orbs(points,orbs,bfs,doshow=False,
               title="Line plot of pyquante orbital"):
    bfmesh = bfs.mesh(points)
    orbmesh = np.dot(bfmesh,orbs)
    ng,norb = orbmesh.shape
    for i in range(norb):
        pl.plot(orbmesh[:,i])
    pl.title(title)
    if doshow:
        pl.show()
    return

def lineplot_bfs(points,bfs,doshow=False,
                 title="Line plot of pyquante bfs"):
    bfmesh = bfs.mesh(points)
    ng,nbf = bfmesh.shape
    for i in range(nbf):
        pl.plot(bfmesh[:,i])
    pl.title(title)
    if doshow:
        pl.show()
    return

def line(origin,destination,npts=50):
    """
    Create a 1d line from the 3-tuple origin to the 3-tuple destination
    Yet another solution from @unutbu at stackoverflow.
    http://stackoverflow.com/questions/17396164/numpythonic-way-to-make-3d-meshes-for-line-plotting
    """
    return np.column_stack((np.linspace(o,d,npts) for o,d in zip(origin,destination)))

def test_plot_bfs():
    from pyquante2 import basisset,h2
    bfs = basisset(h2,'sto3g')
    points = line((0,0,-5),(0,0,5))
    lineplot_bfs(points,bfs,True)
    return

def test_plot_orbs():
    from pyquante2 import basisset,h2
    bfs = basisset(h2,'sto3g')
    orbs = np.array([[1.0,1.0],
                     [1.0,-1.0]],'d')
    
    zvals = np.linspace(-5,5)
    points = line((0,0,-5),(0,0,5))
    lineplot_orbs(points,orbs,bfs,True)
    return

if __name__ == '__main__':
    test_plot_orbs()
    #test_plot_bfs()

        
