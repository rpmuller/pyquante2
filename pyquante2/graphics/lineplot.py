import numpy as np
import pylab as pl

def lineplot_orbs(points,orbs,bfs,doshow=False,
               title="Line plot of pyquante orbital"):
    zvals = [z for (x,y,z) in points]
    for orb in orbs.T: # Transpose makes interations work
        pl.plot(zvals,orb_on_points(points,orb,bfs))
    #pl.savefig('pyq_orb.png')
    pl.title(title)
    if doshow:
        pl.show()
    return

def bf_on_points(points,bf):
    return np.array([bf(point) for point in points],'d')

def orb_on_points(points,orb,bfs):
    f = np.zeros(len(points),'d')
    for c,bf in zip(orb,bfs):
        f = f + c*bf_on_points(points,bf)
    return f

def lineplot_bfs(points,bfs,doshow=False,
               title="Line plot of pyquante bfs"):
    zvals = [z for (x,y,z) in points]
    for bf in bfs:
        pl.plot(zvals,bf_on_points(points,bf))
    pl.title(title)
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
    orbs = np.array([[1.0,1.0],
                     [1.0,-1.0]],'d')
    
    zvals = np.linspace(-5,5)
    points = [(0,0,z) for z in zvals]
    lineplot_orbs(points,orbs,bfs,True)
    return

def test_bf_eval():
    from pyquante2.utils import norm2
    from pyquante2 import basisset,h2
    bfs = basisset(h2,'sto3g')
    zvals = np.linspace(-5,5,10)
    points = np.array([(0,0,z) for z in zvals])
    #bfamps = bf_on_points(points,bfs[0])
    #print bfamps
    d = points - np.array((1,2,3),'d')
    print d
    print [np.dot(d[i,:],d[i,:]) for i in range(d.shape[0])]
    print np.dot(d.T,d)
    print np.diag(np.dot(d,d.T))
    print np.sum(d*d,axis=1)
    xyz = np.array((1,2,3),'d')
    print np.sum(xyz*xyz)
    xyz.dims
    #bfs[0](points)
    

if __name__ == '__main__':
    #test_plot_orbs()
    test_bf_eval()

        
