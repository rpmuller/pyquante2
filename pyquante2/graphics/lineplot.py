import numpy as np
import pylab as pl

def lineplot(points,orbs,bfs,doshow=False):
    for orb in orbs:
        f = np.zeros(len(points),'d')
        for c in orb:
            for bf in bfs:
                for i,(x,y,z) in enumerate(points):
                    f[i] += c*bf(x,y,z)
        print orb
        pl.plot(points,f,label=str(tuple(orb))) # appears to be called 3x per orb???
    pl.legend()
    pl.savefig('pyq_orb.png')
    
    if doshow:
        pl.show()
    return

def plot_fake():
    from pyquante2 import basisset,h2
    bfs = basisset(h2,'sto3g')
    orbs = np.array([[1.0,0.0],
                     [0.0,-1.0]],'d')
    
    zvals = np.linspace(-5,5)
    points = [(0,0,z) for z in zvals]
    lineplot(points,orbs,bfs,True)
    return

if __name__ == '__main__':
    plot_fake()

        
