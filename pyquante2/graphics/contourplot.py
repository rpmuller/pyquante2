#!/usr/bin/env python

import numpy as np
import pylab as pl

# Couldn't figure out a cute way to do the general case,
#  so I did the special case. 
def contour_xy(atoms,orb,bfs,z=0,npts=50,doshow=False):
    bbox = atoms.bbox()
    x = np.linspace(bbox[0],bbox[1],npts)
    y = np.linspace(bbox[2],bbox[3],npts)
    f = np.zeros((npts,npts),'d')
    for i,xi in enumerate(x):
        for j,yj in enumerate(y):
            for c,bf in zip(orb,bfs):
                f[i,j] += c*bf(xi,yj,z)
    X,Y = np.meshgrid(x,y)
    pl.contour(X,Y,f)
    if doshow: pl.show()
    return

def contour_yz(atoms,orb,bfs,x=0,npts=50,doshow=False):
    bbox = atoms.bbox()
    y = np.linspace(bbox[2],bbox[3],npts)
    z = np.linspace(bbox[4],bbox[5],npts)
    f = np.zeros((npts,npts),'d')
    for i,yi in enumerate(y):
        for j,zj in enumerate(z):
            for c,bf in zip(orb,bfs):
                f[i,j] += c*bf(x,yi,zj)
    Y,Z = np.meshgrid(y,z)
    #pl.contour(Y,Z,f) # colored
    cp = pl.contour(Y,Z,f,colors='k') # b&w
    pl.clabel(cp,inline=1,fontsize=8)
    if doshow: pl.show()
    return

def test_contour(doshow=False):
    from pyquante2 import h2,basisset
    contour_yz(h2,np.array([1,-1],'d'),basisset(h2),doshow=doshow)
    return

if __name__ == '__main__':
    import doctest; doctest.testmod()
    test_contour(True)
