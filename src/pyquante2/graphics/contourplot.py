#!/usr/bin/env python

import numpy as np
try:
    import pylab as pl
except ImportError:
    pass

def contourplot(plane,atoms,orb,bfs,val=0,npts=50,doshow=False,
                title="Contour plot of pyquante orbital"):
    if plane == 'xy':
        return contour_xy(atoms,orb,bfs,z=val,npts=npts,doshow=doshow,title=title)
    elif plane == 'yz':
        return contour_yz(atoms,orb,bfs,x=val,npts=npts,doshow=doshow,title=title)
    elif plane == 'xz':
        return contour_xz(atoms,orb,bfs,y=val,npts=npts,doshow=doshow,title=title)
    else:
        raise Exception("Unknown contour plot plane %s" % plane)
    return
        

def contour_xy(atoms,orb,bfs,z=0,npts=50,doshow=False,
               title="Contour plot of pyquante orbital"):
    bbox = atoms.bbox()
    x = np.linspace(bbox[0],bbox[1],npts)
    y = np.linspace(bbox[2],bbox[3],npts)
    f = np.zeros((npts,npts),'d')
    X,Y = np.meshgrid(x,y)
    for c,bf in zip(orb,bfs):
        f += c*bf(X,Y,z) # Can mesh bfs with arrays of points
    cp=pl.contour(X,Y,f,colors='k')
    pl.title(title)
    pl.clabel(cp,inline=1,fontsize=8)
    if doshow: pl.show()
    return

def contour_yz(atoms,orb,bfs,x=0,npts=50,doshow=False,
               title="Contour plot of pyquante orbital"):
    bbox = atoms.bbox()
    y = np.linspace(bbox[2],bbox[3],npts)
    z = np.linspace(bbox[4],bbox[5],npts)
    f = np.zeros((npts,npts),'d')
    Y,Z = np.meshgrid(y,z)
    for c,bf in zip(orb,bfs):
        f += c*bf(x,Y,Z) # Can mesh bfs with arrays of points
    #pl.contour(Y,Z,f) # colored
    cp = pl.contour(Y,Z,f,colors='k') # b&w
    pl.title(title)
    pl.clabel(cp,inline=1,fontsize=8)
    if doshow: pl.show()
    return

def contour_xz(atoms,orb,bfs,y=0,npts=50,doshow=False,
               title="Contour plot of pyquante orbital"):
    bbox = atoms.bbox()
    x = np.linspace(bbox[0],bbox[1],npts)
    z = np.linspace(bbox[4],bbox[5],npts)
    f = np.zeros((npts,npts),'d')
    X,Z = np.meshgrid(x,z)
    for c,bf in zip(orb,bfs):
        f += c*bf(X,y,Z) # Can mesh bfs with arrays of points
    cp = pl.contour(X,Z,f,colors='k') # b&w
    pl.title(title)
    pl.clabel(cp,inline=1,fontsize=8)
    if doshow: pl.show()
    return

def test_contour(doshow=False):
    from pyquante2 import h2,basisset
    contourplot('xy',h2,np.array([1,-1],'d'),basisset(h2),val=1,doshow=doshow)
    return

if __name__ == '__main__':
    import doctest; doctest.testmod()
    test_contour(True)
