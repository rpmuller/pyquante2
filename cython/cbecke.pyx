# cython: boundscheck=False, wraparound=False, cdivision=True
import numpy as np
cimport numpy as np
from pyquante2.grid.data import Bragg
from libc.math cimport sqrt

def becke_reweight_atoms(atoms,agrids,**kwargs):
    cdef int do_becke_hetero = kwargs.get('do_becke_hetero',True)

    cdef double[:] cBragg = np.empty((len(Bragg)),'d')
    cBragg[0] = 0
    for idx in xrange(1,len(Bragg)):
        cBragg[idx] = Bragg[idx]

    cdef long nat = len(atoms)
    cdef double[:,:] rs = np.empty((nat,3))
    for i in xrange(nat):
        for j in xrange(3):
            rs[i,j] = atoms[i].r[j]

    # Precompute norms outside of inner loop
    cdef double[:,:] rijs = np.empty((nat,nat),'d')
    cdef double diff 
    for i in xrange(nat):
        for j in xrange(nat):
            rijs[i,j] = 0.0
            for k in xrange(3):
                diff = (rs[i][k]-rs[j][k])
                rijs[i,j] += diff*diff
            rijs[i,j] = sqrt(rijs[i,j])

    npts = sum(agrid.npts for agrid in agrids)

    # Need list of all grid points of all atoms
    cdef long g_abs = 0
    cdef double[:,:] all_points = np.empty((npts,3))
    for agrid in agrids:
        for g in xrange(agrid.npts):
            for k in xrange(3):
                all_points[g_abs,k] = agrid.points[g,k]
            g_abs += 1

    cdef double[:,:] rjps = np.empty((nat,npts),'d',order="C")
    for at in xrange(nat):
        for g in xrange(npts):
            rjps[at,g] = 0.0
            for k in xrange(3):
                diff = rs[at,k] - all_points[g,k]
                rjps[at,g] += diff*diff
            rjps[at,g] = sqrt(rjps[at,g])

    cdef int[:] atnos = np.empty(nat,dtype=np.intc) 
    for iat in xrange(nat):
        atnos[iat] = atoms[iat].atno

    Ps = np.empty(nat,'d')
    g_abs = 0
    for iat,agrid in enumerate(agrids):
        for g in xrange(agrid.npts):
            for jat in xrange(nat):
                Ps[jat] = becke_atomic_grid_p(jat,g_abs,rijs,rjps,atnos,do_becke_hetero,cBragg)
            P = Ps[iat]/sum(Ps)
            agrid.points[g,3] = P*agrid.points[g,3]
            g_abs += 1
    return

cdef inline double becke_atomic_grid_p(long jat,long g,double[:,:] rijs,double[:,:] rjps,int[:] atnos,int do_becke_hetero,double[:] cBragg):
    cdef double sprod = 1
    cdef double mu
    cdef double chi
    cdef double u
    cdef double a
    for kat in xrange(rijs.shape[0]):
        if jat == kat: continue
        mu = (rjps[jat,g]-rjps[kat,g])/rijs[jat,kat]
        if do_becke_hetero != 0 and atnos[jat] != atnos[kat]:
            chi = cBragg[atnos[jat]]/cBragg[atnos[kat]]
            u = (chi-1.)/(chi+1.)
            a = u/(u*u-1)
            a = min(a,0.5)
            a = max(a,-0.5)
            mu += a*(1-mu*mu)
        sprod *= sbecke(mu)
    return sprod

cdef inline double fbecke(double x,int n=3):
    for i in xrange(n): x = pbecke(x)
    return x
cdef inline double pbecke(double x): 
    cdef double c1 = 1.5
    cdef double c2 = 0.5
    return c1*x-c2*x*x*x
cdef inline double sbecke(double x,int n=3): 
    cdef double c1 = 0.5
    cdef double c2 = 1
    return c1*(c2-fbecke(x,n))
