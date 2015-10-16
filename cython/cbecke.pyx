# cython: boundscheck=False, wraparound=False, cdivision=True
import numpy as np
cimport numpy as np
from pyquante2.grid.data import Bragg
import time
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
            diff = (rs[i][0]-rs[j][0])
            rijs[i,j] += diff*diff
            diff = (rs[i][1]-rs[j][1])
            rijs[i,j] += diff*diff
            diff = (rs[i][2]-rs[j][2])
            rijs[i,j] += diff*diff
            rijs[i,j] = sqrt(rijs[i,j])

    npts = sum(agrid.npts for agrid in agrids)

    t0 = time.time()
    # Need list of all grid points of all atoms
    cdef long g
    cdef long g_abs = 0
    cdef double[:,:] all_points = np.empty((npts,3))
    for agrid in agrids:
        for g in xrange(agrid.npts):
            all_points[g_abs,0] = agrid.points[g,0]
            all_points[g_abs,1] = agrid.points[g,1]
            all_points[g_abs,2] = agrid.points[g,2]
            g_abs += 1
    t1 = time.time()
    print "0:"
    print t1 - t0

    cdef double[:,:] rjps = np.empty((nat,npts),'d',order="C")
    for at in xrange(nat):
        for g in xrange(npts):
            rjps[at,g] = 0.0
            for k in xrange(3):
                diff = rs[at,k] - all_points[g,k]
                rjps[at,g] += diff*diff
            rjps[at,g] = sqrt(rjps[at,g])
    t2 = time.time()            
    print "1: "
    print t2 - t1

    cdef int[:] atnos = np.empty(nat,dtype=np.intc) 
    for iat in xrange(nat):
        atnos[iat] = atoms[iat].atno

    Ps = np.empty(nat,'d')
    cdef long jat
    g_abs = 0
    for iat,agrid in enumerate(agrids):
        for g from 0 <= g < agrid.npts:
            for jat from 0 <= jat < nat:
                Ps[jat] = becke_atomic_grid_p(jat,g_abs,rijs,rjps,atnos,do_becke_hetero,cBragg)
            P = Ps[iat]/sum(Ps)
            agrid.points[g,3] = P*agrid.points[g,3]
            g_abs += 1
    return

cdef inline double becke_atomic_grid_p(long jat,long g,double[:,:] rijs,double[:,:] rjps,int[:] atnos,int do_becke_hetero,double[:] cBragg):
    cdef double sprod = 1
    cdef long kat
    cdef double mu
    cdef double chi
    cdef double u
    cdef double a
    for kat from 0 <= kat < rijs.shape[0]:
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
    cdef int i
    for i from 0 <= i < n: x = pbecke(x)
    return x
cdef inline double pbecke(double x): 
    cdef double c1 = 1.5
    cdef double c2 = 0.5
    return c1*x-c2*x*x*x
cdef inline double sbecke(double x,int n=3): 
    cdef double c1 = 0.5
    cdef double c2 = 1
    return c1*(c2-fbecke(x,n))
