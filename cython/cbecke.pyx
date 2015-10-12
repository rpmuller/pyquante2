# cython: boundscheck=False, wraparound=False, cdivision=True
import numpy as np
cimport numpy as np
from cpython cimport array
import array
from pyquante2.grid.data import Bragg

np.import_array()

def becke_reweight_atoms(atoms,agrids,**kwargs):
    cdef int do_becke_hetero = kwargs.get('do_becke_hetero',True)

    cdef double[:] cBragg = array.array('d',(0 for _ in xrange(len(Bragg))))
    cBragg[0] = 0
    for idx in xrange(1,len(Bragg)):
        cBragg[idx] = Bragg[idx]

    # Precompute norms outside of inner loop
    cdef double[:,:] rijs = np.ndarray((len(atoms),len(atoms)),'d',order="C")
    for i in xrange(len(atoms)):
        for j in xrange(len(atoms)):
            rijs[i,j] = np.linalg.norm(atoms[i].r-atoms[j].r)

    max_agrid_npts = max(agrid.npts for agrid in agrids)
    cdef double[:,:] rjps = np.ndarray((len(atoms),max_agrid_npts),'d',order="C")
    for at,agrid in enumerate(agrids):
        for p in xrange(agrid.npts):
            xyzp = agrid.points[p,:3]
            rjps[at,p] = np.linalg.norm(atoms[at].r-xyzp)

    #cdef int[:] atnos = array.array('i') 
    cdef int[:] atnos = np.empty(len(atoms),dtype=np.intc) 
    for iat in xrange(0,len(atoms)):
        atnos[iat] = atoms[iat].atno
        #atnos.append(atoms[iat].atno)

    Ps = np.empty(len(atoms),'d')
    cdef long g
    cdef long jat

    for iat,agrid in enumerate(agrids):
        for g from 0 <= g < agrid.npts:
            xyzp = agrid.points[i,:3]
            for jat from 0 <= jat < len(atoms):
                Ps[jat] = becke_atomic_grid_p(jat,g,rijs,rjps,atnos,do_becke_hetero,cBragg)
            P = Ps[iat]/sum(Ps)
            agrid.points[g,3] = P*agrid.points[g,3]
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
