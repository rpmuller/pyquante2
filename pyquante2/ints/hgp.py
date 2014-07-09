"""\
 Implementation of Head-Gordon & Pople's scheme for electron repulsion
  integrals (ref), which, in turn, derives from Saika and Obarra's scheme.

 Routines:
 hrr performs the horizontal recursion relationships
 vrr performs the vertical recursion relationship

 The routines in the accompanying chgp module have the same functions, but
 are written in C to be faster.

 This program is part of the PyQuante quantum chemistry program suite.
"""
from numpy import sqrt,exp,pi,isclose
from pyquante2.utils import Fgamma
from pyquante2.ints.one import gaussian_product_center

def ERI_hgp(a,b,c,d):
    """
    >>> from pyquante2.basis.pgbf import pgbf
    >>> s = pgbf(1)
    >>> isclose(ERI_hgp(s,s,s,s),1.128379)
    True
    >>> from pyquante2.basis.cgbf import cgbf
    >>> s = cgbf(exps=[1],coefs=[1])
    >>> isclose(ERI_hgp(s,s,s,s),1.128379)
    True
    >>> s2 = cgbf((0,0,1),(0,0,0),[1],[1])
    >>> isclose(ERI_hgp(s,s,s2,s2),0.842701)
    True
    """ 
    # This should be faster if I can get it to work, but having trouble passing
    # in to double *anorms, etc. 
    if a.contracted and b.contracted and c.contracted and d.contracted:
        acoefs,anorms,aexps = a.cne_list()
        bcoefs,bnorms,bexps = b.cne_list()
        ccoefs,cnorms,cexps = c.cne_list()
        dcoefs,dnorms,dexps = d.cne_list()
        return contr_hrr(
            a.origin,anorms,a.powers,aexps,acoefs,
            b.origin,bnorms,b.powers,bexps,bcoefs,
            c.origin,cnorms,c.powers,cexps,ccoefs,
            d.origin,dnorms,d.powers,dexps,dcoefs)
    if d.contracted:
        return sum(cd*ERI_hgp(pd,c,a,b) for (cd,pd) in d)
    return hrr(a.origin,a.norm,a.powers,a.exponent,
               b.origin,b.norm,b.powers,b.exponent,
               c.origin,c.norm,c.powers,c.exponent,
               d.origin,d.norm,d.powers,d.exponent)

def contr_hrr(xyza,norma,lmna,aexps,acoefs,
              xyzb,normb,lmnb,bexps,bcoefs,
              xyzc,normc,lmnc,cexps,ccoefs,
              xyzd,normd,lmnd,dexps,dcoefs):
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc
    ld,md,nd = lmnd
    xa,ya,za = xyza
    xb,yb,zb = xyzb
    xc,yc,zc = xyzc
    xd,yd,zd = xyzd
    if lb > 0:
        return (contr_hrr(xyza,norma,(la+1,ma,na),aexps,acoefs,
                    xyzb,normb,(lb-1,mb,nb),bexps,bcoefs,
                    xyzc,normc,(lc,mc,nc),cexps,ccoefs,
                    xyzd,normd,(ld,md,nd),dexps,dcoefs)
                + (xa-xb)*contr_hrr(xyza,norma,(la,ma,na),aexps,acoefs,
                              xyzb,normb,(lb-1,mb,nb),bexps,bcoefs,
                              xyzc,normc,(lc,mc,nc),cexps,ccoefs,
                              xyzd,normd,(ld,md,nd),dexps,dcoefs)
                )
    elif mb > 0:
        return (contr_hrr(xyza,norma,(la,ma+1,na),aexps,acoefs,
                    xyzb,normb,(lb,mb-1,nb),bexps,bcoefs,
                    xyzc,normc,(lc,mc,nc),cexps,ccoefs,
                    xyzd,normd,(ld,md,nd),dexps,dcoefs)
                + (ya-yb)*contr_hrr(xyza,norma,(la,ma,na),aexps,acoefs,
                              xyzb,normb,(lb,mb-1,nb),bexps,bcoefs,
                              xyzc,normc,(lc,mc,nc),cexps,ccoefs,
                              xyzd,normd,(ld,md,nd),dexps,dcoefs)
                )
    elif nb > 0:
        return (contr_hrr(xyza,norma,(la,ma,na+1),aexps,acoefs,
                    xyzb,normb,(lb,mb,nb-1),bexps,bcoefs,
                    xyzc,normc,(lc,mc,nc),cexps,ccoefs,
                    xyzd,normd,(ld,md,nd),dexps,dcoefs)
                + (za-zb)*contr_hrr(xyza,norma,(la,ma,na),aexps,acoefs,
                              xyzb,normb,(lb,mb,nb-1),bexps,bcoefs,
                              xyzc,normc,(lc,mc,nc),cexps,ccoefs,
                              xyzd,normd,(ld,md,nd),dexps,dcoefs)
                )
    elif ld > 0:
        return (contr_hrr(xyza,norma,(la,ma,na),aexps,acoefs,
                    xyzb,normb,(lb,mb,nb),bexps,bcoefs,
                    xyzc,normc,(lc+1,mc,nc),cexps,ccoefs,
                    xyzd,normd,(ld-1,md,nd),dexps,dcoefs)
                + (xc-xd)*contr_hrr(xyza,norma,(la,ma,na),aexps,acoefs,
                              xyzb,normb,(lb,mb,nb),bexps,bcoefs,
                              xyzc,normc,(lc,mc,nc),cexps,ccoefs,
                              xyzd,normd,(ld-1,md,nd),dexps,dcoefs)
                )
    elif md > 0:
        return (contr_hrr(xyza,norma,(la,ma,na),aexps,acoefs,
                    xyzb,normb,(lb,mb,nb),bexps,bcoefs,
                    xyzc,normc,(lc,mc+1,nc),cexps,ccoefs,
                    xyzd,normd,(ld,md-1,nd),dexps,dcoefs)
                + (yc-yd)*contr_hrr(xyza,norma,(la,ma,na),aexps,acoefs,
                              xyzb,normb,(lb,mb,nb),bexps,bcoefs,
                              xyzc,normc,(lc,mc,nc),cexps,ccoefs,
                              xyzd,normd,(ld,md-1,nd),dexps,dcoefs)
                )
    elif nd > 0:
        return (contr_hrr(xyza,norma,(la,ma,na),aexps,acoefs,
                    xyzb,normb,(lb,mb,nb),bexps,bcoefs,
                    xyzc,normc,(lc,mc,nc+1),cexps,ccoefs,
                    xyzd,normd,(ld,md,nd-1),dexps,dcoefs)
                + (zc-zd)*contr_hrr(xyza,norma,(la,ma,na),aexps,acoefs,
                              xyzb,normb,(lb,mb,nb),bexps,bcoefs,
                              xyzc,normc,(lc,mc,nc),cexps,ccoefs,
                              xyzd,normd,(ld,md,nd-1),dexps,dcoefs)
                )
    return contr_vrr(xyza,norma,(la,ma,na),aexps,acoefs,
                     xyzb,normb,bexps,bcoefs,
                     xyzc,normc,(lc,mc,nc),cexps,ccoefs,
                     xyzd,normd,dexps,dcoefs)

def contr_vrr(xyza,norma,lmna,aexps,acoefs,
              xyzb,normb,bexps,bcoefs,
              xyzc,normc,lmnc,cexps,ccoefs,
              xyzd,normd,dexps,dcoefs):
    la,ma,na = lmna
    lc,mc,nc = lmnc
    val = 0.
    for i in range(len(aexps)):
        for j in range(len(bexps)):
            for k in range(len(cexps)):
                for l in range(len(dexps)):
                    val = val + acoefs[i]*bcoefs[j]*ccoefs[k]*dcoefs[l]\
                          *vrr(xyza,norma[i],(la,ma,na),aexps[i],
                               xyzb,normb[j],bexps[j],
                               xyzc,normc[k],(lc,mc,nc),cexps[k],
                               xyzd,normd[l],dexps[l],0)
    return val

def hrr(xyza,norma,lmna,alphaa,
        xyzb,normb,lmnb,alphab,
        xyzc,normc,lmnc,alphac,
        xyzd,normd,lmnd,alphad):

    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc
    ld,md,nd = lmnd
    xa,ya,za = xyza
    xb,yb,zb = xyzb
    xc,yc,zc = xyzc
    xd,yd,zd = xyzd
    if lb > 0:
        return (hrr(xyza,norma,(la+1,ma,na),alphaa,
                    xyzb,normb,(lb-1,mb,nb),alphab,
                    xyzc,normc,(lc,mc,nc),alphac,
                    xyzd,normd,(ld,md,nd),alphad)
                + (xa-xb)*hrr(xyza,norma,(la,ma,na),alphaa,
                              xyzb,normb,(lb-1,mb,nb),alphab,
                              xyzc,normc,(lc,mc,nc),alphac,
                              xyzd,normd,(ld,md,nd),alphad)
                )
    elif mb > 0:
        return (hrr(xyza,norma,(la,ma+1,na),alphaa,
                    xyzb,normb,(lb,mb-1,nb),alphab,
                    xyzc,normc,(lc,mc,nc),alphac,
                    xyzd,normd,(ld,md,nd),alphad)
                + (ya-yb)*hrr(xyza,norma,(la,ma,na),alphaa,
                              xyzb,normb,(lb,mb-1,nb),alphab,
                              xyzc,normc,(lc,mc,nc),alphac,
                              xyzd,normd,(ld,md,nd),alphad)
                )
    elif nb > 0:
        return (hrr(xyza,norma,(la,ma,na+1),alphaa,
                    xyzb,normb,(lb,mb,nb-1),alphab,
                    xyzc,normc,(lc,mc,nc),alphac,
                    xyzd,normd,(ld,md,nd),alphad)
                + (za-zb)*hrr(xyza,norma,(la,ma,na),alphaa,
                              xyzb,normb,(lb,mb,nb-1),alphab,
                              xyzc,normc,(lc,mc,nc),alphac,
                              xyzd,normd,(ld,md,nd),alphad)
                )
    elif ld > 0:
        return (hrr(xyza,norma,(la,ma,na),alphaa,
                    xyzb,normb,(lb,mb,nb),alphab,
                    xyzc,normc,(lc+1,mc,nc),alphac,
                    xyzd,normd,(ld-1,md,nd),alphad)
                + (xc-xd)*hrr(xyza,norma,(la,ma,na),alphaa,
                              xyzb,normb,(lb,mb,nb),alphab,
                              xyzc,normc,(lc,mc,nc),alphac,
                              xyzd,normd,(ld-1,md,nd),alphad)
                )
    elif md > 0:
        return (hrr(xyza,norma,(la,ma,na),alphaa,
                    xyzb,normb,(lb,mb,nb),alphab,
                    xyzc,normc,(lc,mc+1,nc),alphac,
                    xyzd,normd,(ld,md-1,nd),alphad)
                + (yc-yd)*hrr(xyza,norma,(la,ma,na),alphaa,
                              xyzb,normb,(lb,mb,nb),alphab,
                              xyzc,normc,(lc,mc,nc),alphac,
                              xyzd,normd,(ld,md-1,nd),alphad)
                )
    elif nd > 0:
        return (hrr(xyza,norma,(la,ma,na),alphaa,
                    xyzb,normb,(lb,mb,nb),alphab,
                    xyzc,normc,(lc,mc,nc+1),alphac,
                    xyzd,normd,(ld,md,nd-1),alphad)
                + (zc-zd)*hrr(xyza,norma,(la,ma,na),alphaa,
                              xyzb,normb,(lb,mb,nb),alphab,
                              xyzc,normc,(lc,mc,nc),alphac,
                              xyzd,normd,(ld,md,nd-1),alphad)
                )
    return vrr(xyza,norma,(la,ma,na),alphaa,
               xyzb,normb,alphab,
               xyzc,normc,(lc,mc,nc),alphac,
               xyzd,normd,alphad,0)

def vrr(xyza,norma,lmna,alphaa,
        xyzb,normb,alphab,
        xyzc,normc,lmnc,alphac,
        xyzd,normd,alphad,M):

    la,ma,na = lmna
    lc,mc,nc = lmnc
    xa,ya,za = xyza
    xb,yb,zb = xyzb
    xc,yc,zc = xyzc
    xd,yd,zd = xyzd

    px,py,pz = xyzp = gaussian_product_center(alphaa,xyza,alphab,xyzb)
    qx,qy,qz = xyzq = gaussian_product_center(alphac,xyzc,alphad,xyzd)
    zeta,eta = float(alphaa+alphab),float(alphac+alphad)
    wx,wy,wz = xyzw = gaussian_product_center(zeta,xyzp,eta,xyzq)

    rab2 = pow(xa-xb,2) + pow(ya-yb,2) + pow(za-zb,2)
    Kab = sqrt(2)*pow(pi,1.25)/(alphaa+alphab)\
          *exp(-alphaa*alphab/(alphaa+alphab)*rab2)
    rcd2 = pow(xc-xd,2) + pow(yc-yd,2) + pow(zc-zd,2)
    Kcd = sqrt(2)*pow(pi,1.25)/(alphac+alphad)\
          *exp(-alphac*alphad/(alphac+alphad)*rcd2)
    rpq2 = pow(px-qx,2) + pow(py-qy,2) + pow(pz-qz,2)
    T = zeta*eta/(zeta+eta)*rpq2

    mtot = la+ma+na+lc+mc+nc+M

    Fgterms = [0]*(mtot+1)
    Fgterms[mtot] = Fgamma(mtot,T)
    for im in range(mtot-1,-1,-1):
        Fgterms[im]=(2.*T*Fgterms[im+1]+exp(-T))/(2.*im+1)

    # Todo: setup this as a regular array

    # Store the vrr values as a 7 dimensional array
    # vrr_terms[la,ma,na,lc,mc,nc,m]
    vrr_terms = {}
    for im in range(mtot+1):
        vrr_terms[0,0,0,0,0,0,im] = (
            norma*normb*normc*normd*Kab*Kcd/sqrt(zeta+eta)*Fgterms[im]
            )

    # Todo: use itertools.product() for the nested for loops
    for i in range(la):
        for im in range(mtot-i):
            vrr_terms[i+1,0,0, 0,0,0, im] = (
                (px-xa)*vrr_terms[i,0,0, 0,0,0, im]
                + (wx-px)*vrr_terms[i,0,0, 0,0,0, im+1]
                )
            if i:
                vrr_terms[i+1,0,0, 0,0,0, im] += (
                    i/2./zeta*( vrr_terms[i-1,0,0, 0,0,0, im]
                               - eta/(zeta+eta)*vrr_terms[i-1,0,0, 0,0,0, im+1]
                               ))

    for j in range(ma):
        for i in range(la+1):
            for im in range(mtot-i-j):
                vrr_terms[i,j+1,0, 0,0,0, im] = (
                    (py-ya)*vrr_terms[i,j,0, 0,0,0, im]
                    + (wy-py)*vrr_terms[i,j,0, 0,0,0, im+1]
                    )
                if j:
                    vrr_terms[i,j+1,0, 0,0,0, im] += (
                        j/2./zeta*(vrr_terms[i,j-1,0, 0,0,0, im]
                                  - eta/(zeta+eta)
                                  *vrr_terms[i,j-1,0, 0,0,0, im+1]
                                  ))


    for k in range(na):
        for j in range(ma+1):
            for i in range(la+1):
                for im in range(mtot-i-j-k):
                    vrr_terms[i,j,k+1, 0,0,0, im] = (
                        (pz-za)*vrr_terms[i,j,k, 0,0,0, im]
                        + (wz-pz)*vrr_terms[i,j,k, 0,0,0, im+1]
                        )
                    if k:
                        vrr_terms[i,j,k+1, 0,0,0, im] += (
                            k/2./zeta*(vrr_terms[i,j,k-1, 0,0,0, im]
                                      - eta/(zeta+eta)
                                      *vrr_terms[i,j,k-1, 0,0,0, im+1]
                                      ))

    for q in range(lc):
        for k in range(na+1):
            for j in range(ma+1):
                for i in range(la+1):
                    for im in range(mtot-i-j-k-q):
                        vrr_terms[i,j,k, q+1,0,0, im] = (
                            (qx-xc)*vrr_terms[i,j,k, q,0,0, im]
                            + (wx-qx)*vrr_terms[i,j,k, q,0,0, im+1]
                            )
                        if q:
                            vrr_terms[i,j,k, q+1,0,0, im] += (
                                q/2./eta*(vrr_terms[i,j,k, q-1,0,0, im]
                                         - zeta/(zeta+eta)
                                         *vrr_terms[i,j,k, q-1,0,0, im+1]
                                         ))
                        if i:
                            vrr_terms[i,j,k, q+1,0,0, im] += (
                                i/2./(zeta+eta)*vrr_terms[i-1,j,k, q,0,0, im+1]
                                )

    for r in range(mc):
        for q in range(lc+1):
            for k in range(na+1):
                for j in range(ma+1):
                    for i in range(la+1):
                        for im in range(mtot-i-j-k-q-r):
                            vrr_terms[i,j,k, q,r+1,0, im] = (
                                (qy-yc)*vrr_terms[i,j,k, q,r,0, im]
                                + (wy-qy)*vrr_terms[i,j,k, q,r,0, im+1]
                                )
                            if r:
                                vrr_terms[i,j,k, q,r+1,0, im] += (
                                    r/2./eta*(vrr_terms[i,j,k, q,r-1,0, im]
                                             - zeta/(zeta+eta)
                                             *vrr_terms[i,j,k, q,r-1,0, im+1]
                                             ))
                            if j:
                                vrr_terms[i,j,k, q,r+1,0, im] += (
                                    j/2./(zeta+eta)*vrr_terms[i,j-1,k,q,r,0,im+1]
                                    )

    for s in range(nc):
        for r in range(mc+1):
            for q in range(lc+1):
                for k in range(na+1):
                    for j in range(ma+1):
                        for i in range(la+1):
                            for im in range(mtot-i-j-k-q-r-s):
                                vrr_terms[i,j,k,q,r,s+1,im] = (
                                    (qz-zc)*vrr_terms[i,j,k,q,r,s,im]
                                    + (wz-qz)*vrr_terms[i,j,k,q,r,s,im+1]
                                    )
                                if s:
                                    vrr_terms[i,j,k,q,r,s+1,im] += (
                                        s/2./eta*(vrr_terms[i,j,k,q,r,s-1,im]
                                                 - zeta/(zeta+eta)
                                                 *vrr_terms[i,j,k,q,r,s-1,im+1]
                                                 ))
                                if k:
                                    vrr_terms[i,j,k,q,r,s+1,im] += (
                                        k/2./(zeta+eta)*vrr_terms[i,j,k-1,q,r,s,im+1]
                                        )
    return vrr_terms[la,ma,na,lc,mc,nc,M]

if __name__ == '__main__':
    import doctest; doctest.testmod()
