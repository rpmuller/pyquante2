"""\
 Implementation of Head-Gordon & Pople's scheme for electron repulsion
  integrals (ref), which, in turn, derives from Saika and Obarra's scheme.

 Routines:
 hrr performs the horizontal recursion relationships
 vrr performs the vertical recursion relationship

 This program is part of the PyQuante quantum chemistry program suite.
"""
from numpy import sqrt,exp,pi,isclose,array,zeros
from itertools import product
from pyquante2.utils import Fgamma,norm2
from pyquante2.ints.one import gaussian_product_center
import time

def vrr(xyza,norma,lmna,alphaa,
        xyzb,normb,alphab,
        xyzc,normc,lmnc,alphac,
        xyzd,normd,alphad):
    return norma*normb*normc*normd*vrr_core(xyza,lmna,alphaa,xyzb,alphab,
                                            xyzc,lmnc,alphac,xyzd,alphad)

def vrr_core(xyza,lmna,alphaa,xyzb,alphab,
             xyzc,lmnc,alphac,xyzd,alphad):

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
    zpe = zeta+eta

    rab2 = norm2(xyza-xyzb)
    Kab = sqrt(2)*pow(pi,1.25)/(alphaa+alphab)\
          *exp(-alphaa*alphab/(alphaa+alphab)*rab2)
    rcd2 = norm2(xyzc-xyzd)
    Kcd = sqrt(2)*pow(pi,1.25)/(alphac+alphad)\
          *exp(-alphac*alphad/(alphac+alphad)*rcd2)
    rpq2 = norm2(xyzp-xyzq)
    T = zeta*eta/zpe*rpq2

    mtot = sum(lmna)+sum(lmnc)

    vrr_terms = zeros((la+1,ma+1,na+1,lc+1,mc+1,nc+1,mtot+1),'d')
    for m in xrange(mtot+1):
        vrr_terms[0,0,0,0,0,0,m] = Fgamma(m,T)*Kab*Kcd/sqrt(zpe)

    for i in xrange(la):
        for m in xrange(mtot-i):
            vrr_terms[i+1,0,0, 0,0,0, m] = (px-xa)*vrr_terms[i,0,0, 0,0,0, m] +\
                                            (wx-px)*vrr_terms[i,0,0, 0,0,0, m+1]
            if i>0:
                vrr_terms[i+1,0,0, 0,0,0, m] += i/2./zeta*(
                    vrr_terms[i-1,0,0, 0,0,0, m]
                    - eta/zpe*vrr_terms[i-1,0,0, 0,0,0, m+1] )

    for i,j in product(xrange(la+1),xrange(ma)):
        for m in range(mtot-i-j):
            vrr_terms[i,j+1,0, 0,0,0, m] = (py-ya)*vrr_terms[i,j,0, 0,0,0, m] +\
                                            (wy-py)*vrr_terms[i,j,0, 0,0,0, m+1]
            if j>0:
                vrr_terms[i,j+1,0, 0,0,0, m] += j/2./zeta*(
                    vrr_terms[i,j-1,0, 0,0,0, m]
                    - eta/zpe*vrr_terms[i,j-1,0, 0,0,0, m+1] )

    for i,j,k in product(xrange(la+1),xrange(ma+1),xrange(na)):
        for m in range(mtot-i-j-k):
            vrr_terms[i,j,k+1, 0,0,0, m] = (pz-za)*vrr_terms[i,j,k, 0,0,0, m] +\
                                            (wz-pz)*vrr_terms[i,j,k, 0,0,0, m+1]
            if k>0:
                vrr_terms[i,j,k+1, 0,0,0, m] += k/2./zeta*(
                    vrr_terms[i,j,k-1, 0,0,0, m]
                    - eta/zpe*vrr_terms[i,j,k-1, 0,0,0, m+1])

    for i,j,k, q in product(xrange(la+1),xrange(ma+1),xrange(na+1),xrange(lc)):
        for m in range(mtot-i-j-k-q):
            vrr_terms[i,j,k, q+1,0,0, m] = (qx-xc)*vrr_terms[i,j,k, q,0,0, m] +\
                                            (wx-qx)*vrr_terms[i,j,k, q,0,0, m+1]
            if q>0:
                vrr_terms[i,j,k, q+1,0,0, m] += q/2./eta*(
                    vrr_terms[i,j,k, q-1,0,0, m]
                    - zeta/zpe*vrr_terms[i,j,k, q-1,0,0, m+1])
            if i>0:
                vrr_terms[i,j,k, q+1,0,0, m] += i/2./zpe*vrr_terms[i-1,j,k, q,0,0, m+1]

    for i,j,k, q,r in product(xrange(la+1),xrange(ma+1),xrange(na+1),
                              xrange(lc+1),xrange(mc)):
        for m in range(mtot-i-j-k-q-r):
            vrr_terms[i,j,k, q,r+1,0, m] = (qy-yc)*vrr_terms[i,j,k, q,r,0, m] +\
                                            (wy-qy)*vrr_terms[i,j,k, q,r,0, m+1]
            if r>0:
                vrr_terms[i,j,k, q,r+1,0, m] += r/2./eta*(
                    vrr_terms[i,j,k, q,r-1,0, m]
                    - zeta/zpe*vrr_terms[i,j,k, q,r-1,0, m+1])
            if j>0:
                vrr_terms[i,j,k, q,r+1,0, m] += j/2./zpe*vrr_terms[i,j-1,k,q,r,0,m+1]

    for i,j,k, q,r,s in product(xrange(la+1),xrange(ma+1),xrange(na+1),
                                xrange(lc+1),xrange(mc+1),xrange(nc)):
        for m in range(mtot-i-j-k-q-r-s):
            vrr_terms[i,j,k,q,r,s+1,m] = (qz-zc)*vrr_terms[i,j,k,q,r,s,m] +\
                                          (wz-qz)*vrr_terms[i,j,k,q,r,s,m+1]
            if s>0:
                vrr_terms[i,j,k,q,r,s+1,m] += s/2./eta*(
                    vrr_terms[i,j,k,q,r,s-1,m]
                    - zeta/zpe*vrr_terms[i,j,k,q,r,s-1,m+1])
            if k>0:
                vrr_terms[i,j,k,q,r,s+1,m] += k/2./zpe*vrr_terms[i,j,k-1,q,r,s,m+1]
    return vrr_terms[la,ma,na,lc,mc,nc,0]

def shell_iterator(am): # From the libint manual
    for i in range(am+1):
        for j in range(i+1):
            yield am-i,i-j,j
    return

def term0(aexpn,xyza,bexpn,xyzb,cexpn,xyzc,dexpn,xyzd,mtot):
    zeta,eta = float(aexpn+bexpn),float(cexpn+dexpn)
    zpe = zeta+eta

    rab2 = norm2(xyza-xyzb)
    Kab = sqrt(2)*pow(pi,1.25)/(aexpn+bexpn)\
          *exp(-aexpn*bexpn/(aexpn+bexpn)*rab2)
    rcd2 = norm2(xyzc-xyzd)
    Kcd = sqrt(2)*pow(pi,1.25)/(cexpn+dexpn)\
          *exp(-cexpn*dexpn/(cexpn+dexpn)*rcd2)
    rpq2 = norm2(xyzp-xyzq)
    T = zeta*eta/zpe*rpq2
    return [Fgamma(m,T)*Kab*Kcd/sqrt(zpe) for m in xrange(mtot+1)]

def copy_and_decrement(input,*args):
    import copy
    l = copy.copy(input)
    for arg in args:
        l[arg] -= 1
    if min(l) < 0:
        # Return an invalid array to trigger errors
        return [None]*len(input) 
    return l
    

def vrr_shell(aexpn,xyza,bexpn,xyzb,cexpn,xyzc,dexpn,xyzd,maxa,maxc):
    import operator,copy
    mtot = 3*maxa+3*maxb
    # Don't need all of these terms:
    #terms = zeros((maxa+1,maxa+1,maxa+1,maxc+1,maxc+1,maxc+1,mtot+1),'d')
    terms = {}
    # Do ama=amc=0 term first:
    ama,amc = 0,0
    for m,t in enumerate(term0(aexpn,xyza,bexpn,xyzb,cexpn,xyzc,mtot)):
        terms[0,0,0, 0,0,0, m] = t

    #ama>0, amc=0
    for ama in xrange(1,maxa+1):
        for aI,aJ,aK in shell_iterator(ama):
            for m in xrange(mtot-aI-aJ-aK):
                l = [aI,aJ,aK, 0,0,0]
                ind,val = max(enumerate(l),key=operator.itermgetter(1))
                reference = copy.copy(l)
                reference[ind] -= 1
                term[aI,aJ,aK, 0,0,0, m] = None
            #end m
        #end cijk
    #end amc
                
    # ama>0, amc>0
    for ama in xrange(1,maxa+1):
        for amc in xrange(1,maxc+1):
            pqac = (xp-xa,yp-ya,zp-za,xq-xc,yq-yc,zq-zc)
            wpq = (xw-xp,yw-yp,zw-zp,xw-xq,yw-yq-zw-zq)
            ze = (eta,eta,eta,zeta,zeta,zeta)
            zezpe = [zei/zpe for zei in ze]
            
            for aI,aJ,aK in shell_iterator(ama):
                for cI,cJ,cK in shell_iterator(amc):
                    l = [aI,aJ,aK, cI,cJ,cK]
                    ind,val = max(enumerate(l),key=operator.itermgetter(1))
                    cind,cval = conj_ind_val(l,ind,val)
                    r1,r2,r3,r4,r5,r6 = copy_and_decrement(l,ind)
                    s1,s2,s3,s4,s5,s6 = copy_and_decrement(l,ind,ind)
                    t1,t2,t3,t4,t5,t6 = copy_and_decrement(l,ind,cind)
                    for m in xrange(mtot-aI-aJ-aK-cI-cJ-cK):
                        terms[aI,aJ,aK, cI,cJ,cK, m] = pqac[ind]*terms[r1,r2,r3, r4,r5,r6,m] +\
                                                       wpq[ind]*terms[r1,r2,r3, r4,r5,r6,m+1]
                        if val>0:
                            terms[aI,aJ,aK, cI,cJ,cK, m] += 0.5*val/ze[ind]*(
                                terms[s1,s2,s3, s4,s5,s6,m] -\
                                zezpe[ind]*terms[s1,s2,s3, s4,s5,s6,m+1])
                        if cval>0:
                            terms[aI,aJ,aK, cI,cJ,cK, m] += 0.5*cval/ze[cind]*terms[t1,t2,t3,t4,t5,t6,m+1]
                    #end m
                #end bijk
            #end aijk
        #end amc
    #end ama
    return terms # and unpack later
    # Here's how you repack here. Seems silly to do this, if we're just going
    # to directly unpack into the integral array:
    #smaller_terms = {}
    #for ama in xrange(1,maxa+1):
    #    for amc in xrange(1,maxc+1):
    #        for aI,aJ,aK in shell_iterator(ama):
    #            for cI,cJ,cK in shell_iterator(amc):
    #                smaller_terms[aI,aJ,aK,cI,cJ,cK] = terms[aI,aJ,aK,cI,cJ,cK,0]
    #return smaller_terms
    
    
def test_vrr():
    ax=ay=az=bx=by=bz=cx=cy=cz=dx=dy=dz=0.0
    aexpn=bexpn=cexpn=dexpn=1.0
    aI=aJ=aK=0
    cI=cJ=cK=0
    M=0

    for (ax,ay,az, aI,aJ,aK, cI,cJ,cK, result) in [
            (0.,0.,0., 0,0,0, 0,0,0, 4.37335456733),
            (0.,0.,0., 1,0,0, 1,0,0, 0.182223107579),
            (0.,0.,0., 0,1,0, 0,1,0, 0.182223107579),
            (0.,0.,0., 0,0,1, 0,0,1, 0.182223107579),

            (0.,0.,0., 2,0,0, 2,0,0,  0.223223306785),
            (0.,0.,0., 0,2,0, 0,2,0,  0.223223306785),
            (0.,0.,0., 0,0,2, 0,0,2,  0.223223306785),

            (1.,2.,3., 1,0,0, 1,0,0, -5.63387712455e-06),
            (1.,2.,3., 0,1,0, 0,1,0, -0.000116463120359),
            (1.,2.,3., 0,0,1, 0,0,1, -0.000301178525749),

            (1.,2.,3., 2,0,0, 2,0,0, 0.000225033081978),
            (1.,2.,3., 0,2,0, 0,2,0, 0.000610247078796),
            (1.,2.,3., 0,0,2, 0,0,2, 0.00134278307956),

            (0.,0.,0., 1,1,0, 1,1,0, 0.0136667330685),
            (0.,0.,0., 0,1,1, 0,1,1, 0.0136667330685),
            (0.,0.,0., 1,0,1, 1,0,1, 0.0136667330685),

            (3.,2.,1., 1,1,0, 1,1,0, 5.97677147819e-05),
            (3.,2.,1., 0,1,1, 0,1,1, 1.57429039496e-06),
            (3.,2.,1., 1,0,1, 1,0,1, 4.00292836291e-06)
        ]:


        t0 = time.time()
        val1 = vrr(array((ax,ay,az)),1,(aI,aJ,aK),aexpn,
                   array((bx,by,bz)),1,bexpn,
                   array((cx,cy,cz)),1,(cI,cJ,cK),cexpn,
                   array((dx,dy,dz)),1,dexpn)
        t1 = time.time()
        val2 = vrr(array((cx,cy,cz)),1,(cI,cJ,cK),cexpn,
                   array((dx,dy,dz)),1,dexpn,
                   array((ax,ay,az)),1,(aI,aJ,aK),aexpn,
                   array((bx,by,bz)),1,bexpn)
        t2 = time.time()
        assert isclose(val1,val2)
        assert isclose(val1,result)
        print t1-t0,t2-t1
    return


if __name__ == '__main__':
    test_vrr()
    test_iteration()
