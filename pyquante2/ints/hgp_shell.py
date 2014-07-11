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

def vrr_shell(aexpn,xyza,bexpn,xyzb,cexpn,xyzc,dexpn,xyzd,maxa,maxc):
    # This is wasteful. For example, for d shells, maxa=2. We would
    # compute integrals through (2,2,2) when we don't need them.
    # The shell_iterator function does a better job of looping through the
    # requisite terms, but the code required for it is completely different,
    # and thus I don't see a simple path to getting it working.
    terms = vrr_array(xyza,(maxa,maxa,maxa),aexpn,xyzb,bexpn,
                      xyzc,(maxc,maxc,maxc),cexpn,xyzd,dexpn)

    # We can reduce the size of the values returned using one of the pack_
    # functions. However, it's good to minimize this, since it gives a good idea
    # of just how wasteful the current scheme is.
    #return pack_full(terms,maxa,maxc)
    return pack_full(terms,maxa,maxc)

# Different methods of reducing the size of the integral record:
def pack_full(d,maxa,maxc,tol=False):
    newd = {}
    #for ama in xrange(maxa+1):
    #    for amc in xrange(maxc+1):
    ama,amc = maxa,maxc
    for aI,aJ,aK in shell_iterator(ama):
        for cI,cJ,cK in shell_iterator(amc):
            t = d[aI,aJ,aK,cI,cJ,cK,0]
            if not tol or abs(t) > tol:
                newd[aI,aJ,aK,cI,cJ,cK] = t
    return newd

def pack_nonzero(d,tol=1e-10):
    newd = {}
    for key,value in d.items():
        if abs(value) > tol:
            newd[key] = value
    #return {key: value for (key,value) in d if abs(value)>tol}
    return newd

def pack_m(d):
    newd = {}
    for key,value in d.items():
        if key[-1] == 0:
            newd[key[:-1]] = value
    return newd

def vrr_core(xyza,lmna,alphaa,xyzb,alphab,
             xyzc,lmnc,alphac,xyzd,alphad):
    terms = vrr_array(xyza,lmna,alphaa,xyzb,alphab,
                      xyzc,lmnc,alphac,xyzd,alphad)
    return terms[lmna[0],lmna[1],lmna[2],lmnc[0],lmnc[1],lmnc[2],0]

def vrr_array(xyza,lmna,alphaa,xyzb,alphab,
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

    vrr_terms = {}#zeros((la+1,ma+1,na+1,lc+1,mc+1,nc+1,mtot+1),'d')
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
    return vrr_terms

def shell_iterator(am): # From the libint manual
    for i in range(am+1):
        for j in range(i+1):
            yield am-i,i-j,j
    return

def term0(aexpn,xyza,bexpn,xyzb,cexpn,xyzc,dexpn,xyzd,mtot):
    zeta,eta = float(aexpn+bexpn),float(cexpn+dexpn)
    zpe = zeta+eta

    xyzp = gaussian_product_center(aexpn,xyza,bexpn,xyzb)
    xyzq = gaussian_product_center(cexpn,xyzc,dexpn,xyzd)
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

def indmax(l):
    import operator
    ind,val=max(enumerate(l),key=operator.itemgetter(1))
    return ind
# This one takes the later of equal values
# def indmax(l):
#     maxval = -1e10
#     maxind = None
#     for i,li in enumerate(l):
#         if li >= maxval:
#             maxval = li
#             maxind = i
#     return maxind

def vrr_shell_2(aexpn,xyza,bexpn,xyzb,cexpn,xyzc,dexpn,xyzd,maxa,maxc):
    # This uses a more intelligent method of looping over the am's to
    # generate the powers required
    import copy
    xyzp = gaussian_product_center(aexpn,xyza,bexpn,xyzb)
    xyzq = gaussian_product_center(cexpn,xyzc,dexpn,xyzd)
    zeta,eta = float(aexpn+bexpn),float(cexpn+dexpn)
    xyzw = gaussian_product_center(zeta,xyzp,eta,xyzq)
    zpe = zeta+eta
    mtot = 3*maxa+3*maxc
    # Don't need all of these terms:
    #terms = zeros((maxa+1,maxa+1,maxa+1,maxc+1,maxc+1,maxc+1,mtot+1),'d')
    terms = {}
    # Do ama=amc=0 term first:
    ama,amc = 0,0
    for m,t in enumerate(term0(aexpn,xyza,bexpn,xyzb,cexpn,xyzc,dexpn,xyzd,mtot)):
        terms[0,0,0, 0,0,0, m] = t
    #end m

    pqac = (xyzp[0]-xyza[0],xyzp[1]-xyza[1],xyzp[2]-xyza[2],
            xyzq[0]-xyzc[0],xyzq[1]-xyzc[1],xyzq[2]-xyzc[2])
    wpq = (xyzw[0]-xyzp[0],xyzw[1]-xyzp[1],xyzw[2]-xyzp[2],
           xyzw[0]-xyzq[0],xyzw[1]-xyzq[1],xyzw[2]-xyzq[2])
    ze = (eta,eta,eta,zeta,zeta,zeta)
    zezpe = [zei/zpe for zei in ze]
            
    #ama>0, amc=0
    for ama in xrange(1,maxa+1):
        for aI,aJ,aK in shell_iterator(ama):
            l = [aI,aJ,aK, 0,0,0]
            ind = indmax(l)
            val = l[ind]
            r1,r2,r3,r4,r5,r6 = copy_and_decrement(l,ind)
            s1,s2,s3,s4,s5,s6 = copy_and_decrement(l,ind,ind)
            for m in xrange(mtot-aI-aJ-aK):
                terms[aI,aJ,aK, 0,0,0, m] = pqac[ind]*terms[r1,r2,r3, r4,r5,r6,m] +\
                                               wpq[ind]*terms[r1,r2,r3, r4,r5,r6,m+1]
                if val>1:
                    terms[aI,aJ,aK, 0,0,0, m] += 0.5*(val-1)/ze[ind]*(
                        terms[s1,s2,s3, s4,s5,s6,m] -\
                        zezpe[ind]*terms[s1,s2,s3, s4,s5,s6,m+1])
            #end m
        #end cijk
    #end ama

    # ama=0, amc>0 here:
    for amc in xrange(1,maxc+1):
        for cI,cJ,cK in shell_iterator(amc):
            l = [0,0,0, cI,cJ,cK]
            ind = indmax(l)
            val = l[ind]
            r1,r2,r3,r4,r5,r6 = copy_and_decrement(l,ind)
            s1,s2,s3,s4,s5,s6 = copy_and_decrement(l,ind,ind)
            for m in xrange(mtot-cI-cJ-cK):
                terms[0,0,0, cI,cJ,cK, m] = pqac[ind]*terms[r1,r2,r3, r4,r5,r6,m] +\
                                            wpq[ind]*terms[r1,r2,r3, r4,r5,r6,m+1]
                if val>1:
                    terms[0,0,0, cI,cJ,cK, m] += 0.5*(val-1)/ze[ind]*(
                        terms[s1,s2,s3, s4,s5,s6,m] -\
                        zezpe[ind]*terms[s1,s2,s3, s4,s5,s6,m+1])
            #end m
        #end cIJK
    #end amc
                
    # ama>0, amc>0
    for ama in xrange(1,maxa+1):
        for amc in xrange(1,maxc+1):
            for aI,aJ,aK in shell_iterator(ama):
                for cI,cJ,cK in shell_iterator(amc):
                    l = [aI,aJ,aK, cI,cJ,cK]
                    ind = indmax(l)
                    val = l[ind]
                    cind = (ind+3)%6
                    cval = l[cind]
                    r1,r2,r3,r4,r5,r6 = copy_and_decrement(l,ind)
                    s1,s2,s3,s4,s5,s6 = copy_and_decrement(l,ind,ind)
                    t1,t2,t3,t4,t5,t6 = copy_and_decrement(l,ind,cind)
                    for m in xrange(mtot-aI-aJ-aK-cI-cJ-cK):
                        terms[aI,aJ,aK, cI,cJ,cK, m] = pqac[ind]*terms[r1,r2,r3, r4,r5,r6,m] +\
                                                       wpq[ind]*terms[r1,r2,r3, r4,r5,r6,m+1]
                        # should this be ...0.5*(val-1)/... ?
                        if val>1:
                            terms[aI,aJ,aK, cI,cJ,cK, m] += 0.5*(val-1)/ze[ind]*(
                                terms[s1,s2,s3, s4,s5,s6,m] -\
                                zezpe[ind]*terms[s1,s2,s3, s4,s5,s6,m+1])
                        #print 'c',(aI,aJ,aK,cI,cJ,cK,m),terms[aI,aJ,aK,cI,cJ,cK,m]
                        if cval>0:
                            terms[aI,aJ,aK, cI,cJ,cK, m] += 0.5*cval/zpe*terms[t1,t2,t3,t4,t5,t6,m+1]
                    #end m
                #end bijk
            #end aijk
        #end amc
    #end ama
    return pack_full(terms,maxa,maxc) # and unpack later
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
        xyza = array((ax,ay,az))
        aIJK = aI,aJ,aK
        xyzb = array((bx,by,bz))
        xyzc = array((cx,cy,cz))
        cIJK = cI,cJ,cK
        xyzd = array((dx,dy,dz))
        val1 = vrr(xyza,1,aIJK,aexpn,xyzb,1,bexpn,
                   xyzc,1,cIJK,cexpn,xyzd,1,dexpn)
        t1 = time.time()
        val2 = vrr(xyzc,1,cIJK,cexpn,xyzd,1,dexpn,
                   xyza,1,aIJK,aexpn,xyzb,1,bexpn)
        t2 = time.time()
        assert isclose(val1,val2)
        assert isclose(val1,result)
        #print (t1-t0,t2-t1)
        maxa = sum(aIJK)
        maxc = sum(cIJK)
        terms2 = vrr_shell_2(aexpn,xyza,bexpn,xyzb,cexpn,xyzc,dexpn,xyzd,maxa,maxc)
        terms = vrr_shell(aexpn,xyza,bexpn,xyzb,cexpn,xyzc,dexpn,xyzd,maxa,maxc)
        val2 = terms2[aI,aJ,aK,cI,cJ,cK]
        val3 = terms[aI,aJ,aK,cI,cJ,cK]
        print (aIJK,cIJK,val1,val2,val3)
        if abs(val1-val2)>1e-7:
            print (terms2)
        print ("")

    return


if __name__ == '__main__':
    test_vrr()
    #print list(shell_iterator(1))
    
    
