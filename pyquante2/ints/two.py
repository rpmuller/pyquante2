from math import factorial,pi,sqrt,exp
from pyquante2.utils import dist2, Fgamma
from pyquante2.ints.one import gaussian_product_center, binomial_prefactor


def contr_coulomb(aexps,acoefs,anorms,xyza,powa,
                  bexps,bcoefs,bnorms,xyzb,powb,
                  cexps,ccoefs,cnorms,xyzc,powc,
                  dexps,dcoefs,dnorms,xyzd,powd):

    """
    Return the coulomb repulsion as in the coulomb_repulsion routine, but
    allow lists of exponents and normalization constants.
    >>> p1 = (0.,0.,0.)
    >>> p2 = (0.,0.,1.)
    >>> lmn = (0,0,0)
    >>> contr_coulomb([1.],[1.],[1.],p1,lmn,[1.],[1.],[1.],p1,lmn,[1.],[1.],[1.],p1,lmn,[1.],[1.],[1.],p1,lmn)
    4.373354581906223
    >>> contr_coulomb([1.],[1.],[1.],p1,lmn,[1.],[1.],[1.],p1,lmn,[1.],[1.],[1.],p2,lmn,[1.],[1.],[1.],p2,lmn)
    3.2661267317941505
    """
    Jij = 0.
    for i in xrange(len(aexps)):
        for j in xrange(len(bexps)):
            for k in xrange(len(cexps)):
                for l in xrange(len(dexps)):
                    incr = coulomb_repulsion(xyza,anorms[i],powa,aexps[i],
                                             xyzb,bnorms[j],powb,bexps[j],
                                             xyzc,cnorms[k],powc,cexps[k],
                                             xyzd,dnorms[l],powd,dexps[l])
                    Jij = Jij + acoefs[i]*bcoefs[j]*ccoefs[k]*dcoefs[l]*incr
    return Jij

def contract(f,a,b,c,d):
    """
    A simpler interface to a contracted coulomb integral
    >>> from pyquante2.basis.cgbf import cgbf
    >>> s = cgbf(exps=[1],coefs=[1])
    >>> contract(ERI,s,s,s,s)
    1.128379167095515
    """
    return sum(ca*cb*cc*cd*f(pa,pb,pc,pd) for (ca,pa) in a
               for (cb,pb) in b for (cc,pc) in c for (cd,pd) in d)

def ERI(a,b,c,d):
    """
    >>> from pyquante2.basis.pgbf import pgbf
    >>> s = pgbf(1)
    >>> ERI(s,s,s,s)
    1.128379167095514
    """ 
    return coulomb_repulsion(a.origin,a.norm,a.powers,a.exponent,
                             b.origin,b.norm,b.powers,b.exponent,
                             c.origin,c.norm,c.powers,c.exponent,
                             d.origin,d.norm,d.powers,d.exponent)

def coulomb_repulsion((xa,ya,za),norma,(la,ma,na),alphaa,
                      (xb,yb,zb),normb,(lb,mb,nb),alphab,
                      (xc,yc,zc),normc,(lc,mc,nc),alphac,
                      (xd,yd,zd),normd,(ld,md,nd),alphad):
    """
    Return the coulomb repulsion between four primitive gaussians a,b,c,d with the given origin
    x,y,z, normalization constants norm, angular momena l,m,n, and exponent alpha.
    >>> p1 = (0.,0.,0.)
    >>> p2 = (0.,0.,1.)
    >>> lmn = (0,0,0)
    >>> coulomb_repulsion(p1,1.,lmn,1.,p1,1.,lmn,1.,p1,1.,lmn,1.,p1,1.,lmn,1.)
    4.373354581906223
    >>> coulomb_repulsion(p1,1.,lmn,1.,p1,1.,lmn,1.,p2,1.,lmn,1.,p2,1.,lmn,1.)
    3.2661267317941505
    """

    rab2 = dist2((xa,ya,za),(xb,yb,zb))
    rcd2 = dist2((xc,yc,zc),(xd,yd,zd))
    xp,yp,zp = gaussian_product_center(alphaa,(xa,ya,za),alphab,(xb,yb,zb))
    xq,yq,zq = gaussian_product_center(alphac,(xc,yc,zc),alphad,(xd,yd,zd))
    rpq2 = dist2((xp,yp,zp),(xq,yq,zq))
    gamma1 = alphaa+alphab
    gamma2 = alphac+alphad
    delta = 0.25*(1/gamma1+1/gamma2)

    Bx = B_array(la,lb,lc,ld,xp,xa,xb,xq,xc,xd,gamma1,gamma2,delta)
    By = B_array(ma,mb,mc,md,yp,ya,yb,yq,yc,yd,gamma1,gamma2,delta)
    Bz = B_array(na,nb,nc,nd,zp,za,zb,zq,zc,zd,gamma1,gamma2,delta)

    sum = 0.
    for I in xrange(la+lb+lc+ld+1):
        for J in xrange(ma+mb+mc+md+1):
            for K in xrange(na+nb+nc+nd+1):
                sum = sum + Bx[I]*By[J]*Bz[K]*Fgamma(I+J+K,0.25*rpq2/delta)

    return 2*pow(pi,2.5)/(gamma1*gamma2*sqrt(gamma1+gamma2)) \
           *exp(-alphaa*alphab*rab2/gamma1) \
           *exp(-alphac*alphad*rcd2/gamma2)*sum*norma*normb*normc*normd

def B_term(i1,i2,r1,r2,u,l1,l2,l3,l4,Px,Ax,Bx,Qx,Cx,Dx,gamma1,gamma2,delta):
    "THO eq. 2.22"
    return fB(i1,l1,l2,Px,Ax,Bx,r1,gamma1) \
           *pow(-1,i2)*fB(i2,l3,l4,Qx,Cx,Dx,r2,gamma2) \
           *pow(-1,u)*fact_ratio2(i1+i2-2*(r1+r2),u) \
           *pow(Qx-Px,i1+i2-2*(r1+r2)-2*u) \
           /pow(delta,i1+i2-2*(r1+r2)-u)

def B_array(l1,l2,l3,l4,p,a,b,q,c,d,g1,g2,delta):
    Imax = l1+l2+l3+l4+1
    B = [0]*Imax
    for i1 in xrange(l1+l2+1):
        for i2 in xrange(l3+l4+1):
            for r1 in xrange(i1/2+1):
                for r2 in xrange(i2/2+1):
                    for u in xrange((i1+i2)/2-r1-r2+1):
                        I = i1+i2-2*(r1+r2)-u
                        B[I] = B[I] + B_term(i1,i2,r1,r2,u,l1,l2,l3,l4,
                                             p,a,b,q,c,d,g1,g2,delta)
    return B

def fB(i,l1,l2,P,A,B,r,g): return binomial_prefactor(i,l1,l2,P-A,P-B)*B0(i,r,g)
def B0(i,r,g): return fact_ratio2(i,r)*pow(4*g,r-i)
def fact_ratio2(a,b): return factorial(a)/factorial(b)/factorial(a-2*b)

def method(**kwargs):
    """
    method returns a two-electron integral method based on either kwargs or
    defaults.
    """
    return ERI # currently only one choice

if __name__ == '__main__':
    import doctest; doctest.testmod()
