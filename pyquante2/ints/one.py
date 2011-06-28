from math import pi,exp,floor,factorial
from pyquante2.utils import dist2, binomial, fact2, Fgamma

# Notes:
# The versions S,T,V include the normalization constants
# The version overlap,kinetic,nuclear_attraction do not.
# This is so, for example, the kinetic routines can call the potential routines
#  without the normalization constants getting in the way.

def contract(f,a,b):
    """
    Can be used to evaluate S,T,V over contracted basis functions.
    """
    return sum(ca*cb*S(pa,pb) for (ca,pa) in a for (cb,pb) in b)

def S(a,b):
    """
    The simple interface to the overlap function, using only primitive basis functions as the arguments.
    >>> from pyquante2.basis.pgbf import pgbf
    >>> s = pgbf(1)
    >>> round(S(s,s),10)
    1.0
    """
    return a.norm*b.norm*overlap(a.exponent,a.powers,a.origin,b.exponent,b.powers,b.origin)

def T(a,b):
    """
    Simple interface to the kinetic function.
    >>> from pyquante2.basis.pgbf import pgbf
    >>> s = pgbf(1)
    >>> round(T(s,s),10)
    1.5
    """
    return a.norm*b.norm*kinetic(a.exponent,a.powers,a.origin,b.exponent,b.powers,b.origin)

def V(a,b,C):
    """
    Simple interface to the nuclear attraction function.
    >>> from pyquante2.basis.pgbf import pgbf
    >>> s = pgbf(1)
    >>> round(V(s,s,(0,0,0)),10)
    -1.5957691216
    """
    return a.norm*b.norm*nuclear_attraction(a.exponent,a.powers,a.origin,
                                            b.exponent,b.powers,b.origin,C)

def overlap(alpha1,(l1,m1,n1),A,alpha2,(l2,m2,n2),B):
    """
    Full form of the overlap integral. Taken from THO eq. 2.12
    >>> round(overlap(1,(0,0,0),(0,0,0),1,(0,0,0),(0,0,0)),10)
    1.9687012432
    """
    rab2 = dist2(A,B)
    gamma = alpha1+alpha2
    P = gaussian_product_center(alpha1,A,alpha2,B)

    pre = pow(pi/gamma,1.5)*exp(-alpha1*alpha2*rab2/gamma)

    wx = overlap1d(l1,l2,P[0]-A[0],P[0]-B[0],gamma)
    wy = overlap1d(m1,m2,P[1]-A[1],P[1]-B[1],gamma)
    wz = overlap1d(n1,n2,P[2]-A[2],P[2]-B[2],gamma)
    return pre*wx*wy*wz

def overlap1d(l1,l2,PAx,PBx,gamma):
    """
    The one-dimensional component of the overlap integral. Taken from THO eq. 2.12
    >>> round(overlap1d(0,0,0,0,1),10)
    1.0
    """
    total = 0
    for i in xrange(1+int(floor(0.5*(l1+l2)))):
        total += binomial_prefactor(2*i,l1,l2,PAx,PBx)* \
                 fact2(2*i-1)/pow(2*gamma,i)
    return total

def gaussian_product_center(alpha1,A,alpha2,B):
    """
    The center of the Gaussian resulting from the product of two Gaussians:
    >>> gaussian_product_center(1,(0,0,0),1,(0,0,0))
    (0, 0, 0)
    """
    gamma = alpha1+alpha2
    return (alpha1*A[0]+alpha2*B[0])/gamma,\
           (alpha1*A[1]+alpha2*B[1])/gamma,\
           (alpha1*A[2]+alpha2*B[2])/gamma

def binomial_prefactor(s,ia,ib,xpa,xpb):
    """
    The integral prefactor containing the binomial coefficients from Augspurger and Dykstra.
    >>> binomial_prefactor(0,0,0,0,0)
    1
    """
    total= 0
    for t in xrange(s+1):
        if s-ia <= t <= ib:
            total +=  binomial(ia,s-t)*binomial(ib,t)* \
                     pow(xpa,ia-s+t)*pow(xpb,ib-t)
    return total

def kinetic(alpha1,(l1,m1,n1),A,alpha2,(l2,m2,n2),B):
    """
    The full form of the kinetic energy integral
    >>> round(kinetic(1,(0,0,0),(0,0,0),1,(0,0,0),(0,0,0)),10)
    5.9061037296
    """
    term0 = alpha2*(2*(l2+m2+n2)+3)*\
            overlap(alpha1,(l1,m1,n1),A,\
                           alpha2,(l2,m2,n2),B)
    term1 = -2*pow(alpha2,2)*\
            (overlap(alpha1,(l1,m1,n1),A,
                            alpha2,(l2+2,m2,n2),B)
             + overlap(alpha1,(l1,m1,n1),A,
                              alpha2,(l2,m2+2,n2),B)
             + overlap(alpha1,(l1,m1,n1),A,
                              alpha2,(l2,m2,n2+2),B))
    term2 = -0.5*(l2*(l2-1)*overlap(alpha1,(l1,m1,n1),A,
                                           alpha2,(l2-2,m2,n2),B) +
                  m2*(m2-1)*overlap(alpha1,(l1,m1,n1),A,
                                           alpha2,(l2,m2-2,n2),B) +
                  n2*(n2-1)*overlap(alpha1,(l1,m1,n1),A,
                                           alpha2,(l2,m2,n2-2),B))
    return term0+term1+term2

def nuclear_attraction(alpha1,(l1,m1,n1),A,alpha2,(l2,m2,n2),B,C):
    """
    Full form of the nuclear attraction integral
    >>> round(nuclear_attraction(1,(0,0,0),(0,0,0),1,(0,0,0),(0,0,0),(0,0,0)),10)
    -3.1415926536
    """
    gamma = alpha1+alpha2

    P = gaussian_product_center(alpha1,A,alpha2,B)
    rab2 = dist2(A,B)
    rcp2 = dist2(C,P)

    Ax = A_array(l1,l2,P[0]-A[0],P[0]-B[0],P[0]-C[0],gamma)
    Ay = A_array(m1,m2,P[1]-A[1],P[1]-B[1],P[1]-C[1],gamma)
    Az = A_array(n1,n2,P[2]-A[2],P[2]-B[2],P[2]-C[2],gamma)

    total = 0.
    for I in xrange(l1+l2+1):
        for J in xrange(m1+m2+1):
            for K in xrange(n1+n2+1):
                total += Ax[I]*Ay[J]*Az[K]*Fgamma(I+J+K,rcp2*gamma)
                
    return -2*pi/gamma*exp(-alpha1*alpha2*rab2/gamma)*total

def A_term(i,r,u,l1,l2,PAx,PBx,CPx,gamma):
    "THO eq. 2.18"
    return pow(-1,i)*binomial_prefactor(i,l1,l2,PAx,PBx)*\
           pow(-1,u)*factorial(i)*pow(CPx,i-2*r-2*u)*\
           pow(0.25/gamma,r+u)/factorial(r)/factorial(u)/factorial(i-2*r-2*u)

def A_array(l1,l2,PA,PB,CP,g):
    "THO eq. 2.18 and 3.1"
    Imax = l1+l2+1
    A = [0]*Imax
    for i in xrange(Imax):
        for r in xrange(int(floor(i/2)+1)):
            for u in xrange(int(floor((i-2*r)/2)+1)):
                I = i-2*r-u
                A[I] = A[I] + A_term(i,r,u,l1,l2,PA,PB,CP,g)
    return A

if __name__ == '__main__':
    import doctest
    doctest.testmod()

